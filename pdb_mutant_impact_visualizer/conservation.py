
"""EBI BLAST+ and Clustal Omega conservation analysis with explicit coordinates.

This reports alignment-column observations (coverage, residue frequencies and
Shannon entropy). It intentionally does not relabel a home-made score as a
ConSurf grade or as a validated functional-impact prediction.
"""

from __future__ import annotations

import math
import re
import time
from collections import Counter
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from typing import Any

import requests

from .models import Mutation

EBI_BASE = "https://www.ebi.ac.uk/Tools/services/rest"
SIFTS_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
DBFETCH_URL = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch"

GET_MAX_ATTEMPTS = 5
GET_MAX_BACKOFF_SECONDS = 30.0
GET_MAX_RETRY_AFTER_SECONDS = 60.0


class ConservationError(RuntimeError):
    """Conservation failure with any submitted EBI job IDs retained."""

    def __init__(
        self,
        message: str,
        *,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.metadata: dict[str, Any] = dict(metadata or {})

    def __str__(self) -> str:
        message = super().__str__()

        job_labels = (
            ("blast_job_id", "BLAST"),
            ("clustalo_job_id", "Clustal Omega"),
        )
        submitted = [
            f"{label}={self.metadata[key]}"
            for key, label in job_labels
            if self.metadata.get(key)
        ]

        if submitted:
            return (
                f"{message} Submitted EBI job IDs: "
                + "; ".join(submitted)
            )

        return message


def _retry_delay(
    response: requests.Response | None,
    attempt: int,
) -> float:
    """Calculate a bounded delay, respecting Retry-After when available."""
    if response is not None:
        retry_after = getattr(response, "headers", {}).get("Retry-After")

        if retry_after:
            try:
                return min(
                    GET_MAX_RETRY_AFTER_SECONDS,
                    max(0.0, float(retry_after)),
                )
            except (TypeError, ValueError):
                try:
                    retry_time = parsedate_to_datetime(retry_after)

                    if retry_time.tzinfo is None:
                        retry_time = retry_time.replace(tzinfo=timezone.utc)

                    delay = (
                        retry_time - datetime.now(timezone.utc)
                    ).total_seconds()

                    return min(
                        GET_MAX_RETRY_AFTER_SECONDS,
                        max(0.0, delay),
                    )
                except (TypeError, ValueError, OverflowError):
                    pass

    return min(GET_MAX_BACKOFF_SECONDS, 2.0 ** attempt)


def _is_nested_rate_limit(response: requests.Response) -> bool:
    """Detect an upstream rate-limit error wrapped in another HTTP response."""
    body = getattr(response, "text", "")[:5000].lower()

    indicators = (
        "rate_limit_exceeded",
        "http 429",
        "status code 429",
        "maximum 10 requests/second",
        "too many requests",
    )

    return any(indicator in body for indicator in indicators)


def _get_with_retry(
    url: str,
    session: requests.Session,
    *,
    params: dict[str, Any] | None = None,
    timeout: tuple[int, int] = (15, 60),
    max_attempts: int = GET_MAX_ATTEMPTS,
) -> requests.Response:
    """Perform a GET with bounded retries for transient service failures.

    Only GET requests use this helper. Job submissions use POST and must
    never be automatically repeated after an ambiguous submission failure.
    """
    if max_attempts < 1:
        raise ValueError("max_attempts must be at least 1.")

    for attempt in range(max_attempts):
        response: requests.Response | None = None

        try:
            if params is None:
                response = session.get(url, timeout=timeout)
            else:
                response = session.get(
                    url,
                    params=params,
                    timeout=timeout,
                )

            status_code = getattr(response, "status_code", 200)
            response_text = getattr(response, "text", "")

            nested_rate_limit = _is_nested_rate_limit(response)

            retryable = (
                status_code == 429
                or status_code >= 500
                or nested_rate_limit
            )

            if not retryable:
                response.raise_for_status()
                return response

            if attempt == max_attempts - 1:
                detail = response_text[:1500]
                raise requests.HTTPError(
                    f"GET request failed after {max_attempts} attempts. "
                    f"HTTP {status_code}: {detail}",
                    response=response,
                )

        except (requests.Timeout, requests.ConnectionError):
            if attempt == max_attempts - 1:
                raise

        if attempt < max_attempts - 1:
            time.sleep(_retry_delay(response, attempt))

    raise ConservationError(
        f"GET request failed after {max_attempts} attempts: {url}"
    )


def _get_json(url: str, session: requests.Session) -> dict[str, Any]:
    try:
        response = _get_with_retry(url, session)
        payload = response.json()

        if not isinstance(payload, dict):
            raise ConservationError(
                f"Expected a JSON object from {url}."
            )

        return payload

    except (requests.RequestException, ValueError) as exc:
        raise ConservationError(
            f"Could not retrieve mapping data from {url}: {exc}"
        ) from exc


def map_pdb_residue_to_uniprot(
    pdb_id: str,
    chain: str,
    residue_number: int,
    insertion_code: str = "",
    session: requests.Session | None = None,
) -> tuple[str, int]:
    """Resolve PDB author numbering through PDBe SIFTS residue mapping."""
    client = session or requests.Session()

    data = _get_json(
        SIFTS_URL.format(pdb_id=pdb_id.lower()),
        client,
    )

    entry = data.get(pdb_id.lower(), {})
    matches: list[tuple[str, int]] = []

    for accession, record in entry.get("UniProt", {}).items():
        for mapping in record.get("mappings", []):
            if str(mapping.get("chain_id")) != str(chain):
                continue

            start = mapping.get("start", {})
            end = mapping.get("end", {})

            try:
                start_num = int(
                    start.get(
                        "author_residue_number",
                        start.get("residue_number"),
                    )
                )
                end_num = int(
                    end.get(
                        "author_residue_number",
                        end.get("residue_number"),
                    )
                )
                unp_start = int(mapping["unp_start"])
                unp_end = int(mapping["unp_end"])

            except (TypeError, ValueError, KeyError):
                continue

            if not (start_num <= residue_number <= end_num):
                continue

            if insertion_code:
                if (
                    str(start.get("author_insertion_code", "")).strip()
                    == insertion_code
                ):
                    if residue_number != start_num:
                        continue

                elif (
                    str(end.get("author_insertion_code", "")).strip()
                    == insertion_code
                ):
                    if residue_number != end_num:
                        continue

                else:
                    continue

            if end_num - start_num != unp_end - unp_start:
                continue

            matches.append(
                (
                    accession,
                    unp_start + residue_number - start_num,
                )
            )

    matches = list(dict.fromkeys(matches))

    if len(matches) != 1:
        raise ConservationError(
            f"Could not map {pdb_id} chain {chain} residue "
            f"{residue_number}{insertion_code} to one unambiguous "
            "UniProt position through SIFTS. Supply an explicit "
            "UniProt accession and position if you have a curated mapping."
        )

    return matches[0]


def fetch_uniprot_sequence(
    accession: str,
    session: requests.Session | None = None,
) -> tuple[str, str]:
    if not re.fullmatch(r"[A-Z0-9_-]+", accession, flags=re.I):
        raise ConservationError("Invalid UniProt accession format.")

    client = session or requests.Session()

    try:
        response = _get_with_retry(
            UNIPROT_FASTA.format(accession=accession),
            client,
        )
    except requests.RequestException as exc:
        raise ConservationError(
            f"Could not fetch UniProt sequence {accession}: {exc}"
        ) from exc

    entries = parse_fasta(response.text)

    if len(entries) != 1 or not entries[0][1]:
        raise ConservationError(
            f"UniProt did not return one usable sequence for {accession}."
        )

    return entries[0]


def parse_fasta(text: str) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []

    for line in text.splitlines():
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(sequence).upper()))

            header = line[1:]
            sequence = []

        elif header is None:
            raise ConservationError(
                "Malformed FASTA: sequence text appeared before a header."
            )

        else:
            sequence.append(line)

    if header is not None:
        entries.append((header, "".join(sequence).upper()))

    return entries


def _submit(
    tool: str,
    email: str,
    session: requests.Session,
    **params: Any,
) -> str:
    """Submit exactly one job; never automatically retry POST submissions."""
    try:
        response = session.post(
            f"{EBI_BASE}/{tool}/run",
            data={"email": email, **params},
            timeout=(15, 60),
        )
        response.raise_for_status()

    except requests.RequestException as exc:
        error_response = getattr(exc, "response", None)
        detail = (
            error_response.text[:1500]
            if error_response is not None
            else str(exc)
        )
        status = (
            error_response.status_code
            if error_response is not None
            else "unknown"
        )

        raise ConservationError(
            f"Could not submit {tool} job (HTTP {status}): {detail}"
        ) from exc

    job_id = response.text.strip()

    if not job_id or " " in job_id:
        raise ConservationError(
            f"EBI returned an invalid {tool} job ID: {job_id!r}"
        )

    return job_id


def _wait(
    tool: str,
    job_id: str,
    timeout_s: int,
    interval_s: int,
    session: requests.Session,
) -> None:
    """Poll a submitted EBI job with a monotonic timeout."""
    if timeout_s < 0:
        raise ValueError("timeout_s must not be negative.")

    if interval_s < 1:
        raise ValueError("poll_interval_s must be at least 1 second.")

    deadline = time.monotonic() + timeout_s

    while True:
        try:
            response = _get_with_retry(
                f"{EBI_BASE}/{tool}/status/{job_id}",
                session,
                timeout=(15, 60),
            )

        except requests.RequestException as exc:
            raise ConservationError(
                f"Could not poll EBI {tool} job {job_id}: {exc}"
            ) from exc

        status = response.text.strip().upper()

        if status == "FINISHED":
            return

        if status not in {"RUNNING", "QUEUED", "PENDING"}:
            raise ConservationError(
                f"EBI {tool} job {job_id} ended with status {status!r}."
            )

        remaining = deadline - time.monotonic()

        if remaining <= 0:
            break

        time.sleep(min(interval_s, remaining))

    raise ConservationError(
        f"EBI {tool} job {job_id} timed out after {timeout_s} seconds."
    )


def _result(
    tool: str,
    job_id: str,
    kind: str,
    session: requests.Session,
) -> str:
    try:
        response = _get_with_retry(
            f"{EBI_BASE}/{tool}/result/{job_id}/{kind}",
            session,
            timeout=(15, 120),
        )

    except requests.RequestException as exc:
        raise ConservationError(
            f"Could not retrieve {tool} {kind} output: {exc}"
        ) from exc

    return response.text


def _fetch_homologs(
    accessions: list[str],
    session: requests.Session,
) -> str:
    try:
        response = _get_with_retry(
            DBFETCH_URL,
            session,
            params={
                "db": "uniprotkb",
                "id": ",".join(accessions),
                "format": "fasta",
                "style": "raw",
            },
            timeout=(15, 120),
        )

    except requests.RequestException as exc:
        raise ConservationError(
            f"Could not retrieve homolog sequences from EBI Dbfetch: {exc}"
        ) from exc

    return response.text


def alignment_column_metrics(
    aligned_fasta: str,
    query_accession: str,
    query_position: int,
) -> dict[str, Any]:
    entries = parse_fasta(aligned_fasta)

    def _matches(header: str) -> bool:
        tokens = re.split(r"[|\s]", header)
        return (
            query_accession in tokens
            or header.split()[0] == query_accession
        )

    query = next(
        (seq for header, seq in entries if _matches(header)),
        None,
    )

    if query is None:
        raise ConservationError(
            f"Query accession {query_accession} was absent from the MSA."
        )

    if query_position < 1:
        raise ConservationError("UniProt positions start at 1.")

    seen = 0
    column = None

    for index, amino_acid in enumerate(query):
        if amino_acid != "-":
            seen += 1

            if seen == query_position:
                column = index
                break

    if column is None:
        raise ConservationError(
            f"UniProt position {query_position} exceeds "
            f"aligned query length {seen}."
        )

    residues = [
        seq[column]
        for _, seq in entries
        if len(seq) > column and seq[column] not in {"-", "."}
    ]

    if not residues:
        raise ConservationError(
            "The alignment column contains no amino acids."
        )

    counts = Counter(residues)
    n = len(residues)

    entropy = -sum(
        (count / n) * math.log2(count / n)
        for count in counts.values()
    )

    return {
        "query_accession": query_accession,
        "query_position_uniprot": query_position,
        "query_amino_acid": query[column],
        "aligned_sequence_count": len(entries),
        "column_coverage_count": n,
        "residue_counts": dict(sorted(counts.items())),
        "shannon_entropy_bits": round(entropy, 6),
        "method": (
            "Observed amino-acid frequencies and Shannon entropy at "
            "the Clustal Omega MSA column; not a ConSurf grade or "
            "pathogenicity score."
        ),
    }


def analyze_conservation(
    pdb_id: str,
    chain: str,
    mutation: Mutation,
    email: str,
    *,
    uniprot_accession: str | None = None,
    uniprot_position: int | None = None,
    insertion_code: str = "",
    max_hits: int = 30,
    timeout_s: int = 600,
    poll_interval_s: int = 10,
    session: requests.Session | None = None,
) -> dict[str, Any]:
    if not email or "@" not in email:
        raise ConservationError(
            "A contact email is required for programmatic "
            "EMBL-EBI job submission."
        )

    if max_hits < 2 or max_hits > 100:
        raise ValueError("max_hits must be between 2 and 100.")

    client = session or requests.Session()

    if uniprot_accession is None or uniprot_position is None:
        uniprot_accession, uniprot_position = (
            map_pdb_residue_to_uniprot(
                pdb_id,
                chain,
                mutation.position,
                insertion_code,
                session=client,
            )
        )

    if uniprot_position < 1:
        raise ConservationError("UniProt positions start at 1.")

    header, sequence = fetch_uniprot_sequence(
        uniprot_accession,
        client,
    )

    if uniprot_position > len(sequence):
        raise ConservationError(
            f"UniProt position {uniprot_position} exceeds "
            f"{uniprot_accession} length {len(sequence)}."
        )

    if sequence[uniprot_position - 1] != mutation.wild_type:
        raise ConservationError(
            f"Curated UniProt mapping points to "
            f"{sequence[uniprot_position - 1]} at "
            f"{uniprot_accession}:{uniprot_position}, "
            f"not mutation wild-type {mutation.wild_type}. "
            "Check isoform and residue mapping."
        )

    submitted_jobs: dict[str, Any] = {}

    try:
        blast_id = _submit(
            "ncbiblast",
            email,
            client,
            program="blastp",
            database="uniprotkb_swissprot",
            stype="protein",
            sequence=sequence,
            alignments=next(
                limit
                for limit in (5, 10, 20, 50, 100)
                if limit >= max_hits
            ),
        )
        submitted_jobs["blast_job_id"] = blast_id

        _wait(
            "ncbiblast",
            blast_id,
            timeout_s,
            poll_interval_s,
            client,
        )

        accessions = [
            line.strip().split(":")[-1]
            for line in _result(
                "ncbiblast",
                blast_id,
                "accs",
                client,
            ).splitlines()
            if line.strip()
        ]

        accessions = list(
            dict.fromkeys([uniprot_accession, *accessions])
        )[:max_hits]

        fasta = _fetch_homologs(accessions, client)
        fetched = parse_fasta(fasta)

        if len(fetched) < 2:
            raise ConservationError(
                f"Dbfetch returned {len(fetched)} sequences; "
                "at least two are required for alignment."
            )

        clustal_id = _submit(
            "clustalo",
            email,
            client,
            sequence=fasta,
            stype="protein",
            outfmt="fa",
        )
        submitted_jobs["clustalo_job_id"] = clustal_id

        _wait(
            "clustalo",
            clustal_id,
            timeout_s,
            poll_interval_s,
            client,
        )

        aligned = _result(
            "clustalo",
            clustal_id,
            "aln-fasta",
            client,
        )

        metrics = alignment_column_metrics(
            aligned,
            uniprot_accession,
            uniprot_position,
        )

        metrics.update({
            "uniprot_fasta_header": header,
            "blast_job_id": blast_id,
            "clustalo_job_id": clustal_id,
            "requested_homolog_count": len(accessions),
            "actual_sequence_count": len(fetched),
            "msa_fasta": aligned,
        })

        return metrics

    except ConservationError as exc:
        exc.metadata.update(submitted_jobs)
        raise
