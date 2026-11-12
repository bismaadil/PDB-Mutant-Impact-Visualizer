
import requests
import pytest

from pdb_mutant_impact_visualizer.conservation import (
    ConservationError,
    _get_with_retry,
    _submit,
    alignment_column_metrics,
    analyze_conservation,
    map_pdb_residue_to_uniprot,
)
from pdb_mutant_impact_visualizer.models import Mutation


def test_alignment_column_metrics_reports_counts_and_entropy():
    aligned = (
        ">sp|P12345|QUERY\n"
        "AC-D\n"
        ">sp|P2|HIT\n"
        "AT-D\n"
        ">sp|P3|HIT\n"
        "AT-D\n"
    )

    metrics = alignment_column_metrics(aligned, "P12345", 2)

    assert metrics["query_amino_acid"] == "C"
    assert metrics["column_coverage_count"] == 3
    assert metrics["residue_counts"] == {"C": 1, "T": 2}
    assert metrics["shannon_entropy_bits"] == pytest.approx(0.918296)
    assert "not a ConSurf grade" in metrics["method"]


class Response:
    """Small requests.Response-like object for isolated tests."""

    def __init__(
        self,
        payload=None,
        status_code=200,
        text="",
        headers=None,
    ):
        self.payload = payload
        self.status_code = status_code
        self.text = text
        self.headers = headers or {}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(
                f"HTTP {self.status_code}",
                response=self,
            )

    def json(self):
        return self.payload


class MappingSession:
    def get(self, url, timeout):
        return Response({
            "1abc": {
                "UniProt": {
                    "P12345": {
                        "mappings": [
                            {
                                "chain_id": "A",
                                "start": {
                                    "author_residue_number": 10
                                },
                                "end": {
                                    "author_residue_number": 20
                                },
                                "unp_start": 100,
                                "unp_end": 110,
                            }
                        ]
                    }
                }
            }
        })


def test_sifts_mapping_converts_author_position_to_uniprot_position():
    assert map_pdb_residue_to_uniprot(
        "1ABC",
        "A",
        15,
        session=MappingSession(),
    ) == ("P12345", 105)


def test_sifts_ambiguous_or_missing_map_is_rejected():
    class Empty:
        def get(self, url, timeout):
            return Response({"1abc": {"UniProt": {}}})

    with pytest.raises(ConservationError, match="unambiguous"):
        map_pdb_residue_to_uniprot(
            "1ABC",
            "A",
            15,
            session=Empty(),
        )


class SequenceGetSession:
    """Returns predefined responses or raises predefined exceptions."""

    def __init__(self, outcomes):
        self.outcomes = list(outcomes)
        self.calls = 0

    def get(self, url, timeout, **kwargs):
        self.calls += 1
        outcome = self.outcomes.pop(0)

        if isinstance(outcome, Exception):
            raise outcome

        return outcome


def test_get_retries_http_429_and_then_succeeds(monkeypatch):
    session = SequenceGetSession([
        Response(status_code=429, text="Too many requests"),
        Response(text="success"),
    ])

    delays = []
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation.time.sleep",
        lambda delay: delays.append(delay),
    )

    response = _get_with_retry(
        "https://example.test/resource",
        session,
        max_attempts=3,
    )

    assert response.text == "success"
    assert session.calls == 2
    assert len(delays) == 1


def test_get_retries_nested_slurm_rate_limit_inside_http_400(monkeypatch):
    nested_error = (
        '{"error":"rate_limit_exceeded",'
        '"message":"Maximum 10 requests/second with burst of 20"}'
    )

    session = SequenceGetSession([
        Response(status_code=400, text=nested_error),
        Response(text="success"),
    ])

    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation.time.sleep",
        lambda delay: None,
    )

    response = _get_with_retry(
        "https://example.test/ebi",
        session,
        max_attempts=3,
    )

    assert response.text == "success"
    assert session.calls == 2


def test_get_retries_http_503_and_then_succeeds(monkeypatch):
    session = SequenceGetSession([
        Response(status_code=503, text="Service unavailable"),
        Response(text="success"),
    ])

    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation.time.sleep",
        lambda delay: None,
    )

    response = _get_with_retry(
        "https://example.test/resource",
        session,
        max_attempts=3,
    )

    assert response.text == "success"
    assert session.calls == 2


def test_get_respects_retry_after_header(monkeypatch):
    session = SequenceGetSession([
        Response(
            status_code=429,
            text="Rate limited",
            headers={"Retry-After": "3"},
        ),
        Response(text="success"),
    ])

    delays = []
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation.time.sleep",
        lambda delay: delays.append(delay),
    )

    response = _get_with_retry(
        "https://example.test/resource",
        session,
        max_attempts=3,
    )

    assert response.text == "success"
    assert delays == [3.0]
    assert session.calls == 2


def test_get_retries_timeout_until_attempt_limit(monkeypatch):
    session = SequenceGetSession([
        requests.Timeout("first timeout"),
        requests.Timeout("second timeout"),
        requests.Timeout("third timeout"),
    ])

    delays = []
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation.time.sleep",
        lambda delay: delays.append(delay),
    )

    with pytest.raises(requests.Timeout):
        _get_with_retry(
            "https://example.test/resource",
            session,
            max_attempts=3,
        )

    assert session.calls == 3
    assert len(delays) == 2


def test_get_does_not_retry_permanent_http_400(monkeypatch):
    session = SequenceGetSession([
        Response(status_code=400, text="Bad request"),
        Response(text="should not be reached"),
    ])

    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation.time.sleep",
        lambda delay: None,
    )

    with pytest.raises(requests.HTTPError):
        _get_with_retry(
            "https://example.test/resource",
            session,
            max_attempts=3,
        )

    assert session.calls == 1


def test_post_job_submission_is_never_retried():
    class FailedPostSession:
        def __init__(self):
            self.calls = 0

        def post(self, url, data, timeout):
            self.calls += 1
            raise requests.Timeout("Submission timed out")

    session = FailedPostSession()

    with pytest.raises(ConservationError, match="Could not submit"):
        _submit(
            "ncbiblast",
            "test@example.com",
            session,
            program="blastp",
        )

    assert session.calls == 1


def test_blast_job_id_is_preserved_when_blast_polling_fails(monkeypatch):
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation.fetch_uniprot_sequence",
        lambda accession, session: (
            "sp|P12345|QUERY",
            "A" * 100,
        ),
    )

    submitted = []

    def fake_submit(tool, email, session, **params):
        submitted.append(tool)
        return "blast-test-job-001"

    def fake_wait(tool, job_id, timeout_s, interval_s, session):
        raise ConservationError("BLAST polling failed")

    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation._submit",
        fake_submit,
    )
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation._wait",
        fake_wait,
    )

    with pytest.raises(ConservationError, match="BLAST polling failed") as exc:
        analyze_conservation(
            "1ABC",
            "A",
            Mutation.parse("A50V"),
            "test@example.com",
            uniprot_accession="P12345",
            uniprot_position=50,
            session=object(),
        )

    assert submitted == ["ncbiblast"]
    assert exc.value.metadata == {
        "blast_job_id": "blast-test-job-001"
    }
    assert "blast-test-job-001" in str(exc.value)


def test_both_job_ids_are_preserved_when_clustal_polling_fails(monkeypatch):
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation.fetch_uniprot_sequence",
        lambda accession, session: (
            "sp|P12345|QUERY",
            "A" * 100,
        ),
    )

    submitted = []

    def fake_submit(tool, email, session, **params):
        job_id = (
            "blast-test-job-002"
            if tool == "ncbiblast"
            else "clustal-test-job-002"
        )
        submitted.append((tool, job_id))
        return job_id

    def fake_wait(tool, job_id, timeout_s, interval_s, session):
        if tool == "clustalo":
            raise ConservationError("Clustal polling failed")

    def fake_result(tool, job_id, kind, session):
        assert tool == "ncbiblast"
        assert job_id == "blast-test-job-002"
        assert kind == "accs"
        return "P12345\nP23456\n"

    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation._submit",
        fake_submit,
    )
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation._wait",
        fake_wait,
    )
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation._result",
        fake_result,
    )
    monkeypatch.setattr(
        "pdb_mutant_impact_visualizer.conservation._fetch_homologs",
        lambda accessions, session: (
            ">sp|P12345|QUERY\n"
            + "A" * 100
            + "\n>sp|P23456|HOMOLOG\n"
            + "A" * 100
            + "\n"
        ),
    )

    with pytest.raises(ConservationError, match="Clustal polling failed") as exc:
        analyze_conservation(
            "1ABC",
            "A",
            Mutation.parse("A50V"),
            "test@example.com",
            uniprot_accession="P12345",
            uniprot_position=50,
            session=object(),
        )

    assert submitted == [
        ("ncbiblast", "blast-test-job-002"),
        ("clustalo", "clustal-test-job-002"),
    ]
    assert exc.value.metadata == {
        "blast_job_id": "blast-test-job-002",
        "clustalo_job_id": "clustal-test-job-002",
    }
    assert "blast-test-job-002" in str(exc.value)
    assert "clustal-test-job-002" in str(exc.value)
