"""Zero-shot protein language-model score for one substitution (ESM-2).

Method: masked-marginal log-likelihood ratio (Meier et al., NeurIPS 2021): mask
the mutated position, read the model's amino-acid distribution there, and
report ``log p(mutant) - log p(wild type)``. NEGATIVE values mean the model
finds the mutant less plausible than the wild type in this sequence context.
This is an evolutionary-plausibility score, not a stability (ddG) value and not
a pathogenicity classifier; it is reported separately from every other metric.

Heavy dependencies (torch, transformers) are imported lazily and only when
this backend is selected. The scoring maths is separate from the model so it
can be tested without downloading weights.
"""

from __future__ import annotations

from typing import Any, Callable, Mapping

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
ESM_MAX_RESIDUES = 1022  # ESM-2 context minus BOS/EOS
DEFAULT_MODEL = "facebook/esm2_t12_35M_UR50D"
REFERENCE = "https://doi.org/10.1101/2021.07.09.450648"  # Meier et al. 2021 (bioRxiv)
ESM2_REFERENCE = "https://doi.org/10.1126/science.ade2574"  # Lin et al. 2023, Science

LogProbFn = Callable[[str, int], Mapping[str, float]]


class ESMError(RuntimeError):
    pass


def sequence_from_residues(residues: list[tuple[int, str, str]]) -> tuple[str, list[tuple[int, str]], int]:
    """Build a chain sequence from (resnum, icode, one-letter) in file order.

    Returns (sequence, keys, gap_count). ``gap_count`` counts jumps in author
    numbering, i.e. residues absent from the coordinates; the language model
    then sees a sequence with those residues removed, which is reported as a
    warning rather than hidden.
    """
    seq = "".join(aa for _, _, aa in residues)
    keys = [(n, i) for n, i, _ in residues]
    gaps = sum(1 for (a, _), (b, _) in zip(keys, keys[1:], strict=False) if b - a > 1)
    return seq, keys, gaps


def crop_window(sequence: str, index: int, max_len: int = ESM_MAX_RESIDUES) -> tuple[str, int]:
    """Crop to at most ``max_len`` residues centred on ``index``; returns the new index."""
    if not 0 <= index < len(sequence):
        raise ESMError(f"Position index {index} is outside the sequence (length {len(sequence)}).")
    if len(sequence) <= max_len:
        return sequence, index
    start = min(max(0, index - max_len // 2), len(sequence) - max_len)
    return sequence[start:start + max_len], index - start


def masked_marginal_llr(
    log_probs: LogProbFn, sequence: str, index: int, wild_type: str, mutant: str
) -> dict[str, Any]:
    if sequence[index] != wild_type:
        raise ESMError(
            f"Sequence position {index + 1} is {sequence[index]!r}, not the stated wild type {wild_type!r}."
        )
    window, local = crop_window(sequence, index)
    masked = window[:local] + "<mask>" + window[local + 1:]
    lp = log_probs(masked, local)
    missing = [a for a in (wild_type, mutant) if a not in lp]
    if missing:
        raise ESMError(f"Model output lacks amino acid(s): {missing}")
    llr = float(lp[mutant]) - float(lp[wild_type])
    ranked = sorted(AMINO_ACIDS, key=lambda a: -float(lp.get(a, float("-inf"))))
    return {
        "esm_llr": round(llr, 6),
        "log_prob_wild_type": round(float(lp[wild_type]), 6),
        "log_prob_mutant": round(float(lp[mutant]), 6),
        "wild_type_rank_among_20": ranked.index(wild_type) + 1,
        "mutant_rank_among_20": ranked.index(mutant) + 1,
        "context_length": len(window),
        "cropped": len(window) != len(sequence),
        "sign_convention": "log p(mutant) - log p(wild type); negative = model finds the mutant less plausible.",
        "note": "Zero-shot evolutionary plausibility. Not a ddG and not a clinical pathogenicity classifier.",
    }


class HuggingFaceESM:
    """Lazy wrapper around a Hugging Face ESM-2 checkpoint (CPU by default)."""

    def __init__(self, model_name: str = DEFAULT_MODEL, device: str = "cpu"):
        try:
            import torch
            from transformers import AutoModelForMaskedLM, AutoTokenizer
        except ImportError as exc:  # pragma: no cover - exercised only with the extra installed
            raise ESMError(
                "The ESM backend needs the optional 'ml' extra: pip install 'pdb-mutant-impact-visualizer[ml]'"
            ) from exc
        self.model_name = model_name
        self._torch = torch
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = AutoModelForMaskedLM.from_pretrained(model_name).to(device).eval()
        self.device = device
        self._ids = {a: self.tokenizer.convert_tokens_to_ids(a) for a in AMINO_ACIDS}

    def __call__(self, masked_sequence: str, index: int) -> Mapping[str, float]:  # pragma: no cover
        torch = self._torch
        text = masked_sequence.replace("<mask>", self.tokenizer.mask_token)
        enc = self.tokenizer(text, return_tensors="pt").to(self.device)
        with torch.no_grad():
            logits = self.model(**enc).logits[0, index + 1]  # +1 for the BOS/CLS token
        logp = torch.log_softmax(logits, dim=-1)
        return {a: float(logp[i]) for a, i in self._ids.items()}


def score_chain_mutation(
    log_probs: LogProbFn, residues: list[tuple[int, str, str]], resnum: int, icode: str,
    wild_type: str, mutant: str, model_name: str,
) -> dict[str, Any]:
    sequence, keys, gaps = sequence_from_residues(residues)
    try:
        index = keys.index((resnum, icode))
    except ValueError as exc:
        raise ESMError(f"Residue {resnum}{icode} not found in the chain sequence.") from exc
    result = masked_marginal_llr(log_probs, sequence, index, wild_type, mutant)
    result.update({
        "backend": "ESM-2 masked marginal",
        "model": model_name,
        "chain_sequence_length": len(sequence),
        "numbering_gaps_in_coordinates": gaps,
        "references": [REFERENCE, ESM2_REFERENCE],
    })
    if gaps:
        result["warning"] = (
            f"{gaps} gap(s) in the coordinate numbering: residues missing from the structure are absent "
            "from the sequence the model saw."
        )
    return result
