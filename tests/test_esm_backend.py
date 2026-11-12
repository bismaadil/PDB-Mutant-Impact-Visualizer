import pytest

from pdb_mutant_impact_visualizer.esm_backend import (
    AMINO_ACIDS, ESMError, crop_window, masked_marginal_llr, score_chain_mutation, sequence_from_residues,
)


def fake_model(preferred):
    """Deterministic stand-in for a language model: favours `preferred`."""
    def fn(masked, index):
        assert "<mask>" in masked
        return {a: (-0.1 if a == preferred else -3.0) for a in AMINO_ACIDS}
    return fn


def test_llr_sign_and_ranks():
    seq = "MKTAYIAKQR"
    good = masked_marginal_llr(fake_model("A"), seq, 3, "A", "V")
    assert good["esm_llr"] == pytest_approx(-2.9)
    assert good["wild_type_rank_among_20"] == 1
    bad = masked_marginal_llr(fake_model("V"), seq, 3, "A", "V")
    assert bad["esm_llr"] > 0 and bad["mutant_rank_among_20"] == 1


def pytest_approx(v):
    return pytest.approx(v, abs=1e-9)


def test_wild_type_must_match_sequence():
    with pytest.raises(ESMError, match="not the stated wild type"):
        masked_marginal_llr(fake_model("A"), "MKTAY", 1, "A", "V")


def test_crop_window_keeps_position_and_bounds():
    seq = "A" * 3000
    w, i = crop_window(seq, 10, 1022)
    assert len(w) == 1022 and i == 10
    w, i = crop_window(seq, 2990, 1022)
    assert len(w) == 1022 and i == 1022 - 10
    w, i = crop_window(seq, 1500, 1022)
    assert len(w) == 1022 and i == 511
    assert crop_window("ACDE", 2) == ("ACDE", 2)
    with pytest.raises(ESMError):
        crop_window("ACDE", 9)


def test_sequence_gaps_are_counted_and_reported():
    residues = [(1, "", "M"), (2, "", "K"), (5, "", "A"), (6, "", "Y")]
    seq, keys, gaps = sequence_from_residues(residues)
    assert seq == "MKAY" and gaps == 1
    out = score_chain_mutation(fake_model("A"), residues, 5, "", "A", "V", "fake/model")
    assert out["numbering_gaps_in_coordinates"] == 1 and "warning" in out
    with pytest.raises(ESMError, match="not found"):
        score_chain_mutation(fake_model("A"), residues, 99, "", "A", "V", "fake/model")
