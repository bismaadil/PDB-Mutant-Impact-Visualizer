"""Serialization helpers. Missing measurements remain null; no proxy score is generated."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .models import AnalysisResult


def write_json(result: AnalysisResult, output: str | Path) -> Path:
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def render_text(result: AnalysisResult | dict[str, Any]) -> str:
    data = result if isinstance(result, dict) else result.to_dict()
    lines = [
        f"Mutation analysis: {data['structure']['pdb_id'] or data['structure']['file_name']} "
        f"chain {data['chain']} {data['mutation']}",
        f"Structure SHA-256: {data['structure']['sha256']}",
        f"Stability backend: {data['stability']['backend'] if data['stability'] else 'not run'}",
    ]
    if data["stability"]:
        lines.append(f"  ddG (kcal/mol): {data['stability']['ddg_kcal_per_mol']}")
    else:
        lines.append("  ddG (kcal/mol): unavailable (no stability backend selected)")
    cons = data["conservation"]
    lines.append(f"Conservation: {cons['method'] if cons else 'not run'}")
    if cons:
        lines.append(
            f"  UniProt {cons['query_accession']}:{cons['query_position_uniprot']} "
            f"{cons['query_amino_acid']}; coverage {cons['column_coverage_count']}; "
            f"entropy {cons['shannon_entropy_bits']} bits"
        )
    local = data["local_structure"]
    lines.append(
        f"Local geometry: {local['neighbor_residue_count']} residues within "
        f"{local['neighbor_cutoff_angstrom']} Angstrom; site mean B-factor "
        f"{local['site_mean_b_factor_angstrom2']} A^2"
    )
    sasa = local.get("solvent_accessibility")
    if sasa:
        lines.append(f"Solvent accessibility: {sasa['sasa_angstrom2']} A^2 (relative {sasa['relative_sasa']})")
    contacts = local.get("contacts")
    if contacts:
        lines.append(
            f"Contacts: {contacts['side_chain_contact_residue_count']} side-chain contact residues "
            f"(<= {contacts['side_chain_contact_cutoff_angstrom']} A); {contacts['hydrogen_bond_candidate_count']} H-bond candidates; "
            f"{contacts['salt_bridge_candidate_count']} salt-bridge candidates; "
            f"disulfide {'yes' if contacts['disulfide_bonds'] else 'no'}; interface residue {'yes' if contacts['is_interface_residue'] else 'no'}"
        )
    lm = data.get("language_model")
    if lm:
        lines.append(f"Language model ({lm['model']}): LLR {lm['esm_llr']} (negative = mutant less plausible than wild type)")
    mlp = data.get("ml_prediction")
    if mlp:
        lines.append(f"Research ML model: ddG {mlp['ddg_kcal_per_mol']} kcal/mol ({mlp['target_sign_convention']}); grouped-CV Spearman {mlp['cv_spearman_grouped']}")
    sub = data.get("substitution")
    if sub:
        lines.append(
            f"Substitution: BLOSUM62 {sub['blosum62_score']}; hydropathy change {sub['hydropathy_change_kyte_doolittle']}; "
            f"volume change {sub['volume_change_angstrom3']} A^3; charge change {sub['charge_change_ph7_approx']}"
        )
    for note in data.get("warnings", []):
        lines.append(f"WARNING: {note}")
    run = data.get("run")
    if run:
        failed = [st["name"] for st in run["stages"] if st["status"] == "failed"]
        ran = sum(1 for st in run["stages"] if st["status"] == "ok")
        lines.append(f"Pipeline: {ran}/{len(run['stages'])} stages ran; failed: {', '.join(failed) or 'none'}; config SHA-256 {run['config_sha256'][:16]}...")
    lines.append(f"Analyzed (UTC): {data['analyzed_at_utc']}")
    lines.append("Interpretation: independent measurements are reported separately; no composite impact score is calculated.")
    return "\n".join(lines)


REQUIRED_KEYS = (
    "schema_version", "software", "dependency_versions", "analyzed_at_utc", "structure", "chain",
    "mutation", "stability", "conservation", "local_structure", "substitution", "warnings", "run",
)


def validate_report(data: dict[str, Any]) -> list[str]:
    """Structural check of a saved report; returns a list of problems (empty = valid)."""
    problems = [f"missing key: {k}" for k in REQUIRED_KEYS if k not in data]
    if data.get("schema_version") != "2.0":
        problems.append(f"unsupported schema_version {data.get('schema_version')!r}")
    structure = data.get("structure") or {}
    if len(str(structure.get("sha256", ""))) != 64:
        problems.append("structure.sha256 is not a 64-character digest")
    run = data.get("run") or {}
    if run and len(str(run.get("config_sha256", ""))) != 64:
        problems.append("run.config_sha256 is not a 64-character digest")
    if "email" in json.dumps(data.get("run", {})).lower().replace("email_provided", ""):
        problems.append("run manifest appears to contain an email field")
    return problems
