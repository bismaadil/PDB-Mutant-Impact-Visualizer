"""Streamlit application for interactive, real-input mutation analysis."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import streamlit as st

from pdb_mutant_impact_visualizer.fetch_structure import download_pdb, load_structure
from pdb_mutant_impact_visualizer.local_structure import validate_mutation
from pdb_mutant_impact_visualizer.models import Mutation
from pdb_mutant_impact_visualizer.pipeline import analyze_mutation
from pdb_mutant_impact_visualizer.report import render_text
from pdb_mutant_impact_visualizer.visualize import build_structure_figure

st.set_page_config(page_title="PDB Mutation Analysis", page_icon="🧬", layout="wide")
st.title("PDB Mutation Analysis")
st.caption("Run structure-based measurements and optional remote predictors. Missing results are reported as missing, never filled with defaults.")


# Widgets live outside st.form so dependent fields (email, FoldX path) appear immediately.
source = st.radio("Structure source", ["RCSB PDB ID", "Upload PDB file"], horizontal=True)
pdb_id = ""
upload = None
refresh_structure = False
if source == "RCSB PDB ID":
    pdb_id = st.text_input("PDB ID", value="", max_chars=4, placeholder="e.g. 1AKI")
    refresh_structure = st.checkbox("Fetch a fresh coordinate file from RCSB", value=False)
else:
    upload = st.file_uploader("PDB coordinate file", type=["pdb", "ent"])
c1, c2, c3 = st.columns(3)
chain = c1.text_input("Chain", value="A", max_chars=8)
mutation_text = c2.text_input("Mutation", placeholder="e.g. S50V")
insertion_code = c3.text_input("Insertion code (optional)", max_chars=1)
cutoff = st.number_input("Neighbor distance cutoff (Å)", min_value=0.1, max_value=50.0, value=10.0)
stability = st.selectbox("Stability prediction", ["none", "dynamut2", "foldx"], index=0)
# Security: never execute a path typed into a web form. FoldX is taken only from the FOLDX_EXE
# environment variable set by whoever deploys the app.
foldx_exe = ""
if stability == "foldx":
    foldx_exe = os.environ.get("FOLDX_EXE", "")
    if foldx_exe:
        st.info("FoldX will run from the server-side FOLDX_EXE setting.")
    else:
        st.warning("FoldX is not configured on this server (set the FOLDX_EXE environment variable). Choose another backend or use the CLI.")
conservation = st.selectbox("Conservation analysis", ["none", "ebi"])
esm_model: str | None = None
if st.checkbox("Add ESM-2 language-model score (needs the 'ml' extra; downloads model weights on first use)"):
    esm_model = st.selectbox("ESM-2 checkpoint", ["facebook/esm2_t12_35M_UR50D", "facebook/esm2_t6_8M_UR50D", "facebook/esm2_t30_150M_UR50D"])
email = uniprot = uniprot_position = ""
if conservation == "ebi":
    email = st.text_input("Contact email (required for EBI jobs)")
    uniprot = st.text_input("Curated UniProt accession (required for uploads, optional for PDB IDs)")
    uniprot_position = st.text_input("Curated UniProt position (1-based; with accession)")

if stability == "dynamut2":
    st.info("DynaMut2 is a public external service. If selected, structure coordinates and the mutation are submitted to it.")
if conservation == "ebi":
    st.info("EBI analysis submits the UniProt sequence query, contact email, and job parameters to EMBL-EBI services.")
submitted = st.button("Run analysis", type="primary")

if submitted:
    for key in ("analysis_result", "analysis_inputs", "analysis_view_structure"):
        st.session_state.pop(key, None)
    structure_file: Path | None = None
    try:
        mutation = Mutation.parse(mutation_text)
        source_name = None
        if source == "RCSB PDB ID":
            if not pdb_id.strip():
                raise ValueError("Enter a PDB identifier.")
            analysis_pdb_id: str | None = pdb_id.strip().upper()
        else:
            if upload is None:
                raise ValueError("Upload a PDB file.")
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as handle:
                handle.write(upload.getvalue())
                structure_file = Path(handle.name)
            analysis_pdb_id, source_name = None, upload.name
        with st.spinner("Checking structure and running selected analyses…"):
            result = analyze_mutation(
                mutation=mutation, chain=chain, pdb_id=analysis_pdb_id,
                structure_path=structure_file, insertion_code=insertion_code or " ",
                neighbor_cutoff=float(cutoff), stability_backend=stability or "none",
                foldx_executable=foldx_exe or None, conservation_backend=conservation or "none",
                email=email or None, uniprot_accession=uniprot.strip() or None,
                uniprot_position=int(uniprot_position) if uniprot_position.strip() else None,
                refresh_structure=refresh_structure, source_name=source_name,
                esm_model=esm_model or None,
            )
        st.session_state["analysis_result"] = result.to_dict()
        st.session_state["analysis_view_structure"] = (
            load_structure(structure_file, "upload") if structure_file else None
        )
        st.session_state["analysis_inputs"] = (analysis_pdb_id, chain.strip(), mutation, float(cutoff), insertion_code or " ")
    except Exception as exc:
        st.error(str(exc))
    finally:
        if structure_file is not None:
            structure_file.unlink(missing_ok=True)

data = st.session_state.get("analysis_result")
if data:
    left, right = st.columns([1, 1])
    with left:
        st.subheader("Measurements and provenance")
        st.code(render_text(data), language="text")
        st.download_button("Download complete JSON", json.dumps(data, indent=2), file_name="mutation-analysis.json", mime="application/json")
        run_info = data.get("run")
        if run_info:
            st.subheader("Pipeline DAG run")
            st.caption(f"Config SHA-256 {run_info['config_sha256'][:16]}…  |  reproduce: `{run_info['reproduce_command']}`")
            st.dataframe(
                [{"stage": s_["name"], "status": s_["status"], "seconds": s_["seconds"], "remote": s_["remote"], "note": s_["error"] or ""}
                 for s_ in run_info["stages"]],
                use_container_width=True, hide_index=True,
            )
    with right:
        st.subheader("Structure view")
        pdb_for_view, view_chain, view_mutation, view_cutoff, view_ins = st.session_state["analysis_inputs"]
        try:
            parsed_structure = load_structure(download_pdb(pdb_for_view), pdb_for_view) if pdb_for_view else st.session_state["analysis_view_structure"]
            validate_mutation(parsed_structure, view_chain, view_mutation, view_ins)
            st.plotly_chart(build_structure_figure(parsed_structure, view_chain, view_mutation.position, view_cutoff, view_ins), use_container_width=True)
            st.caption("Cα trace colored by deposited B-factor; mutation site is largest. This is a coordinate visualization, not a predicted mutant structure.")
        except Exception as exc:
            st.warning(f"Structure view unavailable: {exc}")
