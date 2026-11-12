# Security and data handling

Report security issues privately to the repository maintainer rather than posting credentials or private structures in a public issue.

The Streamlit application binds according to Streamlit's local defaults. Do not expose it to an untrusted network without adding authentication and deployment controls. Uploaded PDB files are written to a temporary file during analysis and removed after the run. When DynaMut2 is selected, structure coordinates and mutation details are submitted to the external service. When EMBL-EBI conservation is selected, the UniProt sequence query, contact email, and job parameters are sent to EMBL-EBI. Consider data-sharing permissions before using non-public structures.

Never commit credentials, private structures, analysis outputs containing sensitive metadata, FoldX executables, or licensed FoldX files.

**FoldX in the web app:** the executable is taken only from the server-side `FOLDX_EXE` environment variable; a path typed into the browser is never executed. Uploads are capped at 50 MB (`.streamlit/config.toml`). The CLI's `--foldx-exe` runs whatever you point it at on your own machine, so only use it with software you trust.

**Reports** never contain the EBI contact email (only `email_provided: true`).

