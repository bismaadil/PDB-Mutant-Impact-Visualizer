# Test in GitHub Codespaces, then publish to GitHub

Everything below is run in a **Codespace** (a cloud VS Code with a terminal and internet access). Replace `YOUR-USERNAME` with your GitHub username.

## A. Create the repository and open a Codespace
1. On github.com click **+ → New repository**. Name: `pdb-mutant-impact-visualizer`. Choose **Private** for now. Do **not** add a README/licence (the project has them). Create.
2. On the empty repo page click **Code → Codespaces → Create codespace on main**. Wait for the editor to load.
3. Drag `pdb-mutant-impact-visualizer-v2_0_0.zip` from your computer into the **Explorer** panel (left side).
4. In the Codespace terminal (**Terminal → New Terminal**):
```bash
unzip pdb-mutant-impact-visualizer-v2_0_0.zip
shopt -s dotglob && mv pdb-mutant-impact-visualizer/* . && rmdir pdb-mutant-impact-visualizer
rm pdb-mutant-impact-visualizer-v2_0_0.zip
ls -a            # you should see pyproject.toml, Dockerfile, .devcontainer, ...
```
5. Rebuild so the devcontainer installs everything: press **F1 → "Codespaces: Rebuild Container"** (first time only). When the terminal returns, run `pdbmi --version`.

## B. Verify (do these in order; stop and fix when something fails)
```bash
make test                      # offline unit tests
make smoke                     # offline run -> validate -> reproduce on the bundled fixture
make typecheck lint            # mypy + ruff: fix anything reported (see note below)
make live                      # NETWORK: real structure 1AKI through the whole offline pipeline
bash scripts/live_check.sh     # NETWORK: run + reproduce + 3-position scan on 1AKI
```
**Note on mypy/ruff:** the v2.0.0 code was written without access to these tools. If they report issues, fix them (usually one-line annotations), then in `.github/workflows/tests.yml` delete the line `continue-on-error: true` so CI enforces them.

If `make live` shows a mismatch (e.g. a wild-type residue error for S50V), check the residue with `pdbmi run --pdb 1AKI --chain A --mutation <WT><N><MUT>`; use the wild type the error message reports.

## C. Containers
```bash
make lock                      # writes constraints.txt (pinned versions) -> commit it
make docker-build
make docker-smoke              # CLI inside the container, offline fixture
docker compose up app          # Codespaces shows a forwarded port 8501 -> open it
```
Record the image digest with results: `export IMAGE_DIGEST=$(docker image inspect pdbmi:2.0.0 --format '{{.Id}}')`.

## D. Remote services (optional, needs a real email address)
```bash
pdbmi run --pdb 1AKI --chain A --mutation S50V --stability dynamut2 --conservation ebi \
          --email you@your-institute.edu --timeout 1200 --output data/results/1aki-remote.json
```
Never commit your email: it is not written to reports, but do not paste it into files.

## E. ML / research run
```bash
pip install -e ".[ml]"                    # PyTorch + Transformers (large download)
pdbmi scan --pdb 1AKI --chain A --positions 45-55 --esm-model facebook/esm2_t12_35M_UR50D --heatmap-metric esm_llr
```
Then follow `docs/RESEARCH_PROTOCOL.md` with a real dataset. A default Codespace has limited RAM; stay with the 35M or 8M ESM-2 checkpoints. **Stop the Codespace** when idle (Command palette → "Codespaces: Stop Current Codespace") to save your free hours.

## F. Make it yours before publishing
- `LICENSE`: replace `[Your Name]` with your name (keep MIT unless your institution requires otherwise).
- `CITATION.cff`: put your real name/affiliation/ORCID in `authors`.
- `README.md`: no placeholders remain, but add a one-paragraph "Project status" with your results once you have them.
- `docs/RESEARCH_PROTOCOL.md`: fill in the bracketed sections; add your results table to `docs/RESULTS.md`.
- Check `git status` shows no `data/`, structures, emails or FoldX files.

## G. Commit and push
In Codespaces your git login is already set up:
```bash
git config user.name "Your Name" && git config user.email "your-github-noreply@users.noreply.github.com"
git add -A && git status        # review the list
git commit -m "Release v2.0.0: DAG pipeline, containers, research layer"
git branch -M main
git push -u origin main
```
(GitHub → Settings → Emails shows your `…@users.noreply.github.com` address; using it keeps your personal email out of the history.)

## H. After pushing
1. **Actions** tab: wait for `tests` (and `docker`) to turn green. Fix red items before sharing the link.
2. **Release**: `git tag -a v2.0.0 -m "v2.0.0" && git push origin v2.0.0`, then GitHub → Releases → *Draft a new release*.
3. **DOI (optional but impressive)**: log in to zenodo.org with GitHub, enable the repository, then publish the GitHub release; Zenodo mints a DOI. Add it to `CITATION.cff` and the README.
4. Flip the repo to **Public** (Settings → Danger zone) when you are ready, add topics (`bioinformatics`, `protein-stability`, `reproducible-research`), and pin it on your profile.
5. Put the repo link, the DAG diagram and your results table at the top of your application.
