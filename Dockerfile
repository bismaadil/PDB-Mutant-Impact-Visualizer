# syntax=docker/dockerfile:1
# One Dockerfile, three build targets (see compose.yaml):
#   EXTRAS=app            -> Streamlit app + CLI            (default, small)
#   EXTRAS=app,ml,plots   -> + ESM-2 (CPU PyTorch) + heatmaps (research image)
# FoldX is licensed separately and is NEVER copied into an image.
ARG BASE_IMAGE=python:3.12-slim-bookworm
FROM ${BASE_IMAGE}

ARG EXTRAS=app
ARG GIT_COMMIT=unknown
ARG TORCH_INDEX=https://download.pytorch.org/whl/cpu

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PDBMI_IN_CONTAINER=1 \
    GIT_COMMIT=${GIT_COMMIT} \
    HF_HOME=/opt/pdb-mutant-impact-visualizer/data/hf-cache \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

WORKDIR /opt/pdb-mutant-impact-visualizer

# Dependency layer first (cached unless pyproject/constraints change).
# constraints.txt is optional: create it with `make lock` for a fully pinned build.
COPY pyproject.toml README.md LICENSE constraints.txt ./
COPY cli.py app.py ./
COPY .streamlit ./.streamlit
COPY pdb_mutant_impact_visualizer ./pdb_mutant_impact_visualizer

RUN python -m pip install --upgrade pip \
    && CONSTRAINT="-c constraints.txt" \
    && case ",${EXTRAS}," in \
         *,ml,*) python -m pip install ${CONSTRAINT} torch --extra-index-url "${TORCH_INDEX}" ;; \
       esac \
    && python -m pip install ${CONSTRAINT} ".[${EXTRAS}]" \
    && useradd --create-home --uid 10001 --shell /usr/sbin/nologin app \
    && mkdir -p data \
    && chown -R app:app /opt/pdb-mutant-impact-visualizer

USER 10001:10001
EXPOSE 8501
VOLUME ["/opt/pdb-mutant-impact-visualizer/data"]

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
  CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=3)"]

# Default: the web app. For the CLI: docker run IMAGE pdbmi run --pdb 1AKI --chain A --mutation S50V
CMD ["streamlit", "run", "app.py", "--server.address=0.0.0.0", "--server.port=8501"]
