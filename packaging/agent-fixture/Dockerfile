# Agent Master fixture-execution image: Python + pyCSAMT core, nothing else.
# Built by build.sh from tracked files only; see
# docs/source/applications/agent_master/code_validation.rst.
FROM python:3.11-slim

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MPLBACKEND=Agg

COPY src/ /opt/pycsamt-src/
RUN pip install /opt/pycsamt-src \
    && rm -rf /opt/pycsamt-src \
    && python -c "import pycsamt, numpy, scipy, pandas, matplotlib; print('pycsamt', pycsamt.__version__)"

USER 65534:65534
WORKDIR /work
