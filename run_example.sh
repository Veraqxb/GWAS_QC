#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN=${PYTHON:-python}

"$PYTHON_BIN" examples/create_smoke_data.py smoke_data

"$PYTHON_BIN" scripts/gwas_mlma_qc_fast.py \
  --manifest smoke_data/manifest.tsv \
  --outdir qc_example \
  --threads 2 \
  --mode metrics
