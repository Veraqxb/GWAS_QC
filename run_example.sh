#!/usr/bin/env bash
set -euo pipefail

Rscript examples/create_smoke_data.R smoke_data

Rscript scripts/gwas_mlma_qc.R \
  --manifest smoke_data/manifest.tsv \
  --outdir qc_example \
  --threads 2
