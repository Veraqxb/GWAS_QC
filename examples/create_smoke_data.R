#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

outdir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  commandArgs(trailingOnly = TRUE)[[1L]]
} else {
  "smoke_data"
}

set.seed(1)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n <- 12000L
d <- data.table(
  Chr = rep(1:2, each = n / 2),
  SNP = paste0("rs", seq_len(n)),
  bp = rep(seq_len(n / 2), 2),
  A1 = "A",
  A2 = "G",
  Freq = runif(n, 0.05, 0.95),
  b = rnorm(n),
  se = runif(n, 0.05, 0.2),
  p = runif(n)
)

mlma_file <- file.path(outdir, "smoke.mlma")
manifest_file <- file.path(outdir, "manifest.tsv")
fwrite(d, mlma_file, sep = "\t")
fwrite(
  data.table(
    population = "pop200",
    trait = "smoke",
    file = normalizePath(mlma_file),
    sample_size = 200
  ),
  manifest_file,
  sep = "\t"
)

message("Wrote ", normalizePath(manifest_file))
