# GWAS MLMA Batch QC

This repository provides a server-ready workflow for batch quality control of GCTA/GEMMA-style GWAS `.mlma` summary result files across multiple populations and traits.

It is designed for use cases like:

- thousands of `.mlma` result files
- automated Manhattan plot, QQ plot, lambda GC, P-value sanity checks, and QC flagging

## What The Pipeline Produces

For every GWAS result file, the pipeline writes:

- combined Manhattan + QQ plot
- optional compressed cleaned summary table
- per-trait QC metrics

For the whole batch, it writes:

- `qc_summary.tsv`: one row per GWAS result
- `qc_fail.tsv`: results requiring review
- `qc_warn.tsv`: borderline results
- `qc_pass.tsv`: results passing automated checks
- `qc_recommendations.tsv`: likely problem pattern and suggested adjustment

## Recommended Directory Layout On Server

```text
project/
  gwas-results/
    pop1/
    pop2/
    pop3/
  manifest.tsv
  gwas-mlma-batch-qc/
```

The manifest controls which files are analyzed.

## Manifest Format

Create a tab-separated `manifest.tsv` with these columns:

```text
population	trait	file
pop200	height	/path/to/gwas-results/pop1/height.mlma
pop500	height	/path/to/gwas-results/pop2/height.mlma
pop540	height	/path/to/gwas-results/pop3/height.mlma
```

Optional columns are allowed and preserved in the output, for example `sample_size`, `batch`, or `note`.

An example file is provided at `examples/manifest.example.tsv`.

## Install R Dependencies

On most Linux servers:

```bash
Rscript -e 'install.packages(c("data.table", "ggplot2"), repos="https://cloud.r-project.org")'
```

No other R packages are required.

## Run QC

```bash
Rscript scripts/gwas_mlma_qc.R \
  --manifest manifest.tsv \
  --outdir qc_output \
  --threads 8
```

For small test runs:

```bash
Rscript examples/create_smoke_data.R smoke_data
Rscript scripts/gwas_mlma_qc.R \
  --manifest smoke_data/manifest.tsv \
  --outdir qc_test \
  --threads 2
```

## Input Column Detection

The script automatically recognizes common `.mlma` columns:

| Meaning | Accepted names |
|---|---|
| chromosome | `Chr`, `CHR`, `chrom`, `chromosome` |
| position | `bp`, `BP`, `pos`, `position` |
| SNP ID | `SNP`, `rs`, `ID`, `marker`, `variant` |
| P value | `p`, `P`, `pval`, `PVAL`, `p_value` |
| effect | `b`, `BETA`, `beta`, `effect` |
| standard error | `se`, `SE`, `stderr` |
| allele frequency | `Freq`, `freq`, `AF`, `maf`, `MAF` |

If your files use unusual column names, use the explicit options:

```bash
Rscript scripts/gwas_mlma_qc.R \
  --manifest manifest.tsv \
  --outdir qc_output \
  --chr-col Chr \
  --bp-col bp \
  --snp-col SNP \
  --p-col p
```

## Default QC Rules

The defaults are intentionally conservative for sample sizes around 200-540.

| QC status | Rule |
|---|---|
| `PASS` | no major automated problem detected |
| `WARN` | borderline lambda GC, few SNPs, many missing P values, or unusual significant count |
| `FAIL` | severe lambda GC deviation, invalid P values dominate, missing required columns, unreadable file |

Default thresholds:

| Parameter | Default |
|---|---:|
| lambda GC fail high | `1.20` |
| lambda GC warn high | `1.10` |
| lambda GC warn low | `0.95` |
| lambda GC fail low | `0.90` |
| minimum SNP count fail | `10000` |
| minimum SNP count warn | `50000` |
| genome-wide significance | `5e-8` |

You can override thresholds:

```bash
Rscript scripts/gwas_mlma_qc.R \
  --manifest manifest.tsv \
  --outdir qc_output \
  --lambda-warn-high 1.10 \
  --lambda-fail-high 1.20 \
  --lambda-warn-low 0.95 \
  --lambda-fail-low 0.90 \
  --min-snps-warn 50000 \
  --min-snps-fail 10000
```

## How To Interpret Problem Patterns

| Pattern | Common cause | Recommended adjustment |
|---|---|---|
| QQ plot globally above diagonal, lambda GC high | population structure, relatedness, batch effect, phenotype confounding | add/top up PCA covariates; check GRM/kinship; include batch/sex/age; rerun by population then meta-analyze |
| QQ plot below diagonal, lambda GC low | overcorrection, too many covariates, phenotype residualization too strong | reduce redundant covariates; check kinship and residualization strategy |
| Many significant SNPs across many chromosomes | phenotype or batch artifact | test phenotype against batch/plate/family/sex; inspect outliers; winsorize or transform trait if justified |
| One chromosome or region has broad elevation | map/build mismatch, local genotyping artifact, structural region, low-quality imputation | check genome build, allele coding, INFO/MAF/missingness; rerun excluding problematic low-quality SNPs |
| Very few SNPs or many invalid P values | file truncation, ID mismatch, overly strict filtering, model convergence failures | check input genotype filters, sample ID matching, model logs, and per-file line counts |
| Only small sample population is noisy | low power and unstable estimates | use as sensitivity analysis; prioritize larger groups or meta-analysis |

## Suggested Review Workflow

1. Run batch QC for all files.
2. Open `qc_summary.tsv` and sort by `status`, `lambda_gc`, `n_significant`, and `n_valid_p`.
3. Review only `FAIL` and `WARN` plots first.
4. Fix upstream model/filtering problems.
5. Rerun the failed subset using a smaller manifest.
6. Preserve the first and rerun QC outputs for audit.

## GitHub Usage

Initial upload from a server:

```bash
git init
git add README.md scripts config examples .gitignore
git commit -m "Add batch GWAS MLMA QC workflow"
git branch -M main
git remote add origin git@github.com:YOUR_ORG/gwas-mlma-batch-qc.git
git push -u origin main
```

After updates:

```bash
git add README.md scripts config examples .gitignore
git commit -m "Update GWAS QC workflow"
git push
```

## Notes

- This pipeline performs summary-statistics QC. It does not replace genotype-level QC.
- For publication-grade interpretation, inspect cohort design, kinship, phenotype transformation, covariates, and genotype filtering logs.
- If ancestry or population structure differs across the 3 groups, analyze each group separately and combine with meta-analysis rather than forcing a single pooled model.
