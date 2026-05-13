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
pop1	height	/path/to/gwas-results/pop1/height.mlma
pop2	height	/path/to/gwas-results/pop2/height.mlma
pop3	height	/path/to/gwas-results/pop3/height.mlma
```

Optional columns are allowed and preserved in the output, for example `sample_size`, `batch`, or `note`.

An example file is provided at `examples/manifest.example.tsv`.

## Install R Dependencies

On most Linux servers:

```bash
Rscript -e 'install.packages(c("data.table", "ggplot2"), repos="https://cloud.r-project.org")'
```

No other R packages are required.

## Three-Step QC Workflow

For large batches, use the three-step workflow below. It is much faster than drawing Manhattan plots for every `.mlma` file.

### Step 1: Numeric QC Only

This reads each `.mlma` file, calculates QC metrics, and does not draw plots.

```bash
Rscript scripts/gwas_mlma_qc.R \
  --manifest manifest.tsv \
  --outdir qc_step1_metrics \
  --threads 8 \
  --mode metrics
```

Main outputs:

- `qc_summary.tsv`
- `qc_pass.tsv`
- `qc_warn.tsv`
- `qc_fail.tsv`
- `qc_quality_solutions.tsv`

### Step 2: QQ Plots For All Results

This draws QQ plots for all manifest rows and classifies the QQ shape.

```bash
Rscript scripts/gwas_mlma_qc.R \
  --manifest manifest.tsv \
  --outdir qc_step2_qq \
  --threads 4 \
  --mode qq
```

Main outputs:

- `qq_plots/`
- `qq_shape_summary.tsv`
- `qc_quality_solutions.tsv`

QQ shape classes:

| QQ class | Meaning |
|---|---|
| `QQ_IDEAL` | QQ body follows the diagonal; acceptable for downstream review |
| `QQ_GOOD_WITH_TAIL_SIGNAL` | QQ body is acceptable, with strong tail deviation compatible with true loci |
| `QQ_INFLATED` | global upward deviation; likely structure, relatedness, batch effect, or confounding |
| `QQ_DEFLATED` | global downward deviation; likely overcorrection or overfitted model |
| `QQ_NOISY` | irregular curve; check sample size, missingness, convergence, or low-count variants |

### Step 3: Full Manhattan + QQ For One Group

After reviewing Step 1 and Step 2, draw complete Manhattan + QQ plots for one population/group.

```bash
Rscript scripts/gwas_mlma_qc.R \
  --manifest manifest.tsv \
  --outdir qc_step3_full_C6 \
  --threads 4 \
  --mode full \
  --only-pop C6
```

You can also draw one specific population-trait pair:

```bash
Rscript scripts/gwas_mlma_qc.R \
  --manifest manifest.tsv \
  --outdir qc_step3_full_C6_height \
  --threads 2 \
  --mode full \
  --only-pop C6 \
  --only-trait height
```

Main outputs:

- `full_plots/`
- `qc_summary.tsv`
- `qc_quality_solutions.tsv`

## Run QC Legacy One-Step Example

```bash
Rscript scripts/gwas_mlma_qc.R \
  --manifest manifest.tsv \
  --outdir qc_output \
  --threads 8 \
  --mode full
```

For small test runs:

```bash
Rscript examples/create_smoke_data.R smoke_data
Rscript scripts/gwas_mlma_qc.R \
  --manifest smoke_data/manifest.tsv \
  --outdir qc_test \
  --threads 2 \
  --mode metrics
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
  --mode metrics \
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

## Final Quality Classes

The file `qc_quality_solutions.tsv` summarizes the combined quality decision and recommended action.

| Quality class | Interpretation | Main action |
|---|---|---|
| `PASS_IDEAL` | Numeric QC passes; no major QQ problem detected | Keep for downstream analysis |
| `PASS_WITH_TAIL_SIGNAL` | QQ body is good with tail signal | Draw full Manhattan + QQ and inspect whether peaks are localized |
| `WARN_EXCESS_SIGNAL` | Too many significant SNPs | Check phenotype-batch associations and outliers |
| `WARN_NOISY_QQ` | QQ is irregular | Check sample size, missingness, low MAC variants, and convergence |
| `FAIL_OR_WARN_INFLATION` | Lambda/QQ inflation | Add/check PCs, kinship, batch covariates, and population-specific analysis |
| `FAIL_OR_WARN_DEFLATION` | Lambda/QQ deflation | Check overcorrection, residualization, and overfitted mixed model |
| `WARN_DATA_COMPLETENESS` | Fewer SNPs or more invalid P values than expected | Check filters and logs; may still be usable after review |
| `FAIL_DATA_COMPLETENESS` | Too few valid SNPs or many invalid P values | Check file completeness, SNP filters, imputation/MAF thresholds, and model logs |
| `ERROR_INPUT` | File or column problem | Check path, delimiter, and required MLMA columns |
| `REVIEW` | Borderline mixed signals | Review QQ plus full Manhattan plot |

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
