#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

usage <- function() {
  cat("
Batch QC for GWAS .mlma result files.

Required:
  --manifest FILE              TSV with columns: population, trait, file
  --outdir DIR                 Output directory

Optional:
  --mode MODE                   metrics, qq, or full; default metrics
  --threads N                  Parallel workers, default 1
  --only-pop VALUE             Optional population filter for plotting one group
  --only-trait VALUE           Optional trait filter for plotting one trait
  --chr-col NAME               Chromosome column override
  --bp-col NAME                Base-pair position column override
  --snp-col NAME               SNP ID column override
  --p-col NAME                 P-value column override
  --lambda-warn-high NUM       Default 1.10
  --lambda-fail-high NUM       Default 1.20
  --lambda-warn-low NUM        Default 0.95
  --lambda-fail-low NUM        Default 0.90
  --min-snps-warn N            Default 50000
  --min-snps-fail N            Default 10000
  --sig-threshold NUM          Default 5e-8
  --save-cleaned TRUE/FALSE    Save cleaned valid rows, default FALSE
  --help                       Show help

Example:
  Rscript scripts/gwas_mlma_qc.R --manifest manifest.tsv --outdir qc_output --threads 8 --mode metrics
")
}

parse_args <- function(args) {
  out <- list(
    mode = "metrics",
    threads = 1L,
    only_pop = NA_character_,
    only_trait = NA_character_,
    chr_col = NA_character_,
    bp_col = NA_character_,
    snp_col = NA_character_,
    p_col = NA_character_,
    lambda_warn_high = 1.10,
    lambda_fail_high = 1.20,
    lambda_warn_low = 0.95,
    lambda_fail_low = 0.90,
    min_snps_warn = 50000L,
    min_snps_fail = 10000L,
    sig_threshold = 5e-8,
    save_cleaned = FALSE
  )
  if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
    usage()
    quit(status = 0)
  }
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for ", key)
    }
    value <- args[[i + 1L]]
    key_clean <- gsub("-", "_", sub("^--", "", key))
    out[[key_clean]] <- value
    i <- i + 2L
  }
  numeric_keys <- c(
    "lambda_warn_high", "lambda_fail_high", "lambda_warn_low",
    "lambda_fail_low", "sig_threshold"
  )
  integer_keys <- c("threads", "min_snps_warn", "min_snps_fail")
  for (key in numeric_keys) {
    out[[key]] <- as.numeric(out[[key]])
  }
  for (key in integer_keys) {
    out[[key]] <- as.integer(out[[key]])
  }
  out$save_cleaned <- tolower(as.character(out$save_cleaned)) %in% c("true", "t", "1", "yes", "y")
  if (!out$mode %in% c("metrics", "qq", "full")) {
    stop("--mode must be one of: metrics, qq, full")
  }
  out
}

find_col <- function(dt, override, candidates, label) {
  names_original <- names(dt)
  names_lower <- tolower(names_original)
  if (!is.na(override) && nzchar(override)) {
    if (!override %in% names_original) {
      stop("Requested ", label, " column not found: ", override)
    }
    return(override)
  }
  hit <- match(tolower(candidates), names_lower, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) == 0L) {
    stop("Could not detect ", label, " column. Tried: ", paste(candidates, collapse = ", "))
  }
  names_original[[hit[[1L]]]]
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  gsub("^_+|_+$", "", x)
}

expected_qq <- function(n) {
  -log10(ppoints(n))
}

observed_qq <- function(p) {
  -log10(sort(p))
}

qq_metrics <- function(p) {
  p <- p[!is.na(p) & is.finite(p) & p > 0 & p <= 1]
  n <- length(p)
  if (n < 10L) {
    return(list(
      qq_body_median_delta = NA_real_,
      qq_body_rmse = NA_real_,
      qq_tail_lift = NA_real_,
      qq_shape = "UNCLASSIFIED",
      qq_interpretation = "Too few valid P values for QQ curve classification."
    ))
  }
  qq_dt <- data.table(expected = expected_qq(n), observed = observed_qq(p))
  qq_dt[, delta := observed - expected]
  body_dt <- qq_dt[expected >= 0.1 & expected <= 2]
  if (nrow(body_dt) == 0L) {
    body_dt <- qq_dt[seq_len(max(1L, floor(n * 0.9)))]
  }
  tail_dt <- qq_dt[expected > 2]
  body_median_delta <- median(body_dt$delta, na.rm = TRUE)
  body_rmse <- sqrt(mean(body_dt$delta^2, na.rm = TRUE))
  tail_lift <- if (nrow(tail_dt) > 0L) max(tail_dt$delta, na.rm = TRUE) else NA_real_

  shape <- "QQ_IDEAL"
  interpretation <- "QQ body follows the null expectation; tail deviation, if present, is compatible with possible association signal."
  if (is.finite(body_median_delta) && body_median_delta > 0.25) {
    shape <- "QQ_INFLATED"
    interpretation <- "QQ body is globally above the diagonal, suggesting population structure, relatedness, batch effects, or phenotype confounding."
  } else if (is.finite(body_median_delta) && body_median_delta < -0.25) {
    shape <- "QQ_DEFLATED"
    interpretation <- "QQ body is below the diagonal, suggesting overcorrection, overfitted kinship, or overly aggressive phenotype residualization."
  } else if (is.finite(body_rmse) && body_rmse > 0.35) {
    shape <- "QQ_NOISY"
    interpretation <- "QQ curve is irregular; review sample size, model convergence, and missingness."
  } else if (is.finite(tail_lift) && tail_lift > 1.0) {
    shape <- "QQ_GOOD_WITH_TAIL_SIGNAL"
    interpretation <- "QQ body is acceptable with strong tail deviation; this can be a real association signal if Manhattan plot is localized."
  }

  list(
    qq_body_median_delta = body_median_delta,
    qq_body_rmse = body_rmse,
    qq_tail_lift = tail_lift,
    qq_shape = shape,
    qq_interpretation = interpretation
  )
}

make_qq_plot <- function(p, title, out_file) {
  p <- p[!is.na(p) & is.finite(p) & p > 0 & p <= 1]
  qq_dt <- data.table(expected = expected_qq(length(p)), observed = observed_qq(p))
  max_axis <- max(qq_dt$expected, qq_dt$observed, na.rm = TRUE)
  p1 <- ggplot(qq_dt, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, color = "#777777", linewidth = 0.4) +
    geom_point(size = 0.45, alpha = 0.75, color = "#3B6EA8") +
    coord_equal(xlim = c(0, max_axis), ylim = c(0, max_axis)) +
    labs(title = title, x = expression(Expected~~-log[10](P)), y = expression(Observed~~-log[10](P))) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11))

  png(out_file, width = 1000, height = 950, res = 170)
  print(p1)
  dev.off()
}

make_full_plot <- function(d, chr_col, bp_col, p_col, title, out_file, sig_threshold) {
  plot_dt <- d[, .(
    chr_raw = get(chr_col),
    bp = as.numeric(get(bp_col)),
    p = as.numeric(get(p_col))
  )]
  plot_dt <- plot_dt[!is.na(chr_raw) & !is.na(bp) & !is.na(p) & p > 0 & p <= 1]
  plot_dt[, chr := suppressWarnings(as.integer(gsub("^chr", "", as.character(chr_raw), ignore.case = TRUE)))]
  plot_dt <- plot_dt[!is.na(chr)]
  setorder(plot_dt, chr, bp)

  chr_offsets <- plot_dt[, .(chr_len = max(bp, na.rm = TRUE)), by = chr]
  chr_offsets[, offset := shift(cumsum(chr_len), fill = 0)]
  plot_dt <- merge(plot_dt, chr_offsets[, .(chr, offset)], by = "chr", all.x = TRUE)
  plot_dt[, pos_cum := bp + offset]
  axis_dt <- plot_dt[, .(center = (min(pos_cum) + max(pos_cum)) / 2), by = chr]

  p1 <- ggplot(plot_dt, aes(x = pos_cum, y = -log10(p), color = factor(chr %% 2))) +
    geom_point(size = 0.35, alpha = 0.75) +
    geom_hline(yintercept = -log10(sig_threshold), color = "#B83232", linewidth = 0.35) +
    scale_x_continuous(label = axis_dt$chr, breaks = axis_dt$center, expand = expansion(mult = c(0.005, 0.005))) +
    scale_color_manual(values = c("#3B6EA8", "#1E1E1E"), guide = "none") +
    labs(title = title, x = "Chromosome", y = expression(-log[10](P))) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.title = element_text(size = 11)
    )

  qq_p <- plot_dt$p
  n_qq <- length(qq_p)
  qq_dt <- data.table(expected = expected_qq(n_qq), observed = observed_qq(qq_p))
  max_axis <- max(qq_dt$expected, qq_dt$observed, na.rm = TRUE)
  p2 <- ggplot(qq_dt, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, color = "#777777", linewidth = 0.4) +
    geom_point(size = 0.45, alpha = 0.75, color = "#3B6EA8") +
    coord_equal(xlim = c(0, max_axis), ylim = c(0, max_axis)) +
    labs(title = "QQ plot", x = expression(Expected~~-log[10](P)), y = expression(Observed~~-log[10](P))) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11))

  png(out_file, width = 2200, height = 950, res = 170)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  dev.off()
}

recommendation_for <- function(status, reasons) {
  if (status == "ERROR") {
    return("Check whether the file exists, is readable, and contains required MLMA columns.")
  }
  if (grepl("lambda_gc_high", reasons)) {
    return("Inspect population structure, GRM/kinship, batch covariates, phenotype outliers, and sample ID matching. Rerun with appropriate PCs/covariates or analyze populations separately before meta-analysis.")
  }
  if (grepl("lambda_gc_low", reasons)) {
    return("Check for overcorrection from redundant covariates, overly aggressive residualization, or an overfitted kinship model.")
  }
  if (grepl("too_few_valid_snps|too_few_snps", reasons)) {
    return("Check genotype filtering, imputation INFO/MAF thresholds, file truncation, chromosome inclusion, and model convergence logs.")
  }
  if (grepl("many_invalid_p", reasons)) {
    return("Check model convergence, phenotype missingness, sample ID matching, and malformed P-value output.")
  }
  if (grepl("many_significant", reasons)) {
    return("Check phenotype-batch association, duplicated/related samples, and trait outliers. Confirm signals are not spread uniformly across the genome.")
  }
  if (status == "WARN") {
    return("Review the Manhattan/QQ plot and compare against other populations for the same trait.")
  }
  "No automated issue detected. Keep genotype-level QC and model logs for audit."
}

solution_for <- function(status, reasons, qq_shape) {
  if (status == "ERROR") {
    return("File or column problem: verify file path, required MLMA columns, delimiters, and whether the GWAS job finished successfully.")
  }
  if (grepl("too_few_valid_snps|many_invalid_p", reasons)) {
    return("Data completeness problem: check genotype filtering, imputation INFO/MAF thresholds, sample ID matching, chromosome inclusion, and model convergence logs; rerun affected GWAS after fixing input completeness.")
  }
  if (grepl("lambda_gc_high", reasons) || qq_shape == "QQ_INFLATED") {
    return("Inflation problem: add/check PCA covariates, batch/plate/field covariates, sex/age if relevant, GRM/kinship construction, duplicate/close relatives, and phenotype outliers; rerun by population and meta-analyze if structure differs.")
  }
  if (grepl("lambda_gc_low", reasons) || qq_shape == "QQ_DEFLATED") {
    return("Deflation problem: reduce redundant covariates, inspect overfitted kinship or overly aggressive residualization, and verify phenotype transformation did not remove true genetic signal.")
  }
  if (grepl("many_significant", reasons)) {
    return("Excess signal problem: test phenotype association with batch/plate/family/source variables, inspect outliers, and use full Manhattan plots to confirm whether peaks are localized rather than genome-wide artifacts.")
  }
  if (qq_shape == "QQ_GOOD_WITH_TAIL_SIGNAL") {
    return("Likely acceptable with candidate signal: draw full Manhattan+QQ for this population/trait, check whether signals are localized, then annotate candidate loci.")
  }
  if (qq_shape == "QQ_NOISY") {
    return("Noisy QQ problem: check sample size, missingness, convergence, low MAC variants, and rerun with stricter variant filters if needed.")
  }
  if (status == "PASS") {
    return("Accept for downstream review; keep model logs and genotype-level QC records.")
  }
  "Review QQ plot together with numeric QC metrics; rerun full Manhattan+QQ for this group/trait before deciding."
}

quality_class_for <- function(status, reasons, qq_shape) {
  if (status == "ERROR") return("ERROR_INPUT")
  if (grepl("too_few_valid_snps_fail", reasons) || (status == "FAIL" && grepl("many_invalid_p", reasons))) {
    return("FAIL_DATA_COMPLETENESS")
  }
  if (grepl("too_few_valid_snps_warn|many_invalid_p", reasons)) return("WARN_DATA_COMPLETENESS")
  if (grepl("lambda_gc_high", reasons) || qq_shape == "QQ_INFLATED") return("FAIL_OR_WARN_INFLATION")
  if (grepl("lambda_gc_low", reasons) || qq_shape == "QQ_DEFLATED") return("FAIL_OR_WARN_DEFLATION")
  if (grepl("many_significant", reasons)) return("WARN_EXCESS_SIGNAL")
  if (qq_shape == "QQ_GOOD_WITH_TAIL_SIGNAL") return("PASS_WITH_TAIL_SIGNAL")
  if (qq_shape == "QQ_NOISY") return("WARN_NOISY_QQ")
  if (status == "PASS") return("PASS_IDEAL")
  "REVIEW"
}

qc_one <- function(row, args) {
  population <- as.character(row$population)
  trait <- as.character(row$trait)
  file <- as.character(row$file)
  id <- paste(safe_name(population), safe_name(trait), sep = "__")
  qq_dir <- file.path(args$outdir, "qq_plots", safe_name(population))
  full_dir <- file.path(args$outdir, "full_plots", safe_name(population))
  cleaned_dir <- file.path(args$outdir, "cleaned", safe_name(population))
  if (args$mode == "qq") dir.create(qq_dir, recursive = TRUE, showWarnings = FALSE)
  if (args$mode == "full") dir.create(full_dir, recursive = TRUE, showWarnings = FALSE)
  if (args$save_cleaned) {
    dir.create(cleaned_dir, recursive = TRUE, showWarnings = FALSE)
  }

  base_result <- data.table(
    population = population,
    trait = trait,
    file = file,
    status = "ERROR",
    reasons = NA_character_,
    n_rows = NA_integer_,
    n_valid_p = NA_integer_,
    n_invalid_p = NA_integer_,
    min_p = NA_real_,
    n_significant = NA_integer_,
    lambda_gc = NA_real_,
    median_chisq = NA_real_,
    qq_body_median_delta = NA_real_,
    qq_body_rmse = NA_real_,
    qq_tail_lift = NA_real_,
    qq_shape = NA_character_,
    qq_interpretation = NA_character_,
    qq_plot_file = if (args$mode == "qq") file.path(qq_dir, paste0(id, ".qq.png")) else NA_character_,
    full_plot_file = if (args$mode == "full") file.path(full_dir, paste0(id, ".mqq.png")) else NA_character_,
    recommendation = NA_character_,
    quality_class = NA_character_,
    solution = NA_character_
  )

  tryCatch({
    if (!file.exists(file)) {
      stop("File does not exist")
    }
    d <- fread(file, showProgress = FALSE)
    if (args$mode == "full") {
      chr_col <- find_col(d, args$chr_col, c("Chr", "CHR", "chrom", "chromosome", "chromosome_name"), "chromosome")
      bp_col <- find_col(d, args$bp_col, c("bp", "BP", "pos", "POS", "position"), "position")
      snp_col <- find_col(d, args$snp_col, c("SNP", "snp", "rs", "ID", "id", "marker", "variant"), "SNP")
    }
    p_col <- find_col(d, args$p_col, c("p", "P", "pval", "PVAL", "p_value", "p.value"), "P-value")

    total_rows <- nrow(d)
    p <- suppressWarnings(as.numeric(d[[p_col]]))
    valid <- !is.na(p) & is.finite(p) & p > 0 & p <= 1
    d_valid <- d[valid]
    p_valid <- p[valid]
    n_valid <- length(p_valid)
    n_invalid <- total_rows - n_valid

    reasons <- character()
    qc_status <- "PASS"
    if (n_valid < args$min_snps_fail) {
      qc_status <- "FAIL"
      reasons <- c(reasons, "too_few_valid_snps_fail")
    } else if (n_valid < args$min_snps_warn) {
      qc_status <- "WARN"
      reasons <- c(reasons, "too_few_valid_snps_warn")
    }
    if (total_rows > 0 && n_invalid / total_rows > 0.05) {
      qc_status <- if (n_invalid / total_rows > 0.20) "FAIL" else if (qc_status == "PASS") "WARN" else qc_status
      reasons <- c(reasons, "many_invalid_p")
    }

    chisq <- qchisq(1 - p_valid, df = 1)
    med_chisq <- median(chisq, na.rm = TRUE)
    gc_lambda <- med_chisq / qchisq(0.5, df = 1)
    if (is.finite(gc_lambda)) {
      if (gc_lambda >= args$lambda_fail_high) {
        qc_status <- "FAIL"
        reasons <- c(reasons, "lambda_gc_high_fail")
      } else if (gc_lambda >= args$lambda_warn_high && qc_status == "PASS") {
        qc_status <- "WARN"
        reasons <- c(reasons, "lambda_gc_high_warn")
      } else if (gc_lambda <= args$lambda_fail_low) {
        qc_status <- "FAIL"
        reasons <- c(reasons, "lambda_gc_low_fail")
      } else if (gc_lambda <= args$lambda_warn_low && qc_status == "PASS") {
        qc_status <- "WARN"
        reasons <- c(reasons, "lambda_gc_low_warn")
      }
    }

    n_sig <- sum(p_valid < args$sig_threshold, na.rm = TRUE)
    if (n_sig > max(50L, ceiling(n_valid * 0.001))) {
      if (qc_status == "PASS") {
        qc_status <- "WARN"
      }
      reasons <- c(reasons, "many_significant")
    }

    qqm <- if (args$mode %in% c("qq", "full")) qq_metrics(p_valid) else list(
      qq_body_median_delta = NA_real_,
      qq_body_rmse = NA_real_,
      qq_tail_lift = NA_real_,
      qq_shape = NA_character_,
      qq_interpretation = NA_character_
    )

    if (args$mode == "qq") {
      make_qq_plot(
        p_valid,
        title = paste(population, trait, qqm$qq_shape),
        out_file = base_result$qq_plot_file
      )
    }

    if (args$mode == "full") {
      make_full_plot(
        d_valid,
        chr_col = chr_col,
        bp_col = bp_col,
        p_col = p_col,
        title = paste(population, trait),
        out_file = base_result$full_plot_file,
        sig_threshold = args$sig_threshold
      )
    }

    if (args$save_cleaned) {
      fwrite(d_valid, file.path(cleaned_dir, paste0(id, ".valid.tsv.gz")), sep = "\t")
    }

    if (length(reasons) == 0L) {
      reasons_text <- "none"
    } else {
      reasons_text <- paste(unique(reasons), collapse = ";")
    }
    qc_recommendation <- recommendation_for(qc_status, reasons_text)
    q_class <- quality_class_for(qc_status, reasons_text, qqm$qq_shape)
    q_solution <- solution_for(qc_status, reasons_text, qqm$qq_shape)
    base_result[, `:=`(
      status = qc_status,
      reasons = reasons_text,
      n_rows = total_rows,
      n_valid_p = n_valid,
      n_invalid_p = n_invalid,
      min_p = min(p_valid, na.rm = TRUE),
      n_significant = n_sig,
      lambda_gc = gc_lambda,
      median_chisq = med_chisq,
      qq_body_median_delta = qqm$qq_body_median_delta,
      qq_body_rmse = qqm$qq_body_rmse,
      qq_tail_lift = qqm$qq_tail_lift,
      qq_shape = qqm$qq_shape,
      qq_interpretation = qqm$qq_interpretation,
      recommendation = qc_recommendation,
      quality_class = q_class,
      solution = q_solution
    )]
    base_result
  }, error = function(e) {
    base_result[, `:=`(
      reasons = conditionMessage(e),
      recommendation = recommendation_for("ERROR", conditionMessage(e)),
      quality_class = "ERROR_INPUT",
      solution = solution_for("ERROR", conditionMessage(e), NA_character_)
    )]
    base_result
  })
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$manifest) || is.null(args$outdir)) {
  usage()
  stop("--manifest and --outdir are required")
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
manifest <- fread(args$manifest)
required <- c("population", "trait", "file")
missing_cols <- setdiff(required, names(manifest))
if (length(missing_cols) > 0L) {
  stop("Manifest is missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!is.na(args$only_pop) && nzchar(args$only_pop)) {
  manifest <- manifest[population == args$only_pop]
}
if (!is.na(args$only_trait) && nzchar(args$only_trait)) {
  manifest <- manifest[trait == args$only_trait]
}
if (nrow(manifest) == 0L) {
  stop("No manifest rows remain after --only-pop/--only-trait filtering.")
}

message("Loaded manifest rows: ", nrow(manifest))
message("Output directory: ", normalizePath(args$outdir, mustWork = FALSE))
message("Run mode: ", args$mode)

rows <- split(manifest, seq_len(nrow(manifest)))
if (args$threads > 1L) {
  results <- parallel::mclapply(rows, qc_one, args = args, mc.cores = args$threads)
} else {
  results <- lapply(rows, qc_one, args = args)
}
summary_dt <- rbindlist(results, fill = TRUE)
summary_dt <- cbind(manifest, summary_dt[, !names(summary_dt) %in% names(manifest), with = FALSE])

fwrite(summary_dt, file.path(args$outdir, "qc_summary.tsv"), sep = "\t")
fwrite(summary_dt[status == "FAIL"], file.path(args$outdir, "qc_fail.tsv"), sep = "\t")
fwrite(summary_dt[status == "WARN"], file.path(args$outdir, "qc_warn.tsv"), sep = "\t")
fwrite(summary_dt[status == "PASS"], file.path(args$outdir, "qc_pass.tsv"), sep = "\t")
fwrite(
  summary_dt[, .(population, trait, file, status, reasons, recommendation, qq_plot_file, full_plot_file)],
  file.path(args$outdir, "qc_recommendations.tsv"),
  sep = "\t"
)
fwrite(
  summary_dt[, .(
    population, trait, file, status, reasons, quality_class,
    lambda_gc, n_valid_p, n_invalid_p, min_p, n_significant,
    qq_shape, qq_body_median_delta, qq_body_rmse, qq_tail_lift,
    qq_interpretation, solution, qq_plot_file, full_plot_file
  )],
  file.path(args$outdir, "qc_quality_solutions.tsv"),
  sep = "\t"
)
if (args$mode %in% c("qq", "full")) {
  fwrite(
    summary_dt[, .(
      population, trait, file, status, quality_class,
      qq_shape, qq_body_median_delta, qq_body_rmse, qq_tail_lift,
      qq_interpretation, qq_plot_file, full_plot_file
    )],
    file.path(args$outdir, "qq_shape_summary.tsv"),
    sep = "\t"
  )
}

message("Done.")
message("PASS: ", nrow(summary_dt[status == "PASS"]))
message("WARN: ", nrow(summary_dt[status == "WARN"]))
message("FAIL: ", nrow(summary_dt[status == "FAIL"]))
message("ERROR: ", nrow(summary_dt[status == "ERROR"]))
