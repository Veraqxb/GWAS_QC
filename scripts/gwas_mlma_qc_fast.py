#!/usr/bin/env python3
"""Fast batch QC for GWAS .mlma files.

This script mirrors the R workflow modes but uses Python streaming/Polars-friendly
data access and matplotlib plotting. It is intended for server-side screening
when R plotting is too slow.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Missing dependency: numpy. Install with: conda install -c conda-forge numpy") from exc

try:
    from scipy.stats import chi2
except ImportError:  # pragma: no cover
    chi2 = None

try:
    import polars as pl
except ImportError:  # pragma: no cover
    pl = None

try:
    import cupy as cp
except ImportError:  # pragma: no cover
    cp = None


P_CANDIDATES = ["p", "P", "pval", "PVAL", "p_value", "p.value"]
CHR_CANDIDATES = ["Chr", "CHR", "chrom", "chromosome", "chromosome_name"]
BP_CANDIDATES = ["bp", "BP", "pos", "POS", "position"]
SNP_CANDIDATES = ["SNP", "snp", "rs", "ID", "id", "marker", "variant"]
SUMMARY_COLUMNS = [
    "population",
    "trait",
    "file",
    "status",
    "reasons",
    "n_rows",
    "n_valid_p",
    "n_invalid_p",
    "min_p",
    "n_significant",
    "lambda_gc",
    "median_chisq",
    "qq_body_median_delta",
    "qq_body_rmse",
    "qq_tail_lift",
    "qq_shape",
    "qq_interpretation",
    "qq_plot_file",
    "manhattan_plot_file",
    "full_plot_file",
    "recommendation",
    "quality_class",
    "solution",
    "assessment_mode",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fast batch QC for GWAS .mlma files")
    parser.add_argument("--manifest", required=True, help="TSV with population, trait, file columns")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--mode", choices=["metrics", "qq", "manhattan", "full"], default="metrics")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--only-pop", default=None)
    parser.add_argument("--only-trait", default=None)
    parser.add_argument("--limit", type=int, default=None, help="Process only the first N manifest rows after filters.")
    parser.add_argument("--p-col", default=None)
    parser.add_argument("--chr-col", default=None)
    parser.add_argument("--bp-col", default=None)
    parser.add_argument("--snp-col", default=None)
    parser.add_argument("--lambda-warn-high", type=float, default=1.10)
    parser.add_argument("--lambda-fail-high", type=float, default=1.20)
    parser.add_argument("--lambda-warn-low", type=float, default=0.95)
    parser.add_argument("--lambda-fail-low", type=float, default=0.90)
    parser.add_argument("--min-snps-warn", type=int, default=50000)
    parser.add_argument("--min-snps-fail", type=int, default=10000)
    parser.add_argument("--sig-threshold", type=float, default=5e-8)
    parser.add_argument("--max-plot-points", type=int, default=200000)
    parser.add_argument("--plot-keep-p", type=float, default=1e-4)
    parser.add_argument("--qq-max-points", type=int, default=100000, help="Maximum points drawn in each QQ plot.")
    parser.add_argument("--qq-ci-points", type=int, default=2000, help="Number of points used for QQ confidence interval ribbons.")
    parser.add_argument("--qq-confidence", type=float, default=0.95, help="QQ confidence interval level.")
    parser.add_argument("--qq-engine", choices=["auto", "cpu", "gpu"], default="auto", help="Use CuPy GPU sorting for QQ when available.")
    parser.add_argument("--plot-dpi", type=int, default=180)
    parser.add_argument(
        "--assessment-mode",
        choices=["full", "sample", "calculate"],
        default="full",
        help="full applies all QC thresholds; sample skips full-file SNP-count/excess-hit thresholds; calculate only reports metrics.",
    )
    return parser.parse_args()


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value))
    return re.sub(r"^_+|_+$", "", value)


def detect_col(header: List[str], override: Optional[str], candidates: List[str], label: str) -> str:
    if override:
        if override not in header:
            raise ValueError(f"Requested {label} column not found: {override}")
        return override
    lower_to_name = {name.lower(): name for name in header}
    for candidate in candidates:
        if candidate.lower() in lower_to_name:
            return lower_to_name[candidate.lower()]
    raise ValueError(f"Could not detect {label} column. Tried: {', '.join(candidates)}")


def sniff_delimiter(path: str) -> Optional[str]:
    with open(path, "r", newline="") as handle:
        first = handle.readline()
    if "\t" in first:
        return "\t"
    if "," in first:
        return ","
    return None


def read_header(path: str) -> Tuple[List[str], Optional[str]]:
    delimiter = sniff_delimiter(path)
    with open(path, "r", newline="") as handle:
        first = handle.readline().strip()
    if delimiter is None:
        return first.split(), delimiter
    return next(csv.reader([first], delimiter=delimiter)), delimiter


def iter_rows_by_header(path: str):
    header, delimiter = read_header(path)
    with open(path, "r", newline="") as handle:
        handle.readline()
        if delimiter is None:
            for line in handle:
                parts = line.strip().split()
                if not parts:
                    continue
                yield header, dict(zip(header, parts))
        else:
            reader = csv.DictReader(handle, fieldnames=header, delimiter=delimiter)
            for row in reader:
                yield header, row


def read_manifest(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    required = {"population", "trait", "file"}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"Manifest is missing required columns: {', '.join(sorted(missing))}")
    return rows


def write_tsv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "" if row.get(key) is None else row.get(key) for key in fieldnames})


def load_pyplot():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Missing dependency: matplotlib. Install with: conda install -c conda-forge matplotlib") from exc


def p_to_chisq(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    if chi2 is not None:
        return chi2.isf(np.clip(p_values, np.nextafter(0, 1), 1.0), 1)
    # Wilson-Hilferty approximation fallback for environments without scipy.
    from statistics import NormalDist

    nd = NormalDist()
    clipped = np.clip(p_values, np.nextafter(0, 1), 1.0 - 1e-16)
    z = np.array([nd.inv_cdf(1.0 - float(p) / 2.0) for p in clipped])
    return z * z


def p_to_chisq_scalar(p_value: float) -> float:
    clipped = min(max(float(p_value), np.nextafter(0, 1)), 1.0)
    if chi2 is not None:
        return float(chi2.isf(clipped, 1))
    return float(p_to_chisq(np.asarray([clipped]))[0])


def read_p_values(path: str, p_override: Optional[str]) -> Tuple[np.ndarray, int, str]:
    header, delimiter = read_header(path)
    if pl is not None and delimiter is not None:
        frame = pl.scan_csv(
            path,
            separator=delimiter,
            infer_schema_length=100,
            ignore_errors=True,
        ).head(0).collect()
        p_col = detect_col(frame.columns, p_override, P_CANDIDATES, "P-value")
        p_series = (
            pl.scan_csv(
                path,
                separator=delimiter,
                infer_schema_length=100,
                ignore_errors=True,
            )
            .select(pl.col(p_col).cast(pl.Float64, strict=False).alias("p"))
            .collect()
            .get_column("p")
        )
        return p_series.to_numpy(), len(p_series), p_col

    p_col = detect_col(header, p_override, P_CANDIDATES, "P-value")
    if delimiter is None:
        p_idx = header.index(p_col)
        try:
            values = np.loadtxt(path, skiprows=1, usecols=(p_idx,), dtype=float)
            values = np.atleast_1d(values)
            return values, int(values.shape[0]), p_col
        except Exception:
            pass

    values = []
    n_rows = 0
    for _, row in iter_rows_by_header(path):
        n_rows += 1
        try:
            values.append(float(row.get(p_col, "nan")))
        except (TypeError, ValueError):
            values.append(float("nan"))
    return np.asarray(values, dtype=float), n_rows, p_col


def read_full_plot_data(path: str, args: argparse.Namespace) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    header, delimiter = read_header(path)
    if pl is not None and delimiter is not None:
        header = pl.scan_csv(path, separator=delimiter, infer_schema_length=100, ignore_errors=True).head(0).collect().columns
        chr_col = detect_col(header, args.chr_col, CHR_CANDIDATES, "chromosome")
        bp_col = detect_col(header, args.bp_col, BP_CANDIDATES, "position")
        p_col = detect_col(header, args.p_col, P_CANDIDATES, "P-value")
        df = (
            pl.scan_csv(path, separator=delimiter, infer_schema_length=100, ignore_errors=True)
            .select(
                pl.col(chr_col).cast(pl.Utf8, strict=False).alias("chr"),
                pl.col(bp_col).cast(pl.Float64, strict=False).alias("bp"),
                pl.col(p_col).cast(pl.Float64, strict=False).alias("p"),
            )
            .collect()
        )
        return (
            df.get_column("chr").to_numpy(),
            df.get_column("bp").to_numpy(),
            df.get_column("p").to_numpy(),
        )

    chr_col = detect_col(header, args.chr_col, CHR_CANDIDATES, "chromosome")
    bp_col = detect_col(header, args.bp_col, BP_CANDIDATES, "position")
    p_col = detect_col(header, args.p_col, P_CANDIDATES, "P-value")
    if delimiter is None:
        chr_idx = header.index(chr_col)
        bp_idx = header.index(bp_col)
        p_idx = header.index(p_col)
        try:
            raw = np.loadtxt(path, skiprows=1, usecols=(chr_idx, bp_idx, p_idx), dtype=str)
            raw = np.atleast_2d(raw)
            return raw[:, 0], raw[:, 1].astype(float), raw[:, 2].astype(float)
        except Exception:
            pass

    chrs, bps, ps = [], [], []
    for _, row in iter_rows_by_header(path):
        chrs.append(row.get(chr_col, ""))
        try:
            bps.append(float(row.get(bp_col, "nan")))
        except (TypeError, ValueError):
            bps.append(float("nan"))
        try:
            ps.append(float(row.get(p_col, "nan")))
        except (TypeError, ValueError):
            ps.append(float("nan"))
    return np.asarray(chrs), np.asarray(bps, dtype=float), np.asarray(ps, dtype=float)


def expected_qq(n: int) -> np.ndarray:
    ranks = np.arange(1, n + 1, dtype=float)
    return -np.log10((ranks - 0.5) / n)


def observed_qq(p_values: np.ndarray) -> np.ndarray:
    return -np.log10(np.sort(p_values))


def valid_p_values(p_values: np.ndarray) -> np.ndarray:
    return p_values[np.isfinite(p_values) & (p_values > 0) & (p_values <= 1)]


def sort_p_values(p_values: np.ndarray, args: argparse.Namespace) -> np.ndarray:
    p = valid_p_values(p_values)
    if len(p) == 0:
        return p
    use_gpu = args.qq_engine == "gpu" or (args.qq_engine == "auto" and cp is not None and len(p) >= 1_000_000)
    if use_gpu:
        if cp is None:
            raise RuntimeError("Requested --qq-engine gpu, but CuPy is not installed.")
        return cp.asnumpy(cp.sort(cp.asarray(p)))
    return np.sort(p)


def qq_plot_indices(n: int, max_points: int) -> np.ndarray:
    if n <= max_points:
        return np.arange(n)
    tail_n = min(max_points // 2, 50000, n)
    body_n = max_points - tail_n
    body = np.unique(np.linspace(0, n - tail_n - 1, max(1, body_n), dtype=int))
    tail = np.arange(max(0, n - tail_n), n, dtype=int)
    return np.unique(np.concatenate([body, tail]))


def qq_ci_band(n: int, idx: np.ndarray, confidence: float) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    if chi2 is None:
        return None, None
    try:
        from scipy.stats import beta
    except ImportError:
        return None, None
    alpha = 1.0 - confidence
    ranks = idx + 1
    lower_p = beta.ppf(alpha / 2.0, ranks, n + 1 - ranks)
    upper_p = beta.ppf(1.0 - alpha / 2.0, ranks, n + 1 - ranks)
    lower_p = np.clip(lower_p, np.nextafter(0, 1), 1.0)
    upper_p = np.clip(upper_p, np.nextafter(0, 1), 1.0)
    # P-value bounds invert on the -log10 scale.
    return -np.log10(upper_p), -np.log10(lower_p)


def qq_metrics_from_sorted(sorted_p: np.ndarray) -> Dict[str, object]:
    p = sorted_p
    n = len(p)
    if n < 10:
        return {
            "qq_body_median_delta": None,
            "qq_body_rmse": None,
            "qq_tail_lift": None,
            "qq_shape": "UNCLASSIFIED",
            "qq_interpretation": "Too few valid P values for QQ curve classification.",
        }
    exp = expected_qq(n)
    obs = -np.log10(p)
    delta = obs - exp
    body_mask = (exp >= 0.1) & (exp <= 2.0)
    if not np.any(body_mask):
        body_mask = np.arange(n) < max(1, int(n * 0.9))
    tail_mask = exp > 2.0
    body_delta = delta[body_mask]
    body_median_delta = float(np.nanmedian(body_delta))
    body_rmse = float(np.sqrt(np.nanmean(body_delta * body_delta)))
    qq_tail_lift = float(np.nanmax(delta[tail_mask])) if np.any(tail_mask) else None

    shape = "QQ_IDEAL"
    interp = "QQ body follows the null expectation; tail deviation, if present, is compatible with possible association signal."
    if math.isfinite(body_median_delta) and body_median_delta > 0.25:
        shape = "QQ_INFLATED"
        interp = "QQ body is globally above the diagonal, suggesting population structure, relatedness, batch effects, or phenotype confounding."
    elif math.isfinite(body_median_delta) and body_median_delta < -0.25:
        shape = "QQ_DEFLATED"
        interp = "QQ body is below the diagonal, suggesting overcorrection, overfitted kinship, or overly aggressive phenotype residualization."
    elif math.isfinite(body_rmse) and body_rmse > 0.35:
        shape = "QQ_NOISY"
        interp = "QQ curve is irregular; review sample size, model convergence, and missingness."
    elif qq_tail_lift is not None and math.isfinite(qq_tail_lift) and qq_tail_lift > 1.0:
        shape = "QQ_GOOD_WITH_TAIL_SIGNAL"
        interp = "QQ body is acceptable with strong tail deviation; this can be a real association signal if Manhattan plot is localized."

    return {
        "qq_body_median_delta": body_median_delta,
        "qq_body_rmse": body_rmse,
        "qq_tail_lift": qq_tail_lift,
        "qq_shape": shape,
        "qq_interpretation": interp,
    }


def qq_metrics(p_values: np.ndarray, args: argparse.Namespace) -> Dict[str, object]:
    return qq_metrics_from_sorted(sort_p_values(p_values, args))


def classify_metrics(p_values: np.ndarray, n_rows: int, args: argparse.Namespace) -> Dict[str, object]:
    valid = np.isfinite(p_values) & (p_values > 0) & (p_values <= 1)
    p_valid = p_values[valid]
    n_valid = int(len(p_valid))
    n_invalid = int(n_rows - n_valid)
    reasons: List[str] = []
    status = "PASS"

    if args.assessment_mode == "calculate":
        status = "CALCULATED"
        reasons.append("calculate_only_no_threshold_classification")
    elif args.assessment_mode == "sample":
        reasons.append("sample_mode_skip_count_thresholds")
    elif n_valid < args.min_snps_fail:
        status = "FAIL"
        reasons.append("too_few_valid_snps_fail")
    elif n_valid < args.min_snps_warn:
        status = "WARN"
        reasons.append("too_few_valid_snps_warn")

    if args.assessment_mode != "calculate" and n_rows > 0 and n_invalid / n_rows > 0.05:
        status = "FAIL" if n_invalid / n_rows > 0.20 else ("WARN" if status == "PASS" else status)
        reasons.append("many_invalid_p")

    if n_valid:
        median_p = float(np.nanmedian(p_valid))
        med_chisq = p_to_chisq_scalar(median_p)
        lambda_gc = med_chisq / 0.454936423119572
        min_p = float(np.nanmin(p_valid))
        n_sig = int(np.sum(p_valid < args.sig_threshold))
    else:
        med_chisq = None
        lambda_gc = None
        min_p = None
        n_sig = 0

    if args.assessment_mode != "calculate" and lambda_gc is not None and math.isfinite(lambda_gc):
        if lambda_gc >= args.lambda_fail_high:
            status = "FAIL"
            reasons.append("lambda_gc_high_fail")
        elif lambda_gc >= args.lambda_warn_high and status == "PASS":
            status = "WARN"
            reasons.append("lambda_gc_high_warn")
        elif lambda_gc <= args.lambda_fail_low:
            status = "FAIL"
            reasons.append("lambda_gc_low_fail")
        elif lambda_gc <= args.lambda_warn_low and status == "PASS":
            status = "WARN"
            reasons.append("lambda_gc_low_warn")

    if args.assessment_mode == "calculate":
        pass
    elif args.assessment_mode == "sample":
        pass
    elif n_valid and n_sig > max(50, math.ceil(n_valid * 0.001)):
        if status == "PASS":
            status = "WARN"
        reasons.append("many_significant")

    reason_text = ";".join(dict.fromkeys(reasons)) if reasons else "none"
    return {
        "status": status,
        "reasons": reason_text,
        "n_rows": int(n_rows),
        "n_valid_p": n_valid,
        "n_invalid_p": n_invalid,
        "min_p": min_p,
        "n_significant": n_sig,
        "lambda_gc": lambda_gc,
        "median_chisq": med_chisq,
    }


def recommendation_for(status: str, reasons: str) -> str:
    if status == "CALCULATED":
        return "Metrics were calculated without threshold-based classification. Use full-file mode for final QC classification."
    if "sample_mode_skip_count_thresholds" in reasons:
        return "Sample-mode metrics were calculated without full-file SNP-count or excess-hit thresholds. Use this for pilot checks, then rerun full files for final QC."
    if status == "ERROR":
        return "Check whether the file exists, is readable, and contains required MLMA columns."
    if "lambda_gc_high" in reasons:
        return "Inspect population structure, GRM/kinship, batch covariates, phenotype outliers, and sample ID matching. Rerun with appropriate PCs/covariates or analyze populations separately before meta-analysis."
    if "lambda_gc_low" in reasons:
        return "Check for overcorrection from redundant covariates, overly aggressive residualization, or an overfitted kinship model."
    if "too_few_valid_snps" in reasons or "too_few_snps" in reasons:
        return "Check genotype filtering, imputation INFO/MAF thresholds, file truncation, chromosome inclusion, and model convergence logs."
    if "many_invalid_p" in reasons:
        return "Check model convergence, phenotype missingness, sample ID matching, and malformed P-value output."
    if "many_significant" in reasons:
        return "Check phenotype-batch association, duplicated/related samples, and trait outliers. Confirm signals are not spread uniformly across the genome."
    if status == "WARN":
        return "Review the QQ plot and compare against other populations for the same trait."
    return "No automated issue detected. Keep genotype-level QC and model logs for audit."


def quality_class_for(status: str, reasons: str, qq_shape: Optional[str]) -> str:
    if status == "CALCULATED":
        return "METRICS_ONLY"
    if status == "ERROR":
        return "ERROR_INPUT"
    if "too_few_valid_snps_fail" in reasons or (status == "FAIL" and "many_invalid_p" in reasons):
        return "FAIL_DATA_COMPLETENESS"
    if "too_few_valid_snps_warn" in reasons or "many_invalid_p" in reasons:
        return "WARN_DATA_COMPLETENESS"
    if "lambda_gc_high" in reasons or qq_shape == "QQ_INFLATED":
        return "FAIL_OR_WARN_INFLATION"
    if "lambda_gc_low" in reasons or qq_shape == "QQ_DEFLATED":
        return "FAIL_OR_WARN_DEFLATION"
    if "many_significant" in reasons:
        return "WARN_EXCESS_SIGNAL"
    if qq_shape == "QQ_GOOD_WITH_TAIL_SIGNAL":
        return "PASS_WITH_TAIL_SIGNAL"
    if qq_shape == "QQ_NOISY":
        return "WARN_NOISY_QQ"
    if status == "PASS":
        return "PASS_IDEAL"
    return "REVIEW"


def solution_for(status: str, reasons: str, qq_shape: Optional[str]) -> str:
    if status == "CALCULATED":
        return "No action assigned because this run only calculated metrics. Rerun with --assessment-mode sample for pilot classification or full for final classification."
    if status == "ERROR":
        return "File or column problem: verify file path, required MLMA columns, delimiters, and whether the GWAS job finished successfully."
    if "too_few_valid_snps" in reasons or "many_invalid_p" in reasons:
        return "Data completeness problem: check genotype filtering, imputation INFO/MAF thresholds, sample ID matching, chromosome inclusion, and model convergence logs; rerun affected GWAS after fixing input completeness."
    if "lambda_gc_high" in reasons or qq_shape == "QQ_INFLATED":
        return "Inflation problem: add/check PCA covariates, batch/plate/field covariates, sex/age if relevant, GRM/kinship construction, duplicate/close relatives, and phenotype outliers; rerun by population and meta-analyze if structure differs."
    if "lambda_gc_low" in reasons or qq_shape == "QQ_DEFLATED":
        return "Deflation problem: reduce redundant covariates, inspect overfitted kinship or overly aggressive residualization, and verify phenotype transformation did not remove true genetic signal."
    if "many_significant" in reasons:
        return "Excess signal problem: test phenotype association with batch/plate/family/source variables, inspect outliers, and use full Manhattan plots to confirm whether peaks are localized rather than genome-wide artifacts."
    if qq_shape == "QQ_GOOD_WITH_TAIL_SIGNAL":
        return "Likely acceptable with candidate signal: draw full Manhattan+QQ for this population/trait, check whether signals are localized, then annotate candidate loci."
    if qq_shape == "QQ_NOISY":
        return "Noisy QQ problem: check sample size, missingness, convergence, low MAC variants, and rerun with stricter variant filters if needed."
    if status == "PASS":
        return "Accept for downstream review; keep model logs and genotype-level QC records."
    return "Review QQ plot together with numeric QC metrics; rerun full Manhattan+QQ for this group/trait before deciding."


def make_qq_plot(sorted_p: np.ndarray, title: str, output: str, args: argparse.Namespace, lambda_gc: Optional[float]) -> None:
    plt = load_pyplot()
    p = sorted_p
    n = len(p)
    if n == 0:
        raise ValueError("No valid P values available for QQ plot.")
    idx = qq_plot_indices(n, args.qq_max_points)
    exp_all = expected_qq(n)
    exp = exp_all[idx]
    obs = -np.log10(p[idx])
    axis_max = float(np.nanmax([np.nanmax(exp), np.nanmax(obs)]))
    ci_idx = qq_plot_indices(n, min(args.qq_ci_points, args.qq_max_points))
    ci_x = exp_all[ci_idx]
    ci_lower, ci_upper = qq_ci_band(n, ci_idx, args.qq_confidence)

    plt.figure(figsize=(6.2, 6.0), dpi=args.plot_dpi)
    if ci_lower is not None and ci_upper is not None:
        plt.fill_between(ci_x, ci_lower, ci_upper, color="#D8DEE9", alpha=0.72, linewidth=0, label=f"{int(args.qq_confidence * 100)}% CI")
    plt.scatter(exp, obs, s=3.2, alpha=0.58, color="#0072B2", linewidths=0, label="Observed")
    plt.plot([0, axis_max], [0, axis_max], color="#444444", linewidth=0.9, label="Expected")
    plt.xlim(0, axis_max)
    plt.ylim(0, axis_max)
    plt.gca().set_aspect("equal", adjustable="box")
    lambda_text = "NA" if lambda_gc is None else f"{lambda_gc:.3f}"
    plt.title(title, fontsize=10, pad=10)
    plt.text(
        0.04,
        0.96,
        f"lambdaGC = {lambda_text}\nN = {n:,}",
        transform=plt.gca().transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#BBBBBB", "alpha": 0.86},
    )
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.legend(loc="lower right", frameon=False, fontsize=8)
    plt.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output)
    plt.close()


def normalize_chr(values: np.ndarray) -> np.ndarray:
    out = []
    for value in values:
        text = str(value).strip().lower().replace("chr", "")
        try:
            out.append(int(float(text)))
        except ValueError:
            out.append(-1)
    return np.asarray(out, dtype=int)


def downsample_for_plot(chr_arr: np.ndarray, bp: np.ndarray, p: np.ndarray, args: argparse.Namespace) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    valid = (chr_arr > 0) & np.isfinite(bp) & np.isfinite(p) & (p > 0) & (p <= 1)
    chr_arr, bp, p = chr_arr[valid], bp[valid], p[valid]
    n = len(p)
    if n <= args.max_plot_points:
        return chr_arr, bp, p
    keep = p <= args.plot_keep_p
    remaining = np.where(~keep)[0]
    n_extra = max(0, args.max_plot_points - int(np.sum(keep)))
    if n_extra < len(remaining):
        rng = np.random.default_rng(1)
        chosen = rng.choice(remaining, size=n_extra, replace=False)
        idx = np.concatenate([np.where(keep)[0], chosen])
    else:
        idx = np.arange(n)
    return chr_arr[idx], bp[idx], p[idx]


def prepare_manhattan_data(path: str, args: argparse.Namespace) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[float], List[str]]:
    chr_raw, bp, p = read_full_plot_data(path, args)
    chr_arr = normalize_chr(chr_raw)
    chr_arr, bp, p = downsample_for_plot(chr_arr, bp, p, args)
    order = np.lexsort((bp, chr_arr))
    chr_arr, bp, p = chr_arr[order], bp[order], p[order]

    offsets = {}
    cursor = 0.0
    centers = []
    labels = []
    for chrom in sorted(set(chr_arr.tolist())):
        mask = chr_arr == chrom
        if not np.any(mask):
            continue
        offsets[chrom] = cursor
        chr_len = float(np.nanmax(bp[mask]))
        centers.append(cursor + chr_len / 2.0)
        labels.append(str(chrom))
        cursor += chr_len
    pos = np.asarray([bp_i + offsets.get(chr_i, 0.0) for chr_i, bp_i in zip(chr_arr, bp)])
    return chr_arr, pos, p, centers, labels


def make_manhattan_plot(path: str, title: str, output: str, args: argparse.Namespace) -> None:
    plt = load_pyplot()
    chr_arr, pos, p, centers, labels = prepare_manhattan_data(path, args)
    y = -np.log10(p)
    plt.figure(figsize=(10.8, 4.8), dpi=args.plot_dpi)
    colors = np.where(chr_arr % 2 == 0, "#0072B2", "#E69F00")
    plt.scatter(pos, y, s=2.0, c=colors, alpha=0.68, linewidths=0)
    plt.axhline(-math.log10(args.sig_threshold), color="#B2182B", linewidth=0.9)
    plt.title(title, fontsize=10, pad=8)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(P)")
    plt.xticks(centers, labels, fontsize=7)
    plt.grid(axis="y", color="#E6E6E6", linewidth=0.5)
    plt.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output)
    plt.close()


def make_full_plot(path: str, title: str, output: str, args: argparse.Namespace) -> None:
    plt = load_pyplot()
    chr_arr, pos, p, centers, labels = prepare_manhattan_data(path, args)
    y = -np.log10(p)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), dpi=args.plot_dpi)
    colors = np.where(chr_arr % 2 == 0, "#0072B2", "#E69F00")
    axes[0].scatter(pos, y, s=1.8, c=colors, alpha=0.65, linewidths=0)
    axes[0].axhline(-math.log10(args.sig_threshold), color="#B2182B", linewidth=0.8)
    axes[0].set_title(title, fontsize=9)
    axes[0].set_xlabel("Chromosome")
    axes[0].set_ylabel("-log10(P)")
    axes[0].set_xticks(centers)
    axes[0].set_xticklabels(labels, fontsize=7)

    p_valid = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    exp = expected_qq(len(p_valid))
    obs = observed_qq(p_valid)
    axis_max = float(np.nanmax([np.nanmax(exp), np.nanmax(obs)]))
    axes[1].scatter(exp, obs, s=2, alpha=0.6, color="#3B6EA8", linewidths=0)
    axes[1].plot([0, axis_max], [0, axis_max], color="#777777", linewidth=0.8)
    axes[1].set_xlim(0, axis_max)
    axes[1].set_ylim(0, axis_max)
    axes[1].set_aspect("equal", adjustable="box")
    axes[1].set_title("QQ plot", fontsize=9)
    axes[1].set_xlabel("Expected -log10(P)")
    axes[1].set_ylabel("Observed -log10(P)")
    plt.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output)
    plt.close()


def process_one(row: Dict[str, str], args: argparse.Namespace) -> Dict[str, object]:
    population = row["population"]
    trait = row["trait"]
    file_path = row["file"]
    ident = f"{safe_name(population)}__{safe_name(trait)}"
    qq_file = str(Path(args.outdir) / "qq_plots" / safe_name(population) / f"{ident}.qq.png") if args.mode == "qq" else None
    manhattan_file = str(Path(args.outdir) / "manhattan_plots" / safe_name(population) / f"{ident}.manhattan.png") if args.mode == "manhattan" else None
    full_file = str(Path(args.outdir) / "full_plots" / safe_name(population) / f"{ident}.mqq.png") if args.mode == "full" else None
    base = {
        "population": population,
        "trait": trait,
        "file": file_path,
        "status": "ERROR",
        "reasons": None,
        "qq_plot_file": qq_file,
        "manhattan_plot_file": manhattan_file,
        "full_plot_file": full_file,
        "recommendation": None,
        "quality_class": "ERROR_INPUT",
        "solution": None,
        "assessment_mode": args.assessment_mode,
    }
    try:
        if not os.path.exists(file_path):
            raise FileNotFoundError("File does not exist")
        p_values, n_rows, _ = read_p_values(file_path, args.p_col)
        metrics = classify_metrics(p_values, n_rows, args)
        base.update(metrics)
        sorted_p = sort_p_values(p_values, args) if args.mode in {"qq", "full"} else None
        qqm = qq_metrics_from_sorted(sorted_p) if sorted_p is not None else {
            "qq_body_median_delta": None,
            "qq_body_rmse": None,
            "qq_tail_lift": None,
            "qq_shape": None,
            "qq_interpretation": None,
        }
        base.update(qqm)
        if args.mode == "qq":
            make_qq_plot(sorted_p, f"{population} {trait} {qqm['qq_shape']}", qq_file, args, base.get("lambda_gc"))
        if args.mode == "manhattan":
            make_manhattan_plot(file_path, f"{population} {trait}", manhattan_file, args)
        if args.mode == "full":
            make_full_plot(file_path, f"{population} {trait}", full_file, args)
        base["recommendation"] = recommendation_for(str(base["status"]), str(base["reasons"]))
        base["quality_class"] = quality_class_for(str(base["status"]), str(base["reasons"]), base.get("qq_shape"))
        base["solution"] = solution_for(str(base["status"]), str(base["reasons"]), base.get("qq_shape"))
        return base
    except Exception as exc:  # pragma: no cover
        base["reasons"] = str(exc)
        base["recommendation"] = recommendation_for("ERROR", str(exc))
        base["quality_class"] = "ERROR_INPUT"
        base["solution"] = solution_for("ERROR", str(exc), None)
        return base


def print_result_preview(rows: List[Dict[str, object]]) -> None:
    print("", file=sys.stderr)
    print("Metric preview:", file=sys.stderr)
    header = ["population", "trait", "status", "n_valid_p", "lambda_gc", "min_p", "n_significant", "reasons"]
    print("\t".join(header), file=sys.stderr)
    for row in rows[:20]:
        lambda_gc = row.get("lambda_gc")
        min_p = row.get("min_p")
        print(
            "\t".join(
                [
                    str(row.get("population", "")),
                    str(row.get("trait", "")),
                    str(row.get("status", "")),
                    str(row.get("n_valid_p", "")),
                    "" if lambda_gc in (None, "") else f"{float(lambda_gc):.4g}",
                    "" if min_p in (None, "") else f"{float(min_p):.4g}",
                    str(row.get("n_significant", "")),
                    str(row.get("reasons", "")),
                ]
            ),
            file=sys.stderr,
        )
    if len(rows) > 20:
        print(f"... {len(rows) - 20} more rows written to qc_summary.tsv", file=sys.stderr)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = read_manifest(args.manifest)
    if args.only_pop:
        rows = [row for row in rows if row["population"] == args.only_pop]
    if args.only_trait:
        rows = [row for row in rows if row["trait"] == args.only_trait]
    if args.limit is not None:
        rows = rows[: args.limit]
    if not rows:
        raise SystemExit("No manifest rows remain after filtering.")

    print(f"Loaded manifest rows: {len(rows)}", file=sys.stderr)
    print(f"Output directory: {outdir.resolve()}", file=sys.stderr)
    print(f"Run mode: {args.mode}", file=sys.stderr)
    print(f"Assessment mode: {args.assessment_mode}", file=sys.stderr)
    print(f"Polars enabled: {'yes' if pl is not None else 'no, using Python csv fallback'}", file=sys.stderr)
    print(f"CuPy GPU QQ sort enabled: {'available' if cp is not None else 'not available'}", file=sys.stderr)
    if args.mode in {"qq", "manhattan", "full"}:
        try:
            load_pyplot()
        except RuntimeError as exc:
            raise SystemExit(str(exc)) from exc

    if args.threads > 1 and len(rows) > 1:
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            future_to_index = {executor.submit(process_one, row, args): i for i, row in enumerate(rows)}
            results = [None] * len(rows)
            for future in as_completed(future_to_index):
                results[future_to_index[future]] = future.result()
        results = [row for row in results if row is not None]
    else:
        results = [process_one(row, args) for row in rows]

    manifest_extra = [name for name in rows[0].keys() if name not in SUMMARY_COLUMNS]
    fieldnames = manifest_extra + [name for name in SUMMARY_COLUMNS if name not in manifest_extra]
    merged = []
    for original, result in zip(rows, results):
        merged_row = dict(original)
        merged_row.update(result)
        merged.append(merged_row)

    write_tsv(outdir / "qc_summary.tsv", merged, fieldnames)
    write_tsv(outdir / "qc_pass.tsv", [row for row in merged if row["status"] == "PASS"], fieldnames)
    write_tsv(outdir / "qc_warn.tsv", [row for row in merged if row["status"] == "WARN"], fieldnames)
    write_tsv(outdir / "qc_fail.tsv", [row for row in merged if row["status"] == "FAIL"], fieldnames)
    write_tsv(outdir / "qc_calculated.tsv", [row for row in merged if row["status"] == "CALCULATED"], fieldnames)
    write_tsv(
        outdir / "qc_recommendations.tsv",
        merged,
        ["population", "trait", "file", "status", "reasons", "recommendation", "qq_plot_file", "manhattan_plot_file", "full_plot_file"],
    )
    write_tsv(
        outdir / "qc_quality_solutions.tsv",
        merged,
        [
            "population",
            "trait",
            "file",
            "status",
            "reasons",
            "quality_class",
            "lambda_gc",
            "n_valid_p",
            "n_invalid_p",
            "min_p",
            "n_significant",
            "qq_shape",
            "qq_body_median_delta",
            "qq_body_rmse",
            "qq_tail_lift",
            "qq_interpretation",
            "solution",
            "qq_plot_file",
            "manhattan_plot_file",
            "full_plot_file",
        ],
    )
    if args.mode in {"qq", "full"}:
        write_tsv(
            outdir / "qq_shape_summary.tsv",
            merged,
            [
                "population",
                "trait",
                "file",
                "status",
                "quality_class",
                "qq_shape",
                "qq_body_median_delta",
                "qq_body_rmse",
                "qq_tail_lift",
                "qq_interpretation",
                "qq_plot_file",
                "manhattan_plot_file",
                "full_plot_file",
            ],
        )
    if args.mode == "manhattan":
        write_tsv(
            outdir / "manhattan_plot_summary.tsv",
            merged,
            [
                "population",
                "trait",
                "file",
                "status",
                "quality_class",
                "lambda_gc",
                "min_p",
                "n_significant",
                "manhattan_plot_file",
            ],
        )
    print("Done.", file=sys.stderr)
    print_result_preview(merged)
    for status in ["PASS", "WARN", "FAIL", "ERROR", "CALCULATED"]:
        print(f"{status}: {sum(row['status'] == status for row in merged)}", file=sys.stderr)


if __name__ == "__main__":
    main()
