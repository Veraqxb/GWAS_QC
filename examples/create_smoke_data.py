#!/usr/bin/env python3
"""Create a tiny .mlma dataset for fast QC smoke tests."""

from __future__ import annotations

import csv
import random
import sys
from pathlib import Path


def main() -> None:
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("smoke_data")
    outdir.mkdir(parents=True, exist_ok=True)
    mlma = outdir / "smoke.mlma"
    manifest = outdir / "manifest.tsv"

    random.seed(1)
    n = 12000
    with open(mlma, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Chr", "SNP", "bp", "A1", "A2", "Freq", "b", "se", "p"])
        for i in range(1, n + 1):
            chrom = 1 if i <= n // 2 else 2
            bp = i if chrom == 1 else i - n // 2
            writer.writerow([
                chrom,
                f"rs{i}",
                bp,
                "A",
                "G",
                f"{random.uniform(0.05, 0.95):.6g}",
                f"{random.gauss(0, 1):.6g}",
                f"{random.uniform(0.05, 0.2):.6g}",
                f"{random.random():.12g}",
            ])

    with open(manifest, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["population", "trait", "file", "sample_size"])
        writer.writerow(["pop200", "smoke", str(mlma.resolve()), 200])

    print(f"Wrote {manifest.resolve()}")


if __name__ == "__main__":
    main()
