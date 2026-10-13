#!/usr/bin/env python3
"""
Re-scores committed Mutect2 calls against the (GRCh37-lifted) COLO829 SNV/indel truth set
without re-running the pipeline, and writes a hap.py-style PASS summary CSV.

Why this exists: the SNV/indel benchmark originally compared GRCh37 calls against a
GRCh38 truth set that had only been renamed, not lifted over (fixed in
modules/bin/convert_colo829_truth_to_vcf.py). Variant calling itself was unaffected, so
the calls already in data/processed/ can be scored again against the corrected truth set.

It runs RTG vcfeval, the engine hap.py uses in modules/benchmark.nf, on PASS calls of the
tumor sample, then counts the result the way hap.py's summary does: an MNP counts as its
individual SNPs, and every non-length-preserving record counts as one indel.

Usage:
    rescore_snv_indel.py --rtg <rtg> --sdf <grch37.sdf> --truth <truth_snv_indel.vcf.gz> \
        --calls <filtered.vcf.gz> --bed <regions.bed> --out-dir <vcfeval_dir> \
        --summary <out.csv> [--tumor-sample COLO829T]
"""
import argparse
import csv
import gzip
import shutil
import subprocess
from pathlib import Path


def tally(vcf_gz):
    """Counts SNPs (MNPs split into their SNPs) and indels the way hap.py reports them."""
    snp = indel = 0
    with gzip.open(vcf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            ref, alts = fields[3], fields[4].split(",")
            for alt in alts:
                if len(alt) == len(ref):
                    snp += sum(1 for r, a in zip(ref, alt) if r != a)
                else:
                    indel += 1
    return {"SNP": snp, "INDEL": indel}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--rtg", required=True, help="path to the rtg executable (RTG Tools 3.13)")
    ap.add_argument("--sdf", required=True, help="GRCh37 reference in RTG SDF format (rtg format genome.fa)")
    ap.add_argument("--truth", required=True, help="lifted truth_snv_indel.vcf.gz (sample TRUTH)")
    ap.add_argument("--calls", required=True, help="Mutect2 *.filtered.vcf.gz")
    ap.add_argument("--bed", required=True, help="comparison regions (exome BED or whole chr9+chr20)")
    ap.add_argument("--out-dir", required=True, help="vcfeval output directory (overwritten)")
    ap.add_argument("--summary", required=True, help="summary CSV to write")
    ap.add_argument("--tumor-sample", default="COLO829T")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)
    # --squash-ploidy matches presence/absence rather than exact zygosity, since the truth
    # set gives every record a placeholder 0/1 genotype. vcfeval only scores PASS calls
    # unless --all-records is given.
    subprocess.run(
        [args.rtg, "vcfeval", "-b", args.truth, "-c", args.calls, "-t", args.sdf,
         "--bed-regions", args.bed, "--sample", f"TRUTH,{args.tumor_sample}",
         "--squash-ploidy", "--ref-overlap", "-o", str(out_dir)],
        check=True,
    )

    tp_base = tally(out_dir / "tp-baseline.vcf.gz")
    fn = tally(out_dir / "fn.vcf.gz")
    tp_call = tally(out_dir / "tp.vcf.gz")
    fp = tally(out_dir / "fp.vcf.gz")

    rows = []
    for vtype in ("INDEL", "SNP"):
        truth_total = tp_base[vtype] + fn[vtype]
        query_total = tp_call[vtype] + fp[vtype]
        recall = tp_base[vtype] / truth_total if truth_total else float("nan")
        precision = tp_call[vtype] / query_total if query_total else float("nan")
        f1 = 2 * recall * precision / (recall + precision) if recall + precision else 0.0
        rows.append({
            "Type": vtype, "Filter": "PASS",
            "TRUTH.TOTAL": truth_total, "TRUTH.TP": tp_base[vtype], "TRUTH.FN": fn[vtype],
            "QUERY.TOTAL": query_total, "QUERY.TP": tp_call[vtype], "QUERY.FP": fp[vtype],
            "METRIC.Recall": f"{recall:.6f}", "METRIC.Precision": f"{precision:.6f}",
            "METRIC.F1_Score": f"{f1:.6f}",
        })
    with open(args.summary, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['Type']:5s} recall={float(row['METRIC.Recall']):.1%} "
              f"precision={float(row['METRIC.Precision']):.1%} F1={float(row['METRIC.F1_Score']):.1%} "
              f"(TP {row['TRUTH.TP']} / FP {row['QUERY.FP']} / FN {row['TRUTH.FN']}, truth {row['TRUTH.TOTAL']})")


if __name__ == "__main__":
    main()
