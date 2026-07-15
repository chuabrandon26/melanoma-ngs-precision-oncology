#!/usr/bin/env python3
"""
Converts a CNVkit .call.cns file to a BED of non-neutral (gain/loss) segments, looking up
the "cn" column by name rather than assuming a fixed column position. CNVkit's exact
column layout is not perfectly stable across versions and call options (e.g. whether
--purity/--ploidy were set, or whether allele-specific cn1/cn2 columns are present), so
a hardcoded column index is a real, easy-to-miss source of silently wrong output, the
same class of bug this project already hit once with the SV truth VCF's END field.

Usage: cns_to_nonneutral_bed.py <input.call.cns> <output.bed> <expected_ploidy>
"""
import csv
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit(f"usage: {sys.argv[0]} <input.call.cns> <output.bed> <expected_ploidy>")

    in_path, out_path, ploidy = sys.argv[1], sys.argv[2], int(sys.argv[3])
    with open(in_path) as fh, open(out_path, "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "cn" not in reader.fieldnames:
            sys.exit(
                f"no 'cn' column found in {in_path} (found columns: {reader.fieldnames}); "
                "this usually means `cnvkit.py call` was not run, or its output format changed"
            )
        n_nonneutral = 0
        n_total = 0
        for row in reader:
            n_total += 1
            cn = int(row["cn"])
            if cn == ploidy:
                continue
            out.write(f"{row['chromosome']}\t{row['start']}\t{row['end']}\t{cn}\n")
            n_nonneutral += 1

    print(
        f"{in_path}: {n_nonneutral}/{n_total} segments are non-neutral (cn != {ploidy}), "
        f"written to {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
