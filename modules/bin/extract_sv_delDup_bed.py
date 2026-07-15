#!/usr/bin/env python3
"""
Extracts DEL/DUP entries from the COLO829 SV truth VCF as a BED file, above a minimum
size threshold. Used as an approximate copy-number-change reference in Module 2's CNV
comparison, since no independently validated CNV truth set exists for COLO829 the way
one does for SNVs, indels, and SVs. A true deletion or tandem duplication above CNVkit's
detection floor should also show up as a copy-number change, so this is a reasonable,
if imperfect, proxy. It is not a substitute for a purpose-built CNV truth set: a balanced
rearrangement would be missed by copy-number callers entirely and correctly so, while a
DEL/DUP call here says nothing about the precise copy-number state CNVkit should report,
only that something changed in that region.

Usage: extract_sv_delDup_bed.py <truth_sv.vcf> <output.bed> <min_size_bp>
"""
import re
import sys


def parse_info(info_str):
    info = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            info[k] = v
        else:
            info[field] = True
    return info


def main():
    if len(sys.argv) != 4:
        sys.exit(f"usage: {sys.argv[0]} <truth_sv.vcf> <output.bed> <min_size_bp>")

    vcf_path, bed_path, min_size = sys.argv[1], sys.argv[2], int(sys.argv[3])
    n_written = 0
    with open(vcf_path) as fh, open(bed_path, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, variant_id, info_str = fields[0], int(fields[1]), fields[2], fields[7]
            info = parse_info(info_str)
            svtype = info.get("SVTYPE", "")
            if svtype not in ("DEL", "DUP"):
                continue
            # Every DEL/DUP in this truth set is written as a pair of breakend (BND-style)
            # mate records sharing one MATEID, e.g. truthset_36_1 / truthset_36_2, both
            # carrying the same SVLEN/SVTYPE. Processing both would double-count each event
            # and, worse, the second mate's own POS + SVLEN overshoots far past the real
            # end (verified against this file: mate "_2" of a 123,393bp DEL sits at the
            # true end coordinate, so anchoring a second interval there and adding SVLEN
            # again lands 123kb beyond the actual deletion). Empirically, across all 45
            # DEL/DUP pairs in this file, the "_1" mate always has the smaller coordinate,
            # so keeping only "_1" records gives exactly one, correctly-anchored interval
            # per event.
            if not variant_id.endswith("_1"):
                continue
            try:
                svlen = abs(int(info.get("SVLEN", "0")))
            except ValueError:
                svlen = 0
            if svlen < min_size:
                continue
            # This truth set encodes SVs as breakend (BND-style) mate pairs, e.g.
            # ALT=C[9:28157692[, not as a symbolic <DEL> with an INFO/END field, so END
            # is usually absent. Falling back to "END defaults to POS" (as if the event
            # were 1bp wide) silently produces near-empty regions for exactly the large
            # events this filter is meant to keep, verified against this exact file: a
            # 123,393bp DEL on chr9 has no INFO/END at all and would collapse to 1bp
            # without this reconstruction. When END truly is absent, reconstruct it from
            # SVLEN instead, which is present on every DEL/DUP record in this file.
            if "END" in info:
                end = int(info["END"])
            else:
                end = pos + svlen
            out.write(f"{chrom}\t{pos - 1}\t{end}\t{svtype}\n")
            n_written += 1

    print(f"{vcf_path}: wrote {n_written} DEL/DUP regions >= {min_size}bp to {bed_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
