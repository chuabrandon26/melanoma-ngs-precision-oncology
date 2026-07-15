#!/usr/bin/env python3
"""
Converts the parklab SMaHT COLO829BLT50/HAPMAP truth-set text files into minimal, sorted
VCFs, restricted to a chosen chromosome list, with the "chr" prefix stripped.

Why this exists: the source files (Truthset_SNV_COLO829.txt, Truthset_Indel_COLO829.txt)
encode each variant as "chr9:130613:G:A" plus a wide set of per-platform VAF columns, not
as a standard VCF, and they use "chr9"/"chr20" contig names. Every other file in this
project (the COLO829 BAMs, the reference FASTA, the SV truth set) uses GRCh37's bare
numeric convention ("9", "20"). Comparing a query VCF against this truth set with hap.py
requires matching contig names, so the rename happens once, here, rather than being
handled ad hoc by each benchmarking step.

Usage:
    convert_colo829_truth_to_vcf.py <input.txt> <output.vcf> <SOURCE_LABEL> <chr9,chr20>
"""
import sys

# GRCh37 chromosome lengths, bare-numeric contig names, karyotypic order. These are the
# exact values read directly from the COLO829 BAM header (@SQ SN:<n> LN:<length>) during
# this project's initial data-sizing pass, not generic GRCh37 assumptions, so they are
# guaranteed to match this project's reference FASTA and BAMs exactly.
GRCH37_LENGTHS = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
    "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
    "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
    "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560, "Y": 59373566,
    "MT": 16569,
}


def main():
    if len(sys.argv) != 5:
        sys.exit(f"usage: {sys.argv[0]} <input.txt> <output.vcf> <SOURCE_LABEL> <chrom1,chrom2,...>")

    in_path, out_path, source_label, chrom_arg = sys.argv[1:5]
    # chrom_arg arrives as bare numeric ("9,20") to match this project's convention
    # everywhere else; the source file itself uses "chr9"/"chr20", so add the prefix
    # back on just for matching against the file's own variant column.
    bare_chroms = [c.strip() for c in chrom_arg.split(",")]
    wanted_chroms = {f"chr{c}" for c in bare_chroms}

    unknown = [c for c in bare_chroms if c not in GRCH37_LENGTHS]
    if unknown:
        sys.exit(
            f"no known GRCh37 length for chromosome(s) {unknown}; add them to "
            f"GRCH37_LENGTHS in this script if you extend this project past chr9/chr20"
        )

    records = []
    skipped_malformed = 0
    with open(in_path) as fh:
        header = fh.readline()
        if not header.startswith("#variant"):
            sys.exit(f"unexpected header in {in_path}: {header!r}")
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            variant = fields[0]
            parts = variant.split(":")
            if len(parts) != 4:
                skipped_malformed += 1
                continue
            chrom, pos, ref, alt = parts
            if chrom not in wanted_chroms:
                continue
            records.append((chrom[3:], int(pos), ref, alt))  # chrom[3:] strips "chr"

    records.sort(key=lambda r: (r[0], r[1]))

    with open(out_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=parklab_SMaHT_COLO829BLT50_HAPMAP_truthset:{source_label}\n")
        # Required, not optional: bcftools (and most other VCF tools) will refuse to sort
        # or index a file that references a contig in its data rows without that contig
        # having been declared in the header first. Found by actually running this
        # against bcftools sort, not caught by testing the script's own row output alone.
        # Karyotypic order (from GRCH37_LENGTHS' own definition order), not whatever
        # order chrom_arg happened to list chromosomes in.
        for chrom in GRCH37_LENGTHS:
            if chrom in bare_chroms:
                out.write(f"##contig=<ID={chrom},length={GRCH37_LENGTHS[chrom]}>\n")
        out.write(
            '##INFO=<ID=SOURCE,Number=1,Type=String,'
            f'Description="Truth-set variant class, {source_label}">\n'
        )
        # hap.py's comparison engines (vcfeval in particular) require a genotyped sample
        # column to operate on, not just a bare list of positions/alleles: a truth VCF
        # with no FORMAT/sample column fails with "Input file has no samples" before any
        # actual comparison happens. Found by actually running hap.py against this
        # file's earlier, sample-less version, not anticipated in advance. This truth set
        # does not give per-position zygosity, so 0/1 is used uniformly for every
        # record, a standard, defensible simplification for a presence/absence truth
        # set: hap.py's primary TP/FP/FN accounting only needs "this variant is present
        # in the truth sample," not an exact zygosity match.
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n")
        for chrom, pos, ref, alt in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tSOURCE={source_label}\tGT\t0/1\n")

    print(
        f"{in_path}: wrote {len(records)} {source_label} records restricted to "
        f"{sorted(wanted_chroms)} ({skipped_malformed} malformed lines skipped)",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
