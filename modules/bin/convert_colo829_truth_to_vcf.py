#!/usr/bin/env python3
"""
Converts the parklab SMaHT COLO829BLT50/HAPMAP truth-set text files into minimal, sorted
VCFs on this project's build (GRCh37, bare numeric contigs), restricted to a chosen
chromosome list.

Why this exists: the source files (Truthset_SNV_COLO829.txt, Truthset_Indel_COLO829.txt)
encode each variant as "chr9:130613:G:A" plus a wide set of per-platform VAF columns, not
as a standard VCF. They are also on a different genome build from everything else here:
the parklab VCF of the same set declares ##reference=.../GRCh38_no_alt/hg38_no_alt.fa
(chr9 length 138394717, chr20 length 64444167), while the COLO829 BAMs, the reference
FASTA and the SV truth set are GRCh37. Stripping "chr" from the contig names is therefore
not enough. An earlier version of this script did only that, which left the truth set on
GRCh38 coordinates and collapsed SNV recall to about 5%, since only the stretch of chr9
where the two builds happen to share coordinates could ever match.

This version lifts every record from GRCh38 to GRCh37 with a chain file, drops records
that cannot be lifted cleanly, and checks every lifted REF allele against the GRCh37
FASTA. More than 1% REF mismatches stops the run, so a build mismatch cannot pass
silently again.

Usage:
    convert_colo829_truth_to_vcf.py <input.txt> <output.vcf> <SOURCE_LABEL> <9,20> \
        <grch38_to_grch37.chain.gz> <grch37_genome.fa>
"""
import sys
from collections import Counter

from pyliftover import LiftOver

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
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
# More than 1% of lifted REF alleles disagreeing with the GRCh37 FASTA means the truth
# set, the chain file or the reference is on the wrong build.
MAX_REF_MISMATCH_FRACTION = 0.01


def read_fasta(path, wanted):
    """Reads only the wanted contigs of an uncompressed FASTA into memory."""
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in wanted:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in wanted:
                buf.append(line.strip())
    if name in wanted:
        seqs[name] = "".join(buf).upper()
    return seqs


def main():
    if len(sys.argv) != 7:
        sys.exit(f"usage: {sys.argv[0]} <input.txt> <output.vcf> <SOURCE_LABEL> "
                 f"<chrom1,chrom2,...> <grch38_to_grch37.chain.gz> <grch37_genome.fa>")

    in_path, out_path, source_label, chrom_arg, chain_path, fasta_path = sys.argv[1:7]
    # chrom_arg arrives as bare numeric ("9,20") to match this project's convention
    # everywhere else, and it names the GRCh37 chromosomes to keep after liftover.
    bare_chroms = [c.strip() for c in chrom_arg.split(",")]
    unknown = [c for c in bare_chroms if c not in GRCH37_LENGTHS]
    if unknown:
        sys.exit(
            f"no known GRCh37 length for chromosome(s) {unknown}; add them to "
            f"GRCH37_LENGTHS in this script if you extend this project past chr9/chr20"
        )

    lifter = LiftOver(chain_path)
    ref37 = read_fasta(fasta_path, set(bare_chroms))
    missing = [c for c in bare_chroms if c not in ref37]
    if missing:
        sys.exit(f"chromosome(s) {missing} not found in {fasta_path}")

    stats = Counter()
    records = []
    with open(in_path) as fh:
        header = fh.readline()
        if not header.startswith("#variant"):
            sys.exit(f"unexpected header in {in_path}: {header!r}")
        for line in fh:
            if not line.strip():
                continue
            parts = line.split("\t", 1)[0].split(":")
            if len(parts) != 4:
                stats["malformed"] += 1
                continue
            chrom38, pos38, ref, alt = parts[0], int(parts[1]), parts[2].upper(), parts[3].upper()
            # Every record is lifted, not only those on chr9/chr20 in GRCh38: what matters
            # is where it lands on GRCh37. pyliftover works with 0-based coordinates.
            start = lifter.convert_coordinate(chrom38, pos38 - 1)
            end = lifter.convert_coordinate(chrom38, pos38 - 1 + len(ref) - 1)
            if not start or not end or len(start) > 1 or len(end) > 1:
                stats["unliftable"] += 1
                continue
            (c_start, p_start, strand_start, _), (c_end, p_end, strand_end, _) = start[0], end[0]
            chrom37 = c_start[3:] if c_start.startswith("chr") else c_start
            if chrom37 not in bare_chroms:
                continue
            if c_start != c_end or strand_start != strand_end:
                stats["split_by_liftover"] += 1
                continue
            if strand_start == "-":
                # A single base on the minus strand only needs complementing; an indel
                # would also need re-anchoring on the other side, so it is dropped.
                if len(ref) != 1 or len(alt) != 1:
                    stats["minus_strand_indel_dropped"] += 1
                    continue
                ref, alt = ref.translate(COMPLEMENT), alt.translate(COMPLEMENT)
            elif p_end - p_start != len(ref) - 1:
                stats["span_changed_by_liftover"] += 1
                continue
            pos37 = p_start + 1
            stats["lifted"] += 1
            if ref37[chrom37][pos37 - 1:pos37 - 1 + len(ref)] != ref:
                stats["ref_mismatch_after_liftover"] += 1
                continue
            records.append((chrom37, pos37, ref, alt))

    lifted = stats["lifted"]
    if lifted == 0 or stats["ref_mismatch_after_liftover"] / lifted > MAX_REF_MISMATCH_FRACTION:
        sys.exit(f"REF check failed after liftover ({dict(stats)}): the truth set and the "
                 f"reference are on different builds, or the chain file is wrong")

    order = {c: i for i, c in enumerate(GRCH37_LENGTHS)}
    records = sorted(set(records), key=lambda r: (order[r[0]], r[1], r[2], r[3]))

    with open(out_path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##source=parklab_SMaHT_COLO829BLT50_HAPMAP_truthset:{source_label}"
                  f";lifted_GRCh38_to_GRCh37\n")
        # Required, not optional: bcftools (and most other VCF tools) will refuse to sort
        # or index a file that references a contig in its data rows without that contig
        # having been declared in the header first.
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
        # actual comparison happens. This truth set does not give per-position zygosity,
        # so 0/1 is used uniformly for every record, a standard, defensible
        # simplification for a presence/absence truth set: hap.py's primary TP/FP/FN
        # accounting only needs "this variant is present in the truth sample," not an
        # exact zygosity match.
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n")
        for chrom, pos, ref, alt in records:
            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tSOURCE={source_label}\tGT\t0/1\n")

    print(
        f"{in_path}: wrote {len(records)} {source_label} records on GRCh37 "
        f"{bare_chroms}; {dict(stats)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
