#!/usr/bin/env bash
# Downloads the small, static reference and truth-set files used by all three pipelines.
# Every file pulled here was size-checked before being added to this script, and all of
# them are well under the 2GB threshold that would otherwise need a sizing check first:
#   - chr9 + chr20 FASTA (GRCh37):     ~36.6MB + ~18.1MB
#   - GENCODE v19 GTF (GRCh37):        ~38.0MB
#   - SV truth set VCF:                 24.6KB
#   - SNV truth set:                    10.1MB
#   - Indel truth set:                  302KB
#   - BICseq2 copy-number reference:    15.0KB (zipped)
# None of this needs more than a couple of minutes even on a slow connection.
#
# This script does NOT touch ERR2752450 / ERR2752449 (the COLO829T/COLO829BL BAMs).
# Those are multi-GB and are fetched separately with an explicit remote-slicing step
# that the user runs themselves, see pipelines/wgs/main.nf and the top-level README.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="$REPO_ROOT/data/reference"
mkdir -p "$REF_DIR/truth_sets" "$REF_DIR/cnv_reference"

echo "== Reference genome: GRCh37, chr9 and chr20 only =="
# Why GRCh37, and why bare numeric contig names ("9", "20") rather than UCSC-style
# ("chr9", "chr20"): the COLO829T/COLO829BL BAMs on ENA (PRJEB27698) were aligned to
# GRCh37 with exactly this contig naming, confirmed by reading the BAM header directly
# (@SQ SN:20 LN:63025520, which matches GRCh37's chr20 length exactly, not GRCh38's).
# All three COLO829 truth sets (SV, SNV, indel) are also GRCh37-based. Using Ensembl's
# GRCh37 archive with its native numeric naming avoids a contig-renaming or liftover
# step for the reference genome itself. Source: Ensembl GRCh37 archive (release is
# whatever "current" resolves to on the GRCh37 archive site, which only receives
# occasional annotation updates, not reference sequence changes).
curl -L --fail -o "$REF_DIR/Homo_sapiens.GRCh37.dna.chromosome.9.fa.gz" \
  "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.9.fa.gz"
curl -L --fail -o "$REF_DIR/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz" \
  "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.20.fa.gz"

echo "== Gene annotation: GENCODE v19 (GRCh37) =="
# Used to build the exome-intersection BED for Module 1 (we derive an "exonic regions"
# BED from GENCODE exon features rather than a specific commercial capture kit, since
# kits like Agilent SureSelect require a vendor account; see README for the tradeoff)
# and to build the transcript-to-gene map used by tximport in the RNA-seq notebook.
curl -L --fail -o "$REF_DIR/gencode.v19.annotation.gtf.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"

echo "== Transcriptome reference for Module 3 (Salmon), genome-wide, not chr9/20-only =="
# Protein-coding + lncRNA transcript sequences only, not every GENCODE biotype, and no
# genome decoy sequence: both are disk-driven simplifications explained in
# modules/quantify.nf, chosen to avoid a multi-GB genome download and a much larger
# Salmon index on top of everything else this project already asks of a tight disk
# budget.
curl -L --fail -o "$REF_DIR/gencode.v19.pc_transcripts.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz"
curl -L --fail -o "$REF_DIR/gencode.v19.lncRNA_transcripts.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz"

echo "== SV truth set (UMCUGenetics / Valle-Inclan et al. 2022, Cell Genomics) =="
# Multi-platform-validated (Illumina, PacBio, Oxford Nanopore, 10x, BioNano). GRCh37,
# bare numeric contigs already, matches the BAMs and the reference genome above directly.
curl -L --fail -o "$REF_DIR/truth_sets/truthset_somaticSVs_COLO829.vcf" \
  "https://raw.githubusercontent.com/UMCUGenetics/COLO829_somaticSV/master/truthset_somaticSVs_COLO829.vcf"

echo "== SNV / indel truth set (parklab SMaHT COLO829BLT50 benchmarking study) =="
# Called from short-read Illumina data, validated against long-read PacBio data, and
# filtered against COLO829BL as a negative control. NOTE: this file uses "chr"-prefixed
# contigs (chr1, chr9, ...), unlike the BAMs/reference/SV-truth-set above. The pipeline's
# prepare_truth_sets process strips this prefix before any comparison, see modules/bin/.
curl -L --fail -o "$REF_DIR/truth_sets/Truthset_SNV_COLO829.txt" \
  "https://raw.githubusercontent.com/parklab/SMaHT_SNV_COLO829BLT50_HAPMAP/main/Resource/Truthset_SNV_COLO829.txt"
curl -L --fail -o "$REF_DIR/truth_sets/Truthset_Indel_COLO829.txt" \
  "https://raw.githubusercontent.com/parklab/SMaHT_SNV_COLO829BLT50_HAPMAP/main/Resource/Truthset_Indel_COLO829.txt"

echo "== Copy-number comparison reference (BICseq2 segmentation, NOT a validated truth set) =="
# There is no independently validated, multi-platform CNV truth set for COLO829 the way
# there is for SNVs/indels/SVs. This BICseq2 segmentation comes from the same truth-set
# study but is itself just one algorithm's output, so Module 2 treats it as a cross-tool
# concordance comparison, not a precision/recall benchmark. Also "chr"-prefixed contigs.
curl -L --fail -o "$REF_DIR/cnv_reference/colo829_somaticSV_copynumber.zip" \
  "https://zenodo.org/records/7515830/files/COLO829_somaticSV_copynumber.zip?download=1"
python3 -m zipfile -e "$REF_DIR/cnv_reference/colo829_somaticSV_copynumber.zip" "$REF_DIR/cnv_reference/"

echo ""
echo "Done. Reference data staged under $REF_DIR"
find "$REF_DIR" -type f | sort
