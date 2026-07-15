#!/usr/bin/env nextflow
/*
 * Module 2: WGS-style somatic analysis, whole chr9 + chr20, COLO829
 *
 * Demonstrates the analysis WES cannot: structural variant and copy-number calling,
 * plus whole-chromosome (not exome-restricted) SNV/indel calling, scoped to chr9+chr20
 * for the resource reasons documented in the README and in Step 0 of the project's
 * original data-sizing pass.
 *
 * Input: reuses the same remote-sliced, duplicate-marked chr9+chr20 BAMs as Module 1,
 * but WITHOUT the exome intersection this time. If pipelines/wes/main.nf has already
 * been run, that slice/dedup work is not repeated here for free (each pipeline currently
 * re-runs SLICE_REMOTE_BAM/MARK_DUPLICATES independently), which matters for this
 * project's tight disk budget: see the README for the recommended sequencing (run one
 * pipeline, confirm its outputs, delete data/raw/*.bam, then run the next).
 *
 * Run from the repository root so nextflow.config is picked up automatically:
 *   nextflow run pipelines/wgs/main.nf -profile conda
 *
 * Output: results/wgs/ (filtered Mutect2 VCF, Manta SV VCF, CNVkit segments/calls) and
 * results/wgs/benchmark/ (hap.py, Truvari, and CNV concordance tables), all read
 * directly by notebooks/02_wgs_sv_cnv_benchmarking.ipynb.
 */

include { SLICE_REMOTE_BAM }        from '../../modules/slice_remote_bam.nf'
include { MARK_DUPLICATES }         from '../../modules/mark_duplicates.nf'
include { PREPARE_REFERENCE; PREPARE_GERMLINE_RESOURCE } from '../../modules/prepare_reference.nf'
include { CALL_SNV_MUTECT2; FILTER_MUTECT_CALLS; EXTRACT_TUMOR_VCF } from '../../modules/call_snv.nf'
include { CALL_SV_MANTA }           from '../../modules/call_sv.nf'
include { CALL_CNV_CNVKIT }         from '../../modules/call_cnv.nf'
include { PREPARE_TRUTH_SETS; PREPARE_CNV_REFERENCE } from '../../modules/prepare_truth_sets.nf'
include { BENCHMARK_SNV_INDEL; BENCHMARK_SV; BENCHMARK_CNV } from '../../modules/benchmark.nf'

workflow {

    log.info """
    ================================================================================
    Module 2: WGS-style analysis (whole COLO829 chr${params.chromosomes.replace(',', ' + chr')})
    ================================================================================
    tumor sample    : ${params.tumor_sample_name}
    normal sample   : ${params.normal_sample_name}
    downsampling    : ${params.downsample_enabled ? "enabled (target ${params.target_tumor_depth}x/${params.target_normal_depth}x vs native ${params.tumor_native_depth}x/${params.normal_native_depth}x)" : 'disabled, full native depth'}
    output dir      : ${params.outdir}/wgs
    ================================================================================
    """.stripIndent()

    // ---- reference genome (chr9+chr20 only), plus Manta/CNVkit region files ----
    ref_fasta_chr9  = file(params.ref_fasta_chr9,  checkIfExists: true)
    ref_fasta_chr20 = file(params.ref_fasta_chr20, checkIfExists: true)
    PREPARE_REFERENCE(ref_fasta_chr9, ref_fasta_chr20)
    PREPARE_GERMLINE_RESOURCE(params.germline_resource_chr9_url, params.germline_resource_chr20_url)

    // ---- COLO829 truth sets and the BICseq2 comparison reference ----
    snv_truth_txt   = file(params.truth_snv_txt,      checkIfExists: true)
    indel_truth_txt = file(params.truth_indel_txt,    checkIfExists: true)
    sv_truth_vcf    = file(params.truth_sv_vcf,       checkIfExists: true)
    bicseq_cna      = file(params.cnv_reference_cna,  checkIfExists: true)
    // This file's own moduleDir is pipelines/wgs/, not modules/, so ../../modules/bin
    // navigates from the repo root back down to the shared helper scripts; passed as
    // real inputs so these processes' -resume cache correctly invalidates if a script is
    // ever edited (see modules/prepare_truth_sets.nf for why this matters).
    convert_script = file("${moduleDir}/../../modules/bin/convert_colo829_truth_to_vcf.py", checkIfExists: true)
    extract_script = file("${moduleDir}/../../modules/bin/extract_sv_delDup_bed.py", checkIfExists: true)
    cns_script      = file("${moduleDir}/../../modules/bin/cns_to_nonneutral_bed.py", checkIfExists: true)
    PREPARE_TRUTH_SETS(snv_truth_txt, indel_truth_txt, sv_truth_vcf, params.chromosomes, convert_script, extract_script)
    PREPARE_CNV_REFERENCE(bicseq_cna, params.chromosomes)

    // ---- pull only chr9+chr20 from the two remote, indexed COLO829 BAMs ----
    // Coverage downsampling fractions: see nextflow.config for the full rationale
    // (real disk pressure during development). min(1.0, ...) means a target depth at or
    // above native depth is simply a no-op (full native depth), not an upsample attempt.
    def tumor_fraction  = params.downsample_enabled ? Math.min(1.0d, params.target_tumor_depth  / params.tumor_native_depth)  : 1.0d
    def normal_fraction = params.downsample_enabled ? Math.min(1.0d, params.target_normal_depth / params.normal_native_depth) : 1.0d

    bam_sources = Channel.of(
        [[id: params.tumor_sample_name,  status: 'tumor',  module: 'wgs', subsample_fraction: tumor_fraction],  params.tumor_bam_url],
        [[id: params.normal_sample_name, status: 'normal', module: 'wgs', subsample_fraction: normal_fraction], params.normal_bam_url]
    )
    SLICE_REMOTE_BAM(bam_sources, params.chromosomes)
    MARK_DUPLICATES(SLICE_REMOTE_BAM.out.bam)
    // No exome intersection here, unlike Module 1: this is the whole-chromosome path.

    // ---- pair tumor and normal (see pipelines/wes/main.nf for why .combine() with no
    // join key is safe here: exactly one tumor and one normal element ever exist) ----
    tumor_ch  = MARK_DUPLICATES.out.bam.filter { meta, bam, bai -> meta.status == 'tumor' }
    normal_ch = MARK_DUPLICATES.out.bam.filter { meta, bam, bai -> meta.status == 'normal' }
    paired_ch = tumor_ch.combine(normal_ch).map { tmeta, tbam, tbai, nmeta, nbam, nbai ->
        [[id: "${params.tumor_sample_name}_vs_${params.normal_sample_name}", module: 'wgs'], tbam, tbai, nbam, nbai]
    }

    // ---- somatic SNV/indel calling, whole chr9+chr20 (genome.access.bed, not the exome BED) ----
    CALL_SNV_MUTECT2(
        paired_ch,
        PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict,
        PREPARE_REFERENCE.out.access_bed,
        PREPARE_GERMLINE_RESOURCE.out.resource
    )
    FILTER_MUTECT_CALLS(
        CALL_SNV_MUTECT2.out.vcf,
        PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict
    )
    // hap.py's comparison engine needs a single-sample VCF; Mutect2's output has both
    // tumor and normal genotypes together, so the tumor sample is split out first.
    EXTRACT_TUMOR_VCF(FILTER_MUTECT_CALLS.out.vcf)
    BENCHMARK_SNV_INDEL(
        EXTRACT_TUMOR_VCF.out.vcf,
        PREPARE_TRUTH_SETS.out.snv_indel_combined,
        PREPARE_REFERENCE.out.access_bed,
        PREPARE_REFERENCE.out.fasta,
        PREPARE_REFERENCE.out.fai
    )

    // ---- somatic structural variants ----
    CALL_SV_MANTA(
        paired_ch,
        PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
        PREPARE_REFERENCE.out.call_regions_bed, PREPARE_REFERENCE.out.call_regions_bed_tbi
    )
    BENCHMARK_SV(
        CALL_SV_MANTA.out.vcf,
        PREPARE_TRUTH_SETS.out.sv,
        PREPARE_REFERENCE.out.access_bed,
        PREPARE_REFERENCE.out.fasta,
        PREPARE_REFERENCE.out.fai
    )

    // ---- somatic copy number ----
    CALL_CNV_CNVKIT(
        paired_ch,
        PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai,
        PREPARE_REFERENCE.out.access_bed
    )
    BENCHMARK_CNV(
        CALL_CNV_CNVKIT.out.calls,
        PREPARE_CNV_REFERENCE.out.bed,
        PREPARE_TRUTH_SETS.out.sv_deldup_bed,
        cns_script
    )
}
