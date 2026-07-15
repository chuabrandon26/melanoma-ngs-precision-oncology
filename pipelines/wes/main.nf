#!/usr/bin/env nextflow
/*
 * Module 1: WES, exome-intersected COLO829 (chr9 + chr20)
 *
 * Rebuilds the same conceptual pipeline as the author's existing Snakemake WES project
 * (BWA-MEM alignment, duplicate marking, GATK Mutect2 somatic calling) in Nextflow DSL2,
 * applied to the chr9+chr20 slice of COLO829T/COLO829BL, restricted to an approximate
 * exome capture region built from GENCODE exons.
 *
 * Input: nothing local to stage beforehand except the small reference/truth-set files
 * under data/reference/ (see scripts/fetch_reference_data.sh). The COLO829 BAMs
 * themselves are pulled directly from ENA over HTTPS range requests by SLICE_REMOTE_BAM,
 * chr9+chr20 only, never downloaded in full.
 *
 * Run from the repository root so nextflow.config is picked up automatically:
 *   nextflow run pipelines/wes/main.nf -profile conda
 *
 * Output: results/wes/ (filtered Mutect2 VCF) and results/wes/benchmark/ (hap.py
 * precision/recall/F1 tables), both read directly by notebooks/01_wes_benchmarking.ipynb.
 */

include { SLICE_REMOTE_BAM }                    from '../../modules/slice_remote_bam.nf'
include { MARK_DUPLICATES }                     from '../../modules/mark_duplicates.nf'
include { PREPARE_REFERENCE; BUILD_EXOME_BED; PREPARE_GERMLINE_RESOURCE } from '../../modules/prepare_reference.nf'
include { INTERSECT_EXOME_BED }                 from '../../modules/intersect_exome.nf'
include { CALL_SNV_MUTECT2; FILTER_MUTECT_CALLS; EXTRACT_TUMOR_VCF } from '../../modules/call_snv.nf'
include { PREPARE_TRUTH_SETS }                  from '../../modules/prepare_truth_sets.nf'
include { BENCHMARK_SNV_INDEL }                 from '../../modules/benchmark.nf'

workflow {

    log.info """
    ================================================================================
    Module 1: WES (exome-intersected COLO829, chr${params.chromosomes.replace(',', ' + chr')})
    ================================================================================
    tumor sample    : ${params.tumor_sample_name}
    normal sample   : ${params.normal_sample_name}
    downsampling    : ${params.downsample_enabled ? "enabled (target ${params.target_tumor_depth}x/${params.target_normal_depth}x vs native ${params.tumor_native_depth}x/${params.normal_native_depth}x)" : 'disabled, full native depth'}
    output dir      : ${params.outdir}/wes
    ================================================================================
    """.stripIndent()

    // ---- reference genome (chr9+chr20 only) and an approximate exome BED ----
    ref_fasta_chr9  = file(params.ref_fasta_chr9,  checkIfExists: true)
    ref_fasta_chr20 = file(params.ref_fasta_chr20, checkIfExists: true)
    gtf             = file(params.gtf,             checkIfExists: true)

    PREPARE_REFERENCE(ref_fasta_chr9, ref_fasta_chr20)
    BUILD_EXOME_BED(gtf, params.chromosomes)
    PREPARE_GERMLINE_RESOURCE(params.germline_resource_chr9_url, params.germline_resource_chr20_url)

    // ---- COLO829 truth sets, prepared once and reused by the benchmarking step ----
    snv_truth_txt   = file(params.truth_snv_txt,   checkIfExists: true)
    indel_truth_txt = file(params.truth_indel_txt, checkIfExists: true)
    sv_truth_vcf    = file(params.truth_sv_vcf,    checkIfExists: true)
    // This file's own moduleDir is pipelines/wes/, not modules/, so ../../modules/bin
    // navigates from the repo root back down to the shared helper scripts; passed as
    // real inputs so PREPARE_TRUTH_SETS's -resume cache correctly invalidates if either
    // script is ever edited (see modules/prepare_truth_sets.nf for why this matters).
    convert_script = file("${moduleDir}/../../modules/bin/convert_colo829_truth_to_vcf.py", checkIfExists: true)
    extract_script = file("${moduleDir}/../../modules/bin/extract_sv_delDup_bed.py", checkIfExists: true)
    PREPARE_TRUTH_SETS(snv_truth_txt, indel_truth_txt, sv_truth_vcf, params.chromosomes, convert_script, extract_script)

    // ---- pull only chr9+chr20 from the two remote, indexed COLO829 BAMs ----
    // Coverage downsampling fractions: see nextflow.config for the full rationale
    // (real disk pressure during development). min(1.0, ...) means a target depth at or
    // above native depth is simply a no-op (full native depth), not an upsample attempt.
    def tumor_fraction  = params.downsample_enabled ? Math.min(1.0d, params.target_tumor_depth  / params.tumor_native_depth)  : 1.0d
    def normal_fraction = params.downsample_enabled ? Math.min(1.0d, params.target_normal_depth / params.normal_native_depth) : 1.0d

    bam_sources = Channel.of(
        [[id: params.tumor_sample_name,  status: 'tumor',  module: 'wes', subsample_fraction: tumor_fraction],  params.tumor_bam_url],
        [[id: params.normal_sample_name, status: 'normal', module: 'wes', subsample_fraction: normal_fraction], params.normal_bam_url]
    )
    SLICE_REMOTE_BAM(bam_sources, params.chromosomes)
    MARK_DUPLICATES(SLICE_REMOTE_BAM.out.bam)
    INTERSECT_EXOME_BED(MARK_DUPLICATES.out.bam, BUILD_EXOME_BED.out.bed)

    // ---- pair the (now exome-intersected) tumor and normal BAMs for somatic calling ----
    // Exactly one tumor and one normal element ever exist in this project (a fixed
    // COLO829T-vs-COLO829BL design, not a general multi-sample cohort pipeline), so a
    // plain .combine() with no join key gives exactly one correctly paired tuple.
    tumor_ch  = INTERSECT_EXOME_BED.out.bam.filter { meta, bam, bai -> meta.status == 'tumor' }
    normal_ch = INTERSECT_EXOME_BED.out.bam.filter { meta, bam, bai -> meta.status == 'normal' }
    paired_ch = tumor_ch.combine(normal_ch).map { tmeta, tbam, tbai, nmeta, nbam, nbai ->
        [[id: "${params.tumor_sample_name}_vs_${params.normal_sample_name}", module: 'wes'], tbam, tbai, nbam, nbai]
    }

    CALL_SNV_MUTECT2(
        paired_ch,
        PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict,
        BUILD_EXOME_BED.out.bed,
        PREPARE_GERMLINE_RESOURCE.out.resource
    )
    FILTER_MUTECT_CALLS(
        CALL_SNV_MUTECT2.out.vcf,
        PREPARE_REFERENCE.out.fasta, PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict
    )
    // hap.py's comparison engine needs a single-sample VCF; Mutect2's output has both
    // tumor and normal genotypes together, so the tumor sample is split out first.
    EXTRACT_TUMOR_VCF(FILTER_MUTECT_CALLS.out.vcf)

    // ---- benchmark against the COLO829 SNV/indel truth set ----
    // The exome BED is passed as hap.py's target-regions argument (-T), which restricts
    // the comparison to that region for both the truth and query VCFs at once, so a
    // variant caller scoped to exome regions is never penalized for "missing" a truth
    // variant that was never in its capture space to begin with.
    BENCHMARK_SNV_INDEL(
        EXTRACT_TUMOR_VCF.out.vcf,
        PREPARE_TRUTH_SETS.out.snv_indel_combined,
        BUILD_EXOME_BED.out.bed,
        PREPARE_REFERENCE.out.fasta,
        PREPARE_REFERENCE.out.fai
    )
}
