#!/usr/bin/env nextflow
/*
 * Module 3: bulk RNA-seq, melanoma
 *
 * Two independent data paths feed one shared downstream analysis (differential
 * expression, GSEA/GO, and the TMB-vs-expression correlation), all of which happens in
 * notebooks/03_rnaseq_de_and_tmb_correlation.ipynb rather than in this pipeline:
 *
 *   1. TCGA-SKCM open-tier gene counts (runs by default): pre-quantified, open-access
 *      STAR gene counts for up to 2 x params.tcga_max_per_group patients from the GDC
 *      API, no login needed. This is the statistically powered path: a real patient
 *      cohort, a real DE comparison (Primary Tumor vs Metastatic), and real TMB data
 *      available for a genuine multi-patient TMB-vs-expression correlation.
 *   2. A small GSE78220 FASTQ subset (opt-in via --run_fastq_demo true): trimming with
 *      fastp and quantification with Salmon, to demonstrate the raw-FASTQ pipeline
 *      mechanics the project brief originally called for (with Salmon standing in for
 *      STAR, see modules/quantify.nf for why). Two samples only, the smallest available
 *      in that series, since RNA-seq FASTQ cannot be region-sliced the way the COLO829
 *      BAMs can, every file here is several GB. This path is a mechanics demonstration,
 *      not a powered study; the TCGA path above carries the actual statistical analysis.
 *
 * Run from the repository root so nextflow.config is picked up automatically:
 *   nextflow run pipelines/rnaseq/main.nf -profile conda
 *   nextflow run pipelines/rnaseq/main.nf -profile conda --run_fastq_demo true
 *
 * Output: data/raw/tcga_skcm/ (TCGA count files + sample sheet) and results/rnaseq/
 * (trimmed-read QC reports and, if the FASTQ demo is enabled, Salmon quant.sf per
 * sample), all read directly by notebooks/03_rnaseq_de_and_tmb_correlation.ipynb, which
 * does all of the differential expression, GSEA/GO, and TMB-correlation work itself.
 */

include { FETCH_TCGA_SKCM_COUNTS }     from '../../modules/fetch_tcga_counts.nf'
include { DOWNLOAD_FASTQ }             from '../../modules/download_fastq.nf'
include { FASTP_TRIM }                 from '../../modules/trim_reads.nf'
include { SALMON_INDEX; SALMON_QUANT } from '../../modules/quantify.nf'

workflow {

    log.info """
    ================================================================================
    Module 3: bulk RNA-seq (melanoma)
    ================================================================================
    TCGA-SKCM open counts   : always runs, up to ${params.tcga_max_per_group} samples/group
    GSE78220 FASTQ demo     : ${params.run_fastq_demo ? "enabled (${params.fastq_demo_samples.size()} samples)" : 'disabled, use --run_fastq_demo true to enable'}
    output dir              : ${params.outdir}/rnaseq
    ================================================================================
    """.stripIndent()

    // This file's own moduleDir is pipelines/rnaseq/, not modules/, so ../../modules/bin
    // navigates from the repo root back down to the shared helper script; passed as a
    // real input so this process's -resume cache correctly invalidates if the script is
    // ever edited (see modules/prepare_truth_sets.nf for why this matters).
    fetch_script = file("${moduleDir}/../../modules/bin/fetch_tcga_skcm_counts.py", checkIfExists: true)
    FETCH_TCGA_SKCM_COUNTS(params.tcga_max_per_group, fetch_script)

    if (params.run_fastq_demo) {
        fastq_sources = Channel.fromList(params.fastq_demo_samples)
            .map { s -> [[id: s.id], s.r1, s.r2] }

        DOWNLOAD_FASTQ(fastq_sources)
        FASTP_TRIM(DOWNLOAD_FASTQ.out.reads)

        pc_transcripts     = file(params.pc_transcripts_fa,     checkIfExists: true)
        lncrna_transcripts = file(params.lncrna_transcripts_fa, checkIfExists: true)
        SALMON_INDEX(pc_transcripts, lncrna_transcripts)
        SALMON_QUANT(FASTP_TRIM.out.reads, SALMON_INDEX.out.index)
    }
}
