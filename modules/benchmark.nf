/*
 * Benchmarking against the COLO829 truth sets. Each process here does the mechanical
 * comparison (running the standard tool for that variant type and writing structured
 * output) and stops there; the precision/recall/F1 arithmetic, plotting, and narrative
 * interpretation happen in the notebooks, which read these outputs directly. That split
 * mirrors how the whole project divides labour: the pipeline is the source of truth for
 * calling and comparison mechanics, the notebook is the source of truth for what the
 * numbers mean.
 */

/*
 * SNV/indel benchmarking with hap.py, the GA4GH/precisionFDA-standard tool for comparing
 * a query VCF against a truth VCF over a defined comparison region, with genotype-aware
 * matching (it accounts for representation differences, e.g. an indel written
 * differently but describing the same edit) rather than requiring identical REF/ALT
 * strings, and separate SNV vs indel accounting.
 *
 * Caveat carried through to the notebook: unlike GIAB reference materials, the COLO829
 * SNV/indel truth set does not ship its own independently-defined "confident regions"
 * BED. This process uses the full comparison region (whole chr9+chr20 for Module 2, the
 * exome-intersected subset for Module 1) as the target space instead, which can inflate
 * apparent false positives/negatives in regions the truth-set authors might otherwise
 * have excluded as unreliable. Treat these precision/recall numbers as good-faith
 * estimates, not GIAB-grade rigor.
 */
process BENCHMARK_SNV_INDEL {
    tag "${meta.id}"
    label 'process_low'

    // Closure required, see modules/call_snv.nf for why a plain interpolated string
    // referencing `meta` does not work in a publishDir directive.
    publishDir(path: { "${params.outdir}/${meta.module}/benchmark" }, mode: 'copy')

    // rtg-tools is a separate package from hap.py itself, needed because --engine=vcfeval
    // below shells out to RTG Tools' own "rtg" binary to actually do the comparison;
    // hap.py's own package does not bundle it. Missing at first (a real failure hit while
    // running this: "Error running rtg tools ... rtg: not found", exit 127), now pinned
    // explicitly. The stock hap.py biocontainer likewise does not include rtg-tools, so
    // -profile docker/singularity for this specific process would need a combined image;
    // -profile conda resolves both packages side by side cleanly.
    conda "bioconda::hap.py=0.3.15 bioconda::rtg-tools=3.13"
    container "quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0"
    // bcftools is deliberately NOT added here: hap.py 0.3.15 is a legacy Python 2.7
    // tool whose pysam dependency needs zlib <1.3.0a0, while a modern bcftools/htslib
    // needs zlib >=1.3.1, directly incompatible ranges that mamba cannot solve together
    // in one environment (a real "Could not solve for environment specs" failure hit
    // while running this). The tumor-sample extraction that used to live here has been
    // moved to its own process (EXTRACT_TUMOR_VCF, modules/call_snv.nf) with its own
    // clean bcftools-only environment instead.

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    tuple path(truth_vcf), path(truth_tbi)
    path comparison_bed
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}.happy.summary.csv"), emit: summary
    tuple val(meta), path("${meta.id}.happy.extended.csv"), emit: extended
    path "versions.yml", emit: versions

    script:
    """
    hap.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        -r ${fasta} \\
        -T ${comparison_bed} \\
        -o ${meta.id}.happy \\
        --engine=vcfeval \\
        --threads ${task.cpus}
    # --engine=vcfeval uses RealTimeGenomics's vcfeval for the actual comparison, hap.py's
    # own documented recommendation over its older internal comparison engine.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        happy: 0.3.15
    END_VERSIONS
    """
}

/*
 * Structural variant benchmarking with Truvari, the current community-standard tool for
 * SV comparison. Truvari matches calls using breakpoint distance plus size/sequence
 * similarity rather than exact position matching, which naive VCF comparison (e.g.
 * bcftools isec) cannot do correctly for SVs, since two callers essentially never report
 * byte-identical breakpoints for the same real event.
 */
process BENCHMARK_SV {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/wgs/benchmark", mode: 'copy'

    conda "bioconda::truvari=5.4.0"
    container "quay.io/biocontainers/truvari:5.4.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(query_vcf), path(query_tbi)
    tuple path(truth_vcf), path(truth_tbi)
    path include_bed
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}_truvari/summary.json"), emit: summary
    path "${meta.id}_truvari", emit: all_outputs
    path "versions.yml", emit: versions

    script:
    """
    truvari bench \\
        -b ${truth_vcf} \\
        -c ${query_vcf} \\
        -f ${fasta} \\
        --includebed ${include_bed} \\
        -o ${meta.id}_truvari \\
        --pick multi \\
        --refdist 500 \\
        --pctseq 0.7 \\
        --pctsize 0.7
    # All four thresholds are Truvari's own documented defaults, not tuned for COLO829:
    #   --refdist 500  max distance (bp) between call and truth breakpoints to even be
    #                  considered a possible match, before size/sequence checks apply
    #   --pctsize 0.7  minimum reciprocal size similarity, min(len)/max(len)
    #   --pctseq 0.7   minimum sequence similarity, for insertion sequences specifically
    #   --pick multi   let a truth call match multiple query calls and vice versa, then
    #                  keep only the single best pairing; Truvari's recommended mode for
    #                  most benchmarking use over the stricter 'single' setting
    # These are reasonable, widely-used starting points. Whether they are optimal for this
    # exact truth set (built from short- plus long-read consensus SV calls) has not been
    # independently verified here.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(truvari version | sed 's/Truvari v//')
    END_VERSIONS
    """
}

/*
 * Copy-number "benchmarking", in scare quotes deliberately: there is no independently
 * validated, multi-platform CNV truth set for COLO829 the way there is for SNVs, indels,
 * and SVs (see README). This process instead prepares two concordance comparisons,
 * neither of which is a precision/recall-against-gold-standard metric:
 *   1. Overlap between CNVkit's non-neutral segments and a BICseq2 segmentation computed
 *      independently in the original COLO829 SV truth-set study (a cross-tool
 *      comparison: both are "callers" here, neither is truth)
 *   2. Overlap between CNVkit's non-neutral segments and the DEL/DUP entries in the
 *      COLO829 SV truth set (which WAS multi-platform validated), on the reasoning that a
 *      genuine deletion or tandem duplication above CNVkit's detection floor should also
 *      show up as a copy-number change
 * This process only produces the raw overlap tables; the notebook computes concordance
 * rates and visualizes them, and is explicit that this is not a true benchmark.
 */
process BENCHMARK_CNV {
    tag "${meta.id}"
    label 'process_single'

    publishDir "${params.outdir}/wgs/benchmark", mode: 'copy'

    conda "bioconda::bedtools=2.31.1 conda-forge::python=3.11"
    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"
    // Same container/profile caveat as PREPARE_TRUTH_SETS: the stock bedtools
    // biocontainer does not bundle python3, needed here for cns_to_nonneutral_bed.py.
    // Conda resolves this cleanly; -profile docker needs a combined image for this step.

    input:
    tuple val(meta), path(cnvkit_calls)
    path bicseq_cna_bed
    path sv_truth_delDup_bed
    // Declared as a real path input (not just referenced by a moduleDir string inside the
    // script block), so a future edit to this script correctly invalidates -resume's
    // cache instead of silently reusing stale output; see modules/prepare_truth_sets.nf
    // for the real instance of this that happened during development.
    path cns_script

    output:
    tuple val(meta), path("${meta.id}.cnvkit_vs_bicseq.overlap.tsv"), emit: bicseq_overlap
    tuple val(meta), path("${meta.id}.cnvkit_vs_svtruth_delDup.overlap.tsv"), emit: svtruth_overlap
    path "versions.yml", emit: versions

    script:
    // COLO829 is a diploid-derived cell line (no whole-genome duplication reported for
    // it), so ploidy 2 is used as the "neutral" baseline, matching CNVkit's own default
    // assumption when no --ploidy override is given to `cnvkit.py call`.
    """
    python3 ${cns_script} ${cnvkit_calls} cnvkit_nonneutral.bed 2
    sort -k1,1 -k2,2n cnvkit_nonneutral.bed > ${meta.id}.cnvkit_nonneutral.sorted.bed

    sort -k1,1 -k2,2n ${bicseq_cna_bed} > bicseq.sorted.bed
    sort -k1,1 -k2,2n ${sv_truth_delDup_bed} > svtruth_delDup.sorted.bed

    bedtools intersect -a ${meta.id}.cnvkit_nonneutral.sorted.bed -b bicseq.sorted.bed -wao \\
        > ${meta.id}.cnvkit_vs_bicseq.overlap.tsv

    bedtools intersect -a ${meta.id}.cnvkit_nonneutral.sorted.bed -b svtruth_delDup.sorted.bed -wao \\
        > ${meta.id}.cnvkit_vs_svtruth_delDup.overlap.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}
