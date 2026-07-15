/*
 * Somatic copy-number calling with CNVkit, in whole-genome mode: `--method wgs` bins
 * evenly across the accessible genome instead of building the exon-target/antitarget bin
 * design CNVkit normally uses for capture data, which fits Module 2's whole-chr9+chr20
 * scope (this process is not used by the WES pipeline).
 */
process CALL_CNV_CNVKIT {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.processed_dir}/wgs", mode: 'copy'

    conda "bioconda::cnvkit=0.9.13"
    container "quay.io/biocontainers/cnvkit:0.9.13--pyhdfd78af_0"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_fai
    path access_bed

    output:
    tuple val(meta), path("${meta.id}.cns"), emit: segments
    tuple val(meta), path("${meta.id}.cnr"), emit: ratios
    tuple val(meta), path("${meta.id}.call.cns"), emit: calls
    path "versions.yml", emit: versions

    script:
    """
    cnvkit.py batch ${tumor_bam} \\
        --normal ${normal_bam} \\
        --method wgs \\
        --fasta ${fasta} \\
        --access ${access_bed} \\
        --output-dir . \\
        -p ${task.cpus}
    # --access restricts CNVkit's binning to chr9+chr20 (see prepare_reference.nf's
    # genome.access.bed); without it CNVkit would try to bin the rest of the genome too,
    # which was never sliced from the source BAMs and so isn't present to bin against.

    TUMOR_PREFIX=\$(basename ${tumor_bam} .bam)
    cp \${TUMOR_PREFIX}.cns ${meta.id}.cns
    cp \${TUMOR_PREFIX}.cnr ${meta.id}.cnr

    cnvkit.py call ${meta.id}.cns -o ${meta.id}.call.cns
    # `cnvkit call` converts continuous log2 copy-ratio segments into discrete calls using
    # CNVkit's own documented default thresholds for a diploid sample: log2 < -1.1 is
    # called a homozygous loss, < -0.4 a single-copy loss, > 0.3 a gain, > 0.7 an
    # amplification. These are the tool's stock defaults, tuned by the CNVkit authors for
    # a reasonably pure tumor sample; COLO829's actual tumor purity has not been
    # independently re-verified for this specific chr9+chr20 slice, so treat the discrete
    # calls as a standard-default interpretation of the log2 ratios, not a
    # purity-corrected one.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version)
    END_VERSIONS
    """
}
