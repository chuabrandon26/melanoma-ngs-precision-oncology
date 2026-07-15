/*
 * Adapter trimming and QC with fastp, chosen over Trim Galore for this project because
 * it does trimming, adapter detection, and QC reporting in a single C++ binary with no
 * Perl dependency chain, and runs noticeably faster on the same input, a real
 * consideration on a resource-constrained laptop. Trim Galore (a Cutadapt+FastQC
 * wrapper) is an equally standard, widely used choice for this step; fastp was picked
 * for speed and a smaller dependency footprint, not because Trim Galore would be wrong.
 */
process FASTP_TRIM {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/rnaseq/qc", mode: 'copy', pattern: "*.{json,html}"

    conda "bioconda::fastp=1.3.6"
    container "quay.io/biocontainers/fastp:1.3.6--h43da1c4_0"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"), path("${meta.id}_R2.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"), emit: json
    path "${meta.id}.fastp.html", emit: html
    path "versions.yml", emit: versions

    script:
    """
    fastp \\
        -i ${reads1} -I ${reads2} \\
        -o ${meta.id}_R1.trimmed.fastq.gz -O ${meta.id}_R2.trimmed.fastq.gz \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe
    # --detect_adapter_for_pe: fastp auto-detects the adapter sequence per read pair by
    # overlap analysis instead of needing a hardcoded adapter FASTA, fastp's documented
    # recommendation for paired-end data, and appropriate here since the exact library
    # prep kit used for the reused GSE78220 samples was not independently confirmed.
    # All other trimming thresholds (quality cutoff, minimum length, etc.) are left at
    # fastp's own defaults.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
    END_VERSIONS
    """
}
