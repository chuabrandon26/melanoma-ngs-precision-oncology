/*
 * Downloads specific public RNA-seq FASTQ files directly over HTTPS from ENA.
 *
 * Unlike the COLO829 BAMs (modules/slice_remote_bam.nf), RNA-seq reads cannot be
 * "sliced" to a small region over the network the same way: expression can come from
 * any gene anywhere in the transcriptome, so there is no small genomic window to seek
 * to, which means downloading a full FASTQ file is unavoidable if raw-read alignment is
 * wanted at all. Every file pulled here is several GB, well over the 2GB threshold this
 * project's hard constraints require explicit confirmation for before downloading. The
 * exact two samples this project uses (the smallest available in GSE78220) were sized
 * against ENA's filereport API and confirmed with the user before being wired in, which
 * is also why this whole path sits behind `--run_fastq_demo true` rather than running
 * by default, see pipelines/rnaseq/main.nf.
 */
process DOWNLOAD_FASTQ {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.raw_dir}/rnaseq_fastq", mode: 'copy'

    input:
    tuple val(meta), val(fastq1_url), val(fastq2_url)

    output:
    tuple val(meta), path("${meta.id}_1.fastq.gz"), path("${meta.id}_2.fastq.gz"), emit: reads

    script:
    """
    curl -L --fail -o ${meta.id}_1.fastq.gz ${fastq1_url}
    curl -L --fail -o ${meta.id}_2.fastq.gz ${fastq2_url}
    """
}
