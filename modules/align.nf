/*
 * BWA-MEM alignment. Included for completeness, since alignment is part of the same
 * conceptual pipeline this project was rebuilt from, but NOT part of this project's
 * default workflow: both pipelines start from the already-aligned, chromosome-sliced
 * COLO829 BAMs (see slice_remote_bam.nf) instead, since re-aligning from scratch would
 * mean downloading the ~125GB/~50GB full FASTQ this project deliberately avoids (hard
 * constraint: no full-genome COLO829 FASTQ or BAM downloads).
 *
 * This module exists, and is kept in working order, so the alignment step is genuinely
 * demonstrated in code rather than only described, and so this pipeline could be pointed
 * at a different, smaller FASTQ dataset later without writing this from scratch. It is
 * not wired into pipelines/wes/main.nf or pipelines/wgs/main.nf by default.
 */
process BWA_MEM_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::bwa=0.7.17 bioconda::samtools=1.19.2"
    // No container pinned here deliberately: this module isn't part of the default
    // workflow, so an exact multi-tool "mulled" container tag has not been verified
    // against the biocontainers registry the way every other container in this project
    // was (see README). Use the conda profile for this specific process, or verify a
    // bwa+samtools container tag on quay.io/biocontainers before using -profile docker
    // with this module.

    input:
    tuple val(meta), path(reads1), path(reads2)
    path fasta
    path bwa_index_files  // output of `bwa index`, staged alongside the fasta

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    // Read group is required by GATK/Mutect2 downstream (SM must match the sample name
    // used later in the Mutect2 -tumor/-normal arguments).
    def read_group = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA"
    """
    bwa mem -t ${task.cpus} -R "${read_group}" ${fasta} ${reads1} ${reads2} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam -
    samtools index ${meta.id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -i Version | sed 's/.*Version: //')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

process BWA_INDEX {
    tag "${fasta}"
    label 'process_medium'

    conda "bioconda::bwa=0.7.17"

    input:
    path fasta

    output:
    path "${fasta}.*", emit: index
    path "versions.yml", emit: versions

    script:
    """
    bwa index ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -i Version | sed 's/.*Version: //')
    END_VERSIONS
    """
}
