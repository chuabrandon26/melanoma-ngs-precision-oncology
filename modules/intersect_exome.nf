/*
 * Restricts a BAM to reads overlapping the exome BED: the step that turns Module 2's
 * whole-chr9+chr20 data into Module 1's exome-intersected WES subset. Runs on the
 * already chromosome-sliced, duplicate-marked BAM, so this never triggers a separate
 * download, it is a fast, local, free filter on data the WGS module already needed.
 */
process INTERSECT_EXOME_BED {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.processed_dir}/wes", mode: 'copy'

    conda "bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(bam), path(bai)
    path exome_bed

    output:
    tuple val(meta), path("${meta.id}.exome.bam"), path("${meta.id}.exome.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    """
    samtools view -h -b -L ${exome_bed} --threads ${task.cpus} -o ${meta.id}.exome.bam ${bam}
    samtools index ${meta.id}.exome.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
