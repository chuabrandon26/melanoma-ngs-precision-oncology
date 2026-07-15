/*
 * Somatic structural variant calling with Manta, tumor vs matched normal.
 *
 * Manta, not GRIDSS: Manta installs from a single bioconda package with no JVM or
 * external-aligner dependency chain, which matters on a resource-constrained laptop.
 * GRIDSS is arguably more sensitive for complex/multi-breakpoint rearrangements, but
 * needs its own BWA index build and a heavier JVM-based assembly stage on top of that.
 * The COLO829 SV truth set used for benchmarking here was itself built from a consensus
 * across several callers and platforms, so Manta's calls are being compared to a
 * multi-caller consensus, not to GRIDSS specifically, which keeps this a fair comparison
 * even though only one caller is run. This is a resource-driven tool choice, not a claim
 * that Manta is unconditionally the better SV caller.
 */
process CALL_SV_MANTA {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.processed_dir}/wgs", mode: 'copy'

    conda "bioconda::manta=1.6.0"
    container "quay.io/biocontainers/manta:1.6.0--py27h9948957_6"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_fai
    path call_regions_bed
    path call_regions_bed_tbi

    output:
    tuple val(meta), path("${meta.id}.manta.somaticSV.vcf.gz"), path("${meta.id}.manta.somaticSV.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    configManta.py \\
        --tumorBam ${tumor_bam} \\
        --normalBam ${normal_bam} \\
        --referenceFasta ${fasta} \\
        --callRegions ${call_regions_bed} \\
        --runDir manta_work

    ./manta_work/runWorkflow.py -m local -j ${task.cpus}
    # runWorkflow.py is executable and carries its own shebang pointing at Manta's bundled
    # python2, which is how Manta's own documentation says to invoke it, rather than
    # relying on whatever "python2" happens to resolve to on PATH.

    cp manta_work/results/variants/somaticSV.vcf.gz ${meta.id}.manta.somaticSV.vcf.gz
    cp manta_work/results/variants/somaticSV.vcf.gz.tbi ${meta.id}.manta.somaticSV.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version)
    END_VERSIONS
    """
}
