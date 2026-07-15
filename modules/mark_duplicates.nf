/*
 * Marks duplicate reads with GATK4 MarkDuplicates.
 *
 * The COLO829T/COLO829BL BAMs from ENA are already named "*_dedup.realigned.bam": the
 * original submitters already ran duplicate marking (and GATK3-era indel realignment, a
 * step GATK4's Best Practices no longer call for since HaplotypeCaller/Mutect2 do local
 * reassembly instead). Running MarkDuplicates again here should be close to a no-op on
 * this data (confirming existing duplicate flags rather than finding large numbers of new
 * duplicates), but the step is kept for two reasons: it matches the same conceptual
 * pipeline as the WES project this was rebuilt from, and it means this exact process
 * would do real work if ever pointed at a fresh, non-deduplicated BAM instead.
 */
process MARK_DUPLICATES {
    tag "${meta.id}"
    label 'process_medium'
    // GATK's JVM overhead plus sorting a chr9+chr20-sized BAM; process_medium's 8GB/1h
    // ceiling is comfortably above what a ~13GB (tumor) or ~5GB (normal) sliced BAM needs,
    // since MarkDuplicates streams through coordinate-sorted input rather than loading it
    // all into memory at once.

    publishDir "${params.processed_dir}/dedup", mode: 'copy', pattern: "*.metrics.txt"

    conda "bioconda::gatk4=4.5.0.0 bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.dedup.bam"), path("${meta.id}.dedup.bam.bai"), emit: bam
    path "${meta.id}.dedup.metrics.txt", emit: metrics
    path "versions.yml", emit: versions

    script:
    // --REMOVE_DUPLICATES false keeps duplicate reads in the file but flags them (GATK
    // Best Practices default), so downstream tools that respect the flag (Mutect2 does)
    // exclude them without the reads being discarded outright, in case anyone wants to
    // inspect them later.
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" MarkDuplicates \\
        --INPUT ${bam} \\
        --OUTPUT ${meta.id}.dedup.bam \\
        --METRICS_FILE ${meta.id}.dedup.metrics.txt \\
        --REMOVE_DUPLICATES false \\
        --CREATE_INDEX false \\
        --VALIDATION_STRINGENCY LENIENT

    samtools index ${meta.id}.dedup.bam
    # Indexing via samtools rather than Picard's own --CREATE_INDEX keeps the index
    # filename convention consistent across the whole pipeline ("<file>.bam.bai"); Picard's
    # own indexer names it "<file>.bai" instead, which would silently break every
    # downstream process here that looks for "<bam>.bai" specifically.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/.*(GATK) v//')
    END_VERSIONS
    """
}
