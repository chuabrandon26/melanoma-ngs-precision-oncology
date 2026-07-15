/*
 * Somatic SNV/indel calling with GATK Mutect2, tumor vs matched normal, followed by
 * GATK's standard FilterMutectCalls post-processing. Shared by both pipelines: the WES
 * pipeline passes an exome-intersected BED as `intervals`, the WGS pipeline passes the
 * whole chr9+chr20 BED, everything else is identical.
 */
process CALL_SNV_MUTECT2 {
    tag "${meta.id}"
    label 'process_high'
    // Mutect2's local reassembly (it builds a De Bruijn graph per active region) is the
    // most memory-hungry step in this whole project. Even scoped to two chromosomes it
    // can spike, so this is one of the few processes pinned at the full 12GB ceiling
    // rather than process_medium's 8GB.

    // publishDir needs a closure here, not a plain interpolated string: directives are
    // otherwise evaluated once at process-definition time, before any task's `meta`
    // input exists yet. `tag "${meta.id}"` above works without one only because `tag`
    // is specifically documented to support lazy per-task evaluation of a plain string;
    // `publishDir` needs the closure to get the same per-task, dynamic behavior. Found
    // by an actual `nextflow run -preview` failure ("No such variable: meta"), not
    // assumed.
    publishDir(path: { "${params.processed_dir}/${meta.module}" }, mode: 'copy')

    conda "bioconda::gatk4=4.5.0.0"
    container "quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_fai
    path fasta_dict
    path intervals  // exome-intersected BED for WES, whole chr9+chr20 BED for WGS
    tuple path(germline_resource), path(germline_resource_tbi)

    output:
    tuple val(meta), path("${meta.id}.mutect2.vcf.gz"), path("${meta.id}.mutect2.vcf.gz.tbi"), path("${meta.id}.mutect2.vcf.gz.stats"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 \\
        -R ${fasta} \\
        -I ${tumor_bam} \\
        -tumor ${params.tumor_sample_name} \\
        -I ${normal_bam} \\
        -normal ${params.normal_sample_name} \\
        -L ${intervals} \\
        --germline-resource ${germline_resource} \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O ${meta.id}.mutect2.vcf.gz
    # --germline-resource: a 1000 Genomes Phase 3 chr9+chr20 population-allele-frequency
    # VCF (see modules/prepare_reference.nf's PREPARE_GERMLINE_RESOURCE for the real
    # finding that motivated adding this: an initial WES run with no germline resource at
    # all showed most "somatic" PASS calls had no relationship to the truth set).
    #
    # Measured, not assumed, and worth being honest about: adding this made only a small
    # difference in practice (checked directly: POPAF, the population-frequency
    # annotation this resource feeds, came back as Mutect2's generic "no match" default
    # of 6.00 for 1187 of ~1223 candidates here). 1000 Genomes has only 2,504 samples and
    # dates to 2013; a real clinical pipeline would use gnomAD (hundreds of thousands of
    # samples) for this, which was ruled out earlier in this project specifically because
    # its own af-only resource is a multi-GB, genome-scale download even before
    # subsetting, the same disk-budget reasoning that shaped most other choices here. So
    # this remains a real, deliberate simplification versus clinical-grade Mutect2 usage,
    # not a solved problem: expect WES precision/recall against the truth set to be
    # noticeably weaker than a clinical-grade run would show, on top of the noise already
    # expected from benchmarking against only 7.4Mb of exonic sequence. No panel-of-
    # normals is supplied either, for the same reason: building one needs many unrelated
    # normal samples run through this same pipeline, which a single-pair portfolio
    # project cannot construct.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/.*(GATK) v//')
    END_VERSIONS
    """
}

process FILTER_MUTECT_CALLS {
    tag "${meta.id}"
    label 'process_medium'

    // publishDir needs a closure here, not a plain interpolated string: directives are
    // otherwise evaluated once at process-definition time, before any task's `meta`
    // input exists yet. `tag "${meta.id}"` above works without one only because `tag`
    // is specifically documented to support lazy per-task evaluation of a plain string;
    // `publishDir` needs the closure to get the same per-task, dynamic behavior. Found
    // by an actual `nextflow run -preview` failure ("No such variable: meta"), not
    // assumed.
    publishDir(path: { "${params.processed_dir}/${meta.module}" }, mode: 'copy')

    conda "bioconda::gatk4=4.5.0.0"
    container "quay.io/biocontainers/gatk4:4.5.0.0--py36hdfd78af_0"

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: vcf
    path "${meta.id}.filtered.vcf.gz.filteringStats.tsv", emit: filtering_stats
    path "versions.yml", emit: versions

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" FilterMutectCalls \\
        -R ${fasta} \\
        -V ${vcf} \\
        -O ${meta.id}.filtered.vcf.gz
    # FilterMutectCalls applies GATK's standard Mutect2 filter suite (a model covering
    # strand bias, base and mapping quality, estimated contamination, read-orientation
    # artifacts, and several others) using the tool's own built-in defaults throughout.
    # These are tuned by the GATK team on real tumor/normal cohorts and are the standard
    # starting point; no individual filter threshold has been re-tuned for this specific
    # dataset. If precision/recall in the benchmarking notebook looks unusually skewed for
    # one filter category in particular (visible in the filteringStats.tsv this step also
    # produces), that threshold is the first place to look, not evidence the pipeline
    # logic itself is wrong.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/.*(GATK) v//')
    END_VERSIONS
    """
}

/*
 * Extracts just the tumor sample's genotype column from a tumor+normal VCF.
 *
 * Mutect2's output carries both samples' genotypes together in one file, which is
 * correct for the VCF itself, but RTG's vcfeval (hap.py's comparison engine, used in
 * BENCHMARK_SNV_INDEL) refuses to evaluate a multi-sample "calls" file without being
 * told which sample to score ("No sample name provided but calls is a multi-sample
 * VCF", a real failure hit while running this). The tumor sample is what is actually
 * being benchmarked against the truth set, so it is split out here into its own small,
 * single-purpose process, kept separate from BENCHMARK_SNV_INDEL specifically because
 * bcftools and hap.py 0.3.15 turned out to be unsolvable in the same conda environment
 * (hap.py's legacy Python 2.7 / pysam stack needs zlib <1.3.0a0, modern bcftools/htslib
 * needs zlib >=1.3.1), not because the two steps are conceptually unrelated.
 */
process EXTRACT_TUMOR_VCF {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::bcftools=1.19"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_0"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.tumor_only.vcf.gz"), path("${meta.id}.tumor_only.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view -s ${params.tumor_sample_name} ${vcf} -O z -o ${meta.id}.tumor_only.vcf.gz
    tabix -p vcf ${meta.id}.tumor_only.vcf.gz
    """
}
