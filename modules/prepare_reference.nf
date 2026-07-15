/*
 * One-time reference preparation shared by the WES and WGS pipelines: builds the
 * combined chr9+chr20 FASTA and its indices, plus the small BED/region files Manta and
 * CNVkit each require in their own specific formats.
 */
process PREPARE_REFERENCE {
    label 'process_low'
    publishDir "${params.reference_dir}/prepared", mode: 'copy'

    conda "bioconda::samtools=1.19.2 bioconda::gatk4=4.5.0.0 bioconda::htslib=1.19.1"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    path fasta_chr9_gz
    path fasta_chr20_gz

    output:
    path "genome.fa", emit: fasta
    path "genome.fa.fai", emit: fai
    path "genome.dict", emit: dict
    path "genome.access.bed", emit: access_bed
    path "genome.callregions.bed.gz", emit: call_regions_bed
    path "genome.callregions.bed.gz.tbi", emit: call_regions_bed_tbi
    path "versions.yml", emit: versions

    script:
    """
    zcat ${fasta_chr9_gz} ${fasta_chr20_gz} > genome.fa
    samtools faidx genome.fa
    gatk CreateSequenceDictionary -R genome.fa -O genome.dict --java-options "-Xmx${task.memory.toGiga()}g"

    # Manta's --callRegions and CNVkit's --access both just need "where in this reference
    # is there sequence to call on", which here is simply the full extent of chr9 and
    # chr20 (we did not slice out any sub-regions of either chromosome).
    awk 'BEGIN{OFS="\\t"} {print \$1, 0, \$2}' genome.fa.fai > genome.access.bed
    sort -k1,1 -k2,2n genome.access.bed | bgzip -c > genome.callregions.bed.gz
    tabix -p bed genome.callregions.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
        gatk4: \$(gatk --version 2>&1 | head -n1 | sed 's/.*(GATK) v//')
    END_VERSIONS
    """
}

/*
 * Builds an approximate exome-capture BED for chr9+chr20 from GENCODE exon features,
 * merging overlapping exons into contiguous intervals.
 *
 * This stands in for a commercial capture kit BED (Agilent SureSelect, Twist, etc.),
 * which typically needs a vendor account to download. A GENCODE-exon BED is a defensible
 * proxy for "which regions an exome kit would target", but it is not identical to any
 * specific kit: real kits add probe padding around exon boundaries and make their own
 * inclusion/exclusion calls for difficult regions (segmental duplications, extreme
 * GC content), so the exact read count landing inside this BED will differ somewhat from
 * what a specific commercial kit would give you. This distinction is repeated in the
 * WES notebook rather than left implicit.
 */
process BUILD_EXOME_BED {
    label 'process_low'
    publishDir "${params.reference_dir}/prepared", mode: 'copy'

    conda "bioconda::bedtools=2.31.1"
    container "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"

    input:
    path gtf_gz
    val chromosomes

    output:
    path "exome_regions.bed", emit: bed
    path "versions.yml", emit: versions

    script:
    // Chromosome filtering happens via awk field comparison ($1=="chr9"), not via a
    // separate grep step with a hand-built "^chr9\t" regex: this project's first attempt
    // at that (in PREPARE_CNV_REFERENCE, same idea) silently matched nothing, because a
    // literal backslash-t inside a single-quoted ERE pattern is not reliably treated as
    // an actual tab character by plain `grep -E`. Folding the chromosome check into the
    // same awk pass that already does the exon filtering avoids the problem entirely.
    def chrom_conditions = chromosomes.tokenize(',').collect { "\$1==\"chr${it}\"" }.join(' || ')
    """
    # GENCODE v19 GTF coordinates are 1-based closed intervals; "\$4 - 1" converts the
    # start to BED's 0-based half-open convention. Filtering happens on the "chr9"/"chr20"
    # names GENCODE actually uses, then "chr" is stripped afterwards to match this
    # project's bare-numeric contig convention (the BAMs, reference FASTA, and SV truth
    # set all use "9"/"20", not "chr9"/"chr20").
    zcat ${gtf_gz} \\
        | awk -v OFS='\\t' '\$3 == "exon" && (${chrom_conditions}) {print \$1, \$4-1, \$5}' \\
        | sed 's/^chr//' \\
        | sort -k1,1 -k2,2n \\
        | bedtools merge -i - \\
        > exome_regions.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}

/*
 * Builds a small, chr9+chr20-only germline population-allele-frequency resource for
 * Mutect2's --germline-resource argument, from 1000 Genomes Phase 3 (GRCh37, confirmed
 * bare-numeric contigs and HTTPS range-request support before this was written, the same
 * check applied to the COLO829 BAMs at the very start of this project).
 *
 * Added after a real finding, not planned from the start: an initial WES run with no
 * germline resource at all produced very poor apparent precision/recall (2.6%/0.3% for
 * SNPs) against the truth set. Cross-checking the actual PASS calls against the full
 * chr9 truth set showed most had no relationship to any known somatic mutation, which,
 * combined with several calls showing a clean-looking but germline-consistent pattern
 * (near-total VAF in "tumor", zero coverage-supported evidence in the downsampled
 * "normal"), pointed at germline sites being misclassified as somatic without a
 * population resource to check them against. This is exactly what --germline-resource is
 * for. The full GATK-recommended resource for this (gnomAD, af-only) is a multi-GB,
 * genome-scale file; 1000 Genomes restricted to just chr9+chr20, with per-sample
 * genotypes dropped (Mutect2 only needs the site-level population AF, not individual
 * genotypes for 2,504 unrelated people), keeps this a comparably small, one-time fetch.
 */
process PREPARE_GERMLINE_RESOURCE {
    label 'process_low'
    publishDir "${params.reference_dir}/prepared", mode: 'copy'

    conda "bioconda::bcftools=1.19 bioconda::htslib=1.19.1"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_0"

    input:
    val chr9_url
    val chr20_url

    output:
    tuple path("germline_resource.vcf.gz"), path("germline_resource.vcf.gz.tbi"), emit: resource
    path "versions.yml", emit: versions

    script:
    """
    # -G / --drop-genotypes: keep only site-level INFO (including AF), discard all
    # 2,504 individual sample genotype columns, which is both unneeded by Mutect2 here
    # and the main reason the source files are as large as they are.
    bcftools view -G ${chr9_url} -O z -o chr9.vcf.gz
    bcftools view -G ${chr20_url} -O z -o chr20.vcf.gz
    tabix -p vcf chr9.vcf.gz
    tabix -p vcf chr20.vcf.gz

    bcftools concat chr9.vcf.gz chr20.vcf.gz -O z -o germline_resource.vcf.gz
    tabix -p vcf germline_resource.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
