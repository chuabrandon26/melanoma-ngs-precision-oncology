/*
 * Prepares all three COLO829 truth sets for use by the benchmarking processes:
 *   - SNV and indel truth: converted from parklab's custom TSV format into minimal VCFs,
 *     restricted to the requested chromosomes, with "chr" stripped from contig names
 *     (see modules/bin/convert_colo829_truth_to_vcf.py: the source files use "chr9"/
 *     "chr20", everything else in this project uses bare "9"/"20")
 *   - SV truth: already GRCh37 bare-numeric, just subset to the requested chromosomes
 *   - A DEL/DUP-only BED derived from the SV truth set, for the CNV concordance check in
 *     BENCHMARK_CNV, since no independent CNV truth set exists for COLO829
 * (see modules/bin/extract_sv_delDup_bed.py for a real bug this project hit and fixed:
 * this VCF encodes SVs as breakend mate pairs with no INFO/END field, not as symbolic
 * <DEL> calls, so END has to be reconstructed from SVLEN and the correct mate, or every
 * region silently collapses to 1bp)
 */
process PREPARE_TRUTH_SETS {
    label 'process_single'
    publishDir "${params.processed_dir}/truth_sets", mode: 'copy'

    // python is needed for the two helper scripts below, not for bcftools itself; adding
    // it explicitly here since the plain bcftools biocontainer does not bundle a Python
    // interpreter. The conda profile resolves this cleanly; if you use -profile docker or
    // singularity instead, this specific process needs a container that has both
    // bcftools and python3, which the stock bcftools biocontainer pinned below does not,
    // so either run this one process with -profile conda or build a small custom image.
    conda "bioconda::bcftools=1.19 bioconda::htslib=1.19.1 conda-forge::python=3.11"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_0"

    input:
    path snv_truth_txt
    path indel_truth_txt
    path sv_truth_vcf
    val chromosomes  // e.g. "9,20"
    // Declared as real path inputs, not just referenced by a moduleDir string inside the
    // script block: Nextflow's -resume cache hash only covers the process's own inline
    // script text plus its declared inputs, not the content of external files a script
    // happens to call by path. This project hit that gap in practice, editing these
    // helper scripts and re-running with -resume kept silently reusing the old, buggy
    // cached output. Declaring them here means a content change is correctly seen as a
    // changed input, invalidating the cache the way it should.
    path convert_script
    path extract_script

    output:
    tuple path("truth_snv.vcf.gz"), path("truth_snv.vcf.gz.tbi"), emit: snv
    tuple path("truth_indel.vcf.gz"), path("truth_indel.vcf.gz.tbi"), emit: indel
    tuple path("truth_snv_indel.vcf.gz"), path("truth_snv_indel.vcf.gz.tbi"), emit: snv_indel_combined
    tuple path("truth_sv.vcf.gz"), path("truth_sv.vcf.gz.tbi"), emit: sv
    path "truth_sv_delDup.bed", emit: sv_deldup_bed
    path "versions.yml", emit: versions

    script:
    """
    python3 ${convert_script} ${snv_truth_txt} truth_snv.unsorted.vcf SNV ${chromosomes}
    python3 ${convert_script} ${indel_truth_txt} truth_indel.unsorted.vcf INDEL ${chromosomes}
    bcftools sort truth_snv.unsorted.vcf -O z -o truth_snv.vcf.gz
    bcftools sort truth_indel.unsorted.vcf -O z -o truth_indel.vcf.gz
    tabix -p vcf truth_snv.vcf.gz
    tabix -p vcf truth_indel.vcf.gz

    # hap.py expects a single truth VCF containing both SNVs and indels together (like
    # GIAB's own truth releases do) and separates them internally for its per-type
    # precision/recall/F1 reporting, rather than being run once per variant type.
    bcftools concat -a truth_snv.vcf.gz truth_indel.vcf.gz -O z -o truth_snv_indel.vcf.gz
    tabix -p vcf truth_snv_indel.vcf.gz

    # The SV truth set is already GRCh37 bare-numeric, it only needs subsetting.
    bgzip -c ${sv_truth_vcf} > sv_truth_full.vcf.gz
    tabix -p vcf sv_truth_full.vcf.gz
    bcftools view -r ${chromosomes} sv_truth_full.vcf.gz -O z -o truth_sv.vcf.gz
    tabix -p vcf truth_sv.vcf.gz

    # 1kb floor: below that, CNVkit's bin resolution on a chr9+chr20-only WGS run at this
    # depth is already close to its own detection floor, so smaller truth events would not
    # be a meaningful comparison point either way.
    python3 ${extract_script} ${sv_truth_vcf} truth_sv_delDup.bed 1000

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}

/*
 * Converts the BICseq2 copy-number reference (COLO829.somatic.bicseq.cna, NOT an
 * independently validated truth set, see BENCHMARK_CNV) into a plain, chr9/chr20-only
 * BED. The source file has a header row, "chr"-prefixed contigs, and several extra
 * columns beyond chrom/start/end (bin count, raw and expected read counts per group,
 * two log2-ratio columns); this keeps only what bedtools needs, columns 1-3 exactly as
 * the file's own header names them (chrom, start, end), so there is no guessing about
 * column order.
 */
process PREPARE_CNV_REFERENCE {
    label 'process_single'
    publishDir "${params.processed_dir}/truth_sets", mode: 'copy'

    conda "conda-forge::gawk=5.3.0"

    input:
    path bicseq_cna
    val chromosomes  // e.g. "9,20", bare numeric

    output:
    path "bicseq_reference.bed", emit: bed

    script:
    // Filtering by awk field comparison ($1=="chr9"), not by piping through grep with a
    // hand-built "^chr9\t" regex: this project tried the grep-with-embedded-tab approach
    // first and it silently matched zero lines, because a literal backslash-t inside a
    // single-quoted ERE pattern is not reliably treated as an actual tab character by
    // plain `grep -E`. awk's own field splitting sidesteps the whole escaping problem.
    def chrom_conditions = chromosomes.tokenize(',').collect { "\$1==\"chr${it}\"" }.join(' || ')
    """
    tail -n +2 ${bicseq_cna} \\
        | awk -v OFS='\\t' '${chrom_conditions} {print \$1, \$2, \$3}' \\
        | sed 's/^chr//' \\
        | sort -k1,1 -k2,2n \\
        > bicseq_reference.bed
    """
}
