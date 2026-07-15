/*
 * Downloads open-access TCGA-SKCM gene-level RNA-seq counts from the GDC API. No login
 * or data access application is needed: gene-level counts are on GDC's open tier (unlike
 * raw reads or germline variant calls, which require dbGaP authorization), confirmed
 * directly against the live API before any pipeline code was written, returning 473
 * open-access "STAR - Counts" files with no authentication.
 *
 * Note on genome builds: these files use GENCODE v36 / GRCh38 gene models (GDC's
 * harmonized pipeline), while the WES/WGS modules in this project use GRCh37. This does
 * not create a coordinate-mismatch problem anywhere in this project, because the only
 * place TCGA expression and the chr9/chr20 variant calls meet is the TMB-vs-expression
 * correlation in the RNA-seq notebook, which operates on one scalar TMB value and one
 * scalar expression-signature score per sample, not on shared genomic coordinates.
 */
process FETCH_TCGA_SKCM_COUNTS {
    label 'process_single'
    publishDir "${params.raw_dir}/tcga_skcm", mode: 'copy'

    conda "conda-forge::python=3.11"

    input:
    val max_per_group  // samples per sample_type group (Primary Tumor, Metastatic); ~4.2MB each
    // Declared as a real path input (not just referenced by a moduleDir string inside the
    // script block), so a future edit to this script correctly invalidates -resume's
    // cache instead of silently reusing stale output; see modules/prepare_truth_sets.nf
    // for the real instance of this that happened during development.
    path fetch_script

    output:
    path "*.rna_seq.augmented_star_gene_counts.tsv", emit: counts
    path "tcga_skcm_sample_sheet.tsv", emit: sample_sheet

    script:
    """
    python3 ${fetch_script} . ${max_per_group}
    """
}
