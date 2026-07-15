/*
 * Gene-level quantification with Salmon, used here instead of STAR (the tool named in
 * the original project brief) for one specific, resource-driven reason: a full human
 * STAR genome index needs roughly 27-30GB of disk regardless of RAM tuning, which alone
 * would consume nearly this project's entire disk budget once stacked on top of the
 * chr9+chr20 BAM slice. Salmon indexes the transcriptome instead of the genome. This is
 * a real substitution, documented here and repeated in the RNA-seq notebook rather than
 * made silently: Salmon does pseudo-alignment (k-mer/quasi-mapping against transcripts)
 * rather than STAR's full spliced, base-level alignment to the genome, so no BAM of
 * aligned reads comes out of this step, only transcript-level quantifications, which is
 * exactly what downstream tximport-based gene-level summarization (done in the notebook)
 * expects as input anyway.
 *
 * A second, related simplification: this index is NOT decoy-aware. Salmon's documented
 * best practice adds the full genome as "decoy" sequence so reads that truly originate
 * from unannotated or intergenic DNA are recognised as such, rather than being
 * force-matched to the most similar transcript. Building a decoy-aware index means
 * downloading the full GRCh37 genome FASTA (roughly 830MB-3GB depending on compression),
 * on top of everything else this project already asks of a tight disk budget. Skipping
 * decoys can mildly inflate a small number of low-abundance transcripts' estimated
 * counts; for a relative differential-expression comparison between conditions, rather
 * than an absolute-quantification use case, this is a reasonable, explicitly-flagged
 * tradeoff, not a silent shortcut. The index is built from GENCODE's own pre-extracted
 * protein-coding and lncRNA transcript FASTAs, which together cover the biotypes that
 * essentially all GSEA/GO gene sets are organized around.
 */
process SALMON_INDEX {
    label 'process_medium'
    publishDir "${params.reference_dir}/prepared", mode: 'copy'

    conda "bioconda::salmon=2.3.1"
    container "quay.io/biocontainers/salmon:2.3.1--hfa8f182_0"

    input:
    path pc_transcripts_fa_gz
    path lncrna_transcripts_fa_gz

    output:
    path "salmon_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    cat ${pc_transcripts_fa_gz} ${lncrna_transcripts_fa_gz} > transcripts.fa.gz

    salmon index \\
        -t transcripts.fa.gz \\
        -i salmon_index \\
        -k 31 \\
        -p ${task.cpus}
    # k=31 is Salmon's own default k-mer size, tuned for typical Illumina read lengths
    # (roughly 75bp or longer); not re-tuned for this dataset specifically.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}

process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/rnaseq/salmon", mode: 'copy'

    conda "bioconda::salmon=2.3.1"
    container "quay.io/biocontainers/salmon:2.3.1--hfa8f182_0"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path salmon_index

    output:
    tuple val(meta), path("${meta.id}"), emit: quant_dir
    tuple val(meta), path("${meta.id}/quant.sf"), emit: quant_sf
    path "versions.yml", emit: versions

    script:
    """
    salmon quant \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${reads1} -2 ${reads2} \\
        --validateMappings \\
        --gcBias \\
        -p ${task.cpus} \\
        -o ${meta.id}
    # -l A: auto-detect library strandedness instead of assuming one, Salmon's documented
    # recommendation whenever the library prep protocol has not been independently
    # confirmed, as is the case reusing GSE78220's public samples here.
    # --gcBias: corrects for fragment GC-content bias introduced during library prep. Off
    # by default for backward compatibility, but Salmon's own documented recommendation
    # for new analyses, so enabled explicitly.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}
