/*
 * Pulls only the requested chromosomes out of a remote, indexed BAM over HTTPS range
 * requests, without ever downloading the full file, optionally downsampling coverage
 * at the same time.
 *
 * This is the single most important module in the project: it is what keeps this
 * pipeline's disk and bandwidth footprint down to a few GB per sample instead of the
 * 267.9GB the two full BAMs add up to. Before this module was written, the exact ENA
 * URLs below were checked directly: a ranged HTTPS request against the 196.8GB tumor
 * BAM returned HTTP 206 Partial Content with a correct Content-Range header, and the
 * matching .bai index resolves over HTTPS too, which is what lets samtools seek straight
 * to chr9 and chr20 instead of streaming the whole file.
 *
 * meta.subsample_fraction (set by the calling pipeline from params.downsample_enabled/
 * target_*_depth/native_*_depth) is applied in the same samtools pass as the region
 * filter, at native depth (fraction 1.0, no -s flag added) when downsampling is
 * disabled. See nextflow.config for why the default changed to downsampled: real disk
 * pressure during development on the machine this was built for.
 */
process SLICE_REMOTE_BAM {
    tag "${meta.id}"
    label 'process_high'
    // Not CPU-bound so much as "give it enough wall-clock and memory headroom for a slow
    // connection and BGZF re-compression of a multi-GB region"; process_high's 2h/12GB
    // ceiling is a safety margin, not an expectation that either number will be hit.

    publishDir "${params.raw_dir}", mode: 'copy'

    conda "bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), val(bam_url)
    val chromosomes   // e.g. "9,20", bare numeric to match this BAM's own contig naming

    output:
    tuple val(meta), path("${meta.id}.sliced.bam"), path("${meta.id}.sliced.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    def regions = chromosomes.tokenize(',').join(' ')
    def max_attempts = 5
    // Fixed seed (100): makes the subsample deterministic and reproducible across runs
    // for the same fraction, rather than a different random subset every re-run.
    def fraction = (meta.subsample_fraction ?: 1.0) as double
    def frac_str = String.format('%.3f', fraction).replaceFirst('^0\\.', '')
    def subsample_flag = fraction < 0.999 ? "-s 100.${frac_str}" : ""
    // Same exact-match, per-field awk approach as BUILD_EXOME_BED in
    // prepare_reference.nf, not a hand-built grep regex (that already caused a real,
    // silent-no-match bug earlier in this project). Scans every field of each @SQ line
    // for an exact "SN:<chrom>" match rather than assuming SN: is always field 2: true
    // in practice for every @SQ line samtools itself emits, but not guaranteed by the
    // SAM spec, so not worth relying on.
    def sq_keep_conditions = chromosomes.tokenize(',').collect { "\$i==\"SN:${it}\"" }.join(' || ')
    // Same idea applied to the two other places a chromosome name can appear in a read
    // beyond its own RNAME: the RNEXT (mate reference) field, and the SA:Z:/XA:Z:
    // optional tags (split-read / alternate-hit alignments). See the two comment blocks
    // further down for why both need fixing up, not just the header.
    def rnext_keep_conditions = chromosomes.tokenize(',').collect { "\$7==\"${it}\"" }.join(' || ')
    def tag_chrom_keep_conditions = chromosomes.tokenize(',').collect { "chrom==\"${it}\"" }.join(' || ')
    """
    # -h keeps the header, -b writes BAM. The region list after the URL restricts the read
    # to only those contigs: samtools uses the remote .bai (expected at "<bam_url>.bai")
    # to seek directly to each one rather than streaming the whole file past them. This
    # requires samtools built with libcurl/HTS remote-file support, which the bioconda
    # and biocontainer builds pinned above both have. -s applies the coverage downsample
    # (a per-read random keep-probability) in the same pass, when enabled.
    #
    # Wrapped in a retry loop: this is a long-running series of HTTPS range requests
    # against a multi-GB remote file, which is exposed to transient network blips a
    # single attempt has no way to recover from. This exact failure showed up in
    # practice (a bgzf_read error partway into the first requested region, the
    # connection apparently dropped mid-stream), not a problem with the range-request
    # approach itself, which was verified working (a correct 206 Partial Content
    # response) before this pipeline was written. Retrying is cheap and safe here: the
    # operation is idempotent and each attempt overwrites rather than appends.
    #
    # Output of this loop is an intermediate (${meta.id}.raw.bam), not the process's
    # declared output: samtools view/-b here only restricts which ALIGNMENT RECORDS come
    # out, it does not prune the @SQ header dictionary, so the result still declares all
    # 25 original GRCh37 contigs even though only ${chromosomes} reads are actually
    # present. See the reheader step below for why that matters and how it's fixed.
    attempt=1
    max_attempts=${max_attempts}
    until samtools view -h -b \\
            --threads ${task.cpus} \\
            ${subsample_flag} \\
            -o ${meta.id}.raw.bam \\
            ${bam_url} ${regions} \\
          && samtools quickcheck ${meta.id}.raw.bam
    do
        echo "SLICE_REMOTE_BAM attempt \${attempt}/\${max_attempts} failed" >&2
        rm -f ${meta.id}.raw.bam
        if [ "\${attempt}" -ge "\${max_attempts}" ]; then
            echo "giving up after \${max_attempts} attempts" >&2
            exit 1
        fi
        attempt=\$((attempt + 1))
        sleep \$((attempt * 15))
    done

    # ${meta.id}.raw.bam's header still lists all 25 original contigs (see above), which
    # MarkDuplicates/Mutect2/CNVkit tolerate (they only validate contigs actually used by
    # reads/intervals against the reference) but Manta's configManta.py does not: it
    # validates every @SQ entry against the reference FASTA regardless of use, and fails
    # ("Reference fasta file is missing a chromosome found in the Normal BAM/CRAM file:
    # '1'") since genome.fa intentionally only contains ${chromosomes}.
    #
    # Fixed with a full SAM-text round-trip, not `samtools reheader`: reheader only swaps
    # the header bytes and leaves each alignment record's binary reference-ID integer
    # untouched, so dropping @SQ entries with it corrupts every record whose refID then
    # points past the truncated dictionary. Confirmed directly against this project's own
    # sliced COLO829R BAM before relying on it: samtools quickcheck and view -H on a
    # reheader'd file both looked fine, but samtools view and samtools index on the same
    # file both failed with "Numerical result out of range". Round-tripping through SAM
    # text is safe because RNAME there is the string "9"/"20", re-resolved by name against
    # whichever header is in effect at encode time, not a positional index.
    #
    # Two more places a now-removed chromosome name can still be lurking in a read besides
    # the header, both fixed in the same awk pass rather than left for samtools/Manta to
    # trip over later:
    #
    # 1. RNEXT (mate reference, field 7): a read on chr9/chr20 whose mate lands outside
    #    ${chromosomes} still has its own FLAG claiming "mate mapped" once encoded against
    #    the trimmed header, while the mate's chromosome no longer resolves -- an
    #    internally inconsistent record (`samtools reheader`-style tools cannot fix this;
    #    `samtools fixmate` cannot either, confirmed directly: it only reconciles metadata
    #    between two mates that are both actually present in the file, and leaves a
    #    singleton whose mate is entirely absent untouched). Manta's SVLocusScanner does
    #    not tolerate this and throws a fatal, unrecoverable error ("SVbreakend has unknown
    #    or invalid chromosome id"), confirmed directly by running Manta against data with
    #    this left unfixed. Fix: whenever RNEXT points outside ${chromosomes}, force it to
    #    "*", zero PNEXT/TLEN, and correct the FLAG (set mate-unmapped 0x8, clear
    #    proper-pair 0x2) so the record is internally self-consistent. Confirmed by an
    #    exact before/after comparison of every read's own QNAME/RNAME/POS (never touched)
    #    plus idxstats (mapped/unmapped counts per contig unchanged, i.e. no reads lost).
    #
    # 2. SA:Z:/XA:Z: (supplementary/alternate-alignment tags): split-read evidence encodes
    #    the other alignment segment's chromosome as plain text in these optional tags,
    #    completely independent of RNEXT/FLAG. Manta parses SA:Z: directly for split-read
    #    SV evidence and throws the same kind of fatal error ("Split alignment segment maps
    #    to an unknown chromosome") when an entry points outside ${chromosomes}, confirmed
    #    directly the same way as (1). XA:Z: (BWA's alternate-hit tag) has the identical
    #    "semicolon-separated chrom,pos,...  entries" shape and is stripped the same way as
    #    a defensive measure, whether or not Manta itself reads it. Fix: for either tag,
    #    keep only the semicolon-separated entries whose chromosome is in ${chromosomes},
    #    dropping the whole tag if none remain. This only ever removes a claim about where
    #    an unmapped-in-our-scope segment sits; it never touches the read's own RNAME/POS/
    #    CIGAR/FLAG-as-primary-record.
    samtools view -h ${meta.id}.raw.bam \\
        | awk '
    BEGIN { OFS = "\\t" }
    \$1 == "@SQ" {
        keep = 0
        for (i = 2; i <= NF; i++) { if (${sq_keep_conditions}) keep = 1 }
        if (keep) print
        next
    }
    /^@/ { print; next }
    {
        if (\$7 != "*" && \$7 != "=" && !(${rnext_keep_conditions})) {
            flag = \$2
            if (int(flag / 8) % 2 == 0) flag += 8
            if (int(flag / 2) % 2 == 1) flag -= 2
            \$2 = flag
            \$7 = "*"
            \$8 = 0
            \$9 = 0
        }
        line = \$1
        for (i = 2; i <= 11; i++) line = line OFS \$i
        for (i = 12; i <= NF; i++) {
            tag = \$i
            if (substr(tag, 1, 5) == "SA:Z:" || substr(tag, 1, 5) == "XA:Z:") {
                prefix = substr(tag, 1, 5)
                n = split(substr(tag, 6), entries, ";")
                newval = ""
                for (j = 1; j <= n; j++) {
                    if (entries[j] == "") continue
                    split(entries[j], parts, ",")
                    chrom = parts[1]
                    if (${tag_chrom_keep_conditions}) newval = newval entries[j] ";"
                }
                if (newval != "") line = line OFS prefix newval
            } else {
                line = line OFS tag
            }
        }
        print line
    }
    ' \\
        | samtools view -b --threads ${task.cpus} -o ${meta.id}.sliced.bam -

    rm -f ${meta.id}.raw.bam
    samtools quickcheck ${meta.id}.sliced.bam
    samtools index ${meta.id}.sliced.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
