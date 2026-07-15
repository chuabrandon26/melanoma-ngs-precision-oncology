#!/usr/bin/env bash
# Runs all three Nextflow pipelines back to back, in the order that makes sense given
# this project's tight disk budget: WES and WGS first (they share the same chr9+chr20
# BAM slice), then RNA-seq. Stops immediately if any pipeline fails, rather than
# continuing on to the next one with an incomplete prior stage.
#
# This does not replace watching each run yourself: check
# results/<module>/pipeline_info/execution_trace.txt (the peak_rss column specifically)
# while these run, that is the actual measured RAM each step used.
#
# Usage: bash scripts/run_all_pipelines.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

echo "================================================================"
echo "Module 1: WES"
echo "================================================================"
nextflow run pipelines/wes/main.nf -profile conda

echo "================================================================"
echo "Module 2: WGS"
echo "================================================================"
nextflow run pipelines/wgs/main.nf -profile conda

echo "================================================================"
echo "Cleaning up the chr9+chr20 BAM slice (WES and WGS are the only"
echo "consumers of it, both are done now; RNA-seq needs none of this)"
echo "================================================================"
# Two places this data actually lives, both need clearing to reclaim the space:
#   1. data/raw/*.bam*: the published copy from each pipeline's publishDir
#   2. .nextflow's work/ cache: WES and WGS each independently pulled their own
#      copy into their own work directories, `nextflow clean` is what actually
#      reclaims that, deleting the published copy alone leaves the bigger chunk
#      of disk usage sitting in work/ untouched.
# Safe to remove: everything either pipeline actually needed to keep has already
# been copied out to results/ and data/processed/ via publishDir (mode: 'copy'),
# work/ is a cache for resumability, not the canonical location for any output.
rm -fv data/raw/*.bam data/raw/*.bam.bai 2>/dev/null || true
nextflow clean -f -q

echo "================================================================"
echo "Module 3: RNA-seq (TCGA-SKCM open counts only, no FASTQ demo)"
echo "================================================================"
nextflow run pipelines/rnaseq/main.nf -profile conda

echo ""
echo "All three pipelines finished. Check results/*/pipeline_info/execution_trace.txt"
echo "for actual RAM/CPU usage per step, then open the notebooks in notebooks/."
