#!/usr/bin/env python3
"""
Queries the GDC API for open-access TCGA-SKCM gene-level RNA-seq counts (the
"STAR - Counts" workflow), selects up to N samples per sample_type group, downloads each
selected file, and writes a sample sheet the RNA-seq notebook uses to build its DESeq2
sample metadata.

No authentication is needed: GDC serves these specific files (gene-level counts, not raw
reads or germline variants) on its open tier over plain HTTPS with no token required,
confirmed during this project's initial data-sizing pass before any code was written.

Usage:
    fetch_tcga_skcm_counts.py <output_dir> <max_per_group>
"""
import json
import os
import sys
import urllib.parse
import urllib.request

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data"

# Primary Tumor vs Metastatic is TCGA-SKCM's own natural, well-populated two-group split
# (roughly 103 primary vs roughly 368 metastatic samples cohort-wide), and a real,
# biologically meaningful contrast (metastatic melanoma has documented expression
# differences from primary lesions), not an arbitrary grouping invented for this project.
WANTED_GROUPS = ["Primary Tumor", "Metastatic"]


def gdc_query(max_results):
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-SKCM"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,cases.case_id,cases.submitter_id,cases.samples.sample_type",
        "format": "JSON",
        "size": str(max_results),
    }
    url = f"{GDC_FILES_URL}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=60) as resp:
        return json.load(resp)


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <output_dir> <max_per_group>")
    out_dir, max_per_group = sys.argv[1], int(sys.argv[2])
    os.makedirs(out_dir, exist_ok=True)

    result = gdc_query(max_results=2000)
    hits = result["data"]["hits"]
    print(f"GDC returned {len(hits)} open-access TCGA-SKCM STAR-count files total", file=sys.stderr)

    by_group = {}
    for hit in hits:
        cases = hit.get("cases", [])
        if not cases:
            continue
        samples = cases[0].get("samples", [])
        if not samples:
            continue
        sample_type = samples[0].get("sample_type", "Unknown")
        by_group.setdefault(sample_type, []).append(hit)

    print("Sample type breakdown:", {k: len(v) for k, v in by_group.items()}, file=sys.stderr)

    selected = []
    for group in WANTED_GROUPS:
        group_hits = by_group.get(group, [])[:max_per_group]
        selected.extend((group, hit) for hit in group_hits)

    sample_sheet_path = os.path.join(out_dir, "tcga_skcm_sample_sheet.tsv")
    with open(sample_sheet_path, "w") as sheet:
        sheet.write("file_id\tfile_name\tcase_id\tsubmitter_id\tsample_type\tlocal_path\n")
        for group, hit in selected:
            file_id = hit["file_id"]
            file_name = hit["file_name"]
            case = hit["cases"][0]
            case_id = case.get("case_id", "NA")
            submitter_id = case.get("submitter_id", "NA")
            local_path = os.path.join(out_dir, file_name)

            data_url = f"{GDC_DATA_URL}/{file_id}"
            urllib.request.urlretrieve(data_url, local_path)

            sheet.write(f"{file_id}\t{file_name}\t{case_id}\t{submitter_id}\t{group}\t{local_path}\n")
            print(f"downloaded {file_name} ({group})", file=sys.stderr)

    print(f"Done: {len(selected)} files, sample sheet at {sample_sheet_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
