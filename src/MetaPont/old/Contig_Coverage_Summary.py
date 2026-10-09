import argparse
import csv
import os
import re
from glob import glob
from statistics import mean, median

try:
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *


def find_samples(root_dir, dir_tags):
    samples = []
    for tag in dir_tags:
        pattern = os.path.join(root_dir, f"{tag}*")
        for path in glob(pattern):
            if os.path.isdir(path):
                samples.append(os.path.abspath(path))
    return sorted(samples)

def find_contig_summary(sample_dir):
    # Common patterns used in the workflow
    candidates = []
    candidates += glob(os.path.join(sample_dir, "*_readmapped", "*contig*summary*.txt"))
    candidates += glob(os.path.join(sample_dir, "*_readmapped", "*idxstats*.txt"))
    candidates += glob(os.path.join(sample_dir, "*_readmapped", "*readmapped_contig_summary*.txt"))
    # fallback: any contig summary in sample dir tree
    candidates += glob(os.path.join(sample_dir, "**", "*contig*summary*.txt"), recursive=True)
    candidates = [c for c in candidates if os.path.isfile(c)]
    return candidates[0] if candidates else None

def parse_contig_summary_file(path):

    contigs = []
    with open(path, 'r') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            parts = re.split(r'\s+', line)
            # detect header like "Contig_Name\tMapped_Reads\tTotal_Reads\tProportion"
            if re.match(r'(?i)contig', parts[0]):
                # skip header line
                continue
            if len(parts) == 5:
                # likely idxstats: name, length, mapped, unmapped
                name = parts[0]
                length = int(parts[1])
                mapped = int(parts[2])
                proportion = float(parts[4])
                contigs.append({"name": name, "length": length, "mapped": mapped, "proportion": proportion})
            # elif len(parts) >= 3 and is_int(parts[1]) and is_float_like(parts[2]):
            #     # e.g. name mapped total/proportion - ambiguous; treat as name,mapped,maybe_total
            #     name = parts[0]
            #     mapped = safe_int(parts[1])
            #     proportion = safe_float(parts[2])
            #     contigs.append({"name": name, "length": None, "mapped": mapped, "proportion": proportion})
            # elif len(parts) >= 2 and is_int(parts[1]):
            #     name = parts[0]
            #     mapped = safe_int(parts[1])
            #     contigs.append({"name": name, "length": None, "mapped": mapped, "proportion": None})
            # else:
            #     # try CSV-like with tabs
            #     cols = line.split('\t')
            #     if len(cols) >= 2 and is_int(cols[1]):
            #         contigs.append({"name": cols[0], "length": None, "mapped": safe_int(cols[1]), "proportion": None})
            #     else:
            #         # skip unparsable line
            #         continue
    return contigs



def compute_stats(contigs, read_length=None):
    mapped_list = [c["mapped"] for c in contigs if c["mapped"] is not None]
    #lengths = [c["length"] for c in contigs if c["length"] and c["length"] > 0]
    rpk_list = []
    coverage_list = []
    prop_list = [c["proportion"] for c in contigs if c.get("proportion") is not None]
    for c in contigs:
        if c.get("length") and c["length"] > 0:
            rpk = c["mapped"] / (c["length"] / 1000.0)
            rpk_list.append(rpk)
            if read_length:
                cov = (c["mapped"] * read_length) / c["length"]
                coverage_list.append(cov)
    stats = {}
    stats["num_contigs"] = len(contigs)
    stats["total_mapped_reads"] = sum(mapped_list) if mapped_list else 0
    stats["avg_mapped_reads"] = mean(mapped_list) if mapped_list else 0
    stats["min_mapped_reads"] = min(mapped_list) if mapped_list else 0
    stats["max_mapped_reads"] = max(mapped_list) if mapped_list else 0
    stats["avg_rpk"] = mean(rpk_list) if rpk_list else None
    stats["min_rpk"] = min(rpk_list) if rpk_list else None
    stats["max_rpk"] = max(rpk_list) if rpk_list else None
    stats["avg_proportion"] = mean(prop_list) if prop_list else None
    stats["median_mapped_reads"] = median(mapped_list) if mapped_list else None
    stats["avg_estimated_coverage"] = mean(coverage_list) if coverage_list else None
    stats["min_estimated_coverage"] = min(coverage_list) if coverage_list else None
    stats["max_estimated_coverage"] = max(coverage_list) if coverage_list else None
    return stats

def write_csv(out_path, rows, fieldnames):
    with open(out_path, 'w', newline='') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

def main():
    parser = argparse.ArgumentParser(description='MetaPont ' + MetaPont_Version + '- Contig-Coverage-Summary: Aggregate readmapping contig summaries and compute overview stats per sample.')
    parser.add_argument("--root_dir", "-d", required=True,
                        help="Root directory containing sample folders (use `root_dir` path).")
    parser.add_argument("--prefix", "-p", required=True,
                        help="Comma-separated directory tags to search for (e.g. E,L,P).")
    parser.add_argument("--read-length", "-r", type=int, default=None,
                        help="Optional average read length to compute estimated coverage.")
    parser.add_argument("--output", "-o",
                        help="Output CSV path (default: `root_dir/readmap_overview.csv`).")
    options = parser.parse_args()

    # determine output path: if user provided one use it, otherwise put default file in root_dir
    if options.output:
        output_path = options.output
    else:
        output_path = os.path.join(options.root_dir, "readmap_overview.csv")

    dir_tags = [t.strip() for t in options.prefix.split(",") if t.strip()]
    samples = find_samples(options.root_dir, dir_tags)
    rows = []
    for s in samples:
        sample_name = os.path.basename(s.rstrip("/"))
        summary_path = find_contig_summary(s)
        if not summary_path:
            # no contig summary found, skip
            continue
        contigs = parse_contig_summary_file(summary_path)
        if not contigs:
            continue
        stats = compute_stats(contigs, read_length=options.read_length)
        row = {
            "sample": sample_name,
            #"summary_path": summary_path,
            "num_contigs": stats["num_contigs"],
            "total_mapped_reads": stats["total_mapped_reads"],
            "avg_mapped_reads": f"{stats['avg_mapped_reads']:.2f}",
            "median_mapped_reads": f"{stats['median_mapped_reads']:.2f}" if stats['median_mapped_reads'] is not None else "",
            "min_mapped_reads": stats["min_mapped_reads"],
            "max_mapped_reads": stats["max_mapped_reads"],
            "avg_rpk": f"{stats['avg_rpk']:.4f}" if stats["avg_rpk"] is not None else "",
            "min_rpk": f"{stats['min_rpk']:.4f}" if stats["min_rpk"] is not None else "",
            "max_rpk": f"{stats['max_rpk']:.4f}" if stats["max_rpk"] is not None else "",
            "avg_proportion": f"{stats['avg_proportion']:.6f}" if stats["avg_proportion"] is not None else "",
            "avg_estimated_coverage": f"{stats['avg_estimated_coverage']:.4f}" if stats["avg_estimated_coverage"] is not None else "",
            "min_estimated_coverage": f"{stats['min_estimated_coverage']:.4f}" if stats["min_estimated_coverage"] is not None else "",
            "max_estimated_coverage": f"{stats['max_estimated_coverage']:.4f}" if stats["max_estimated_coverage"] is not None else "",
        }
        rows.append(row)

    if rows:
        fieldnames = [
            "sample", "num_contigs", "total_mapped_reads",
            "avg_mapped_reads", "median_mapped_reads", "min_mapped_reads", "max_mapped_reads",
            "avg_rpk", "min_rpk", "max_rpk", "avg_proportion",
            "avg_estimated_coverage", "min_estimated_coverage", "max_estimated_coverage"
        ]
        write_csv(output_path, rows, fieldnames)
        print(f"Wrote overview for {len(rows)} samples to {output_path}")
    else:
        print("No samples with contig summaries found.")

if __name__ == "__main__":
    main()

## NOT WORKING