"""Compatibility command entry points retained from MetaPont 0.0.10."""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
from pathlib import Path
from typing import Any

from .analysis import coverage_summary, query_function, query_taxon
from .constants import FEATURE_COLUMNS
from .layout import discover_samples
from .parsers import open_text
from .pipeline import OUTPUT_DIRNAME, run_analysis, write_tsv


def _find_gene_table(path_str: str) -> Path:
    path = Path(path_str).expanduser().resolve()
    if path.is_file():
        return path
    candidates = [
        path / "combined" / "gene_annotations.tsv.gz",
        path / OUTPUT_DIRNAME / "combined" / "gene_annotations.tsv.gz",
        path / "MetaPont_Results" / "combined" / "gene_annotations.tsv.gz",
        path.parent / "combined" / "gene_annotations.tsv.gz",
        path.parent / "MetaPont_Results" / "combined" / "gene_annotations.tsv.gz",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"Integrated gene table was not found near {path}. Run `metapont run -d ROOT -p PREFIX` first, "
        "then query MetaPont_Results/combined/gene_annotations.tsv.gz."
    )


def combine_main() -> int:
    parser = argparse.ArgumentParser(description="Integrate eggNOG, Kraken2 and readmapping output across samples.")
    parser.add_argument("-d", "--parent_directory_path", "--root-dir", dest="root_dir", required=True)
    parser.add_argument("-p", "--prefix", required=True, help="Comma-separated sample prefixes, e.g. 3,4,N,P.")
    parser.add_argument("-o", "--output-dir", default=None)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-heal-flags", action="store_true")
    args = parser.parse_args()
    result = run_analysis(args.root_dir, args.prefix, args.output_dir, heal_flags=not args.no_heal_flags, force=args.force)
    print(f"MetaPont integration written to {result['output_dir']}")
    print(f"Samples: {result['samples']}; samples with warnings: {result.get('samples_with_warnings', 0)}")
    print(f"Completion marker: {result.get('success_flag', '')}")
    return 0


def coverage_main() -> int:
    parser = argparse.ArgumentParser(description="Summarise contig-level readmapping support by sample.")
    parser.add_argument("-d", "--root_dir", "--root-dir", dest="root_dir", required=True)
    parser.add_argument("-p", "--prefix", required=True)
    parser.add_argument("-r", "--read-length", type=int, default=None)
    parser.add_argument("-o", "--output", default=None)
    args = parser.parse_args()
    out, count = coverage_summary(args.root_dir, args.prefix, args.output, args.read_length)
    print(f"Wrote coverage summary for {count} samples: {out}")
    return 0


def _query_by_function_main() -> int:
    parser = argparse.ArgumentParser(description="Identify taxa with support for a functional annotation.")
    parser.add_argument("-d", "--directory", required=True, help="MetaPont output, root directory, or combined gene table path.")
    parser.add_argument("-f", "--function_id", required=True)
    parser.add_argument("-a", "--file_affix", default=None, help="Retained for compatibility; new query uses the integrated gene table.")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-m", "--min_proportion", type=float, default=None)
    parser.add_argument("-top", "--top_taxa", type=int, default=None)
    parser.add_argument("-c", "--category", default=None)
    args = parser.parse_args()
    if args.min_proportion is None and args.top_taxa is None:
        parser.error("specify either --min_proportion or --top_taxa")
    if args.min_proportion is not None and args.top_taxa is not None:
        parser.error("--min_proportion and --top_taxa are mutually exclusive")
    table = _find_gene_table(args.directory)
    # Use a temporary query result if proportion filtering is requested.
    if args.min_proportion is not None:
        with tempfile.TemporaryDirectory(prefix="metapont-query-") as tmp:
            temp_out = Path(tmp) / "query.tsv"
            query_function(table, args.function_id, temp_out, category=args.category, rank="Genus", top_n=None)
            rows = _read_tsv_rows(temp_out)
            threshold = args.min_proportion * 100 if args.min_proportion <= 1 else args.min_proportion
            rows = [row for row in rows if float(row.get("Percent_of_function_gene_read_support", 0) or 0) >= threshold]
            write_tsv(Path(args.output), rows, list(rows[0].keys()) if rows else ["Sample", "Taxonomic_Rank", "Taxon", "Function_ID", "Category", "Annotated_Genes", "Gene_Read_Support", "Percent_of_function_gene_read_support"])
    else:
        query_function(table, args.function_id, args.output, category=args.category, rank="Genus", top_n=args.top_taxa)
    print(f"Wrote results to {args.output}")
    return 0


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract_by_taxa_main() -> int:
    parser = argparse.ArgumentParser(description="Extract top eggNOG functions linked to a taxon.")
    parser.add_argument("-d", "--directory", required=True)
    parser.add_argument("-t", "--taxon", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-func", "--functional_classes", required=True, help="Comma-separated MetaPont category labels, e.g. GO,EC,KEGG_KO.")
    parser.add_argument("-top", "--top_functions", type=int, default=3)
    parser.add_argument("--rank", default="Genus")
    args = parser.parse_args()
    table = _find_gene_table(args.directory)
    categories = [v.strip() for v in args.functional_classes.split(",") if v.strip()]
    all_rows: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="metapont-tax-query-") as tmp:
        for category in categories:
            out = Path(tmp) / f"{category.replace('/', '_')}.tsv"
            query_taxon(table, args.taxon, out, category, args.top_functions, args.rank)
            all_rows.extend(_read_tsv_rows(out))
    all_rows.sort(key=lambda r: (r.get("Sample", ""), -int(r.get("Gene_Read_Support", 0) or 0)))
    write_tsv(Path(args.output), all_rows, ["Sample", "Taxon_Query", "Taxonomic_Rank", "Taxon", "Taxon_Lineage", "Category", "Feature", "Annotated_Genes", "Gene_Read_Support"])
    print(f"Wrote results to {args.output}")
    return 0


def extract_function_by_taxa_main() -> int:
    parser = argparse.ArgumentParser(description="Estimate the share of gene-read support for a function associated with a taxon.")
    parser.add_argument("-d", "--directory", required=True)
    parser.add_argument("-t", "--taxon", required=True)
    parser.add_argument("-f", "--function", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-c", "--category", default=None)
    args = parser.parse_args()
    table = _find_gene_table(args.directory)
    categories = [args.category] if args.category else list(FEATURE_COLUMNS)
    taxon_support: dict[str, int] = {}
    all_support: dict[str, int] = {}
    gene_hits: dict[str, int] = {}
    with open_text(table) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample = row.get("Sample", "")
            if not any(args.function in [x.strip() for x in row.get(cat, "").replace(",", "|").split("|") if x.strip()] for cat in categories):
                continue
            support = int(float(row.get("Gene_Read_Support", 0) or 0))
            all_support[sample] = all_support.get(sample, 0) + support
            gene_hits[sample] = gene_hits.get(sample, 0) + 1
            if args.taxon.lower() in (row.get("Lineage", "") + "|" + row.get("Genus", "")).lower():
                taxon_support[sample] = taxon_support.get(sample, 0) + support
    samples = sorted(set(all_support) | set(taxon_support))
    rows = [{
        "Sample": sample,
        "Taxon_Query": args.taxon,
        "Function_ID": args.function,
        "Category": args.category or "any",
        "Gene_Read_Support_In_Taxon": taxon_support.get(sample, 0),
        "Gene_Read_Support_All_Taxa": all_support.get(sample, 0),
        "Percent_of_function_support_in_taxon": round(100 * taxon_support.get(sample, 0) / all_support[sample], 6) if all_support.get(sample, 0) else 0,
        "Matched_Genes_All_Taxa": gene_hits.get(sample, 0),
    } for sample in samples]
    write_tsv(Path(args.output), rows)
    print(f"Wrote results to {args.output}")
    return 0


def extract_reads_by_function_main() -> int:
    parser = argparse.ArgumentParser(description="Build read-weighted function and taxon summaries from the current sample tree.")
    parser.add_argument("-d", "--directory", required=True)
    parser.add_argument("-p", "--prefix", default="3,4,N,P")
    parser.add_argument("-t", "--ex_taxon", default="", help="Retained for CLI compatibility; no taxa are excluded by default in the integrated analysis.")
    parser.add_argument("-m", "--minlen", type=int, default=2500)
    parser.add_argument("-o", "--outdir", default=None)
    args = parser.parse_args()
    result = run_analysis(args.directory, args.prefix, args.outdir, min_contig_length=args.minlen, heal_flags=True, force=False)
    print(f"Integrated function-by-read-support analysis written to {result['output_dir']}/functional")
    if args.ex_taxon:
        print("NOTE: --ex_taxon is accepted for compatibility but does not filter the new default output. Use the taxon query command to create a filtered view.")
    return 0


def function_by_taxa_main() -> int:
    parser = argparse.ArgumentParser(description="Build function-by-taxon tables from the sample output directory tree.")
    parser.add_argument("directory", nargs="?")
    parser.add_argument("-d", "--root-dir", dest="root_dir", default=None)
    parser.add_argument("-p", "--prefix", default="PN")
    parser.add_argument("-o", "--output-dir", default=None)
    args = parser.parse_args()
    root = args.root_dir or args.directory
    if not root:
        parser.error("provide the sample root as a positional argument or with -d")
    result = run_analysis(root, args.prefix, args.output_dir, heal_flags=True, force=False)
    print(f"Function-by-taxon results written to {result['output_dir']}/functional/taxon_function_abundance.tsv")
    return 0


def lineage_main() -> int:
    parser = argparse.ArgumentParser(description="Report contig-associated taxonomic read support across samples.")
    parser.add_argument("-d", "--dir_path", required=True)
    parser.add_argument("-p", "--prefix", default="3,4,N,P")
    parser.add_argument("-o", "--output", default=None)
    parser.add_argument("-s", "--separate-taxa", default=None)
    parser.add_argument("-r", "--remove-taxa", default=None)
    args = parser.parse_args()
    source = Path(args.dir_path).expanduser().resolve()
    result_file = source / "taxonomy" / "taxonomy_abundance.tsv"
    if not result_file.is_file():
        result = run_analysis(source, args.prefix, None, heal_flags=True, force=False)
        result_file = Path(result["output_dir"]) / "taxonomy" / "taxonomy_abundance.tsv"
    rows = _read_tsv_rows(result_file)
    if args.separate_taxa:
        separate = [x.strip() for x in args.separate_taxa.split(",") if x.strip()]
        rows = [r for r in rows if any(t in (r.get("Taxon", "") + r.get("Rank", "")) for t in separate)]
    if args.remove_taxa:
        remove = [x.strip() for x in args.remove_taxa.split(",") if x.strip()]
        rows = [r for r in rows if not any(t in r.get("Taxon", "") for t in remove)]
    out = Path(args.output) if args.output else source / "taxonomy_lineage_report.tsv"
    write_tsv(out, rows, ["Sample", "Rank", "Taxon", "Taxon_Lineage", "Contig_Count", "Mapped_Reads", "Percent_of_sample_mapped_reads"])
    print(f"Taxonomic report written to {out}")
    return 0
