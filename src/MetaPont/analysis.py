"""High-level analyses over MetaPont's integrated per-gene and per-contig tables."""

from __future__ import annotations

import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Any

from .constants import FEATURE_COLUMNS
from .layout import discover_samples
from .parsers import open_text, parse_contig_summary
from .pipeline import write_tsv


def _read_gene_table(path: str | Path):
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row


def coverage_summary(root_dir: str | Path, prefixes: str, output: str | Path | None = None,
                     read_length: int | None = None) -> tuple[Path, int]:
    root = Path(root_dir).expanduser().resolve()
    out = Path(output).expanduser().resolve() if output else root / "MetaPont_Results" / "summary" / "contig_coverage_summary.tsv"
    samples = discover_samples(root, prefixes)
    rows: list[dict[str, Any]] = []
    for sample in samples:
        if sample.contig_summary is None:
            rows.append({"Sample": sample.name, "Status": "missing_contig_summary"})
            continue
        n_contigs = total_mapped = total_reported = 0
        mapped_counts: list[int] = []
        rpk_values: list[float] = []
        proportions: list[float] = []
        estimated_coverage: list[float] = []
        min_length = None
        max_length = 0
        total_bp = 0
        for rec in parse_contig_summary(sample.contig_summary):
            n_contigs += 1
            length = rec["length"]
            mapped = rec["mapped"]
            total = rec["total"]
            total_bp += length
            min_length = length if min_length is None else min(min_length, length)
            max_length = max(max_length, length)
            total_mapped += mapped
            total_reported = max(total_reported, total)
            mapped_counts.append(mapped)
            proportions.append(float(rec.get("proportion") or 0.0))
            if length > 0:
                rpk_values.append(mapped / (length / 1000.0))
                if read_length:
                    estimated_coverage.append((mapped * read_length) / length)
        rows.append({
            "Sample": sample.name,
            "Status": "ok" if n_contigs else "empty_contig_summary",
            "Contig_Summary": str(sample.contig_summary),
            "Contig_Count": n_contigs,
            "Total_Contig_BP": total_bp,
            "Mean_Contig_Length_BP": round(total_bp / n_contigs, 3) if n_contigs else 0,
            "Min_Contig_Length_BP": min_length or 0,
            "Max_Contig_Length_BP": max_length,
            "Total_Mapped_Read_Counts": total_mapped,
            "Reported_Total_Reads": total_reported,
            "Mean_Mapped_Reads_Per_Contig": round(mean(mapped_counts), 4) if mapped_counts else 0,
            "Median_Mapped_Reads_Per_Contig": median(mapped_counts) if mapped_counts else 0,
            "Min_Mapped_Reads_Per_Contig": min(mapped_counts) if mapped_counts else 0,
            "Max_Mapped_Reads_Per_Contig": max(mapped_counts) if mapped_counts else 0,
            "Mean_Reads_Per_Kb": round(mean(rpk_values), 6) if rpk_values else "",
            "Mean_Mapping_Proportion": round(mean(proportions), 8) if proportions else "",
            "Mean_Estimated_Coverage": round(mean(estimated_coverage), 6) if estimated_coverage else "",
            "Read_Length_Used": read_length or "",
        })
    out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Sample", "Status", "Contig_Summary", "Contig_Count", "Total_Contig_BP", "Mean_Contig_Length_BP",
        "Min_Contig_Length_BP", "Max_Contig_Length_BP", "Total_Mapped_Read_Counts", "Reported_Total_Reads",
        "Mean_Mapped_Reads_Per_Contig", "Median_Mapped_Reads_Per_Contig", "Min_Mapped_Reads_Per_Contig",
        "Max_Mapped_Reads_Per_Contig", "Mean_Reads_Per_Kb", "Mean_Mapping_Proportion",
        "Mean_Estimated_Coverage", "Read_Length_Used",
    ]
    write_tsv(out, rows, fieldnames)
    return out, len(rows)


def _feature_values(cell: str) -> list[str]:
    return [part.strip() for part in cell.replace(",", "|").replace(";", "|").split("|") if part.strip() and part.strip() != "-"]


def query_function(gene_table: str | Path, function_id: str, output: str | Path,
                   category: str | None = None, rank: str = "Genus", top_n: int | None = None) -> int:
    """Rank taxa contributing gene-read support to a functional ID."""
    path = Path(gene_table)
    if not path.is_file():
        raise FileNotFoundError(f"Gene annotation table not found: {path}")
    cat_norm = category.upper().replace("-", "_") if category else None
    cat_aliases = {k.upper(): k for k in FEATURE_COLUMNS}
    if cat_norm and cat_norm not in cat_aliases:
        raise ValueError(f"Unknown category {category!r}. Valid categories: {', '.join(FEATURE_COLUMNS)}")
    chosen_cat = cat_aliases.get(cat_norm) if cat_norm else None
    rank_key = rank.title()
    if rank_key not in {"Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"}:
        raise ValueError("rank must be Domain, Kingdom, Phylum, Class, Order, Family or Species")
    totals: dict[tuple[str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    denominator: dict[str, int] = defaultdict(int)
    categories = [chosen_cat] if chosen_cat else list(FEATURE_COLUMNS)
    for row in _read_gene_table(path):
        found = any(function_id in _feature_values(row.get(cat, "")) for cat in categories)
        if not found:
            continue
        sample = row.get("Sample", "")
        lineage = row.get("Lineage", "")
        support = int(float(row.get("Gene_Read_Support", 0) or 0))
        taxon = row.get(rank_key, "unknown") or "unknown"
        lineage = row.get("Lineage", "") or ""
        if taxon.lower() == "unknown":
            continue
        totals[(sample, taxon, lineage)][0] += 1
        totals[(sample, taxon, lineage)][1] += support
        denominator[sample] += support
    rows = []
    for (sample, taxon, lineage), (genes, support) in totals.items():
        rows.append({
            "Sample": sample, "Taxonomic_Rank": rank_key, "Taxon": taxon, "Taxon_Lineage": lineage,
            "Function_ID": function_id, "Category": chosen_cat or "any", "Annotated_Genes": genes,
            "Gene_Read_Support": support,
            "Percent_of_function_gene_read_support": round(100 * support / denominator[sample], 6) if denominator[sample] else 0,
        })
    rows.sort(key=lambda r: (r["Sample"], -r["Gene_Read_Support"], r["Taxon"], r["Taxon_Lineage"]))
    if top_n:
        kept = []
        by_sample: dict[str, int] = Counter()
        for row in rows:
            if by_sample[row["Sample"]] < top_n:
                kept.append(row)
                by_sample[row["Sample"]] += 1
        rows = kept
    out = Path(output)
    write_tsv(out, rows, ["Sample", "Taxonomic_Rank", "Taxon", "Taxon_Lineage", "Function_ID", "Category", "Annotated_Genes", "Gene_Read_Support", "Percent_of_function_gene_read_support"])
    return len(rows)


def query_taxon(gene_table: str | Path, taxon_query: str, output: str | Path,
                category: str, top_n: int = 20, rank: str = "Genus") -> int:
    """Return top functional annotations linked to a taxon at a selected rank."""
    path = Path(gene_table)
    if not path.is_file():
        raise FileNotFoundError(f"Gene annotation table not found: {path}")
    category_key = category.upper().replace("-", "_")
    aliases = {key.upper(): key for key in FEATURE_COLUMNS}
    if category_key not in aliases:
        raise ValueError(f"Unknown category {category!r}. Valid categories: {', '.join(FEATURE_COLUMNS)}")
    category_name = aliases[category_key]
    rank_key = rank.title()
    if rank_key not in {"Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"}:
        raise ValueError("Invalid taxonomic rank")
    counts: dict[tuple[str, str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    for row in _read_gene_table(path):
        value = row.get(rank_key, "") or ""
        lineage = row.get("Lineage", "") or ""
        if taxon_query.lower() not in value.lower() and taxon_query.lower() not in lineage.lower():
            continue
        sample = row.get("Sample", "")
        support = int(float(row.get("Gene_Read_Support", 0) or 0))
        for feature in _feature_values(row.get(category_name, "")):
            counts[(sample, value, lineage, feature)][0] += 1
            counts[(sample, value, lineage, feature)][1] += support
    per_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for (sample, value, lineage, feature), (genes, support) in counts.items():
        per_sample[sample].append({"Sample": sample, "Taxon_Query": taxon_query, "Taxonomic_Rank": rank_key,
                                   "Taxon": value, "Taxon_Lineage": lineage,
                                   "Category": category_name, "Feature": feature,
                                   "Annotated_Genes": genes, "Gene_Read_Support": support})
    rows = []
    for sample in sorted(per_sample):
        rows.extend(sorted(per_sample[sample], key=lambda r: (-r["Gene_Read_Support"], -r["Annotated_Genes"], r["Feature"]))[:top_n])
    out = Path(output)
    write_tsv(out, rows, ["Sample", "Taxon_Query", "Taxonomic_Rank", "Taxon", "Taxon_Lineage", "Category", "Feature", "Annotated_Genes", "Gene_Read_Support"])
    return len(rows)
