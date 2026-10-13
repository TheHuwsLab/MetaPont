"""Integrated inventory, join and abundance analysis for MetaPont workflows."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TextIO

from . import __version__
from .constants import CONTIG_TABLE_COLUMNS, FEATURE_COLUMNS, LINEAGE_RANKS, UNKNOWN
from .layout import SamplePaths, discover_samples
from .parsers import (
    empty_lineage, gene_to_contig, length_from_contig_name, lineage_string,
    parse_cds_summary, parse_contig_summary, parse_fasta, parse_kraken_classifications,
    parse_kraken_report, read_eggnog_annotations,
)

OUTPUT_DIRNAME = "MetaPont_Results"


def _nonempty(path: Path | None) -> bool:
    return bool(path and path.is_file() and path.stat().st_size > 0)


def _touch_flag(path: Path, log: list[str], label: str) -> bool:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(datetime.now().astimezone().isoformat(timespec="seconds") + "\n", encoding="utf-8")
        log.append(f"HEALED\t{label}\t{path}")
        return True
    return False


def _stage_outputs(sample: SamplePaths) -> dict[str, bool]:
    """Assess output completeness and freshness before a missing flag is healed."""
    pyro_paths = (sample.pyrodigal_gff, sample.pyrodigal_nt, sample.pyrodigal_aa, sample.pyrodigal_qc)
    pyro = bool(sample.contigs and _nonempty(sample.contigs) and all(_nonempty(p) for p in pyro_paths))
    if pyro and sample.contigs:
        source_mtime = sample.contigs.stat().st_mtime_ns
        pyro = all(p.stat().st_mtime_ns >= source_mtime for p in pyro_paths if p is not None)

    eggnog_paths = (
        sample.eggnog_annotations, sample.eggnog_hits, sample.eggnog_seed_orthologs,
        sample.eggnog_decorated_gff, sample.eggnog_xlsx, sample.eggnog_qc,
    )
    eggnog = bool(sample.pyrodigal_aa and _nonempty(sample.pyrodigal_aa) and all(_nonempty(p) for p in eggnog_paths))
    if eggnog and sample.pyrodigal_aa:
        source_mtime = sample.pyrodigal_aa.stat().st_mtime_ns
        eggnog = all(p.stat().st_mtime_ns >= source_mtime for p in eggnog_paths if p is not None)

    kraken_paths = (sample.kraken_contigs, sample.kraken_contig_report, sample.kraken_mpa_report)
    kraken = bool(sample.contigs and all(_nonempty(p) for p in kraken_paths))
    if kraken and sample.contigs:
        source_mtime = sample.contigs.stat().st_mtime_ns
        kraken = all(p.stat().st_mtime_ns >= source_mtime for p in kraken_paths if p is not None)

    readmap_paths = (sample.contig_summary, sample.cds_summary)
    readmap = all(_nonempty(p) for p in (*readmap_paths, sample.bam, sample.bam_index))
    if readmap and sample.bam:
        source_mtime = sample.bam.stat().st_mtime_ns
        readmap = all(p.stat().st_mtime_ns >= source_mtime for p in readmap_paths if p is not None)

    assembly = bool(sample.contigs and _nonempty(sample.contigs) and _nonempty(sample.assembly_qc))
    if assembly and sample.contigs and sample.assembly_qc:
        assembly = sample.assembly_qc.stat().st_mtime_ns >= sample.contigs.stat().st_mtime_ns

    return {
        "assembly": assembly,
        "pyrodigal": pyro,
        "eggnog": eggnog,
        "functional": pyro and eggnog,
        "kraken_contigs": kraken,
        "kraken_reads": bool(_nonempty(sample.kraken_reads) and _nonempty(sample.kraken_reads_report)),
        "readmapping": readmap,
    }



def heal_sample_flags(sample: SamplePaths) -> list[str]:
    """Recreate stage flags only when the corresponding output set is complete."""
    events: list[str] = []
    flags = sample.flag_paths
    stages = _stage_outputs(sample)
    for stage in ("assembly", "pyrodigal", "eggnog", "functional", "kraken_contigs", "readmapping"):
        if not stages[stage]:
            continue
        _touch_flag(flags[stage], events, stage)
    return events


def inventory_rows(samples: list[SamplePaths], heal_flags: bool = False) -> tuple[list[dict[str, Any]], list[str]]:
    rows: list[dict[str, Any]] = []
    events: list[str] = []
    for sample in samples:
        if heal_flags:
            events.extend(heal_sample_flags(sample))
        stage_outputs = _stage_outputs(sample)
        flags = sample.flag_paths
        has_read1 = any(sample.sample_dir.glob("*_R1_001.fastq.gz")) or any(sample.sample_dir.glob("*_R1_*.fastq.gz"))
        has_read2 = any(sample.sample_dir.glob("*_R2_001.fastq.gz")) or any(sample.sample_dir.glob("*_R2_*.fastq.gz"))
        warnings: list[str] = []
        if not stage_outputs["assembly"]:
            warnings.append("assembly_outputs_missing_or_stale")
        if not stage_outputs["pyrodigal"]:
            warnings.append("pyrodigal_outputs_missing_or_stale")
        if not stage_outputs["eggnog"]:
            warnings.append("eggnog_outputs_missing_or_stale")
        if not stage_outputs["kraken_contigs"]:
            warnings.append("kraken_contig_outputs_missing_or_stale")
        if not stage_outputs["readmapping"]:
            warnings.append("readmapping_outputs_missing_or_stale")
        if not _nonempty(sample.contig_summary):
            warnings.append("readmapping_contig_summary_missing")
        if not _nonempty(sample.cds_summary):
            warnings.append("readmapping_cds_summary_missing")
        row = {
            "Sample": sample.name,
            "Sample_Directory": str(sample.sample_dir),
            "Raw_R1_Present": int(has_read1),
            "Raw_R2_Present": int(has_read2),
            "Fastp_Flag": int(flags["fastp"].is_file()),
            "Decontamination_Flag": int(flags["decontamination"].is_file()),
            "Assembly_Complete": int(stage_outputs["assembly"]),
            "Assembly_Flag": int(flags["assembly"].is_file()),
            "Pyrodigal_Complete": int(stage_outputs["pyrodigal"]),
            "Pyrodigal_Flag": int(flags["pyrodigal"].is_file()),
            "eggNOG_Complete": int(stage_outputs["eggnog"]),
            "eggNOG_Flag": int(flags["eggnog"].is_file()),
            "Readmapping_Complete": int(stage_outputs["readmapping"]),
            "Readmapping_Flag": int(flags["readmapping"].is_file()),
            "Kraken2_Contigs_Complete": int(stage_outputs["kraken_contigs"]),
            "Kraken2_Flag": int(flags["kraken_contigs"].is_file()),
            "Kraken2_Reads_Complete": int(stage_outputs["kraken_reads"]),
            "Functional_Complete": int(stage_outputs["functional"]),
            "Functional_Flag": int(flags["functional"].is_file()),
            "Contig_Summary": str(sample.contig_summary or ""),
            "CDS_Summary": str(sample.cds_summary or ""),
            "Kraken_Contig_Output": str(sample.kraken_contigs or ""),
            "Kraken_Report": str(sample.kraken_contig_report or ""),
            "eggNOG_Annotations": str(sample.eggnog_annotations or ""),
            "Status": "ready" if not warnings else "partial",
            "Warnings": ";".join(warnings),
        }
        rows.append(row)
    return rows, events


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_inventory(root_dir: str | Path, prefixes: str, output_dir: str | Path | None = None,
                    heal_flags: bool = False) -> tuple[Path, list[dict[str, Any]], list[str]]:
    root = Path(root_dir).expanduser().resolve()
    output_root = Path(output_dir).expanduser().resolve() if output_dir else root / OUTPUT_DIRNAME
    inventory_dir = output_root / "inventory"
    inventory_dir.mkdir(parents=True, exist_ok=True)
    samples = discover_samples(root, prefixes)
    rows, events = inventory_rows(samples, heal_flags=heal_flags)
    out = inventory_dir / "sample_inventory.tsv"
    fieldnames = list(rows[0].keys()) if rows else ["Sample", "Sample_Directory", "Status", "Warnings"]
    write_tsv(out, rows, fieldnames)
    if events:
        with (inventory_dir / "flag_healing.log").open("a", encoding="utf-8") as handle:
            handle.write("\n".join(events) + "\n")
    return out, rows, events


def _lineage_taxon(lineage: dict[str, str], rank: str) -> str:
    value = (lineage or {}).get(rank, UNKNOWN) or UNKNOWN
    return value


def _iter_base_contigs(sample: SamplePaths, assignments: dict[str, dict[str, Any]], annotated: dict[str, dict[str, Any]]):
    """Yield base contigs from the mapping summary, or stream FASTA as fallback."""
    if sample.contig_summary and sample.contig_summary.is_file():
        yield from parse_contig_summary(sample.contig_summary)
    elif sample.contigs and sample.contigs.is_file():
        for record in parse_fasta(sample.contigs):
            yield {"name": record["name"], "length": record["length"], "mapped": 0, "total": 0, "proportion": 0.0}
    else:
        for contig in sorted(set(assignments) | set(annotated)):
            info = assignments.get(contig, {})
            yield {
                "name": contig,
                "length": info.get("length") or length_from_contig_name(contig),
                "mapped": 0, "total": 0, "proportion": 0.0,
            }


def _annotation_maps(sample: SamplePaths, contig_ids: set[str], assignments: dict[str, dict[str, Any]], gene_support: dict[str, int], gene_writer: csv.DictWriter, include_taxon_functions: bool = False):
    by_contig: dict[str, dict[str, Any]] = {}
    gene_count = 0
    functional: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    taxon_function: dict[tuple[str, str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    for ann in read_eggnog_annotations(sample.eggnog_annotations):
        gene = ann["query"]
        contig = gene_to_contig(gene)
        # If known contig IDs show a more exact parent relation, prefer that match.
        if contig not in contig_ids and contig not in assignments and contig not in by_contig:
            # EggNOG query IDs normally end in _<gene-number>; the suffix remover covers the
            # Pyrodigal naming convention. Keep this fallback explicit for unusual IDs.
            candidate = gene.rsplit("_", 1)[0]
            if candidate in contig_ids or candidate in assignments:
                contig = candidate
        features: dict[str, list[str]] = ann["features"]
        support = gene_support.get(gene, 0)
        all_features = [feature for values in features.values() for feature in values]
        annotated_gene = int(bool(all_features))
        record = by_contig.setdefault(contig, {
            "gene_count": 0, "annotated_gene_count": 0, "gene_read_support": 0,
            "features": defaultdict(set),
        })
        record["gene_count"] += 1
        record["annotated_gene_count"] += annotated_gene
        if annotated_gene:
            record["gene_read_support"] += support
        for category, values in features.items():
            for feature in values:
                record["features"][category].add(feature)
                functional[(category, feature)][0] += 1
                functional[(category, feature)][1] += support
                if include_taxon_functions:
                    lineage = assignments.get(contig, {}).get("lineage", empty_lineage())
                    for rank in LINEAGE_RANKS:
                        taxon = _lineage_taxon(lineage, rank)
                        if taxon != UNKNOWN:
                            lineage_value = lineage_string(lineage)
                            taxon_function[(rank, taxon, lineage_value, f"{category}\t{feature}")][0] += 1
                            taxon_function[(rank, taxon, lineage_value, f"{category}\t{feature}")][1] += support
        tax = assignments.get(contig, {})
        gene_count += 1
        gene_lineage = tax.get("lineage", empty_lineage())
        gene_writer.writerow({
            "Sample": sample.name,
            "Gene": gene,
            "Contig": contig,
            "Gene_Read_Support": support,
            "TaxID": tax.get("taxid", UNKNOWN),
            "Lineage": tax.get("lineage_string", "|".join(f"{p}{UNKNOWN}" for p in ("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"))),
            **{rank: gene_lineage.get(rank, UNKNOWN) for rank in LINEAGE_RANKS},
            "Description": ann.get("description", ""),
            "Preferred_Name": ann.get("preferred_name", ""),
            **{category: "|".join(values) for category, values in features.items()},
        })
    return by_contig, gene_count, functional, taxon_function


def _contig_row(sample_name: str, base: dict[str, Any], tax: dict[str, Any] | None, ann: dict[str, Any] | None) -> dict[str, Any]:
    tax = tax or {}
    lineage = tax.get("lineage", empty_lineage())
    ann = ann or {}
    features = ann.get("features", {})
    out: dict[str, Any] = {
        "Sample": sample_name,
        "Contig": base["name"],
        "Contig_Length": base.get("length", 0),
        "Mapped_Reads": base.get("mapped", 0),
        "Total_Reads": base.get("total", 0),
        "Mapping_Proportion": base.get("proportion", 0.0),
        "Kraken_Classification": tax.get("classification", ""),
        "TaxID": tax.get("taxid", UNKNOWN),
        "Scientific_Name": tax.get("scientific_name", UNKNOWN),
        **{rank: lineage.get(rank, UNKNOWN) for rank in LINEAGE_RANKS},
        "Lineage": tax.get("lineage_string", lineage_string(empty_lineage())),
        "Gene_Count": ann.get("gene_count", 0),
        "Annotated_Gene_Count": ann.get("annotated_gene_count", 0),
        "Gene_Read_Support": ann.get("gene_read_support", 0),
    }
    for category in FEATURE_COLUMNS:
        out[category] = "|".join(sorted(features.get(category, set())))
    return out


def _input_signature(samples: list[SamplePaths], min_contig_length: int,
                    include_taxon_function_abundance: bool) -> dict[str, Any]:
    records = []
    for sample in samples:
        selected = [
            sample.contigs, sample.contig_summary, sample.cds_summary, sample.kraken_contigs,
            sample.kraken_contig_report, sample.kraken_mpa_report, sample.eggnog_annotations,
            sample.pyrodigal_gff, sample.pyrodigal_aa,
        ]
        for path in selected:
            if path and path.is_file():
                st = path.stat()
                records.append({"sample": sample.name, "path": str(path), "size": st.st_size, "mtime_ns": st.st_mtime_ns})
    parameters = {
        "min_contig_length": int(min_contig_length),
        "include_taxon_function_abundance": include_taxon_function_abundance,
    }
    payload = {"inputs": records, "parameters": parameters}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return {"version": __version__, "signature": digest, "parameters": parameters, "inputs": records}


def _valid_cached_outputs(output_root: Path, signature: dict[str, Any], force: bool,
                          samples: list[SamplePaths], include_taxon_function_abundance: bool) -> bool:
    manifest = output_root / "metadata" / "input_signature.json"
    required = [
        output_root / "combined" / "all_contigs.tsv.gz",
        output_root / "combined" / "gene_annotations.tsv.gz",
        output_root / "summary" / "sample_summary.tsv",
        output_root / "summary" / "sample_inventory.tsv",
        output_root / "taxonomy" / "taxonomy_abundance.tsv",
        output_root / "functional" / "functional_abundance.tsv",
        output_root / "summary" / "pipeline_warnings.tsv",
        output_root / "metadata" / "run_manifest.json",
    ]
    if include_taxon_function_abundance:
        required.append(
            output_root / "functional" / "taxon_function_abundance.tsv"
        )
    required.extend(output_root / "per_sample_contigs" / f"{sample.name}_Contigs.tsv.gz" for sample in samples)
    for category in FEATURE_COLUMNS:
        required.append(output_root / "functional" / "matrices" / f"{category}_gene_counts.tsv")
        required.append(output_root / "functional" / "matrices" / f"{category}_gene_read_support.tsv")
    if force or not manifest.is_file() or not all(_nonempty(p) for p in required):
        return False
    try:
        old = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return old.get("signature") == signature.get("signature")


def run_analysis(root_dir: str | Path, prefixes: str, output_dir: str | Path | None = None,
                 min_contig_length: int = 2500, heal_flags: bool = True, force: bool = False, include_taxon_function_abundance: bool = True) -> dict[str, Any]:
    """Build the integrated MetaPont result set; partial samples are included with warnings."""
    root = Path(root_dir).expanduser().resolve()
    output_root = Path(output_dir).expanduser().resolve() if output_dir else root / OUTPUT_DIRNAME
    dirs = {name: output_root / name for name in ("inventory", "per_sample_contigs", "combined", "summary", "taxonomy", "functional", "metadata")}
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)

    samples = discover_samples(root, prefixes)
    if not samples:
        raise RuntimeError(f"No sample directories found under {root} for prefixes {prefixes!r}.")
    healed_events: list[str] = []
    if heal_flags:
        for sample in samples:
            healed_events.extend(heal_sample_flags(sample))
    inv_rows, _ = inventory_rows(samples, heal_flags=False)
    write_tsv(dirs["inventory"] / "sample_inventory.tsv", inv_rows)
    if healed_events:
        with (dirs["inventory"] / "flag_healing.log").open("a", encoding="utf-8") as handle:
            handle.write("\n".join(healed_events) + "\n")

    signature = _input_signature(samples, min_contig_length, include_taxon_function_abundance)
    if _valid_cached_outputs(output_root, signature, force, samples, include_taxon_function_abundance):
        warnings_path = dirs["summary"] / "pipeline_warnings.tsv"
        has_warnings = _nonempty(warnings_path) and sum(1 for _ in warnings_path.open(encoding="utf-8")) > 1
        complete_flag = output_root / "MetaPont_analysis.success"
        warning_flag = output_root / "completed_with_warnings.flag"
        if has_warnings:
            complete_flag.unlink(missing_ok=True)
            _touch_flag(warning_flag, [], "analysis with warnings")
        else:
            warning_flag.unlink(missing_ok=True)
            _touch_flag(complete_flag, [], "analysis")
        return {"output_dir": str(output_root), "samples": len(samples), "cached": True, "warnings": int(has_warnings)}

    # A new/forced run must not leave a success marker from an older analysis visible.
    (output_root / "MetaPont_analysis.success").unlink(missing_ok=True)
    (output_root / "completed_with_warnings.flag").unlink(missing_ok=True)

    contigs_path = dirs["combined"] / "all_contigs.tsv.gz"
    genes_path = dirs["combined"] / "gene_annotations.tsv.gz"
    sample_summaries: list[dict[str, Any]] = []
    warning_rows: list[dict[str, Any]] = []
    functional: dict[tuple[str, str], list[Any]] = defaultdict(lambda: [0, 0, set()])
    per_sample_functional: dict[tuple[str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    taxon_functions: dict[tuple[str, str, str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    taxonomy: dict[tuple[str, str, str, str], list[int]] = defaultdict(lambda: [0, 0])

    with gzip.open(contigs_path, "wt", encoding="utf-8") as all_contigs, gzip.open(genes_path, "wt", encoding="utf-8") as all_genes:
        contig_writer = csv.DictWriter(all_contigs, fieldnames=CONTIG_TABLE_COLUMNS, delimiter="\t", extrasaction="ignore")
        gene_fields = ["Sample", "Gene", "Contig", "Gene_Read_Support", "TaxID", "Lineage", *LINEAGE_RANKS, "Description", "Preferred_Name", *FEATURE_COLUMNS.keys()]
        gene_writer = csv.DictWriter(all_genes, fieldnames=gene_fields, delimiter="\t", extrasaction="ignore")
        contig_writer.writeheader()
        gene_writer.writeheader()

        for sample in samples:
            warnings: list[str] = []
            current_stages = _stage_outputs(sample)
            if not current_stages["assembly"]:
                warnings.append("assembly_outputs_missing_or_stale")
            if not current_stages["pyrodigal"]:
                warnings.append("pyrodigal_outputs_missing_or_stale")
            if not current_stages["eggnog"]:
                warnings.append("eggnog_outputs_missing_or_stale")
            if not current_stages["kraken_contigs"]:
                warnings.append("kraken_contig_outputs_missing_or_stale")
            if not current_stages["readmapping"]:
                warnings.append("readmapping_outputs_missing_or_stale")
            if not _nonempty(sample.contig_summary):
                warnings.append("contig_readmapping_summary_missing")
            if not _nonempty(sample.cds_summary):
                warnings.append("cds_readmapping_summary_missing")

            report_map = parse_kraken_report(sample.kraken_contig_report)
            assignments = parse_kraken_classifications(sample.kraken_contigs, report_map)
            if min_contig_length > 0:
                assignments = {contig: record for contig, record in assignments.items()
                               if not record.get("length") or int(record["length"]) >= min_contig_length}
            gene_support = parse_cds_summary(sample.cds_summary)
            known_contigs = set(assignments)
            annotation_map, gene_count, sample_functional, sample_taxon_functions = _annotation_maps(
                sample, known_contigs, assignments, gene_support, gene_writer, include_taxon_functions=include_taxon_function_abundance,

            )
            for (category, feature), counts in sample_functional.items():
                functional[(category, feature)][0] += counts[0]
                functional[(category, feature)][1] += counts[1]
                functional[(category, feature)][2].add(sample.name)
                per_sample_functional[(sample.name, category, feature)][0] += counts[0]
                per_sample_functional[(sample.name, category, feature)][1] += counts[1]
            for (rank, taxon, taxon_lineage, category_feature), counts in sample_taxon_functions.items():
                category, feature = category_feature.split("\t", 1)
                taxon_functions[(sample.name, rank, taxon, taxon_lineage, category, feature)][0] += counts[0]
                taxon_functions[(sample.name, rank, taxon, taxon_lineage, category, feature)][1] += counts[1]

            # Only assignments/annotations need a set for extra-row handling; all mapping-summary
            # contigs are streamed and are not accumulated in memory.
            remaining = set(assignments) | set(annotation_map)
            sample_summary = {
                "Sample": sample.name,
                "Sample_Directory": str(sample.sample_dir),
                "Contigs_Processed": 0,
                "Classified_Contigs": 0,
                "Annotated_Contigs": 0,
                "Genes_Annotated": gene_count,
                "Functional_Gene_Read_Support": 0,
                "Total_Mapped_Reads": 0,
                "Status": "ready",
                "Warnings": "",
            }
            per_sample_path = dirs["per_sample_contigs"] / f"{sample.name}_Contigs.tsv.gz"
            with gzip.open(per_sample_path, "wt", encoding="utf-8") as per_sample_handle:
                sample_writer = csv.DictWriter(per_sample_handle, fieldnames=CONTIG_TABLE_COLUMNS, delimiter="\t", extrasaction="ignore")
                sample_writer.writeheader()
                for base in _iter_base_contigs(sample, assignments, annotation_map):
                    contig = base["name"]
                    remaining.discard(contig)
                    tax = assignments.get(contig)
                    ann = annotation_map.get(contig)
                    row = _contig_row(sample.name, base, tax, ann)
                    sample_writer.writerow(row)
                    contig_writer.writerow(row)
                    sample_summary["Contigs_Processed"] += 1
                    mapped = int(base.get("mapped", 0) or 0)
                    sample_summary["Total_Mapped_Reads"] += mapped
                    if tax and tax.get("classification") == "C":
                        sample_summary["Classified_Contigs"] += 1
                    if ann and ann.get("annotated_gene_count", 0):
                        sample_summary["Annotated_Contigs"] += 1
                        sample_summary["Functional_Gene_Read_Support"] += ann.get("gene_read_support", 0)
                    if tax:
                        lineage = tax.get("lineage", empty_lineage())
                        for rank in LINEAGE_RANKS:
                            taxon = lineage.get(rank, UNKNOWN)
                            if taxon != UNKNOWN:
                                counts = taxonomy[(sample.name, rank, taxon, lineage_string(lineage))]
                                counts[0] += 1
                                counts[1] += mapped
                # Kraken classifications cover the contigs above the requested size threshold;
                # include relevant rows absent from the readmapping summary as well.
                for contig in sorted(remaining):
                    base = {
                        "name": contig,
                        "length": assignments.get(contig, {}).get("length") or length_from_contig_name(contig),
                        "mapped": 0, "total": 0, "proportion": 0.0,
                    }
                    tax = assignments.get(contig)
                    ann = annotation_map.get(contig)
                    row = _contig_row(sample.name, base, tax, ann)
                    sample_writer.writerow(row)
                    contig_writer.writerow(row)
                    sample_summary["Contigs_Processed"] += 1
                    if tax and tax.get("classification") == "C":
                        sample_summary["Classified_Contigs"] += 1
                    if ann and ann.get("annotated_gene_count", 0):
                        sample_summary["Annotated_Contigs"] += 1
                        sample_summary["Functional_Gene_Read_Support"] += ann.get("gene_read_support", 0)
                    if tax:
                        lineage = tax.get("lineage", empty_lineage())
                        for rank in LINEAGE_RANKS:
                            taxon = lineage.get(rank, UNKNOWN)
                            if taxon != UNKNOWN:
                                counts = taxonomy[(sample.name, rank, taxon, lineage_string(lineage))]
                                counts[0] += 1
                                counts[1] += 0

            if warnings:
                sample_summary["Status"] = "partial"
                sample_summary["Warnings"] = ";".join(warnings)
                warning_rows.append({"Sample": sample.name, "Sample_Directory": str(sample.sample_dir), "Warnings": ";".join(warnings)})
            sample_summaries.append(sample_summary)

    # Long-form output tables can be filtered directly in R, Python or a spreadsheet.
    functional_rows = []
    for (category, feature), (gene_count, support, present_samples) in sorted(functional.items()):
        functional_rows.append({
            "Category": category, "Feature": feature, "Annotated_Gene_Count": gene_count,
            "Gene_Read_Support": support, "Samples_Present": len(present_samples),
        })
    write_tsv(dirs["functional"] / "functional_abundance.tsv", functional_rows,
              ["Category", "Feature", "Annotated_Gene_Count", "Gene_Read_Support", "Samples_Present"])

    # Matrix outputs are convenient for ordination, differential analysis and spreadsheet review.
    matrix_dir = dirs["functional"] / "matrices"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    sample_names = [sample.name for sample in samples]
    for category in FEATURE_COLUMNS:
        features_for_category = sorted({feature for (_, cat, feature) in per_sample_functional if cat == category})
        for metric_idx, metric_suffix in ((0, "gene_counts"), (1, "gene_read_support")):
            matrix_rows = []
            for feature in features_for_category:
                row = {"Feature": feature}
                for sample_name in sample_names:
                    row[sample_name] = per_sample_functional.get((sample_name, category, feature), [0, 0])[metric_idx]
                matrix_rows.append(row)
            write_tsv(matrix_dir / f"{category}_{metric_suffix}.tsv", matrix_rows, ["Feature", *sample_names])

    taxon_function_rows = [
        {"Sample": sample, "Taxonomic_Rank": rank, "Taxon": taxon, "Taxon_Lineage": taxon_lineage,
         "Category": category, "Feature": feature, "Annotated_Gene_Count": count, "Gene_Read_Support": support}
        for (sample, rank, taxon, taxon_lineage, category, feature), (count, support) in sorted(taxon_functions.items())
    ]
    write_tsv(dirs["functional"] / "taxon_function_abundance.tsv", taxon_function_rows,
              ["Sample", "Taxonomic_Rank", "Taxon", "Taxon_Lineage", "Category", "Feature", "Annotated_Gene_Count", "Gene_Read_Support"])

    mapped_by_sample: dict[str, int] = defaultdict(int)
    for summary in sample_summaries:
        mapped_by_sample[summary["Sample"]] = int(summary["Total_Mapped_Reads"])
    taxonomy_rows = []
    for (sample_name, rank, taxon, taxon_lineage), (contigs_count, mapped) in sorted(taxonomy.items()):
        denominator = mapped_by_sample.get(sample_name, 0)
        taxonomy_rows.append({
            "Sample": sample_name, "Rank": rank, "Taxon": taxon, "Taxon_Lineage": taxon_lineage,
            "Contig_Count": contigs_count, "Mapped_Reads": mapped,
            "Percent_of_sample_mapped_reads": round(100.0 * mapped / denominator, 6) if denominator else 0.0,
        })
    write_tsv(dirs["taxonomy"] / "taxonomy_abundance.tsv", taxonomy_rows,
              ["Sample", "Rank", "Taxon", "Taxon_Lineage", "Contig_Count", "Mapped_Reads", "Percent_of_sample_mapped_reads"])

    write_tsv(dirs["summary"] / "sample_summary.tsv", sample_summaries)
    write_tsv(dirs["summary"] / "pipeline_warnings.tsv", warning_rows,
              ["Sample", "Sample_Directory", "Warnings"])
    compact_inventory = [{k: row[k] for k in ("Sample", "Status", "Warnings", "Assembly_Complete", "Pyrodigal_Complete",
                                                "eggNOG_Complete", "Readmapping_Complete", "Kraken2_Contigs_Complete",
                                                "Functional_Complete")} for row in inv_rows]
    write_tsv(dirs["summary"] / "sample_inventory.tsv", compact_inventory)

    (dirs["metadata"] / "input_signature.json").write_text(json.dumps(signature, indent=2) + "\n", encoding="utf-8")
    run_manifest = {
        "metapont_version": __version__,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "root_dir": str(root),
        "prefixes": prefixes,
        "output_dir": str(output_root),
        "sample_count": len(samples),
        "samples_with_warnings": len(warning_rows),
        "min_contig_length_for_taxonomy": min_contig_length,
        "notes": [
            "Gene_Read_Support is the sum of read-to-CDS overlaps from the CDS summary; it is not a count of unique reads across genes.",
            "A read may support multiple genes/functional annotations, so functional support values are not mutually exclusive.",
            "Taxonomic relative abundance uses total mapped-read support as denominator; individual contig mapping counts depend on the readmapping method.",
        ],
    }
    (dirs["metadata"] / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2) + "\n", encoding="utf-8")

    clean_flag = output_root / "MetaPont_analysis.success"
    warning_flag = output_root / "completed_with_warnings.flag"
    if warning_rows:
        clean_flag.unlink(missing_ok=True)
        warning_flag.write_text(datetime.now().astimezone().isoformat(timespec="seconds") + "\n", encoding="utf-8")
    else:
        warning_flag.unlink(missing_ok=True)
        clean_flag.write_text(datetime.now().astimezone().isoformat(timespec="seconds") + "\n", encoding="utf-8")

    return {
        "output_dir": str(output_root),
        "samples": len(samples),
        "samples_with_warnings": len(warning_rows),
        "cached": False,
        "success_flag": str(warning_flag if warning_rows else clean_flag),
    }
