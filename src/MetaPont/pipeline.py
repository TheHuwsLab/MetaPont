"""Integrated inventory, join and abundance analysis for MetaPont workflows."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import sqlite3
import tempfile
import zlib
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .constants import CONTIG_TABLE_COLUMNS, FEATURE_COLUMNS, LINEAGE_RANKS, UNKNOWN
from .layout import SamplePaths, discover_samples
from .parsers import (
    empty_lineage, gene_to_contig, length_from_contig_name, lineage_string, split_features, open_text,
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
    """Assess whether the files MetaPont actually consumes are present and non-empty.

    Modification times are deliberately *not* used. Scratch-retention commands
    commonly touch every file, which destroys the chronological relationship
    between a source file and its derived outputs without changing their
    contents. MetaPont cannot infer freshness from those timestamps.
    """
    assembly = _nonempty(sample.contigs)
    # GFF and protein FASTA are the core Pyrodigal products relevant to the
    # annotation chain. NT FASTA and QC sidecars are useful but not consumed by
    # MetaPont and therefore do not gate its analysis.
    pyrodigal = _nonempty(sample.pyrodigal_gff) and _nonempty(sample.pyrodigal_aa)
    # The annotations file is the only required eggNOG product used here.
    # .hits/.seed_orthologs/.xlsx/decorated.gff are optional side products.
    eggnog = _nonempty(sample.eggnog_annotations)
    # Kraken's per-contig classifications and standard report are sufficient
    # for taxonomic lineage assignment. The MPA report is optional.
    kraken = _nonempty(sample.kraken_contigs) and _nonempty(sample.kraken_contig_report)
    # MetaPont consumes the tabular summaries, not the original BAM/index.
    readmap = _nonempty(sample.contig_summary) and _nonempty(sample.cds_summary)
    return {
        "assembly": assembly,
        "pyrodigal": pyrodigal,
        "eggnog": eggnog,
        "functional": pyrodigal and eggnog,
        "kraken_contigs": kraken,
        "kraken_reads": bool(_nonempty(sample.kraken_reads) and _nonempty(sample.kraken_reads_report)),
        "readmapping": readmap,
    }


def _sample_warning_codes(sample: SamplePaths, stages: dict[str, bool]) -> list[str]:
    """Return precise, actionable warnings for missing MetaPont input products."""
    warnings: list[str] = []
    if not stages["assembly"]:
        warnings.append("assembly_contigs_missing_or_empty")
    if not stages["pyrodigal"]:
        missing = []
        if not _nonempty(sample.pyrodigal_gff):
            missing.append("gff")
        if not _nonempty(sample.pyrodigal_aa):
            missing.append("protein_fasta")
        warnings.append("pyrodigal_core_outputs_missing:" + ",".join(missing))
    if not stages["eggnog"]:
        warnings.append("eggnog_annotations_missing_or_empty")
    if not _nonempty(sample.kraken_contigs):
        warnings.append("kraken_contig_classifications_missing_or_empty")
    if not _nonempty(sample.kraken_contig_report):
        warnings.append("kraken_contig_report_missing_or_empty")
    if not _nonempty(sample.contig_summary):
        warnings.append("readmapping_contig_summary_missing_or_empty")
    if not _nonempty(sample.cds_summary):
        warnings.append("readmapping_cds_summary_missing_or_empty")
    return warnings



def _stage_flag_complete(sample: SamplePaths, stage: str) -> bool:
    """Use the full expected stage product set before recreating upstream flags.

    This is deliberately stricter than MetaPont's own input-availability check:
    an integrated analysis may be possible with a subset of upstream products,
    but that is not sufficient evidence to mark the upstream workflow stage as
    successfully completed.
    """
    required: tuple[Path | None, ...]
    if stage == "assembly":
        required = (sample.contigs, sample.assembly_qc)
    elif stage == "pyrodigal":
        required = (sample.pyrodigal_gff, sample.pyrodigal_nt, sample.pyrodigal_aa, sample.pyrodigal_qc)
    elif stage == "eggnog":
        required = (
            sample.eggnog_annotations, sample.eggnog_hits, sample.eggnog_seed_orthologs,
            sample.eggnog_decorated_gff, sample.eggnog_xlsx, sample.eggnog_qc,
        )
    elif stage == "functional":
        return _stage_flag_complete(sample, "pyrodigal") and _stage_flag_complete(sample, "eggnog")
    elif stage == "kraken_contigs":
        required = (sample.contigs, sample.kraken_contigs, sample.kraken_contig_report, sample.kraken_mpa_report)
    elif stage == "readmapping":
        required = (sample.contig_summary, sample.cds_summary, sample.bam, sample.bam_index)
    else:
        raise ValueError(f"Unknown stage flag: {stage}")
    return all(_nonempty(path) for path in required)


def heal_sample_flags(sample: SamplePaths) -> list[str]:
    """Recreate missing stage flags only when every expected stage product exists."""
    events: list[str] = []
    flags = sample.flag_paths
    for stage in ("assembly", "pyrodigal", "eggnog", "functional", "kraken_contigs", "readmapping"):
        if _stage_flag_complete(sample, stage):
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
        warnings = _sample_warning_codes(sample, stage_outputs)
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


def _iter_base_contigs(sample: SamplePaths, assignments: dict[str, dict[str, Any]], annotated: dict[str, dict[str, Any]]):
    """Yield base contigs from the mapping summary, or stream FASTA as fallback."""
    if _nonempty(sample.contig_summary):
        yield from parse_contig_summary(sample.contig_summary)
    elif _nonempty(sample.contigs):
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


def _annotation_maps(sample: SamplePaths, contig_ids: set[str], assignments: dict[str, dict[str, Any]],
                     gene_support: dict[str, int], gene_writer: csv.DictWriter):
    by_contig: dict[str, dict[str, Any]] = {}
    gene_count = 0
    genes_with_functions = 0
    functional: dict[tuple[str, str], list[int]] = defaultdict(lambda: [0, 0])
    for ann in read_eggnog_annotations(sample.eggnog_annotations):
        gene = ann["query"]
        contig = gene_to_contig(gene)
        # If known contig IDs show a more exact parent relation, prefer that match.
        if contig not in contig_ids and contig not in assignments and contig not in by_contig:
            candidate = gene.rsplit("_", 1)[0]
            if candidate in contig_ids or candidate in assignments:
                contig = candidate
        features: dict[str, list[str]] = ann["features"]
        support = gene_support.get(gene, 0)
        all_features = [feature for values in features.values() for feature in values]
        annotated_gene = int(bool(all_features))
        genes_with_functions += annotated_gene
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
        tax = assignments.get(contig, {})
        gene_lineage = tax.get("lineage", empty_lineage())
        gene_writer.writerow({
            "Sample": sample.name,
            "Gene": gene,
            "Contig": contig,
            "Gene_Read_Support": support,
            "Has_Functional_Annotation": annotated_gene,
            "TaxID": tax.get("taxid", UNKNOWN),
            "Lineage": tax.get("lineage_string", lineage_string(empty_lineage())),
            **{rank: gene_lineage.get(rank, UNKNOWN) for rank in LINEAGE_RANKS},
            "Description": ann.get("description", ""),
            "Preferred_Name": ann.get("preferred_name", ""),
            **{category: "|".join(values) for category, values in features.items()},
        })
        gene_count += 1
    return by_contig, gene_count, genes_with_functions, functional



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


def _file_fingerprint(path: Path, mode: str = "sampled", chunk_size: int = 128 * 1024) -> str:
    """Content fingerprint insensitive to ``touch``; full mode hashes every byte.

    Sampled mode hashes seven evenly spaced 128 KiB windows, including the
    beginning and end (or the complete file when it is small). This is a
    practical change detector for multi-GB workflow inputs; use full mode when
    cryptographic whole-file verification is needed.
    """
    if mode not in {"sampled", "full"}:
        raise ValueError("signature_mode must be 'sampled' or 'full'.")
    size = path.stat().st_size
    digest = hashlib.sha256()
    digest.update(str(size).encode("ascii"))
    with path.open("rb") as handle:
        if mode == "full" or size <= chunk_size * 7:
            while True:
                block = handle.read(1024 * 1024)
                if not block:
                    break
                digest.update(block)
        else:
            max_offset = max(0, size - chunk_size)
            offsets = sorted({(max_offset * index) // 6 for index in range(7)})
            for offset in offsets:
                handle.seek(offset)
                block = handle.read(chunk_size)
                digest.update(offset.to_bytes(8, "big"))
                digest.update(len(block).to_bytes(8, "big"))
                digest.update(block)
    return digest.hexdigest()



def _input_signature(samples: list[SamplePaths], min_contig_length: int,
                     include_taxon_function_abundance: bool, signature_mode: str = "sampled") -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    selected_fields = (
        "contigs", "contig_summary", "cds_summary", "kraken_contigs",
        "kraken_contig_report", "eggnog_annotations", "pyrodigal_gff", "pyrodigal_aa",
    )
    for sample in samples:
        for field in selected_fields:
            path = getattr(sample, field)
            if path and path.is_file():
                st = path.stat()
                records.append({
                    "sample": sample.name,
                    "role": field,
                    "path": str(path),
                    "size": st.st_size,
                    # Kept for diagnostics only; deliberately excluded from the digest.
                    "mtime_ns": st.st_mtime_ns,
                    "content_fingerprint": _file_fingerprint(path, signature_mode),
                })
    parameters = {
        "min_contig_length": int(min_contig_length),
        "include_taxon_function_abundance": bool(include_taxon_function_abundance),
        "signature_mode": signature_mode,
    }
    payload_inputs = [
        {k: record[k] for k in ("sample", "role", "path", "size", "content_fingerprint")}
        for record in records
    ]
    sample_records = [
        {"sample": sample.name, "sample_dir": str(sample.sample_dir)}
        for sample in samples
    ]
    payload = {
        "version": __version__,
        "signature_algorithm": f"sha256_{signature_mode}_v1",
        "samples": sample_records,
        "inputs": payload_inputs,
        "parameters": parameters,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return {
        "version": __version__,
        "signature_algorithm": payload["signature_algorithm"],
        "signature": digest,
        "parameters": parameters,
        "samples": sample_records,
        "inputs": records,
    }


def _gzip_has_header(path: Path) -> bool:
    """Return whether a gzip TSV opens and contains a non-empty first line."""
    try:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return bool(handle.readline().strip())
    except (OSError, EOFError, UnicodeDecodeError, gzip.BadGzipFile, zlib.error):
        return False


def _valid_cached_outputs(output_root: Path, signature: dict[str, Any], force: bool,
                          samples: list[SamplePaths], include_taxon_function_abundance: bool) -> bool:
    manifest = output_root / "metadata" / "input_signature.json"
    in_progress = output_root / "metadata" / "run_in_progress.flag"
    required = [
        output_root / "combined" / "all_contigs.tsv.gz",
        output_root / "combined" / "gene_annotations.tsv.gz",
        output_root / "summary" / "sample_summary.tsv",
        output_root / "summary" / "sample_inventory.tsv",
        output_root / "taxonomy" / "taxonomy_abundance.tsv",
        output_root / "functional" / "functional_abundance.tsv",
        output_root / "summary" / "pipeline_warnings.tsv",
        output_root / "inventory" / "sample_inventory.tsv",
        output_root / "metadata" / "run_manifest.json",
    ]
    if include_taxon_function_abundance:
        required.append(output_root / "functional" / "taxon_function_abundance.tsv")
    required.extend(output_root / "per_sample_contigs" / f"{sample.name}_Contigs.tsv.gz" for sample in samples)
    for category in FEATURE_COLUMNS:
        required.append(output_root / "functional" / "matrices" / f"{category}_gene_counts.tsv")
        required.append(output_root / "functional" / "matrices" / f"{category}_gene_read_support.tsv")
    complete_marker = output_root / "MetaPont_analysis.success"
    warning_marker = output_root / "completed_with_warnings.flag"
    if (force or in_progress.exists() or not (complete_marker.is_file() or warning_marker.is_file())
            or not manifest.is_file() or not all(_nonempty(p) for p in required)):
        return False
    if not _gzip_has_header(required[0]) or not _gzip_has_header(required[1]):
        return False
    for path in required:
        if path.suffix == ".gz" and path.parent.name == "per_sample_contigs" and not _gzip_has_header(path):
            return False
    try:
        old = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return (
        old.get("version") == __version__
        and old.get("signature_algorithm") == signature.get("signature_algorithm")
        and old.get("signature") == signature.get("signature")
    )


def _atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def build_taxon_function_abundance(gene_table: str | Path, output: str | Path) -> int:
    """Aggregate a gene table into rank-specific taxon/function counts on disk.

    SQLite is used as an external aggregation store so the optional expanded
    table does not require retaining every taxon/function combination in RAM.
    """
    source = Path(gene_table).expanduser().resolve()
    target = Path(output).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Gene annotation table not found: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, db_name = tempfile.mkstemp(prefix=".metapont_taxon_function_", suffix=".sqlite3", dir=target.parent)
    os.close(fd)
    conn: sqlite3.Connection | None = None
    temporary_output = target.with_name(target.name + ".tmp")
    headers = ["Sample", "Taxonomic_Rank", "Taxon", "Taxon_Lineage", "Category", "Feature", "Annotated_Gene_Count", "Gene_Read_Support"]
    row_count = 0
    try:
        conn = sqlite3.connect(db_name)
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("PRAGMA temp_store=FILE")
        conn.execute("""CREATE TABLE counts (
            sample TEXT NOT NULL, rank TEXT NOT NULL, taxon TEXT NOT NULL,
            lineage TEXT NOT NULL, category TEXT NOT NULL, feature TEXT NOT NULL,
            gene_count INTEGER NOT NULL, read_support INTEGER NOT NULL,
            PRIMARY KEY (sample, rank, taxon, lineage, category, feature)
        ) WITHOUT ROWID""")
        upsert = """INSERT INTO counts VALUES (?, ?, ?, ?, ?, ?, 1, ?)
            ON CONFLICT(sample, rank, taxon, lineage, category, feature)
            DO UPDATE SET gene_count = gene_count + 1,
                          read_support = read_support + excluded.read_support"""
        pending: list[tuple[Any, ...]] = []
        with open_text(source) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError(f"Gene annotation table has no header: {source}")
            required_columns = {"Sample", "Gene_Read_Support", *LINEAGE_RANKS, *FEATURE_COLUMNS.keys()}
            missing = sorted(required_columns.difference(reader.fieldnames))
            if missing:
                raise ValueError("Gene annotation table is missing required columns: " + ", ".join(missing))
            for row in reader:
                sample = (row.get("Sample") or "").strip()
                if not sample:
                    continue
                try:
                    support = int(float(row.get("Gene_Read_Support") or 0))
                except (TypeError, ValueError):
                    support = 0
                lineage = {rank: (row.get(rank) or UNKNOWN).strip() or UNKNOWN for rank in LINEAGE_RANKS}
                taxon_lineage = lineage_string(lineage)
                for rank in LINEAGE_RANKS:
                    taxon = lineage[rank]
                    if taxon.lower() == UNKNOWN:
                        continue
                    for category, _column in FEATURE_COLUMNS.items():
                        for feature in split_features(row.get(category), category):
                            pending.append((sample, rank, taxon, taxon_lineage, category, feature, support))
                if len(pending) >= 20000:
                    conn.executemany(upsert, pending)
                    conn.commit()
                    pending.clear()
            if pending:
                conn.executemany(upsert, pending)
                conn.commit()
        with temporary_output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(headers)
            cursor = conn.execute("SELECT sample, rank, taxon, lineage, category, feature, gene_count, read_support FROM counts ORDER BY sample, rank, taxon, lineage, category, feature")
            for values in cursor:
                writer.writerow(values)
                row_count += 1
        os.replace(temporary_output, target)
    except sqlite3.Error as error:
        raise RuntimeError(f"SQLite taxon-function aggregation failed: {error}") from error
    finally:
        if conn is not None:
            conn.close()
        Path(db_name).unlink(missing_ok=True)
        Path(db_name + "-journal").unlink(missing_ok=True)
        temporary_output.unlink(missing_ok=True)
    return row_count


def run_analysis(root_dir: str | Path, prefixes: str, output_dir: str | Path | None = None,
                 min_contig_length: int = 2500, heal_flags: bool = True, force: bool = False,
                 include_taxon_function_abundance: bool = False,
                 signature_mode: str = "sampled") -> dict[str, Any]:
    """Build integrated MetaPont results; incomplete samples are retained with warnings.

    ``signature_mode='sampled'`` makes cache validation insensitive to ``touch``
    while checking representative content blocks. ``full`` hashes every byte.
    """
    if min_contig_length < 0:
        raise ValueError("min_contig_length must be zero or a positive integer.")
    if signature_mode not in {"sampled", "full"}:
        raise ValueError("signature_mode must be 'sampled' or 'full'.")

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

    manifest_path = dirs["metadata"] / "input_signature.json"
    taxon_function_path = dirs["functional"] / "taxon_function_abundance.tsv"
    if not include_taxon_function_abundance:
        # An optional product from an earlier opt-in run must not survive as a
        # misleading apparently-current result when this run has it disabled.
        taxon_function_path.unlink(missing_ok=True)

    signature = _input_signature(samples, min_contig_length, include_taxon_function_abundance, signature_mode)
    if _valid_cached_outputs(output_root, signature, force, samples, include_taxon_function_abundance):
        # Keep both inventory views synchronised even when the analytical tables
        # are reused from cache (e.g. a flag was healed without changing inputs).
        compact_inventory = [{k: row[k] for k in (
            "Sample", "Status", "Warnings", "Assembly_Complete", "Pyrodigal_Complete",
            "eggNOG_Complete", "Readmapping_Complete", "Kraken2_Contigs_Complete", "Functional_Complete",
        )} for row in inv_rows]
        write_tsv(dirs["summary"] / "sample_inventory.tsv", compact_inventory)
        # Refresh diagnostic mtimes in the manifest without changing the
        # content-derived signature or rewriting analytical result tables.
        _atomic_write_json(manifest_path, signature)
        warnings_path = dirs["summary"] / "pipeline_warnings.tsv"
        if _nonempty(warnings_path):
            with warnings_path.open(encoding="utf-8", newline="") as handle:
                warning_count = sum(1 for row in csv.DictReader(handle, delimiter="\t") if row.get("Sample"))
        else:
            warning_count = 0
        marker = output_root / ("completed_with_warnings.flag" if warning_count else "MetaPont_analysis.success")
        return {
            "output_dir": str(output_root), "samples": len(samples), "cached": True,
            "warnings": warning_count, "samples_with_warnings": warning_count,
            "success_flag": str(marker),
        }

    # A crash-safe marker prevents a half-written output set from being accepted
    # as a valid cache on the next invocation. Remove the old signature before
    # starting so a previous successful manifest cannot bless partial files.
    clean_flag = output_root / "MetaPont_analysis.success"
    warning_flag = output_root / "completed_with_warnings.flag"
    clean_flag.unlink(missing_ok=True)
    warning_flag.unlink(missing_ok=True)
    manifest_path.unlink(missing_ok=True)
    in_progress = dirs["metadata"] / "run_in_progress.flag"
    in_progress.write_text(datetime.now(timezone.utc).isoformat() + "\n", encoding="utf-8")

    # Remove any previous table before rebuilding; an interrupted opt-in run
    # must not leave an older table looking like the current result.
    taxon_function_path.unlink(missing_ok=True)

    contigs_path = dirs["combined"] / "all_contigs.tsv.gz"
    genes_path = dirs["combined"] / "gene_annotations.tsv.gz"
    sample_summaries: list[dict[str, Any]] = []
    warning_rows: list[dict[str, Any]] = []
    functional: dict[tuple[str, str], list[Any]] = defaultdict(lambda: [0, 0, set()])
    per_sample_functional: dict[tuple[str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    taxonomy: dict[tuple[str, str, str, str], list[int]] = defaultdict(lambda: [0, 0])

    with gzip.open(contigs_path, "wt", encoding="utf-8") as all_contigs, gzip.open(genes_path, "wt", encoding="utf-8") as all_genes:
        contig_writer = csv.DictWriter(all_contigs, fieldnames=CONTIG_TABLE_COLUMNS, delimiter="\t", extrasaction="ignore")
        gene_fields = [
            "Sample", "Gene", "Contig", "Gene_Read_Support", "Has_Functional_Annotation",
            "TaxID", "Lineage", *LINEAGE_RANKS, "Description", "Preferred_Name", *FEATURE_COLUMNS.keys(),
        ]
        gene_writer = csv.DictWriter(all_genes, fieldnames=gene_fields, delimiter="\t", extrasaction="ignore")
        contig_writer.writeheader()
        gene_writer.writeheader()

        for sample in samples:
            current_stages = _stage_outputs(sample)
            warnings = _sample_warning_codes(sample, current_stages)

            report_map = parse_kraken_report(sample.kraken_contig_report)
            raw_assignments = parse_kraken_classifications(sample.kraken_contigs, report_map)
            assignments = raw_assignments
            if min_contig_length > 0:
                assignments = {
                    contig: record for contig, record in raw_assignments.items()
                    if not record.get("length") or int(record["length"]) >= min_contig_length
                }
            gene_support = parse_cds_summary(sample.cds_summary)
            known_contigs = set(assignments)
            annotation_map, gene_count, genes_with_functions, sample_functional = _annotation_maps(
                sample, known_contigs, assignments, gene_support, gene_writer,
            )
            for (category, feature), counts in sample_functional.items():
                functional[(category, feature)][0] += counts[0]
                functional[(category, feature)][1] += counts[1]
                functional[(category, feature)][2].add(sample.name)
                per_sample_functional[(sample.name, category, feature)][0] += counts[0]
                per_sample_functional[(sample.name, category, feature)][1] += counts[1]

            # Mapping-summary contigs are streamed. Only the Kraken and annotated
            # contig IDs (usually much smaller) are retained for extra-row handling.
            remaining = set(assignments) | set(annotation_map)
            sample_summary: dict[str, Any] = {
                "Sample": sample.name,
                "Sample_Directory": str(sample.sample_dir),
                "Contigs_Processed": 0,
                "Classified_Contigs": 0,
                "Unclassified_Contigs": 0,
                "Taxonomy_Filtered_Contigs": max(0, len(raw_assignments) - len(assignments)),
                "Annotated_Contigs": 0,
                # Retained for backwards compatibility; this equals eggNOG rows
                # processed, not necessarily genes with a functional feature.
                "Genes_Annotated": gene_count,
                "Genes_Processed": gene_count,
                "Genes_With_Function": genes_with_functions,
                "Genes_Without_Function": max(0, gene_count - genes_with_functions),
                "Functional_Annotation_Rate_Percent": round(100.0 * genes_with_functions / gene_count, 4) if gene_count else 0.0,
                "Functional_Gene_Read_Support": 0,
                "Total_Mapped_Reads": 0,
                "Mapped_Reads_on_Classified_Contigs": 0,
                "Mapped_Reads_on_Unclassified_Contigs": 0,
                "Mapped_Reads_on_Kraken_Unassigned_Contigs": 0,
                "Mapped_Reads_on_Length_Filtered_Contigs": 0,
                "Percent_Mapped_Reads_on_Classified_Contigs": 0.0,
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
                    raw_tax = raw_assignments.get(contig)
                    tax = assignments.get(contig)
                    ann = annotation_map.get(contig)
                    row = _contig_row(sample.name, base, tax, ann)
                    sample_writer.writerow(row)
                    contig_writer.writerow(row)

                    sample_summary["Contigs_Processed"] += 1
                    mapped = int(base.get("mapped", 0) or 0)
                    sample_summary["Total_Mapped_Reads"] += mapped
                    if raw_tax is not None and tax is None:
                        sample_summary["Mapped_Reads_on_Length_Filtered_Contigs"] += mapped
                    elif tax and tax.get("classification") == "C":
                        sample_summary["Classified_Contigs"] += 1
                        sample_summary["Mapped_Reads_on_Classified_Contigs"] += mapped
                    elif tax and tax.get("classification") == "U":
                        sample_summary["Unclassified_Contigs"] += 1
                        sample_summary["Mapped_Reads_on_Unclassified_Contigs"] += mapped
                    else:
                        sample_summary["Mapped_Reads_on_Kraken_Unassigned_Contigs"] += mapped

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

                # Include Kraken classifications and annotated genes that are not
                # present in the readmapping contig summary, with zero read support.
                for contig in sorted(remaining):
                    tax = assignments.get(contig)
                    ann = annotation_map.get(contig)
                    base = {
                        "name": contig,
                        "length": (tax or {}).get("length") or length_from_contig_name(contig),
                        "mapped": 0, "total": 0, "proportion": 0.0,
                    }
                    row = _contig_row(sample.name, base, tax, ann)
                    sample_writer.writerow(row)
                    contig_writer.writerow(row)
                    sample_summary["Contigs_Processed"] += 1
                    if tax and tax.get("classification") == "C":
                        sample_summary["Classified_Contigs"] += 1
                    elif tax and tax.get("classification") == "U":
                        sample_summary["Unclassified_Contigs"] += 1
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

            total_mapped = int(sample_summary["Total_Mapped_Reads"])
            sample_summary["Percent_Mapped_Reads_on_Classified_Contigs"] = (
                round(100.0 * sample_summary["Mapped_Reads_on_Classified_Contigs"] / total_mapped, 6)
                if total_mapped else 0.0
            )
            if warnings:
                sample_summary["Status"] = "partial"
                sample_summary["Warnings"] = ";".join(warnings)
                warning_rows.append({
                    "Sample": sample.name,
                    "Sample_Directory": str(sample.sample_dir),
                    "Warnings": ";".join(warnings),
                })
            sample_summaries.append(sample_summary)

    # Long-form functional results and per-category matrices.
    functional_rows = []
    for (category, feature), (gene_count, support, present_samples) in sorted(functional.items()):
        functional_rows.append({
            "Category": category, "Feature": feature, "Annotated_Gene_Count": gene_count,
            "Gene_Read_Support": support, "Samples_Present": len(present_samples),
        })
    write_tsv(dirs["functional"] / "functional_abundance.tsv", functional_rows,
              ["Category", "Feature", "Annotated_Gene_Count", "Gene_Read_Support", "Samples_Present"])

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

    # Generate optional taxon/function abundance from the completed gene table.
    # This is disk-backed (SQLite), avoiding a potentially very large Python dict.
    if include_taxon_function_abundance:
        taxon_function_rows = build_taxon_function_abundance(genes_path, taxon_function_path)
    else:
        taxon_function_path.unlink(missing_ok=True)
        taxon_function_rows = 0

    mapped_by_sample = {row["Sample"]: int(row["Total_Mapped_Reads"]) for row in sample_summaries}
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
    compact_inventory = [{k: row[k] for k in (
        "Sample", "Status", "Warnings", "Assembly_Complete", "Pyrodigal_Complete",
        "eggNOG_Complete", "Readmapping_Complete", "Kraken2_Contigs_Complete", "Functional_Complete",
    )} for row in inv_rows]
    write_tsv(dirs["summary"] / "sample_inventory.tsv", compact_inventory)

    run_manifest = {
        "metapont_version": __version__,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "root_dir": str(root),
        "prefixes": prefixes,
        "output_dir": str(output_root),
        "sample_count": len(samples),
        "samples_with_warnings": len(warning_rows),
        "samples_missing_functional_annotations": sum(not _stage_outputs(sample)["eggnog"] for sample in samples),
        "samples_missing_pyrodigal_core_outputs": sum(not _stage_outputs(sample)["pyrodigal"] for sample in samples),
        "samples_missing_taxonomy_inputs": sum(not _stage_outputs(sample)["kraken_contigs"] for sample in samples),
        "min_contig_length_for_taxonomy": min_contig_length,
        "taxon_function_abundance_enabled": bool(include_taxon_function_abundance),
        "taxon_function_abundance_rows": taxon_function_rows,
        "input_signature_algorithm": signature["signature_algorithm"],
        "input_signature": signature["signature"],
        "notes": [
            "Gene_Read_Support is the sum of read-to-CDS overlaps from the CDS summary; it is not a count of unique reads across genes.",
            "A read may support multiple genes/functional annotations, so functional support values are not mutually exclusive.",
            "Taxonomic relative abundance uses total mapped-read support as denominator; individual contig mapping counts depend on the readmapping method.",
            "Stage completeness checks file presence and non-zero size, not timestamp freshness; scratch-retention touch operations make modification times unsuitable for freshness validation.",
            "The sampled input signature hashes seven evenly spaced 128 KiB windows per input (full mode hashes every byte); modification times are recorded for diagnostics but excluded from cache identity.",
            "Genes_Annotated is retained for compatibility and equals eggNOG annotation rows processed; Genes_With_Function counts rows with at least one parsed functional feature.",
        ],
    }
    _atomic_write_json(manifest_path, signature)
    _atomic_write_json(dirs["metadata"] / "run_manifest.json", run_manifest)

    if warning_rows:
        warning_flag.write_text(datetime.now().astimezone().isoformat(timespec="seconds") + "\n", encoding="utf-8")
        success_flag = warning_flag
    else:
        clean_flag.write_text(datetime.now().astimezone().isoformat(timespec="seconds") + "\n", encoding="utf-8")
        success_flag = clean_flag
    # This is the final commit marker. If the process crashes before this line,
    # the manifest cannot be used as a cache on the next run.
    in_progress.unlink(missing_ok=True)

    return {
        "output_dir": str(output_root),
        "samples": len(samples),
        "samples_with_warnings": len(warning_rows),
        "warnings": len(warning_rows),
        "taxon_function_abundance_rows": taxon_function_rows,
        "cached": False,
        "success_flag": str(success_flag),
    }

