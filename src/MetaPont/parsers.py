"""Parsers for Kraken2, readmapping and eggNOG-mapper output files."""

from __future__ import annotations

import csv
import gzip
import re
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from .constants import EMAPPER_FALLBACK_COLUMNS, FEATURE_COLUMNS, LINEAGE_PREFIXES, LINEAGE_RANKS, UNKNOWN

csv.field_size_limit(2**31 - 1)

TAXON_NAME_RE = re.compile(r"\s*\(taxid\s+(\d+)\)\s*$", re.IGNORECASE)
GENE_SUFFIX_RE = re.compile(r"_\d+$")
LENGTH_RE = re.compile(r"(?:^|_)length_(\d+)(?:_|$)")


def open_text(path: str | Path):
    """Open plain text or gzip-compressed text based on the filename suffix."""
    p = Path(path)
    if p.suffix.lower() == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace", newline="")
    return p.open("r", encoding="utf-8", errors="replace", newline="")


def parse_int(value: Any, default: int = 0) -> int:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        try:
            return int(float(str(value).strip()))
        except (TypeError, ValueError):
            return default


def parse_float(value: Any, default: float | None = None) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return default


def parse_fasta(path: str | Path) -> Iterator[dict[str, Any]]:
    """Yield FASTA ID, description and sequence length without retaining sequences."""
    name = None
    description = ""
    length = 0
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield {"name": name, "description": description, "length": length}
                header = line[1:].strip()
                name = header.split(None, 1)[0] if header else ""
                description = header
                length = 0
            elif name is not None:
                length += len("".join(line.split()))
    if name is not None:
        yield {"name": name, "description": description, "length": length}


def parse_contig_summary(path: str | Path) -> Iterator[dict[str, Any]]:
    """Read readmapping's tabular contig summary, with conservative legacy fallback."""
    with open_text(path) as handle:
        first_data = True
        header_map: dict[str, int] | None = None
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if first_data:
                first_data = False
                normalized = [re.sub(r"[^a-z0-9]+", "_", p.strip().lower()).strip("_") for p in parts]
                if normalized and normalized[0] in {"contig_name", "contig", "reference", "name"}:
                    header_map = {name: idx for idx, name in enumerate(normalized)}
                    continue
            if header_map:
                def field(aliases: tuple[str, ...], default: str = "") -> str:
                    for alias in aliases:
                        if alias in header_map and header_map[alias] < len(parts):
                            return parts[header_map[alias]].strip()
                    return default
                name = field(("contig_name", "contig", "reference", "name"))
                length = parse_int(field(("contig_length", "length", "sequence_length")), 0)
                mapped = parse_int(field(("mapped_reads", "mapped", "mapped_count")), 0)
                total = parse_int(field(("total_reads", "reads_total", "total")), 0)
                proportion = parse_float(field(("proportion", "mapped_proportion", "fraction_mapped")), None)
            else:
                # Current readmapping format has exactly five columns:
                # name, length, mapped reads, total reads, mapped / total.
                if len(parts) < 4:
                    continue
                name = parts[0].strip()
                length = parse_int(parts[1], 0)
                mapped = parse_int(parts[2], 0)
                total = parse_int(parts[3], 0)
                proportion = parse_float(parts[4], None) if len(parts) > 4 else (mapped / total if total else 0.0)
            if not name or name.lower() in {"*", "no_feature"}:
                continue
            if proportion is None:
                proportion = mapped / total if total else 0.0
            yield {
                "name": name,
                "length": max(0, length),
                "mapped": max(0, mapped),
                "total": max(0, total),
                "proportion": proportion,
            }


def parse_cds_summary(path: str | Path | None) -> dict[str, int]:
    """Read `sort | uniq -c` CDS support summaries into gene -> overlap-read count."""
    counts: dict[str, int] = {}
    if path is None or not Path(path).is_file():
        return counts
    with open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            bits = line.split(None, 1)
            if len(bits) != 2:
                continue
            try:
                count = int(bits[0])
            except ValueError:
                # Headered alternative: Gene_ID <tab> Read_Count.
                low = line.lower()
                if "gene" in low or "read" in low:
                    continue
                continue
            gene = bits[1].strip()
            if gene.startswith("ID="):
                gene = gene[3:]
            gene = gene.split(";", 1)[0]
            counts[gene] = counts.get(gene, 0) + count
    return counts


def empty_lineage() -> dict[str, str]:
    return {rank: UNKNOWN for rank in LINEAGE_RANKS}


def lineage_string(lineage: dict[str, str]) -> str:
    return "|".join(f"{LINEAGE_PREFIXES[rank]}{lineage.get(rank, UNKNOWN) or UNKNOWN}" for rank in LINEAGE_RANKS)


def parse_kraken_report(path: str | Path | None) -> dict[str, dict[str, Any]]:
    """Map NCBI taxids to standardized lineage using report indentation and rank codes."""
    taxid_map: dict[str, dict[str, Any]] = {}
    if path is None or not Path(path).is_file():
        return taxid_map

    stack: list[tuple[int, dict[str, str]]] = []
    rank_map = {
        "D": "Domain", "K": "Kingdom", "P": "Phylum", "C": "Class",
        "O": "Order", "F": "Family", "G": "Genus", "S": "Species",
    }
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            parts = line.split("\t", 5)
            if len(parts) < 6:
                continue
            rank_code = parts[3].strip().upper()
            taxid = parts[4].strip()
            raw_name = parts[5]
            indent = len(raw_name) - len(raw_name.lstrip(" "))
            name = raw_name.strip()
            if not taxid:
                continue
            while stack and stack[-1][0] >= indent:
                stack.pop()
            lineage = stack[-1][1].copy() if stack else empty_lineage()
            if name.lower() in {"root", "unclassified"}:
                pass
            elif rank_code in rank_map:
                lineage[rank_map[rank_code]] = name
            # For subordinate ranks (e.g. G1, S1), inherit the closest primary rank.
            taxid_map[taxid] = {
                "taxid": taxid,
                "scientific_name": name,
                "lineage": lineage,
                "lineage_string": lineage_string(lineage),
                "rank_code": rank_code,
            }
            stack.append((indent, lineage.copy()))
    return taxid_map


def parse_kraken_classifications(path: str | Path | None, taxid_map: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Read standard Kraken2 per-sequence output and attach lineage from its report."""
    assignments: dict[str, dict[str, Any]] = {}
    if path is None or not Path(path).is_file():
        return assignments
    with open_text(path) as handle:
        for raw in handle:
            parts = raw.rstrip("\r\n").split("\t")
            if len(parts) < 4:
                continue
            status, sequence_id, taxon_field = parts[0].strip(), parts[1].strip(), parts[2].strip()
            length = parse_int(parts[3], 0)
            taxid_match = TAXON_NAME_RE.search(taxon_field)
            if taxid_match:
                taxid = taxid_match.group(1)
                scientific_name = taxon_field[:taxid_match.start()].strip()
            elif taxon_field.isdigit():
                taxid = taxon_field
                scientific_name = taxid_map.get(taxid, {}).get("scientific_name", UNKNOWN)
            else:
                taxid = "0" if status.upper() == "U" else ""
                scientific_name = taxon_field or UNKNOWN
            if status.upper() == "U" or taxid == "0" or taxon_field.lower() == "unclassified":
                taxid = "0"
                scientific_name = UNKNOWN
                lineage = empty_lineage()
                rank_code = "U"
            else:
                info = taxid_map.get(taxid, {})
                lineage = info.get("lineage", empty_lineage()).copy()
                rank_code = info.get("rank_code", "")
                if scientific_name in {"", UNKNOWN}:
                    scientific_name = info.get("scientific_name", UNKNOWN)
            assignments[sequence_id] = {
                "classification": status.upper(), "taxid": taxid or UNKNOWN,
                "scientific_name": scientific_name, "length": length,
                "lineage": lineage, "lineage_string": lineage_string(lineage), "rank_code": rank_code,
            }
    return assignments


def split_features(value: str | None, category: str) -> list[str]:
    """Split eggNOG functional values using the encoding convention for each category."""
    if not value or value.strip() in {"-", "", "NA", "na"}:
        return []
    text = value.strip()
    if category == "COG":
        # COG categories are represented by one-letter functional classes.
        return sorted(set(char for char in text if char.isalpha()))
    values = re.split(r"[,|;]", text)
    return sorted(set(part.strip() for part in values if part.strip() and part.strip() != "-"))


def read_eggnog_annotations(path: str | Path | None) -> Iterator[dict[str, Any]]:
    """Yield eggNOG annotation records using the file's header when available."""
    if path is None or not Path(path).is_file():
        return
    header: list[str] | None = None
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                candidate = line[1:].split("\t")
                normalized = [c.strip() for c in candidate]
                if normalized and normalized[0].lower() in {"query", "#query"} and ("seed_ortholog" in normalized or "GOs" in normalized):
                    header = normalized
                continue
            parts = line.split("\t")
            columns = header or EMAPPER_FALLBACK_COLUMNS
            if len(parts) < 1:
                continue
            row = {columns[i]: parts[i] if i < len(parts) else "" for i in range(len(columns))}
            query = row.get("query", row.get("#query", parts[0])).strip()
            if not query:
                continue
            features = {category: split_features(row.get(column), category) for category, column in FEATURE_COLUMNS.items()}
            yield {
                "query": query,
                "seed_ortholog": row.get("seed_ortholog", ""),
                "description": row.get("Description", ""),
                "preferred_name": row.get("Preferred_name", ""),
                "features": features,
                "raw": row,
            }


def gene_to_contig(gene_id: str) -> str:
    """Convert a typical Pyrodigal gene ID to its parent MetaSPAdes contig ID."""
    value = gene_id.strip()
    if value.startswith("ID="):
        value = value[3:]
    value = value.split(";", 1)[0]
    return GENE_SUFFIX_RE.sub("", value)


def length_from_contig_name(contig_id: str) -> int:
    match = LENGTH_RE.search(contig_id)
    return int(match.group(1)) if match else 0
