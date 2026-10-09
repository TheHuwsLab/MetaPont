"""Discovery of the per-sample output layout used by MetaPont workflows."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(slots=True)
class SamplePaths:
    name: str
    sample_dir: Path
    contigs: Path | None
    assembly_qc: Path | None
    pyrodigal_dir: Path | None
    pyrodigal_gff: Path | None
    pyrodigal_nt: Path | None
    pyrodigal_aa: Path | None
    pyrodigal_qc: Path | None
    eggnog_dir: Path | None
    eggnog_annotations: Path | None
    eggnog_hits: Path | None
    eggnog_seed_orthologs: Path | None
    eggnog_decorated_gff: Path | None
    eggnog_xlsx: Path | None
    eggnog_qc: Path | None
    kraken_dir: Path | None
    kraken_contigs: Path | None
    kraken_contig_report: Path | None
    kraken_mpa_report: Path | None
    kraken_reads: Path | None
    kraken_reads_report: Path | None
    readmapping_dir: Path | None
    contig_summary: Path | None
    cds_summary: Path | None
    bam: Path | None
    bam_index: Path | None

    @property
    def flag_paths(self) -> dict[str, Path]:
        name, root = self.name, self.sample_dir
        return {
            "fastp": root / "fastp.succeed",
            "decontamination": root / f"{name}_decontamination_done.flag",
            "assembly": root / f"{name}_metaspades_assembly_done.flag",
            "pyrodigal": (self.pyrodigal_dir or root / f"{name}_pyrodigal") / f"{name}_pyrodigal_done.flag",
            "eggnog": (self.eggnog_dir or root / f"{name}_eggnog_mapper") / f"{name}_eggnogmapper.success",
            "kraken_contigs": root / f"{name}_kraken2.success",
            "kraken_reads_legacy": (self.kraken_dir or root / f"{name}_kraken2") / f"{name}_kraken2_reads_done.flag",
            "functional": root / f"{name}_functional.success",
            "readmapping": (self.readmapping_dir or root / f"{name}_readmapped") / f"{name}_readmapping.success",
        }


def _first_file(candidates: Iterable[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def _glob_first(base: Path, patterns: Iterable[str]) -> Path | None:
    for pattern in patterns:
        matches = sorted(base.glob(pattern))
        for match in matches:
            if match.is_file():
                return match
    return None


def _first_dir(candidates: Iterable[Path]) -> Path | None:
    return next((p for p in candidates if p.is_dir()), None)


def discover_sample(root_dir: str | Path, sample_dir: str | Path) -> SamplePaths:
    root = Path(root_dir).resolve()
    sample = Path(sample_dir).resolve()
    name = sample.name

    pyrodigal_dir = _first_dir([sample / f"{name}_pyrodigal", sample / "pyrodigal"])
    eggnog_dir = _first_dir([sample / f"{name}_eggnog_mapper", sample / "eggnog_mapper"])
    kraken_dir = _first_dir([sample / f"{name}_kraken2", sample / "kraken2"])

    contigs = _first_file([
        sample / f"{name}_metaspades" / "contigs.fasta",
        sample / f"{name}_unmapped_metaspades" / "contigs.fasta",
        sample / f"{name}_unmapped_trimmed_metaspades" / "contigs.fasta",
        sample / "metaspades" / "contigs.fasta",
    ]) or _glob_first(sample, ["**/contigs.fasta", "**/contigs.fasta.gz"])

    assembly_qc = _first_file([sample / f"{name}_assembly_qc.tsv", sample / "assembly_qc.tsv"])

    pyrodigal_dir = pyrodigal_dir or _first_dir([p.parent for p in [
        _glob_first(sample, [f"**/{name}_pyrodigal.gff"]),
        _glob_first(sample, [f"**/{name}_pyrodigal_aa.fa"]),
    ] if p is not None])
    pyrodigal_gff = _first_file([
        (pyrodigal_dir / f"{name}_pyrodigal.gff") if pyrodigal_dir else sample / f"{name}_pyrodigal/{name}_pyrodigal.gff",
    ]) or _glob_first(sample, [f"**/{name}_pyrodigal.gff"])
    pyrodigal_nt = _first_file([
        (pyrodigal_dir / f"{name}_pyrodigal_nt.fa") if pyrodigal_dir else sample / f"{name}_pyrodigal/{name}_pyrodigal_nt.fa",
        sample / f"{name}_pyrodigal" / f"{name}_pyrodigal.nt.fa",
    ]) or _glob_first(sample, [f"**/{name}_pyrodigal_nt.fa", f"**/{name}_pyrodigal.nt.fa"])
    pyrodigal_aa = _first_file([
        (pyrodigal_dir / f"{name}_pyrodigal_aa.fa") if pyrodigal_dir else sample / f"{name}_pyrodigal/{name}_pyrodigal_aa.fa",
        sample / f"{name}_pyrodigal" / f"{name}_pyrodigal.faa",
    ]) or _glob_first(sample, [f"**/{name}_pyrodigal_aa.fa", f"**/{name}_pyrodigal.faa"])
    pyrodigal_qc = _first_file([
        (pyrodigal_dir / f"{name}_pyrodigal_qc.tsv") if pyrodigal_dir else sample / f"{name}_pyrodigal/{name}_pyrodigal_qc.tsv",
    ]) or _glob_first(sample, [f"**/{name}_pyrodigal_qc.tsv"])

    eggnog_expected = f"{name}_pyrodigal_eggnog_mapped.emapper"
    eggnog_annotations = _first_file([
        (eggnog_dir / f"{eggnog_expected}.annotations") if eggnog_dir else sample / f"{name}_eggnog_mapper/{eggnog_expected}.annotations",
    ]) or _glob_first(sample, ["**/*.emapper.annotations"])
    eggnog_dir = eggnog_dir or (eggnog_annotations.parent if eggnog_annotations else None)
    eggnog_hits = _first_file([
        eggnog_dir / f"{eggnog_expected}.hits" if eggnog_dir else sample / f"{name}_eggnog_mapper/{eggnog_expected}.hits",
    ]) or _glob_first(sample, ["**/*.emapper.hits"])
    eggnog_seed = _first_file([
        eggnog_dir / f"{eggnog_expected}.seed_orthologs" if eggnog_dir else sample / f"{name}_eggnog_mapper/{eggnog_expected}.seed_orthologs",
    ]) or _glob_first(sample, ["**/*.emapper.seed_orthologs"])
    eggnog_gff = _first_file([
        eggnog_dir / f"{eggnog_expected}.decorated.gff" if eggnog_dir else sample / f"{name}_eggnog_mapper/{eggnog_expected}.decorated.gff",
    ]) or _glob_first(sample, ["**/*.emapper.decorated.gff"])
    eggnog_xlsx = _first_file([
        eggnog_dir / f"{eggnog_expected}.annotations.xlsx" if eggnog_dir else sample / f"{name}_eggnog_mapper/{eggnog_expected}.annotations.xlsx",
    ]) or _glob_first(sample, ["**/*.emapper.annotations.xlsx"])
    eggnog_qc = _first_file([
        eggnog_dir / f"{name}_eggnogmapper_qc.tsv" if eggnog_dir else sample / f"{name}_eggnog_mapper/{name}_eggnogmapper_qc.tsv",
    ]) or _glob_first(sample, [f"**/{name}_eggnogmapper_qc.tsv"])

    kraken_contigs = _first_file([
        kraken_dir / f"{name}_kraken2.txt" if kraken_dir else sample / f"{name}_kraken2/{name}_kraken2.txt",
        sample / f"{name}_kraken2.txt",
    ]) or _glob_first(sample, [f"**/{name}_kraken2.txt"])
    kraken_contig_report = _first_file([
        kraken_dir / f"{name}_kraken2_report.txt" if kraken_dir else sample / f"{name}_kraken2/{name}_kraken2_report.txt",
        sample / f"{name}_kraken2_report.txt",
    ]) or _glob_first(sample, [f"**/{name}_kraken2_report.txt"])
    kraken_mpa = _first_file([
        kraken_dir / f"{name}_kraken2_report_mpa.txt" if kraken_dir else sample / f"{name}_kraken2/{name}_kraken2_report_mpa.txt",
        sample / f"{name}_kraken2_report_mpa.txt",
    ]) or _glob_first(sample, [f"**/{name}_kraken2_report_mpa.txt"])
    kraken_reads = _first_file([
        kraken_dir / f"{name}_kraken2_reads.txt" if kraken_dir else sample / f"{name}_kraken2/{name}_kraken2_reads.txt",
    ]) or _glob_first(sample, [f"**/{name}_kraken2_reads.txt"])
    kraken_reads_report = _first_file([
        kraken_dir / f"{name}_kraken2_reads_report.txt" if kraken_dir else sample / f"{name}_kraken2/{name}_kraken2_reads_report.txt",
    ]) or _glob_first(sample, [f"**/{name}_kraken2_reads_report.txt"])
    kraken_dir = kraken_dir or (kraken_contigs.parent if kraken_contigs else None)

    candidate_map_dirs = [
        sample / f"{name}_readmapped",
        root / "readmapping" / name,
        root / "readmapping_v3" / name,
        root / "readmapped" / name,
        sample / "readmapping",
    ]
    readmapping_dir = _first_dir(candidate_map_dirs)
    summary_names = [f"{name}_readmapped_contig_summary.txt", "readmapped_contig_summary.txt"]
    cds_names = [f"{name}_readmapped_cds_summary.txt", "readmapped_cds_summary.txt"]
    contig_summary = _first_file([
        *(d / filename for d in candidate_map_dirs for filename in summary_names),
        sample / f"{name}_readmapped_contig_summary.txt",
    ]) or _glob_first(sample, [f"**/{name}_readmapped_contig_summary.txt", "**/readmapped_contig_summary.txt"])
    cds_summary = _first_file([
        *(d / filename for d in candidate_map_dirs for filename in cds_names),
        sample / f"{name}_readmapped_cds_summary.txt",
    ]) or _glob_first(sample, [f"**/{name}_readmapped_cds_summary.txt", "**/readmapped_cds_summary.txt"])
    if contig_summary:
        readmapping_dir = contig_summary.parent
    elif cds_summary:
        readmapping_dir = cds_summary.parent

    bam = _first_file([
        sample / f"{name}_sorted.bam",
        sample / f"{name}_readmapped_sorted.bam",
        (readmapping_dir / f"{name}_readmapped_sorted.bam") if readmapping_dir else sample / f"{name}.bam",
    ]) or _glob_first(sample, ["**/*sorted.bam", "**/*.bam"])
    bam_index = _first_file([
        Path(f"{bam}.bai") if bam else sample / "missing.bai",
        Path(str(bam).removesuffix(".bam") + ".bai") if bam else sample / "missing.bai",
        Path(f"{bam}.csi") if bam else sample / "missing.csi",
        Path(str(bam).removesuffix(".bam") + ".csi") if bam else sample / "missing.csi",
    ])

    return SamplePaths(
        name=name, sample_dir=sample, contigs=contigs, assembly_qc=assembly_qc,
        pyrodigal_dir=pyrodigal_dir, pyrodigal_gff=pyrodigal_gff, pyrodigal_nt=pyrodigal_nt,
        pyrodigal_aa=pyrodigal_aa, pyrodigal_qc=pyrodigal_qc,
        eggnog_dir=eggnog_dir, eggnog_annotations=eggnog_annotations, eggnog_hits=eggnog_hits,
        eggnog_seed_orthologs=eggnog_seed, eggnog_decorated_gff=eggnog_gff,
        eggnog_xlsx=eggnog_xlsx, eggnog_qc=eggnog_qc,
        kraken_dir=kraken_dir, kraken_contigs=kraken_contigs, kraken_contig_report=kraken_contig_report,
        kraken_mpa_report=kraken_mpa, kraken_reads=kraken_reads, kraken_reads_report=kraken_reads_report,
        readmapping_dir=readmapping_dir, contig_summary=contig_summary, cds_summary=cds_summary,
        bam=bam, bam_index=bam_index,
    )


def parse_prefixes(prefixes: str | Iterable[str]) -> tuple[str, ...]:
    if isinstance(prefixes, str):
        values = prefixes.split(",")
    else:
        values = list(prefixes)
    return tuple(p.strip() for p in values if p and p.strip())


def discover_samples(root_dir: str | Path, prefixes: str | Iterable[str]) -> list[SamplePaths]:
    root = Path(root_dir).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Sample root does not exist or is not a directory: {root}")
    tags = parse_prefixes(prefixes)
    if not tags:
        raise ValueError("At least one sample-directory prefix is required, e.g. '3,4,N,P'.")
    dirs = sorted(
        (path for path in root.iterdir() if path.is_dir() and path.name.startswith(tags)),
        key=lambda path: path.name,
    )
    return [discover_sample(root, path) for path in dirs]
