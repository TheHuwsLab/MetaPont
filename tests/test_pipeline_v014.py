import csv
import gzip
import json
import os
from pathlib import Path

from MetaPont import __version__
from MetaPont import pipeline
from MetaPont.cli import main
from MetaPont.layout import discover_samples


def _write_sample(root: Path) -> Path:
    sample = root / "S1"
    contigs = sample / "S1_metaspades" / "contigs.fasta"
    contigs.parent.mkdir(parents=True)
    contigs.write_text(
        ">ctg_1_length_3000\n" + "A" * 3000 + "\n"
        ">ctg_2_length_3000\n" + "C" * 3000 + "\n",
        encoding="utf-8",
    )

    pyro = sample / "S1_pyrodigal"
    pyro.mkdir()
    (pyro / "S1_pyrodigal.gff").write_text("##gff-version 3\n", encoding="utf-8")
    (pyro / "S1_pyrodigal_aa.fa").write_text(">ctg_1_length_3000_1\nMAAA\n", encoding="utf-8")
    # NT FASTA and QC are intentionally absent: neither is a required input
    # to MetaPont's integration stage.

    eggnog = sample / "S1_eggnog_mapper"
    eggnog.mkdir()
    header = [
        "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
        "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
        "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
        "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
    ]
    values = [
        "ctg_1_length_3000_1", "-", "1e-10", "50", "COG1@1", "-", "AB", "description",
        "geneA", "GO:0001,GO:0002", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-",
    ]
    with (eggnog / "S1_pyrodigal_eggnog_mapped.emapper.annotations").open("w", encoding="utf-8") as handle:
        handle.write("#" + "\t".join(header) + "\n")
        handle.write("\t".join(values) + "\n")

    kraken = sample / "S1_kraken2"
    kraken.mkdir()
    (kraken / "S1_kraken2.txt").write_text(
        "C\tctg_1_length_3000\t1236\t3000\t1:3000\n"
        "U\tctg_2_length_3000\t0\t3000\t0:3000\n",
        encoding="utf-8",
    )
    (kraken / "S1_kraken2_report.txt").write_text(
        "100.00\t1\t0\tR\t1\troot\n"
        "100.00\t1\t0\tD\t2\t  Bacteria\n"
        "100.00\t1\t0\tP\t1224\t    Firmicutes\n"
        "100.00\t1\t0\tG\t1236\t      Bacillus\n",
        encoding="utf-8",
    )

    readmapping = root / "readmapping" / "S1"
    readmapping.mkdir(parents=True)
    (readmapping / "S1_readmapped_contig_summary.txt").write_text(
        "ctg_1_length_3000\t3000\t7\t10\t0.7\n"
        "ctg_2_length_3000\t3000\t3\t10\t0.3\n",
        encoding="utf-8",
    )
    (readmapping / "S1_readmapped_cds_summary.txt").write_text(
        "5 ctg_1_length_3000_1\n", encoding="utf-8"
    )
    return sample


def _read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stage_completeness_ignores_touch_timestamps_and_optional_sidecars(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    sample_dir = _write_sample(root)
    sample = discover_samples(root, "S")[0]
    statuses = pipeline._stage_outputs(sample)
    assert statuses["assembly"]
    assert statuses["pyrodigal"]
    assert statuses["eggnog"]
    assert statuses["kraken_contigs"]
    assert statuses["readmapping"]

    # Artificially make derived products older than their sources. This must
    # not turn complete files into false stale warnings.
    old = 1_000_000_000
    for path in (sample.pyrodigal_gff, sample.pyrodigal_aa, sample.eggnog_annotations):
        os.utime(path, ns=(old, old))
    statuses_after_touch = pipeline._stage_outputs(sample)
    assert statuses_after_touch == statuses


def test_run_is_touch_resilient_and_summarises_functional_and_taxonomic_coverage(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    sample_dir = _write_sample(root)

    first = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=False)
    out = Path(first["output_dir"])
    assert not (out / "functional" / "taxon_function_abundance.tsv").exists()
    assert first["cached"] is False

    summary = _read_tsv(out / "summary" / "sample_summary.tsv")[0]
    assert summary["Genes_Processed"] == "1"
    assert summary["Genes_With_Function"] == "1"
    assert summary["Genes_Without_Function"] == "0"
    assert summary["Mapped_Reads_on_Classified_Contigs"] == "7"
    assert summary["Mapped_Reads_on_Unclassified_Contigs"] == "3"
    assert summary["Percent_Mapped_Reads_on_Classified_Contigs"] == "70.0"

    # Touch every sample file to mimic scratch-retention policy. The cache must
    # remain valid because mtimes are diagnostics, not part of its identity.
    for path in sample_dir.rglob("*"):
        if path.is_file():
            os.utime(path, None)
    second = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=False)
    assert second["cached"] is True

    # A content change to a small input invalidates the cache.
    cds = root / "readmapping" / "S1" / "S1_readmapped_cds_summary.txt"
    cds.write_text("6 ctg_1_length_3000_1\n", encoding="utf-8")
    third = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=False)
    assert third["cached"] is False


def test_taxon_function_can_be_built_posthoc_from_existing_gene_table(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    result = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=False)
    output_root = Path(result["output_dir"])
    gene_table = output_root / "combined" / "gene_annotations.tsv.gz"
    out = output_root / "functional" / "taxon_function_abundance.tsv"

    rows_written = pipeline.build_taxon_function_abundance(gene_table, out)
    rows = _read_tsv(out)
    assert rows_written == len(rows)
    assert rows
    assert {row["Taxon"] for row in rows} == {"Bacteria", "Firmicutes", "Bacillus"}
    assert {row["Feature"] for row in rows} == {"A", "B", "COG1@1", "GO:0001", "GO:0002"}
    assert all(row["Gene_Read_Support"] == "5" for row in rows)


def test_taxon_function_cli_can_use_gene_table_without_root_dir(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    result = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=False)
    gene_table = Path(result["output_dir"]) / "combined" / "gene_annotations.tsv.gz"
    out = tmp_path / "posthoc.tsv"
    rc = main(["taxon-function-abundance", "--gene-table", str(gene_table), "--output", str(out)])
    assert rc == 0
    assert out.is_file()


def test_query_function_cli_accepts_gene_table_without_root_dir_and_handles_cog(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    result = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=False)
    gene_table = Path(result["output_dir"]) / "combined" / "gene_annotations.tsv.gz"
    out = tmp_path / "cog_query.tsv"
    rc = main([
        "query-function", "--gene-table", str(gene_table), "--function-id", "A",
        "--category", "COG", "--output", str(out),
    ])
    assert rc == 0
    rows = _read_tsv(out)
    assert len(rows) == 1
    assert rows[0]["Function_ID"] == "A"


def test_signature_hash_includes_version(tmp_path, monkeypatch):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    samples = discover_samples(root, "S")
    sig1 = pipeline._input_signature(samples, 2500, False, "sampled")
    monkeypatch.setattr(pipeline, "__version__", "99.0.0")
    sig2 = pipeline._input_signature(samples, 2500, False, "sampled")
    assert sig1["signature"] != sig2["signature"]
    assert sig1["version"] == __version__


def test_opt_in_taxon_table_is_created_and_cached(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    result = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=True)
    output_root = Path(result["output_dir"])
    table = output_root / "functional" / "taxon_function_abundance.tsv"
    assert table.is_file()
    assert _read_tsv(table)
    manifest = json.loads((output_root / "metadata" / "run_manifest.json").read_text())
    assert manifest["taxon_function_abundance_enabled"] is True
    assert manifest["taxon_function_abundance_rows"] > 0
    cached = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=True)
    assert cached["cached"] is True


def test_interrupted_run_marker_prevents_cache_reuse(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    result = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=False)
    output_root = Path(result["output_dir"])
    marker = output_root / "metadata" / "run_in_progress.flag"
    marker.write_text("simulated interruption\n", encoding="utf-8")
    rerun = pipeline.run_analysis(root, "S", heal_flags=False, include_taxon_function_abundance=False)
    assert rerun["cached"] is False
    assert not marker.exists()


def test_full_hash_detects_changes_outside_sampled_windows(tmp_path):
    path = tmp_path / "large-input.bin"
    path.write_bytes(b"A" * (16 * 1024 * 1024))
    sampled_before = pipeline._file_fingerprint(path, "sampled")
    full_before = pipeline._file_fingerprint(path, "full")
    # 1 MiB lies between the fixed sampled windows for this 16 MiB input.
    with path.open("r+b") as handle:
        handle.seek(1024 * 1024)
        handle.write(b"B")
    assert pipeline._file_fingerprint(path, "sampled") == sampled_before
    assert pipeline._file_fingerprint(path, "full") != full_before


def test_cached_run_reports_actual_warning_sample_count(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    (root / "S2").mkdir()
    (root / "S3").mkdir()

    first = pipeline.run_analysis(root, "S", heal_flags=False)
    assert first["samples_with_warnings"] == 2
    second = pipeline.run_analysis(root, "S", heal_flags=False)
    assert second["cached"] is True
    assert second["samples_with_warnings"] == 2
    assert second["warnings"] == 2


def test_empty_contig_summary_falls_back_to_assembly_fasta(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    summary = root / "readmapping" / "S1" / "S1_readmapped_contig_summary.txt"
    summary.write_text("", encoding="utf-8")

    result = pipeline.run_analysis(root, "S", heal_flags=False)
    output_root = Path(result["output_dir"])
    row = _read_tsv(output_root / "summary" / "sample_summary.tsv")[0]
    assert row["Contigs_Processed"] == "2"
    assert row["Total_Mapped_Reads"] == "0"
    assert "readmapping_contig_summary_missing_or_empty" in row["Warnings"]


def test_coverage_command_honours_custom_output_dir(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    custom_output = tmp_path / "custom_metapont_results"
    rc = main([
        "coverage", "--root-dir", str(root), "--prefix", "S",
        "--output-dir", str(custom_output),
    ])
    assert rc == 0
    assert (custom_output / "summary" / "contig_coverage_summary.tsv").is_file()


def test_flag_healing_does_not_mark_partial_upstream_stages_complete(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    sample = discover_samples(root, "S")[0]

    events = pipeline.heal_sample_flags(sample)
    assert events == []
    assert not sample.flag_paths["assembly"].exists()
    assert not sample.flag_paths["pyrodigal"].exists()
    assert not sample.flag_paths["eggnog"].exists()
    assert not sample.flag_paths["functional"].exists()


def test_empty_or_corrupt_gzip_is_not_accepted_as_cached_output(tmp_path):
    empty_gzip = tmp_path / "empty.tsv.gz"
    with gzip.open(empty_gzip, "wt", encoding="utf-8"):
        pass
    corrupt_gzip = tmp_path / "corrupt.tsv.gz"
    corrupt_gzip.write_bytes(b"not a gzip stream")
    assert pipeline._gzip_has_header(empty_gzip) is False
    assert pipeline._gzip_has_header(corrupt_gzip) is False


def test_adding_empty_sample_directory_invalidates_cached_run(tmp_path):
    root = tmp_path / "workflow"
    root.mkdir()
    _write_sample(root)
    first = pipeline.run_analysis(root, "S", heal_flags=False)
    assert first["cached"] is False

    (root / "S2").mkdir()
    second = pipeline.run_analysis(root, "S", heal_flags=False)
    assert second["cached"] is False
    assert second["samples"] == 2
    summaries = _read_tsv(Path(second["output_dir"]) / "summary" / "sample_summary.tsv")
    assert {row["Sample"] for row in summaries} == {"S1", "S2"}
