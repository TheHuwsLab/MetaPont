# Changelog

## 0.1.4

- Removed modification-time freshness comparisons from stage inventory and flag healing. Scratch `touch` operations no longer make otherwise present outputs appear stale.
- Kept automatic success-flag healing stricter than MetaPont's input-availability checks: flags are recreated only when the full expected upstream stage output set, including QC/sidecars where applicable, is present.
- Stage completeness now checks the non-empty files MetaPont actually consumes; optional sidecars and original BAM files no longer gate the integration analysis.
- Replaced mtime-dependent cache identity with SHA256 content fingerprints. The default sampled mode hashes seven evenly spaced blocks and is touch-resilient; `--signature-mode full` hashes every byte.
- Included the MetaPont version, signature algorithm and discovered sample set in the cache digest and explicit cache validation, so newly added empty sample directories also invalidate cached summaries.
- Added an in-progress marker and delayed manifest commit so interrupted runs cannot be accepted as a valid cache based on an old manifest and partially overwritten outputs.
- Made `taxon_function_abundance.tsv` genuinely optional: it is not created when disabled, and stale tables from previous opt-in runs are removed.
- Added disk-backed SQLite aggregation and the `taxon-function-abundance` / `taxon-function` command to generate the optional taxon-function table from an existing gene annotation table.
- Added gene-level functional-annotation status and per-sample counts/rates for genes with and without parsed functional features.
- Added sample-level mapped-read support summaries for classified, unclassified, length-filtered and Kraken-unassigned contigs.
- Fixed query commands using `--gene-table` without `--root-dir` and corrected COG feature splitting in query operations.
- Aligned the Python API default with the CLI (`include_taxon_function_abundance=False`); the legacy function-by-taxa command explicitly enables the table.
- Added regression tests for touch-resilient cache validation, stage checks, functional/taxonomic summaries, post-hoc taxon-function aggregation and CLI query paths.

## 0.1.3

- Previous release.
