# MetaPont v0.1.2

MetaPont integrates existing metagenomic workflow outputs into per-sample and
cross-sample tables for contigs, taxonomy, readmapping support, and eggNOG-mapper
functional annotations. It is designed for outputs from the
[HuwsLab Metagenomic Workflow](https://github.com/TheHuwsLab/Metagenomic_Workflow),
but does not run assembly, readmapping, Kraken2, or eggNOG-mapper itself.

## Installation

MetaPont requires Python 3.10 or later and has no required third-party runtime
dependencies.

```bash
pip install MetaPont
```

Check the installed version and available commands with:

```bash
metapont --version
metapont --help
```

## Input layout

Pass MetaPont the directory containing sample directories. Sample directories
are selected when their names start with one of the comma-separated prefixes
provided with `--prefix`.

```text
workflow-output/
├── E1.0/
│   ├── E1.0_metaspades/contigs.fasta
│   ├── E1.0_pyrodigal/
│   ├── E1.0_eggnog_mapper/
│   │   └── E1.0_pyrodigal_eggnog_mapped.emapper.annotations
│   ├── E1.0_kraken2/
│   │   ├── E1.0_kraken2.txt
│   │   └── E1.0_kraken2_report.txt
│   └── E1.0_readmapped/
│       ├── E1.0_readmapped_contig_summary.txt
│       └── E1.0_readmapped_cds_summary.txt
└── E2.0/
    └── ...
```

MetaPont discovers common workflow naming variants and searches within each
sample directory for the relevant files. The main inputs are:

| Input | Used for |
| --- | --- |
| eggNOG-mapper `.emapper.annotations` | Gene annotations and functional features |
| Kraken2 contig classifications and report | Contig taxonomic assignments and lineages |
| Readmapping contig summary | Contig lengths, mapped-read counts, and mapping proportions |
| Readmapping CDS summary | Read support attributed to annotated genes |
| Contig FASTA | Contig inventory and lengths when a readmapping contig summary is absent |

Inputs are discovered independently, so incomplete samples can still be
included. Missing or stale outputs are reported as warnings; MetaPont does not
rerun missing upstream analyses. Kraken2 assignments shorter than the selected
minimum contig length are omitted from taxonomy attribution by default. This
filter does not rerun Kraken2.

## Commands

### Build integrated results

```bash
metapont run \
  --root-dir /path/to/workflow-output \
  --prefix E
```

`--prefix` accepts comma-separated values, for example `E,L,P`. By default,
results are written under `ROOT/MetaPont_Results`. Override that location with
`--output-dir`.

Useful options:

| Option | Default | Description |
| --- | --- | --- |
| `--force` | Off | Regenerate outputs even when the input signature matches a cached run. |
| `--no-heal-flags` | Off | Do not recreate missing upstream success flags when their outputs are complete and current. |
| `--min-contig-length N` | `2500` | Minimum Kraken2-assigned contig length retained for taxonomy attribution; use `0` to disable the filter. |
| `--taxon-function-abundance` | Off | Build the expanded taxon-by-function table. This can use substantial memory and disk space. |
| `--plots` | Off | Write optional quick-look PNG plots. Requires the `plots` installation extra. |

`metapont combine` is an alias for `metapont run`. The `metapont-combine`
command is also available for compatibility with the earlier combine command.

### Audit sample inputs

```bash
metapont inventory --root-dir /path/to/workflow-output --prefix E
```

The inventory reports discovered inputs, expected outputs, success flags, and
warnings for each sample. Add `--heal-flags` to recreate a missing stage flag
only when its corresponding output files are complete and current.

### Summarise readmapping coverage

```bash
metapont coverage \
  --root-dir /path/to/workflow-output \
  --prefix E \
  --read-length 150
```

The optional read length adds an approximate mapped-read coverage estimate.
Without `--output`, the TSV is written to
`MetaPont_Results/summary/contig_coverage_summary.tsv`.

### Query functions and taxa

Find taxa supporting a functional ID:

```bash
metapont query-function \
  --root-dir /path/to/workflow-output \
  --function-id GO:0016597 \
  --category GO \
  --rank Genus \
  --top 10 \
  --output go_taxa.tsv
```

Omit `--category` to search all supported annotation categories. `--top`
limits the results to the highest-support taxa per sample. Alternatively,
provide `--gene-table` to query a particular integrated gene table, or
`--output-dir` to select a non-default results directory.

Find the most supported features associated with a taxon:

```bash
metapont query-taxon \
  --root-dir /path/to/workflow-output \
  --taxon Bacillus \
  --category GO \
  --rank Genus \
  --top 20 \
  --output bacillus_go.tsv
```

Taxon matching is case-insensitive and searches the chosen rank as well as the
full lineage. `--category` is required for this query. Both query commands
accept an explicit `--gene-table` path; plain TSV and gzip-compressed TSV files
are supported.

## Results

The default output tree is:

```text
MetaPont_Results/
├── combined/
│   ├── all_contigs.tsv.gz
│   └── gene_annotations.tsv.gz
├── functional/
│   ├── functional_abundance.tsv
│   ├── matrices/
│   └── taxon_function_abundance.tsv
├── inventory/
│   └── sample_inventory.tsv
├── metadata/
│   ├── input_signature.json
│   └── run_manifest.json
├── per_sample_contigs/
│   └── SAMPLE_Contigs.tsv.gz
├── summary/
│   ├── pipeline_warnings.tsv
│   ├── sample_inventory.tsv
│   └── sample_summary.tsv
└── taxonomy/
    └── taxonomy_abundance.tsv
```

The main tables are tab-delimited; `.tsv.gz` files are gzip-compressed TSVs.
Functional matrices are written separately for every supported annotation
category, with both gene counts and gene read-support values. The optional
taxon-function abundance table contains taxon/category/feature combinations
across taxonomic ranks.

Supported functional categories include `COG`, `eggNOG_OG`, `GO`, `EC`,
`KEGG_KO`, `KEGG_Pathway`, `KEGG_Module`, `KEGG_Reaction`, `KEGG_rclass`,
`BRITE`, `KEGG_TC`, `CAZy`, `BiGG_Reaction`, and `PFAM`.

`Gene_Read_Support` is derived from read-to-CDS overlaps in the CDS summary; it
is not a count of unique reads across genes. A read can support multiple genes
or annotations, so functional support values are not mutually exclusive.
Taxonomic relative abundance uses total mapped-read support as its denominator.

Samples with missing or stale workflow outputs are retained where possible and
listed in `summary/pipeline_warnings.tsv`. A clean run writes
`MetaPont_analysis.success`; a run with warnings writes
`completed_with_warnings.flag`.

## License

MetaPont is distributed under the
[GNU General Public License v3.0 or later](LICENSE).
