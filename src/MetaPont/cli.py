"""Command-line interface for MetaPont."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from . import MetaPont_Version, __version__
from .analysis import coverage_summary, query_function, query_taxon
from .layout import discover_samples
from .pipeline import OUTPUT_DIRNAME, build_taxon_function_abundance, run_analysis, write_inventory


def _add_root_prefix(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-d", "--root-dir", "--directory", dest="root_dir", required=True,
                        help="Root directory containing the main sample directories.")
    parser.add_argument("-p", "--prefix", required=True,
                        help="Comma-separated sample-directory prefixes, e.g. 3,4,N,P or E,L,P.")
    parser.add_argument("-o", "--output-dir", default=None,
                        help=f"Results directory (default: ROOT/{OUTPUT_DIRNAME}).")


def _output_root(root_dir: str | None, output_dir: str | None) -> Path:
    if output_dir:
        return Path(output_dir).expanduser().resolve()
    if not root_dir:
        raise ValueError("Provide --root-dir, --output-dir or --gene-table.")
    root = Path(root_dir).expanduser().resolve()
    return root if root.name == OUTPUT_DIRNAME else root / OUTPUT_DIRNAME


def _gene_table(root_dir: str | None, output_dir: str | None, supplied: str | None) -> Path:
    if supplied:
        return Path(supplied).expanduser().resolve()
    return _output_root(root_dir, output_dir) / "combined" / "gene_annotations.tsv.gz"


def _create_plots(output_root: Path) -> list[str]:
    """Create optional quick-look plots; core analysis does not require matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return ["matplotlib is not installed; plots skipped. Install MetaPont with the [plots] extra."]

    plot_dir = output_root / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []
    summary_file = output_root / "summary" / "sample_summary.tsv"
    if summary_file.is_file():
        with summary_file.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if rows:
            rows.sort(key=lambda r: int(r.get("Total_Mapped_Reads") or 0))
            fig, ax = plt.subplots(figsize=(max(10, len(rows) * 0.25), 6))
            ax.bar(range(len(rows)), [int(r.get("Total_Mapped_Reads") or 0) for r in rows])
            ax.set_xticks(range(len(rows)))
            ax.set_xticklabels([r["Sample"] for r in rows], rotation=90, fontsize=7)
            ax.set_ylabel("Sum of contig mapped-read counts")
            ax.set_title("Readmapping support per sample")
            fig.tight_layout()
            fig.savefig(plot_dir / "mapped_read_support_per_sample.png", dpi=180)
            plt.close(fig)

    func_file = output_root / "functional" / "functional_abundance.tsv"
    if func_file.is_file():
        with func_file.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows.sort(key=lambda r: int(r.get("Gene_Read_Support") or 0), reverse=True)
        top = rows[:20]
        if top:
            fig, ax = plt.subplots(figsize=(11, 7))
            labels = [f"{r['Category']}:{r['Feature']}" for r in reversed(top)]
            vals = [int(r.get("Gene_Read_Support") or 0) for r in reversed(top)]
            ax.barh(range(len(labels)), vals)
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlabel("Gene read-support (CDS overlaps; non-exclusive)")
            ax.set_title("Top functional features by gene read-support")
            fig.tight_layout()
            fig.savefig(plot_dir / "top_functional_features.png", dpi=180)
            plt.close(fig)
    if not list(plot_dir.glob("*.png")):
        warnings.append("No plot-ready data were found.")
    return warnings


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="metapont",
        description="MetaPont " + MetaPont_Version + " Integrates MetaSPAdes, readmapping, Kraken2 and eggNOG-mapper outputs for per-sample and cross-sample analysis.",
    )
    parser.add_argument("--version", action="version", version=f"MetaPont {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run_p = sub.add_parser("run", aliases=["combine"], help="Build the integrated sample inventory, contig table and abundance summaries.")
    _add_root_prefix(run_p)
    run_p.add_argument("--force", action="store_true", help="Regenerate outputs even when the recorded input signature is unchanged.")
    run_p.add_argument("--no-heal-flags", action="store_true", help="Do not recreate missing upstream flags, even when the full expected stage output set is present.")
    run_p.add_argument("--min-contig-length", type=int, default=2500, help="Minimum contig length to retain an existing Kraken2 assignment; classifications are not rerun.")
    run_p.add_argument("--plots", action="store_true", help="Generate optional quick-look PNG plots (requires matplotlib).")
    run_p.add_argument(
        "--signature-mode", choices=("sampled", "full"), default="sampled",
        help="Input change detection: sampled hashes seven evenly spaced blocks (touch-safe, default), or full SHA256 of every input byte (slower).",
    )
    run_p.add_argument(
        "--taxon-function-abundance",
        action="store_true",
        help=(
            "Generate the expanded taxon-function abundance table using a disk-backed aggregation. "
            "This requires additional disk I/O and temporary space."
        ),
    )
    inv_p = sub.add_parser("inventory", help="Audit per-sample inputs, pipeline outputs and success flags.")
    _add_root_prefix(inv_p)
    inv_p.add_argument("--heal-flags", action="store_true", help="Recreate missing upstream flags only when the full expected stage output set is present; timestamps are not used.")

    cov_p = sub.add_parser("coverage", help="Summarise per-contig readmapping statistics by sample.")
    _add_root_prefix(cov_p)
    cov_p.add_argument("--read-length", type=int, default=None, help="Optional read length for an approximate mapped-read coverage estimate.")
    cov_p.add_argument("--output", default=None, help="Output TSV path (defaults under MetaPont_Results/summary).")

    qf_p = sub.add_parser("query-function", help="Find taxa supporting a particular eggNOG functional ID.")
    qf_p.add_argument("--root-dir", "-d", default=None,
                      help="Workflow root used to find the default gene table; not needed with --gene-table or --output-dir.")
    qf_p.add_argument("--gene-table", default=None,
                      help="Gene annotation TSV/TSV.GZ; defaults to the integrated output table.")
    qf_p.add_argument("--output-dir", default=None)
    qf_p.add_argument("--function-id", "-f", required=True, help="Exact function ID, e.g. GO:0016597 or EC:2.7.11.1.")
    qf_p.add_argument("--category", "-c", default=None, help="Optional annotation category, e.g. GO, EC, KEGG_KO, PFAM.")
    qf_p.add_argument("--rank", default="Genus", help="Taxonomic rank for attribution (default: Genus).")
    qf_p.add_argument("--top", type=int, default=None, help="Keep the top N contributing taxa per sample.")
    qf_p.add_argument("--output", "-o", required=True)

    qt_p = sub.add_parser("query-taxon", help="Report the most supported functions for a taxon.")
    qt_p.add_argument("--root-dir", "-d", default=None)
    qt_p.add_argument("--gene-table", default=None,
                      help="Gene annotation TSV/TSV.GZ; defaults to the integrated output table.")
    qt_p.add_argument("--output-dir", default=None)
    qt_p.add_argument("--taxon", "-t", required=True, help="Taxon name or lineage token, e.g. Bacillus or g__Bacillus.")
    qt_p.add_argument("--category", "-c", required=True, help="Annotation category, e.g. GO, EC, KEGG_KO.")
    qt_p.add_argument("--rank", default="Genus")
    qt_p.add_argument("--top", type=int, default=20)
    qt_p.add_argument("--output", "-o", required=True)

    tf_p = sub.add_parser(
        "taxon-function-abundance", aliases=["taxon-function"],
        help="Build the expanded taxon-by-function table from an existing gene annotation table without rerunning the main analysis.",
    )
    tf_p.add_argument("--root-dir", "-d", default=None, help="Workflow root or MetaPont_Results directory.")
    tf_p.add_argument("--gene-table", default=None, help="Existing gene_annotations.tsv or .tsv.gz file.")
    tf_p.add_argument("--output-dir", default=None, help="MetaPont results directory (default derived from --root-dir).")
    tf_p.add_argument("--output", "-o", default=None, help="Output TSV path (default: OUTPUT_DIR/functional/taxon_function_abundance.tsv).")

    args = parser.parse_args(argv)
    try:
        if args.command in {"run", "combine"}:
            result = run_analysis(
                args.root_dir, args.prefix, args.output_dir,
                min_contig_length=args.min_contig_length,
                heal_flags=not args.no_heal_flags,
                force=args.force,
                include_taxon_function_abundance=args.taxon_function_abundance,
                signature_mode=args.signature_mode,
            )
            if args.plots:
                output_root = Path(result["output_dir"])
                plot_warnings = _create_plots(output_root)
                if plot_warnings:
                    (output_root / "summary" / "plot_warnings.txt").write_text("\n".join(plot_warnings) + "\n", encoding="utf-8")
                else:
                    print(f"Plots written to {output_root / 'plots'}")
            print(f"MetaPont {MetaPont_Version}: results written to {result['output_dir']}")
            print(f"Samples: {result['samples']}; samples with warnings: {result.get('samples_with_warnings', result.get('warnings', 0))}")
            print(f"Completion marker: {result.get('success_flag', 'existing output cache validated')}")
            return 0
        if args.command == "inventory":
            out, rows, events = write_inventory(args.root_dir, args.prefix, args.output_dir, heal_flags=args.heal_flags)
            print(f"Inventory written for {len(rows)} samples: {out}")
            if events:
                print(f"Healed {len(events)} flags.")
            return 0
        if args.command == "coverage":
            supplied_root = Path(args.root_dir).expanduser().resolve()
            samples_root = supplied_root.parent if supplied_root.name == OUTPUT_DIRNAME else supplied_root
            default_output_root = (
                Path(args.output_dir).expanduser().resolve() if args.output_dir
                else supplied_root if supplied_root.name == OUTPUT_DIRNAME
                else supplied_root / OUTPUT_DIRNAME
            )
            default_output = default_output_root / "summary" / "contig_coverage_summary.tsv"
            out, count = coverage_summary(
                samples_root, args.prefix, args.output or default_output, args.read_length,
            )
            print(f"Coverage summary written for {count} samples: {out}")
            return 0
        if args.command == "query-function":
            if args.root_dir is None and args.gene_table is None and args.output_dir is None:
                raise ValueError("Provide --root-dir, --output-dir or --gene-table.")

            # A MetaPont_Results directory can be supplied directly as root-dir.
            if args.root_dir and Path(args.root_dir).expanduser().resolve().name == OUTPUT_DIRNAME and not args.output_dir:
                args.output_dir = args.root_dir

            table = _gene_table(args.root_dir, args.output_dir, args.gene_table)
            count = query_function(table, args.function_id, args.output, args.category, args.rank, args.top)
            print(f"Wrote {count} contributing taxon rows to {args.output}")
            return 0

        if args.command == "query-taxon":
            if args.root_dir is None and args.gene_table is None and args.output_dir is None:
                raise ValueError("Provide --root-dir, --output-dir or --gene-table.")

            # A MetaPont_Results directory can be supplied directly as root-dir.
            if args.root_dir and Path(args.root_dir).expanduser().resolve().name == OUTPUT_DIRNAME and not args.output_dir:
                args.output_dir = args.root_dir

            table = _gene_table(args.root_dir, args.output_dir, args.gene_table)
            count = query_taxon(table, args.taxon, args.output, args.category, args.top, args.rank)

            print(f"Wrote {count} taxon-associated function rows to {args.output}")
            return 0
        if args.command in {"taxon-function-abundance", "taxon-function"}:
            if args.gene_table is None and args.root_dir is None and args.output_dir is None:
                raise ValueError("Provide --gene-table, --root-dir or --output-dir.")
            table = _gene_table(args.root_dir, args.output_dir, args.gene_table)
            if args.output:
                output = Path(args.output).expanduser().resolve()
            elif args.output_dir or args.root_dir:
                output = _output_root(args.root_dir, args.output_dir) / "functional" / "taxon_function_abundance.tsv"
            else:
                output = table.parent.parent / "functional" / "taxon_function_abundance.tsv"
            count = build_taxon_function_abundance(table, output)
            print(f"Wrote {count} taxon/function rows to {output}")
            return 0
    except (FileNotFoundError, RuntimeError, ValueError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    parser.error(f"Unhandled command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
