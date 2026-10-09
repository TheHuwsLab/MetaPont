import argparse
import os
import csv
import sys
from collections import defaultdict

try:
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *

# Needed to account for the large CSV/TSV files we are working with
csv.field_size_limit(sys.maxsize)


def process_tsv(file_path, function_id):
    """
    Processes a TSV file to calculate:
    1. Total reads for each taxon (within the function).
    2. Total reads for the specified function.
    3. Total reads across all contigs in the file.
    """
    taxa_reads_function = defaultdict(int)  # Reads for each taxon within the function
    total_reads_function = 0  # Total reads for the specified function
    total_reads_all = 0  # Total reads across all contigs

    with open(file_path, "r") as tsv_file:
        reader = csv.reader(tsv_file, delimiter="\t")
        next(reader)  # Skip the first row (sample name)
        headers = next(reader)  # Read headers from the second row

        # Locate the necessary columns
        taxa_idx = headers.index("Lineage")
        reads_idx = headers.index("Mapped_Reads")

        # Process each row
        for idx, row in enumerate(reader):
            if len(row) < len(headers):
                continue  # Skip malformed rows

            lineage = row[taxa_idx]
            reads = int(row[reads_idx])  # Get the number of reads for this row
            total_reads_all += reads  # Increment total reads for all contigs

            # Check if this row matches the specified function ID
            for cell in row[6:]:
                if cell and any(function_id in part for part in cell.replace(',', '|').split('|')):
                    if "g__" in lineage:
                        lineage = lineage.split("s__")[0]#.split("|")[0]
                        taxa_reads_function[lineage] += reads
                        total_reads_function += reads
                    else:
                        sys.exit("'__g' not in lineage.")
                    break  # Stop checking further functional columns for this row

    return taxa_reads_function, total_reads_function, total_reads_all


def main():
    parser = argparse.ArgumentParser(description='MetaPont ' + MetaPont_Version + ': Extract-By-Function - Identify taxa contributing to a specific function.')
    parser.add_argument(
        "-d", "--directory", required=True,
        help="Directory containing TSV files to analyse."
    )
    parser.add_argument(
        "-f", "--function_id", required=True,
        help="Specific function ID to search for (e.g., 'GO:0016597')."
    )
    parser.add_argument(
        "-a", "--file_affix", required=True,
        help="Input TSV file affix to use when searching contig files."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output file to save results."
    )
    parser.add_argument(
        "-m", "--min_proportion", type=float,
        help="Minimum proportion threshold for taxa to be included in the output."
    )
    parser.add_argument(
        "-top", "--top_taxa", type=int,
        help="Top n taxa to be included in the output."
    )

    options = parser.parse_args()
    if not options.min_proportion and not options.top_taxa:
        sys.exit("Error: Please specify either a minimum proportion or the number of top taxa to include in the output.")
    elif options.min_proportion and options.top_taxa:
        sys.exit("Error: Please specify either a minimum proportion or the number of top taxa to include in the output, not both.")
    print("Running MetaPont: Extract-By-Function " + MetaPont_Version)

    input_path = os.path.abspath(options.directory)
    output_path = os.path.abspath(options.output)

    all_results = {}

    # Process each TSV file in the directory
    for file_name in os.listdir(input_path):
        if file_name.endswith(options.file_affix):
            file_path = os.path.join(options.directory, file_name)
            print(f"Processing file: {file_name}")
            taxa_reads_function, total_reads_function, total_reads_all = process_tsv(file_path, options.function_id)
            all_results[file_name] = (taxa_reads_function, total_reads_function, total_reads_all)

    # Write results to output
    with open(output_path, "w") as out:
        out.write("Function ID: " + options.function_id + "\n")
        out.write("Sample\tFull Lineage\tGenus\tReads Assigned (Function)\tProportion (Function)\tProportion (Total Reads)\n")
        for sample, (taxa_reads_function, total_reads_function, total_reads_all) in sorted(all_results.items()):
            if options.min_proportion:
                for lineage, reads_function in taxa_reads_function.items():
                    genus = lineage.split("g__")[1].split("|")[0]
                    proportion_function = reads_function / total_reads_function if total_reads_function > 0 else 0
                    proportion_total = reads_function / total_reads_all if total_reads_all > 0 else 0
                    if proportion_function >= options.min_proportion:  # Apply minimum proportion filter
                        out.write(f"{sample}\t{lineage}\t{genus}\t{reads_function}\t{proportion_function:.3f}\t{proportion_total:.3f}\n")
            elif options.top_taxa:
                sorted_taxa_reads = sorted(taxa_reads_function.items(), key=lambda x: x[1], reverse=True)
                for i, (lineage, reads_function) in enumerate(sorted_taxa_reads):
                    if i < options.top_taxa:
                        genus = lineage.split("g__")[1].split("|")[0]
                        proportion_function = reads_function / total_reads_function if total_reads_function > 0 else 0
                        proportion_total = reads_function / total_reads_all if total_reads_all > 0 else 0
                        out.write(f"{sample.replace(options.file_affix,'')}\t{lineage}\t{genus}\t{reads_function}\t{proportion_function:.3f}\t{proportion_total:.3f}\n")
                    else:
                        break

    print(f"Results saved to {options.output}")


if __name__ == "__main__":
    main()
