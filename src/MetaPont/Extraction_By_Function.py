import argparse
import os
import csv
import sys
from collections import Counter

from HuwsLab.MetaPont.src.MetaPont.constants import MetaPont_Version

# Needed to account for the large CSV/TSV files we are working with
csv.field_size_limit(sys.maxsize)


def process_tsv(file_path, function_id):
    """
    Processes a TSV file to calculate taxa counts and total matches for a given function ID.
    """
    taxa_counts = Counter()
    total_matches = 0

    with open(file_path, "r") as tsv_file:
        reader = csv.reader(tsv_file, delimiter="\t")
        next(reader)  # Skip the first row (sample name)
        headers = next(reader)  # Read headers from the second row

        # Locate the Lineage column
        taxa_idx = headers.index("Lineage")

        # Process each row to find matches
        for idx, row in enumerate(reader):
            if len(row) < len(headers):
                continue  # Skip malformed rows

            lineage = row[taxa_idx]

            # Search functional columns starting from column 6
            for cell in row[6:]:
                if cell and any(function_id in part for part in cell.replace(',', '|').split('|')):
                    # Extract genus from the Lineage field
                    if "g__" in lineage:
                        genus = lineage.split("g__")[1].split("|")[0]
                        taxa_counts[genus] += 1
                        total_matches += 1
                    break  # Stop checking further functional columns for this row

    return taxa_counts, total_matches


def main():

    parser = argparse.ArgumentParser(description='MetaPont ' + MetaPont_Version + ': Extract-By-Function - Identify taxa contributing to a specific function.')
    parser.add_argument(
        "-d", "--directory", required=True,
        help="Directory containing TSV files to analyse."
    )
    parser.add_argument(
        "-f", "--function_id", required=True,
        help="Specific function ID to search for (e.g., 'GO:0002')."
    )
    parser.add_argument(
        "-o", "--output", default="output_taxa_proportions.tsv",
        help="Output file to save results (default: output_taxa_proportions.tsv)."
    )
    parser.add_argument(
        "-m", "--min_proportion", type=float, default=0.05,
        help="Minimum proportion threshold for taxa to be included in the output (default: 0.05)."
    )

    options = parser.parse_args()
    print("Running MetaPont: Extract-By-Function " + MetaPont_Version)

    all_results = {}

    # Process each TSV file in the directory
    for file_name in os.listdir(options.directory):
        if file_name.endswith("_Final_Contig.tsv"):
            file_path = os.path.join(options.directory, file_name)
            print(f"Processing file: {file_name}")
            taxa_counts, total_matches = process_tsv(file_path, options.function_id)
            all_results[file_name] = (taxa_counts, total_matches)

    # Write results to output
    with open(options.output, "w") as out:
        out.write("Function ID: " + options.function_id + "\n")
        out.write("Sample\tTaxa\tProportion\n")
        for sample, (taxa_counts, total_matches) in all_results.items():
            for taxa, count in taxa_counts.items():
                proportion = count / total_matches if total_matches > 0 else 0
                if proportion >= options.min_proportion:  # Apply minimum proportion filter
                    out.write(f"{sample}\t{taxa}\t{proportion:.6f}\n")

    print(f"Results saved to {options.output}")


if __name__ == "__main__":
    main()
