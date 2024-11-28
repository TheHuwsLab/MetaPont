import argparse
import os
import csv
import sys
from collections import defaultdict

# Adjust constants import for module vs. standalone script usage
try:
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError):
    from constants import *

# Needed to account for the large CSV/TSV files
csv.field_size_limit(sys.maxsize)


def process_tsv_by_taxon(file_path, target_taxon):
    """
    Processes a TSV file to calculate the top functions grouped by taxon.

    Parameters:
        file_path (str): Path to the TSV file.
        target_taxon (str): Taxon to search for in the lineage column.

    Returns:
        taxon_function_reads (dict): A dictionary mapping functions to their total reads for the specified taxon.
        total_reads_taxon (int): Total reads assigned to the specified taxon across all functions.
    """
    taxon_function_presence = defaultdict(int)  # Reads for each function for the specified taxon
    total_reads_taxon = 0  # Total reads for the specified taxon

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

            # Check if this row matches the specified taxon
            if target_taxon in lineage:
                total_reads_taxon += reads  # Increment total reads for the taxon

                # Loop through the function columns (row[6:] onward)
                for col_idx, cell in enumerate(row[6:], start=6):
                    if cell:
                        header = headers[col_idx]
                        if header not in taxon_function_presence:
                            taxon_function_presence[header] = defaultdict(int)
                        for function in cell.replace(',', '|').split('|'):
                            taxon_function_presence[header][function] += 1

    return taxon_function_presence, total_reads_taxon


def main():
    parser = argparse.ArgumentParser(description='MetaPont: Extract Top Functions by Taxon')
    parser.add_argument(
        "-d", "--directory", required=True,
        help="Directory containing TSV files to analyse."
    )
    parser.add_argument(
        "-t", "--taxon", required=True,
        help="Target taxon to search for (e.g., 'g__Escherichia')."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output file to save results."
    )
    parser.add_argument(
        "-func", "--functional_classes", required=True,
        help="Which functional classes to report (e.g. GO,EC,KEGG etc)."
    )
    parser.add_argument(
        "-top", "--top_functions", type=int, default=3,
        help="Top n functions to include in the output for each sample (default: 3)."
    )

    options = parser.parse_args()
    print("Running MetaPont: Extract Top Functions by Taxon")

    input_path = os.path.abspath(options.directory)
    output_path = os.path.abspath(options.output)

    all_results = {}

    # Process each TSV file in the directory
    for file_name in os.listdir(input_path):
        if file_name.endswith("_Final_Contig.tsv"):
            file_path = os.path.join(options.directory, file_name)
            print(f"Processing file: {file_name}")
            taxon_function_presence, total_reads_taxon = process_tsv_by_taxon(file_path, options.taxon)
            all_results[file_name] = (taxon_function_presence, total_reads_taxon)

    # Write results to output
    with open(output_path, "w") as out:
        out.write("Selected Taxon: " + options.taxon + "\n")
        out.write("Sample\tFunction\tNum of Assignments\n")
        for sample, (taxon_function_presence, total_assignments_taxon) in all_results.items():
            for function, assignments_dict in taxon_function_presence.items():
                if any(func_class in function for func_class in options.functional_classes.split(',')):
                    sorted_functions = sorted(assignments_dict.items(), key=lambda x: x[1], reverse=True)
                    for i, (sub_function, assignments) in enumerate(sorted_functions):
                        if i < options.top_functions:
                            #proportion_taxon = assignments / total_reads_taxon if total_reads_taxon > 0 else 0
                            out.write(f"{sample.replace('_Final_Contig.tsv', '')}\t{sub_function}\t{assignments}\n")
                        else:
                            break
        #     sorted_functions = sorted(taxon_function_presence.items(), key=lambda x: x[1], reverse=True)
        #     for i, (function, reads) in enumerate(sorted_functions):
        #         if i < options.top_functions:
        #             proportion_taxon = reads / total_reads_taxon if total_reads_taxon > 0 else 0
        #             out.write(f"{sample.replace('_Final_Contig.tsv', '')}\t{function}\t{reads}\t{proportion_taxon:.3f}\n")
        #         else:
        #             break

    print(f"Results saved to {options.output}")


if __name__ == "__main__":
    main()
