import argparse
import os
import csv
import sys
from collections import defaultdict

# Needed to account for large CSV/TSV files
csv.field_size_limit(sys.maxsize)

def process_tsv_by_taxon_and_function(file_path, target_taxon, target_function):
    """
    Processes a TSV file to calculate the proportion of reads assigned to a taxon and function
    relative to all reads assigned to that function across all taxa.

    Parameters:
        file_path (str): Path to the TSV file.
        target_taxon (str): Taxon to search for in the lineage column.
        target_function (str): Function to search for in the function columns.

    Returns:
        function_proportion (float): Proportion of reads assigned to the function within the taxon
                                     relative to all reads assigned to the function.
        reads_taxon_function (int): Total reads assigned to the taxon and function.
        total_reads_function (int): Total reads assigned to the function across all taxa.
    """
    reads_taxon_function = 0  # Reads assigned to the specified taxon and function
    total_reads_function = 0  # Total reads assigned to the function across all taxa

    with open(file_path, "r") as tsv_file:
        reader = csv.reader(tsv_file, delimiter="\t")
        next(reader)  # Skip the first row (sample name)
        headers = next(reader)  # Read headers from the second row

        # Ensure required columns exist
        try:
            taxa_idx = headers.index("Lineage")
            reads_idx = headers.index("Mapped_Reads")
        except ValueError:
            print(f"Error: Required columns not found in {file_path}. Skipping...")
            return 0, 0, 0

        # Process each row
        for row in reader:
            if len(row) < len(headers):
                continue  # Skip malformed rows

            lineage = row[taxa_idx]
            reads = int(row[reads_idx])  # Get the number of reads for this row

            # Check if this row matches the specified function ID
            for cell in row[6:]:
                if cell and any(target_function in part for part in cell.replace(',', '|').split('|')):
                    if target_taxon in lineage:
                       # lineage = lineage.split("s__")[0]  # .split("|")[0]
                        reads_taxon_function += reads
                        total_reads_function += reads
                    else:
                        total_reads_function += reads
                    break  # Stop checking further functional columns for this row


            # # Check if the function is present in this row
            # function_present = any(target_function in cell for cell in row[6:] if cell)
            #
            # if function_present:
            #     total_reads_function += reads  # Add to total function reads across all taxa
            #
            #     # If the taxon is also present, add to taxon-specific function reads
            #     if target_taxon in lineage:
            #         reads_taxon_function += reads

    # Compute the proportion
    function_proportion = reads_taxon_function / total_reads_function if total_reads_function > 0 else 0

    return function_proportion, reads_taxon_function, total_reads_function


def main():
    parser = argparse.ArgumentParser(description='MetaPont: Extract Reads Proportions for a Specific Taxon and Function')
    parser.add_argument(
        "-d", "--directory", required=True,
        help="Directory containing TSV files to analyse."
    )
    parser.add_argument(
        "-t", "--taxon", required=True,
        help="Target taxon to search for (e.g., 'g__Escherichia')."
    )
    parser.add_argument(
        "-f", "--function", required=True,
        help="Target function to extract (e.g., 'EC:2.7.11.1')."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output file to save results."
    )

    options = parser.parse_args()
    print(f"Running extraction for Taxon: {options.taxon} and Function: {options.function}")

    input_path = os.path.abspath(options.directory)
    output_path = os.path.abspath(options.output)

    # Store results
    results = []

    # Process each TSV file in the directory
    for file_name in os.listdir(input_path):
        if file_name.endswith("_Final_Contig.tsv"):
            file_path = os.path.join(options.directory, file_name)
            print(f"Processing file: {file_name}")

            function_proportion, reads_taxon_function, total_reads_function = process_tsv_by_taxon_and_function(
                file_path, options.taxon, options.function
            )

            # Store the result
            results.append((file_name.replace("_Final_Contig.tsv", ""), function_proportion, reads_taxon_function, total_reads_function))

    # Write results to output file
    with open(output_path, "w") as out:
        out.write(f"Taxon: {options.taxon}, Function: {options.function}\n")
        out.write("Sample\tFunction Proportion\tReads in Taxon & Function\tTotal Reads in Function\n")
        for sample, function_proportion, reads_taxon_function, total_reads_function in results:
            out.write(f"{sample}\t{function_proportion:.6f}\t{reads_taxon_function}\t{total_reads_function}\n")

    print(f"Results saved to {options.output}")


if __name__ == "__main__":
    main()
