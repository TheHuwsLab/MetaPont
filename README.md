# MetaPont
**MetaPont**  - A tool to bridge the gap between the output of metagenomic tools and the analysis of the data

## Features - These are the current aims of this project  - Still under development

- **Targeted Functional Analysis:** Search for specific functional IDs (e.g., GO terms) within the `.tsv` files provided by the HuwsLab Metagenome Workflow (https://github.com/TheHuwsLab/Metagenome_Workflow) .
- **Taxonomic Breakdown:** Extract genus-level taxonomy information and calculate their proportions in the dataset.
- **Batch Processing:** Analyse all `.tsv` files in a specified directory.
- **Customisable Output:** Save results in a format suitable for downstream analysis.

---

## Installation

### Prerequisites

Ensure you have the following installed:

- Python ~3.6 or later
- Required Python libraries: `argparse`, `csv`, and `collections` (standard libs).

### Installation via pip

MetaPont is provided as a pip distribution. 

```bash
pip install MetaPont 
```

---

## Usage

### Command-line Arguments
```Extract-By-Function -h ``` 
```bash
usage: Extract-By-Function [-h] -d DIRECTORY -f FUNCTION_ID [-o OUTPUT] [-m MIN_PROPORTION]

MetaPont v0.0.2: Extract-By-Function - Identify taxa contributing to a specific function.

options:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Directory containing TSV files to analyse.
  -f FUNCTION_ID, --function_id FUNCTION_ID
                        Specific function ID to search for (e.g., 'GO:0002').
  -o OUTPUT, --output OUTPUT
                        Output file to save results (default: output_taxa_details.tsv).
  -m MIN_PROPORTION, --min_proportion MIN_PROPORTION
                        Minimum proportion threshold for taxa to be included in the output (default: 0.05).
```

The `Extract-By-Function` tool provides several command-line options:

| Option                   | Description                                   | Required | Default                       |
|--------------------------|-----------------------------------------------|----------|-------------------------------|
| `-d`, `--directory`      | Directory containing `.tsv` files to analyse. | Yes      | None                          |
| `-f`, `--function_id`    | Functional ID to search for (e.g., `GO:0002`). | Yes      | None                          |
| `-m`, `--min_proportion` | Minimum proportion needed for reporting.      | Yes      | 0.05 (5%)                     |
| `-o`, `--output`         | Output file name to save results.             | No       | `output_taxa_proportions.tsv` |

### Example

To search for the functional ID `GO:0002` in all `.tsv` files within the `data/` directory:

```bash
ExtractByFunction -d .../test_data/Final_contig/ -f GO:0002 -m 0.10 -o .../test_data/Final_Contig/Extract_By_Function_Out/results.tsv
```

---

## Output

The tool generates a tab-delimited output file with the following columns:

1. **Sample:** Name of the processed `.tsv` file.
2. **Taxa:** Genus-level taxonomic assignment extracted from the `Lineage` column.
3. **Proportion:** Proportion of matches to the given functional ID within the sample.

Example output:

```
Function ID: GO:0002
Sample	Taxa	Reads Assigned (Function)	Proportion (Function)	Proportion (Total Reads)
PN0536_0003_S83_Final_Contig.tsv	Gordonibacter	60788	0.075	0.002
PN0536_0003_S83_Final_Contig.tsv	Streptomyces	115671	0.142	0.004
PN0536_0003_S83_Final_Contig.tsv	unknown	80890	0.099	0.003
PN0536_0003_S83_Final_Contig.tsv	Clostridium	51018	0.063	0.002
PN0536_0003_S83_Final_Contig.tsv	Lactobacillus	149909	0.184	0.005
PN0536_0003_S83_Final_Contig.tsv	Limosilactobacillus	79694	0.098	0.003
```

---

## Implementation Details

### Workflow

1. The script reads `.tsv` files from the specified directory.
2. For each file, it searches for occurrences of the given functional ID within specific columns.
3. Matches are associated with genus-level taxonomic information extracted from the `Lineage` column.
4. Taxa proportions are calculated and saved to the output file.

### Large File Handling (Might be a failure point)

The script uses `csv.field_size_limit` to handle exceptionally large `.tsv` files.

---

## Future Plans

- Add support for additional file formats (e.g., `.csv`, `.txt`).
- Expand functionality for more complex taxonomic and functional analyses.
---

