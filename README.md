# MetaPont
**MetaPont**  - A tool to bridge the gap between the output of metagenomic tools and the analysis of the data

## Features - These are the current aims of this project  - Still under development

- **Targeted Functional Analysis:** Search for specific functional IDs (e.g., GO terms) within the `_Final_Contig.tsv` files provided by the HuwsLab Metagenome Workflow (https://github.com/TheHuwsLab/Metagenome_Workflow) .
- **Taxonomic Breakdown:** Extract genus-level taxonomy information and calculate their proportions in the dataset.
- **Batch Processing:** Analyse all `_Final_Contig.tsv` files in a specified directory.
- **Customisable Output:** Save results in a format suitable for downstream analysis.

---

## Installation

### Prerequisites

Ensure you have the following installed:

- Python ~3.6 or later

### Installation via pip

MetaPont is provided as a pip distribution. 

```bash
pip install MetaPont 
```

---

## Usage

--- 
### Extract-By-Function Command-line Arguments
```Extract-By-Function -h ``` 
```bash
usage: Extract_By_Function.py [-h] -d DIRECTORY -f FUNCTION_ID -o OUTPUT
                              [-m MIN_PROPORTION] [-top TOP_TAXA]

MetaPont v0.0.6: Extract-By-Function - Identify taxa contributing to a
specific function.

options:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Directory containing TSV files to analyse.
  -f FUNCTION_ID, --function_id FUNCTION_ID
                        Specific function ID to search for (e.g.,
                        'GO:0016597').
  -o OUTPUT, --output OUTPUT
                        Output file to save results.
  -m MIN_PROPORTION, --min_proportion MIN_PROPORTION
                        Minimum proportion threshold for taxa to be included
                        in the output.
  -top TOP_TAXA, --top_taxa TOP_TAXA
                        Top n taxa to be included in the output.

```

The `Extract-By-Function` tool provides several command-line options: \
Note: Either -m or -top is required.

| Option                   | Description                                                | Required | Default |
|--------------------------|------------------------------------------------------------|----------|---------|
| `-d`, `--directory`      | Directory containing `_Final_Contig.tsv` files to analyse. | Yes      | None    |
| `-f`, `--function_id`    | Functional ID to search for (e.g., `GO:0016597`).          | Yes      | None    |
| `-m`, `--min_proportion` | Minimum proportion needed for reporting.                   | Yes/No   | None    |
| `-top`, `--top_taxa`     | Number of taxa to report.                                  | Yes/No   | None    |
| `-o`, `--output`         | Output file name to save results.                          | Yes      | None    |

### Example

To search for the functional ID `GO:0016597` in all `_Final_Contig.tsv` files within the `test_data/` directory:

```bash
Extract=By-Function -d .../test_data/Final_contig/ -f GO:0016597 -top 3 -o .../test_data/Final_Contig/Extract_By_Function_Out/results.tsv
```

---

## Output

The tool generates a tab-delimited output file with the following columns:

1. **Sample:** Name of the processed Sample.
2. **Taxa:** Genus-level taxonomic assignment extracted from the `Lineage` column.
3. **Reads Assigned (Function):** Number of reads assigned to contigs with the given functional ID.
3. **Proportion:** Proportion of reads assigned to contigs of stated Taxa with the given functional ID within the sample.
4. **Proportion (Total Reads):** Proportion of reads assigned to contigs of stated Taxa with the given functional ID within the total reads of the sample.

Example output:

```
Function ID: GO:0016597
Sample	Taxa	Reads Assigned (Function)	Proportion (Function)	Proportion (Total Reads)
PN0536_0001_S1_Final_Contig.tsv	Lactobacillus	111963	0.602	0.004
PN0536_0003_S83_Final_Contig.tsv	Lactobacillus	20072	0.457	0.001
PN0536_0002_S2_Final_Contig.tsv	Acutalibacter	145222	0.795	0.005
PN0536_0004_S3_Final_Contig.tsv	Lactobacillus	40076	0.404	0.002
```

---



### Workflow - unfinished

1. The script reads `_Final_Contig.tsv` files from the specified directory.
2. For each file, it searches for occurrences of the given functional ID within specific columns.
3. Matches are associated with genus-level taxonomic information extracted from the `Lineage` column.
4. Taxa proportions are calculated and saved to the output file.

---
## Extract-By-Taxa Command-line Arguments
```Extract-By-Taxa -h ``` 
```bash
usage: Extract_By_Taxa.py [-h] -d DIRECTORY -t TAXON -o OUTPUT -func
                          FUNCTIONAL_CLASSES [-top TOP_FUNCTIONS]

MetaPont: Extract Top Functions by Taxon

options:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Directory containing TSV files to analyse.
  -t TAXON, --taxon TAXON
                        Target taxon to search for (e.g., 'g__Bacillus').
  -o OUTPUT, --output OUTPUT
                        Output file to save results.
  -func FUNCTIONAL_CLASSES, --functional_classes FUNCTIONAL_CLASSES
                        Which functional classes to report (e.g. GO,EC,KEGG
                        etc).
  -top TOP_FUNCTIONS, --top_functions TOP_FUNCTIONS
                        Top n functions to include in the output for each
                        sample (default: 3).

```

The `Extract-By-Taxa` tool provides several command-line options: 


| Option                            | Description                                                 | Required | Default |
|-----------------------------------|-------------------------------------------------------------|----------|---------|
| `-d`, `--directory`               | Directory containing `_Fincal_Contig.tsv` files to analyse. | Yes      | None    |
| `-t`, `--taxon`                   | Taxa to search for (e.g., `g__Bacillus`).                   | Yes      | None    |
| `-func`, `--functional_classes`   | Functional classes to report (e.g. GO,EC,KEGG etc).         | Yes      | None |
| `-top`, `--top_taxa`              | Number of functions to report (default 3).                  | No       | None    |
| `-o`, `--output`                  | Output file name to save results.                           | Yes      | None    |

### Example

To search for the top reported functions for taxon `g__Bacillus` in all `_Final_Contig.tsv` files within the `test_data/` directory:

```bash
Extract-By-Taxa -d .../test_data/Final_Contig -t g__Bacillus -o .../test_data/Final_Contig/Extract_By_Taxa/results.tsv  -func GO
```

---
## Output

The tool generates a tab-delimited output file with the following columns:

1. **Sample:** Name of the processed Sample.
2. **Function:** Reported 'top' function.
3. **Num of Assignments (Functions):** Number of times the function has been assigned across all contigs reported as chosen Taxon.

Example output:

```
Selected Taxon: g__Bacillus
Sample	Function	Num of Assignments
PN0536_0001_S1	GO:0008150	296
PN0536_0001_S1	GO:0003674	285
PN0536_0001_S1	GO:0005575	254
PN0536_0003_S83	GO:0005575	45
PN0536_0003_S83	GO:0008150	44
PN0536_0003_S83	GO:0003674	43
PN0536_0002_S2	GO:0005575	5
PN0536_0002_S2	GO:0008150	5
PN0536_0002_S2	GO:0005623	4
PN0536_0004_S3	GO:0008150	4
PN0536_0004_S3	GO:0003674	3
PN0536_0004_S3	GO:0005488	3

```

---

### Large File Handling (Might be a failure point)

The script uses `csv.field_size_limit` to handle exceptionally large `.tsv` files.

---


