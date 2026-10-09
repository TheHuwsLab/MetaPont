import collections
import re
import glob
import argparse
import os

try:
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *

def read_files(files_list):
    all_entries = []
    read_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    total_reads = collections.defaultdict(int)
    for file_path in files_list:
        sample = file_path.split('/')[-1].split('_Contigs')[0] # _Contigs is specific to this dataset
        with open(file_path, 'r') as file:
            for line in file:
                line = line.replace('\n','')
                if line.startswith('NODE'): # Only works for this dataset/metaspades
                    line_data = line.split('\t')
                    mapped_reads = int(line_data[3])
                    total_reads[sample] = int(line_data[2])
                    lineage = line_data[5]
                    genus = lineage.rsplit('|', 1)[0] # We are currently reporting only down to the genus level
                    read_counts[sample][genus] += mapped_reads
                    if genus not in all_entries:
                        all_entries.append(genus)
    return all_entries, read_counts, total_reads


def write_out(output_dir, read_counts, total_reads):
    for substr, counts in read_counts.items():
        safe_substr = substr
        if substr == 'd__unknown|k__unknown|p__unknown|c__unknown|o__unknown|f__unknown|g__unknown':
            safe_substr = 'd__unknown'
        output_file = os.path.join(output_dir, f"{safe_substr}_output.tsv")
        sample_names = sorted(total_reads.keys())
        keys_as_string = '\t'.join(sample_names)
        values_as_string = '\t'.join([str(total_reads[key]) for key in sample_names])

        all_current_taxa = sorted({key for sample_counts in counts.values() for key in sample_counts})

        with open(output_file, 'w') as outfile:
            outfile.write('Lineage_Genus\t' + keys_as_string + '\n')
            outfile.write('Total_Num_Reads' + '\t' + values_as_string + '\n')
            for current_taxa in all_current_taxa:
                outfile.write(current_taxa)
                for sample_key in sample_names:
                    outfile.write('\t' + str(counts[sample_key].get(current_taxa, 0)))
                outfile.write('\n')

    # Combined output across all separated taxa
    combined_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    for counts in read_counts.values():
        for sample, taxa_map in counts.items():
            for taxa, cnt in taxa_map.items():
                combined_counts[sample][taxa] += cnt

    combined_file = os.path.join(output_dir, "combined_separated_taxa_output.tsv")
    sample_names = sorted(total_reads.keys())
    keys_as_string = '\t'.join(sample_names)
    values_as_string = '\t'.join([str(total_reads[key]) for key in sample_names])
    all_combined_taxa = sorted({key for sample_counts in combined_counts.values() for key in sample_counts})

    with open(combined_file, 'w') as outfile:
        outfile.write('Lineage_Genus\t' + keys_as_string + '\n')
        outfile.write('Total_Num_Reads' + '\t' + values_as_string + '\n')
        for taxa in all_combined_taxa:
            outfile.write(taxa)
            for sample_key in sample_names:
                outfile.write('\t' + str(combined_counts[sample_key].get(taxa, 0)))
            outfile.write('\n')


def extract_numeric_suffix(key):
    match = re.search(r'\d+$', key)
    return int(match.group()) if match else 0

def group_by_substrings(read_counts, substrings, remove_substrings):
    grouped_counts = {substr: collections.defaultdict(lambda: collections.defaultdict(int)) for substr in substrings}
    for sample, counts in read_counts.items():
        for key, value in counts.items():
            if any(remove_substr in key for remove_substr in remove_substrings):
                continue
            for substr in substrings:
                if substr in key:
                    grouped_counts[substr][sample][key] += value
    return grouped_counts


def main():

    parser = argparse.ArgumentParser(description='MetaPont ' + MetaPont_Version + '- Reporter-Contig-Lineage:  Report contig lineage read counts across samples, grouping by specified taxa substrings.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-d', action='store', dest='dir_path', required=True,
                        help='Define the directory path containing the files')

    optional = parser.add_argument_group('Optional Arguments')
    parser.add_argument("--output", "-o",
                        help="Output CSV path (default: `root_dir/lineage_overview.csv`).")
    optional.add_argument('-s', '--separate-taxa', dest='separate_taxa', action='store',
                          help='Comma-separated list of taxa to separate (e.g. d__Bacteria,d__Archaea). If omitted, defaults are used.',
                          default=None)
    optional.add_argument('-r', '--remove-taxa', dest='remove_taxa', action='store',
                          help='Comma-separated list of taxa to remove. If omitted, defaults are used.',
                          default=None)

    options = parser.parse_args()

    # determine output path: if user provided one use it, otherwise put default file in root_dir
    if options.output:
        output_path = options.output
    else:
        output_path = options.dir_path # os.path.join(options.dir_path, "lineage_overview.csv")

    # Use glob to find files ending with '_Final_Output.tsv'
    files_list = glob.glob(f"{options.dir_path}/*_Contigs.tsv")
    all_entries, read_counts, total_reads = read_files(files_list)

    # Default lists
    _default_separate_taxa = ['d__Bacteria', 'd__Archaea', 'd__Eukaryota', 'd__Viruses', 'k__Fungi',
                              'd__unknown|k__unknown|p__unknown|c__unknown|o__unknown|f__unknown|g__unknown']
    _default_remove_taxa = [ # Remove mammalian and plant lineages by default
        'c__Mammalia', 'k__Viridiplantae'
        # 'd__Eukaryota|k__unknown|p__Evosea|c__Eumycetozoa|o__Dictyosteliales|f__Dictyosteliaceae|g__Dictyostelium',
        # 'd__Eukaryota|k__unknown|p__Euglenozoa|c__Kinetoplastea|o__Trypanosomatida|f__Trypanosomatidae|g__Leishmania',
        # 'd__Eukaryota|k__unknown|p__Apicomplexa|c__Conoidasida|o__Eucoccidiorida|f__Sarcocystidae|g__Toxoplasma',
        # 'd__Eukaryota|k__unknown|p__Apicomplexa|c__Aconoidasida|o__Haemosporida|f__Plasmodiidae|g__Plasmodium',
        # 'd__Eukaryota|k__unknown|p__Apicomplexa|c__Aconoidasida|o__Piroplasmida|f__Theileriidae|g__Theileria',
        # 'd__Eukaryota|k__unknown|p__Parabasalia|c__unknown|o__Trichomonadida|f__Trichomonadidae|g__Trichomonas',
        # 'd__Eukaryota|k__unknown|p__Apicomplexa|c__Conoidasida|o__Eucoccidiorida|f__Sarcocystidae|g__Besnoitia',
        # 'd__Eukaryota|k__unknown|p__Euglenozoa|c__Kinetoplastea|o__Trypanosomatida|f__Trypanosomatidae|g__Trypanosoma',
        # 'd__Eukaryota|k__unknown|p__Ciliophora|c__Oligohymenophorea|o__Peniculida|f__Parameciidae|g__Paramecium',
        # 'd__Eukaryota|k__Fungi|p__Microsporidia|c__unknown|o__unknown|f__Unikaryonidae|g__Encephalitozoon',
        # 'd__Eukaryota|k__unknown|p__unknown|c__Cryptophyceae|o__Cryptomonadales|f__Cryptomonadaceae|g__Cryptomonas'
    ]

    # Apply user-provided comma-separated values or fall back to defaults
    if options.separate_taxa:
        separate_taxa = [t.strip() for t in options.separate_taxa.split(',') if t.strip()]
    else:
        separate_taxa = _default_separate_taxa

    if options.remove_taxa:
        remove_taxa = [t.strip() for t in options.remove_taxa.split(',') if t.strip()]
    else:
        remove_taxa = _default_remove_taxa

    read_counts = group_by_substrings(read_counts, separate_taxa, remove_taxa)

    read_counts = dict(sorted(read_counts.items(), key=lambda x: extract_numeric_suffix(x[0])))

    write_out(output_path, read_counts, total_reads)




if __name__ == "__main__":
    main()
    print("Complete")



