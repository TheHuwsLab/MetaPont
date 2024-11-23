import collections
import re
import glob
import argparse
import os


def read_files(files_list):
    all_entries = []
    read_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    total_reads = collections.defaultdict(int)
    for file_path in files_list:
        sample = file_path.split('/')[-1].split('_Final')[0] # _Final is specific to this dataset
        with open(file_path, 'r') as file:
            for line in file:
                line = line.replace('\n','')
                if line.startswith('NODE'): # Only works for this dataset/metaspades
                    line_data = line.split('\t')
                    mapped_reads = int(line_data[3])
                    total_reads[sample] = int(line_data[2])
                    lineage = line_data[5]
                    genus = lineage.rsplit('|', 1)[0] # We are reporting only down to the genus level
                    read_counts[sample][genus] += mapped_reads
                    if genus not in all_entries:
                        all_entries.append(genus)
    return all_entries, read_counts, total_reads


def write_out(output_dir, read_counts, total_reads):
    for substr, counts in read_counts.items():
        if substr == 'd__unknown|k__unknown|p__unknown|c__unknown|o__unknown|f__unknown|g__unknown':
            substr = 'd__unknown'
        output_file = os.path.join(output_dir, f"{substr}_output.tsv")
        sample_names = sorted(total_reads.keys())
        keys_as_string = '\t'.join(sample_names)
        values_as_string = '\t'.join([str(total_reads[key]) for key in sample_names])


        all_current_taxa = list({key for sample_counts in counts.values() for key in sample_counts})



        with open(output_file, 'w') as outfile:
            outfile.write('Lineage_Genus\t'+keys_as_string+'\n')
            outfile.write('Total_Num_Reads')
            outfile.write('\t' + str(values_as_string))
            outfile.write('\n')
            for current_taxa in all_current_taxa:
                outfile.write(current_taxa)
                for sample_key in sample_names:
                    outfile.write('\t' + str(counts[sample_key][current_taxa]))
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

    parser = argparse.ArgumentParser(description='....')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-d', action='store', dest='dir_path', required=True,
                        help='Define the directory path containing the files')
    required.add_argument('-o', action='store', dest='output', help='Outdir',
                        required=True)

    options = parser.parse_args()

    # Use glob to find files ending with '_Final_Output.tsv'
    files_list = glob.glob(f"{options.dir_path}/*_Final_*.tsv")
    all_entries, read_counts, total_reads = read_files(files_list)

    separate_taxa = ['d__Bacteria', 'd__Archaea', 'd__Eukaryota', 'd__Viruses','k__Fungi',
                  'd__unknown|k__unknown|p__unknown|c__unknown|o__unknown|f__unknown|g__unknown']
    remove_taxa = ['c__Mammalia','k__Viridiplantae'
                        ,'d__Eukaryota|k__unknown|p__Evosea|c__Eumycetozoa|o__Dictyosteliales|f__Dictyosteliaceae|g__Dictyostelium'
                        ,'d__Eukaryota|k__unknown|p__Euglenozoa|c__Kinetoplastea|o__Trypanosomatida|f__Trypanosomatidae|g__Leishmania'
                        ,'d__Eukaryota|k__unknown|p__Apicomplexa|c__Conoidasida|o__Eucoccidiorida|f__Sarcocystidae|g__Toxoplasma'
                        ,'d__Eukaryota|k__unknown|p__Apicomplexa|c__Aconoidasida|o__Haemosporida|f__Plasmodiidae|g__Plasmodium'
                        ,'d__Eukaryota|k__unknown|p__Apicomplexa|c__Aconoidasida|o__Piroplasmida|f__Theileriidae|g__Theileria'
                        ,'d__Eukaryota|k__unknown|p__Parabasalia|c__unknown|o__Trichomonadida|f__Trichomonadidae|g__Trichomonas'
                        ,'d__Eukaryota|k__unknown|p__Apicomplexa|c__Conoidasida|o__Eucoccidiorida|f__Sarcocystidae|g__Besnoitia'
                        ,'d__Eukaryota|k__unknown|p__Euglenozoa|c__Kinetoplastea|o__Trypanosomatida|f__Trypanosomatidae|g__Trypanosoma'
                        ,'d__Eukaryota|k__unknown|p__Ciliophora|c__Oligohymenophorea|o__Peniculida|f__Parameciidae|g__Paramecium'
                        ,'d__Eukaryota|k__Fungi|p__Microsporidia|c__unknown|o__unknown|f__Unikaryonidae|g__Encephalitozoon'
                        ,'d__Eukaryota|k__unknown|p__unknown|c__Cryptophyceae|o__Cryptomonadales|f__Cryptomonadaceae|g__Cryptomonas']
    #make a overview list of samples and report 0s even if group has none for that sample
    #output total number of reads


    read_counts = group_by_substrings(read_counts, separate_taxa, remove_taxa)

    read_counts = dict(sorted(read_counts.items(), key=lambda x: extract_numeric_suffix(x[0])))

    write_out(options.output, read_counts, total_reads)




if __name__ == "__main__":
    main()
    print("Complete")



