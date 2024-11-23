import collections
import re
import glob


outfile = open('/mnt/Internal/Nextcloud/Collab/Zhenbin/2024_project/Contig_Final_Outputs/Lineage_Genus.tsv', 'w')

KO_reads = collections.defaultdict(lambda: collections.defaultdict(int))

# Define the directory path containing the files
directory_path = '/mnt/Internal/Nextcloud/Collab/Zhenbin/2024_project/Contig_Final_Outputs/'



all_KOs = []

# Iterate through each file found
for file_path in files_list:
    sample = file_path.split('/')[-1].split('_Final')[0]
    # Perform actions on each file (for example, read the file content)
    with open(file_path, 'r') as file:
        for line in file:
            line = line.replace('\n','')
            seen_KOs = []
            if line.startswith('NODE'): # Only works for this dataset
                line_data = line.split('\t')
                mapped_reads = int(line_data[3])
                lineage = line_data[5]
                #for KO in KOs:
                #    if KO not in seen_KOs:
                #        KO_reads[sample][KO] += mapped_reads
                #        seen_KOs.append(KO)
                genus = lineage.rsplit('|', 1)[0]
                KO_reads[sample][genus] += mapped_reads

                if genus not in all_KOs:
                    all_KOs.append(genus)
print("")


def extract_numeric_suffix(key):
    match = re.search(r'\d+$', key)
    return int(match.group()) if match else 0

# Reorder the dictionary based on the last numeric characters of the keys
KO_reads = dict(sorted(KO_reads.items(), key=lambda x: extract_numeric_suffix(x[0])))




keys_as_string = '\t'.join(KO_reads.keys())


outfile.write('Lineage_Genus\t'+keys_as_string+'\n')

for KO in all_KOs:
    outfile.write(KO)
    for sample in KO_reads.keys():
        outfile.write('\t'+str(KO_reads[sample][KO]))
    outfile.write('\n')




















