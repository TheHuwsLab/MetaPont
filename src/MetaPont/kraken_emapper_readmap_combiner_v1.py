import os
from collections import defaultdict
import re
import sys

################

parent_directory_path = sys.argv[1]


for subdir, dirs, files in os.walk(parent_directory_path, topdown=True):
    for dir in dirs:
        try:
            current_dir = os.path.basename(dir)
            parent_of_current_dir = os.path.basename(subdir)
            if current_dir.startswith('PN') and parent_of_current_dir == os.path.basename(parent_directory_path):
                print(current_dir)
                specific_dir_path = os.path.join(subdir, dir)
                first_file_path = os.path.join(specific_dir_path, f"{dir}_eggnog_mapper/{dir}_pyrodigal_eggnog_mapped.emapper.annotations")
                second_file_path = os.path.join(specific_dir_path, f"{dir}_kraken2/{dir}_kraken2_report.txt")
                third_file_path = os.path.join(specific_dir_path, f"{dir}_kraken2/{dir}_kraken2_report_mpa.txt")
                forth_file_path = os.path.join(specific_dir_path, f"{dir}_kraken2/{dir}_kraken2.txt")
                fifth_file_path = os.path.join(specific_dir_path, f"{dir}_readmapped/{dir}_readmapped_contig_summary.txt")

                if os.path.exists(first_file_path) and os.path.exists(second_file_path) and os.path.exists(third_file_path) and os.path.exists(forth_file_path) and os.path.exists(fifth_file_path):
                    print(f"Master directory: {specific_dir_path}")
                    print(f"First file path: {first_file_path}")
                    print(f"Second file path: {second_file_path}")
                    print(f"Third file path: {third_file_path}")
                    print(f"Forth file path: {forth_file_path}")
                    print(f"Fifth file path: {fifth_file_path}")

                contigs = defaultdict(list)
                genes = defaultdict(list)
                lineages = {}
                taxa_ids = []

                with open(first_file_path, 'r') as emapper_in:
                    for line in emapper_in:
                        if not line.startswith('#'):
                            line_data = line.strip().split('\t')
                            contig = line_data[0].rsplit('_', 1)[0]
                            COGs = line_data[6]
                            GOs = line_data[9]
                            EC = line_data[10]
                            KEGG_ko = line_data[11]
                            KEGG_Pathway = line_data[12]
                            KEGG_Module = line_data[13]
                            KEGG_Reaction = line_data[14]
                            KEGG_rclass = line_data[15]
                            BRITE = line_data[16]
                            KEGG_TC = line_data[17]
                            CAZy = line_data[18]
                            BiGG_Reaction = line_data[19]
                            PFAMs = line_data[20]
                            genes[contig].append(
                               [COGs, GOs, EC, KEGG_ko, KEGG_Pathway, KEGG_Module, KEGG_Reaction, KEGG_rclass, BRITE,
                                KEGG_TC, CAZy, BiGG_Reaction, PFAMs])

                with open(second_file_path, 'r') as f:
                    kraken2_lines = f.readlines()


                with open(third_file_path, 'r') as kraken_lineage_in:
                    for line in kraken_lineage_in:
                        lineages_to_check = {'d__':'d__unknown', 'k__':'k__unknown', 'p__':'p__unknown', 'c__':'c__unknown',
                                             'o__':'o__unknown', 'f__':'f__unknown', 'g__':'g__unknown', 's__':'s__unknown'}
                        if line.startswith('d__'):
                            line_data = line.split('|')
                            contig_length = line_data[-1].split('\t')[1].strip()  # check
                            line_data[-1] = line_data[-1].split('\t')[0] #check
                            # Check for missing lineages
                            for data in line_data:
                                first_three_chars = data[:3]
                                if first_three_chars in lineages_to_check:
                                    lineages_to_check[first_three_chars] = data

                            combined = ''
                            for key, value in lineages_to_check.items():
                                #if value is not None and value != '' and '__unknown' not in value:
                                combined+=value+'|'
                            combined = combined[:-1]
                            #print(combined)

                            lineages[line.split('\t')[0].split('_')[-1]] = combined

                counter = 0
                with open(forth_file_path, 'r') as kraken_in:
                    for line in kraken_in:
                        counter +=1
                        line_data = line.strip().split('\t')
                        contig = line_data[1]
                        taxa_id = line_data[2].split(' (taxid')[0]
                        contig_length = int(line_data[3])
                        if contig_length >= 2500:
                            full_lineage = None
                            if taxa_id == 'unclassified' or taxa_id == 'root':
                                full_lineage = 'd__unknown|k__unknown|p__unknown|c__unknown|o__unknown|f__unknown|g__unknown|s__unknown'
                            else:
                                try:
                                    full_lineage = lineages[taxa_id]
                                except KeyError:
                                    while True:
                                        taxa_id = taxa_id.rsplit(' ', 1)[0]
                                        try:
                                            full_lineage = lineages[taxa_id]
                                            break
                                        except KeyError:
                                            tmp_size = taxa_id.split(' ')
                                            if len(tmp_size) == 1:
                                                taxa_id = line_data[2].split(' (taxid')[0]
                                                escaped_pattern = re.escape(taxa_id)
                                                matching_line = next((line for line in kraken2_lines if re.search(escaped_pattern, line)), None)
                                                # Loop through the lines and search for the pattern

                                                index = kraken2_lines.index(matching_line)
                                                line_above = kraken2_lines[index - 1]
                                                while True:
                                                    if not any(char.isdigit() for char in line_above.split('\t')[3]):
                                                        break
                                                    else:
                                                        index = kraken2_lines.index(line_above)
                                                        line_above = kraken2_lines[index - 1]
                                                taxa_id = line_above.strip().split('\t')[-1].lstrip()
                                                try:
                                                    full_lineage = lineages[taxa_id]
                                                    break
                                                except KeyError:
                                                    if taxa_id == 'unclassified' or taxa_id == 'root':
                                                        full_lineage = 'd__unknown|k__unknown|p__unknown|c__unknown|o__unknown|f__unknown|g__unknown|s__unknown'
                                                    break
                                            continue
                                    if full_lineage == None:
                                        full_lineage = 'd__unknown|k__unknown|p__unknown|c__unknown|o__unknown|f__unknown|g__unknown|s__unknown'
                                        print("Contig without lineage. ")
                            contigs[contig] = [line_data[2],line_data[3],full_lineage]
                            #print(counter)
                        else:
                            print("Contig length less than 2,500")
                            break




                print("RM")
                with open(fifth_file_path, 'r') as readmapped_in:
                    for line in readmapped_in:
                        if not line.startswith('Contig'): #sloppy
                            line_data = line.strip().split('\t')
                            contig = line_data[0]
                            if contig in contigs:
                                info = [line_data[1],line_data[2]]
                                info.extend(contigs[contig])
                                contigs[contig] = info

                final_dict = defaultdict(list)

                for contig, data in contigs.items():  # Easy to read but not efficient at all
                    data.extend([[], [], [], [], [], [], [], [], [], [], [], [], []])
                    final_dict[contig] = data
                    #for gene in data[0]:
                    gene_info = genes[contig]
                    for current in gene_info:

                        try:
                            if current[0] != '-':
                                final_dict[contig][5].append(current[0])
                            if current[1] != '-':
                                final_dict[contig][6].append(current[1])
                            if current[2] != '-':
                                final_dict[contig][7].append(current[2])
                            if current[3] != '-':
                                final_dict[contig][8].append(current[3])
                            if current[4] != '-':
                                final_dict[contig][9].append(current[4])
                            if current[5] != '-':
                                final_dict[contig][10].append(current[5])
                            if current[6] != '-':
                                final_dict[contig][11].append(current[6])
                            if current[7] != '-':
                                final_dict[contig][12].append(current[7])
                            if current[8] != '-':
                                final_dict[contig][13].append(current[8])
                            if current[9] != '-':
                                final_dict[contig][14].append(current[9])
                            if current[10] != '-':
                                final_dict[contig][15].append(current[10])
                            if current[11] != '-':
                                final_dict[contig][16].append(current[11])
                            if current[12] != '-':
                                final_dict[contig][17].append(current[12])
                            if current[13] != '-':
                                final_dict[contig][18].append(current[13])
                        except IndexError:
                            continue

                output_path = os.path.join(parent_directory_path, 'Contig_Final_Outputs')
                if not os.path.exists(output_path):
                    os.makedirs(output_path)
                output = open(os.path.join(output_path, dir + '_Final_Contig.tsv'),'w')


                output.write(dir+'\n')
                output.write('Contig\tContig_Length\tTotal_Reads\tMapped_Reads\tScientific_Name\tLineage\tCOG\tGO\tEC\tKEGG_KO\tKEGG_Pathway\tKEGG_Module\t'
                             'KEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n')

                for contig, data in final_dict.items():
                    try:
                        #print(contig)
                        output.write(contig+'\t'+data[3]+'\t'+data[1]+'\t'+data[0]+'\t'+data[2]+'\t'+data[4]+'\t'+'|'.join(data[5])+'\t'+'|'.join(data[6])+'\t'+'|'.join(data[7])+'\t'+'|'.join(data[8])
                                     +'\t'+'|'.join(data[9])+'\t'+'|'.join(data[10])+'\t'+'|'.join(data[11])+'\t'+'|'.join(data[12])+'\t'+'|'.join(data[13])+'\t'+'|'.join(data[14])
                                     +'\t'+'|'.join(data[15])+'\t'+'|'.join(data[16])+'\t'+'|'.join(data[17])+'\n')
                    except TypeError as e:
                        print("Some Error: " + e)

                print("Done")


        except FileNotFoundError as e:
            print(f"Error processing directory {dir}: {e}")


