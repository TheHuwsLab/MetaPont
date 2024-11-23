import os
from collections import defaultdict
import re
import sys


################

parent_directory_path = sys.argv[1]

all_COGs = []
all_GOs = []
all_ECs = []
all_KEGG_kos = []
all_KEGG_Pathways = []
all_KEGG_Modules = []
all_KEGG_Reactions = []
all_KEGG_rclasses = []
all_BRITEs = []
all_KEGG_TCs = []
all_CAZys =[]
all_BiGG_Reactions = []
all_PFAMs = []

functions = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))


for subdir, dirs, files in os.walk(parent_directory_path, topdown=True):
    for dir in dirs:
        try:
            current_dir = os.path.basename(dir)
            parent_of_current_dir = os.path.basename(subdir)
            if current_dir.startswith('PN') and parent_of_current_dir == os.path.basename(parent_directory_path):
                print(current_dir)
                specific_dir_path = os.path.join(subdir, dir)
                first_file_path = os.path.join(specific_dir_path, f"{dir}_readmapped/{dir}_readmapped_cds_summary.txt")
                second_file_path = os.path.join(specific_dir_path, f"{dir}_eggnog_mapper/{dir}_pyrodigal_eggnog_mapped.emapper.annotations")

                if os.path.exists(first_file_path) and os.path.exists(second_file_path):
                    print(f"Master directory: {specific_dir_path}")
                    print(f"First file path: {first_file_path}")
                    print(f"Second file path: {second_file_path}")

                genes = defaultdict(int)
                with open(first_file_path, 'r') as readmap_cds_in:
                    for line in readmap_cds_in:
                        reads = line.split()[0]
                        gene = line.split()[1].replace('ID=','')
                        contig_length = int(gene.split('length_')[1].split('_cov')[0])
                        if contig_length >= 2500:
                            genes.update({gene:reads})




                lineages = {}
                taxa_ids = []

                with open(second_file_path, 'r') as emapper_in:
                    for line in emapper_in:
                        if not line.startswith('#'):
                            line_data = line.strip().split('\t')
                            gene = line_data[0]
                            if gene in genes:
                                gene_read_count = int(genes[gene])

                                COGs = line_data[6]
                                if COGs != '-':
                                    #cog_list = COGs.split(',')
                                    cog_list = [cog for cog in COGs]
                                    for cog in cog_list:
                                        functions[dir]['COG'][cog] += gene_read_count
                                        if cog not in all_COGs:
                                            all_COGs.append(cog)

                                GOs = line_data[9]
                                if GOs != '-':
                                    go_list = GOs.split(',')
                                    for go in go_list:
                                        functions[dir]['GO'][go] += gene_read_count
                                        if go not in all_GOs:
                                            all_GOs.append(go)

                                EC = line_data[10]
                                if EC != '-':
                                    ec_list = EC.split(',')
                                    for ec in ec_list:
                                        functions[dir]['EC'][ec] += gene_read_count
                                        if ec not in all_ECs:
                                            all_ECs.append(ec)

                                KEGG_ko = line_data[11]
                                if KEGG_ko != '-':
                                    kegg_ko_list = KEGG_ko.split(',')
                                    for kegg_ko in kegg_ko_list:
                                        functions[dir]['KEGG_ko'][kegg_ko] += gene_read_count
                                        if kegg_ko not in all_KEGG_kos:
                                            all_KEGG_kos.append(kegg_ko)

                                KEGG_Pathway = line_data[12]
                                if KEGG_Pathway != '-':
                                    kegg_pathway_list = KEGG_Pathway.split(',')
                                    for kegg_pathway in kegg_pathway_list:
                                        functions[dir]['KEGG_Pathway'][kegg_pathway] += gene_read_count
                                        if kegg_pathway not in all_KEGG_Pathways:
                                            all_KEGG_Pathways.append(kegg_pathway)

                                KEGG_Module = line_data[13]
                                if KEGG_Module != '-':
                                    kegg_module_list = KEGG_Module.split(',')
                                    for kegg_module in kegg_module_list:
                                        functions[dir]['KEGG_Module'][kegg_module] += gene_read_count
                                        if kegg_module not in all_KEGG_Modules:
                                            all_KEGG_Modules.append(kegg_module)

                                KEGG_Reaction = line_data[14]
                                if KEGG_Reaction != '-':
                                    kegg_reaction_list = KEGG_Reaction.split(',')
                                    for kegg_reaction in kegg_reaction_list:
                                        functions[dir]['KEGG_Reaction'][kegg_reaction] += gene_read_count
                                        if kegg_reaction not in all_KEGG_Reactions:
                                            all_KEGG_Reactions.append(kegg_reaction)

                                KEGG_rclass = line_data[15]
                                if KEGG_rclass != '-':
                                    kegg_rclass_list = KEGG_rclass.split(',')
                                    for kegg_rclass in kegg_rclass_list:
                                        functions[dir]['KEGG_rclass'][kegg_rclass] += gene_read_count
                                        if kegg_rclass not in all_KEGG_rclasses:
                                            all_KEGG_rclasses.append(kegg_rclass)

                                BRITE = line_data[16]
                                if BRITE != '-':
                                    brite_list = BRITE.split(',')
                                    for brite in brite_list:
                                        functions[dir]['BRITE'][brite] += gene_read_count
                                        if brite not in all_BRITEs:
                                            all_BRITEs.append(brite)

                                KEGG_TC = line_data[17]
                                if KEGG_TC != '-':
                                    kegg_tc_list = KEGG_TC.split(',')
                                    for kegg_tc in kegg_tc_list:
                                        functions[dir]['KEGG_TC'][kegg_tc] += gene_read_count
                                        if kegg_tc not in all_KEGG_TCs:
                                            all_KEGG_TCs.append(kegg_tc)

                                CAZy = line_data[18]
                                if CAZy != '-':
                                    cazy_list = CAZy.split(',')
                                    for cazy in cazy_list:
                                        functions[dir]['CAZy'][cazy] += gene_read_count
                                        if cazy not in all_CAZys:
                                            all_CAZys.append(cazy)

                                BiGG_Reaction = line_data[19]
                                if BiGG_Reaction != '-':
                                    bigg_reaction_list = BiGG_Reaction.split(',')
                                    for bigg_reaction in bigg_reaction_list:
                                        functions[dir]['BiGG_Reaction'][bigg_reaction] += gene_read_count
                                        if bigg_reaction not in all_BiGG_Reactions:
                                            all_BiGG_Reactions.append(bigg_reaction)

                                PFAMs = line_data[20]
                                if PFAMs != '-':
                                    pfams_list = PFAMs.split(',')
                                    for pfam in pfams_list:
                                        functions[dir]['PFAM'][pfam] += gene_read_count
                                        if pfam not in all_PFAMs:
                                            all_PFAMs.append(pfam)


                print("Sample " + dir + " Done")



        except FileNotFoundError as e:
            print(e)

functions = dict(sorted(functions.items()))


def write_transposed_output(output_path, category, all_items, functions):
    output = open(output_path, 'w')
    output.write('\t' + '\t'.join(functions.keys()) + '\n')

    for item in all_items:
        row_values = [item] + [str(functions[sample][category][item]) for sample in functions.keys()]
        output.write('\t'.join(row_values) + '\n')

    output.close()

output_path = os.path.join(parent_directory_path, 'CDS_Final_Outputs')
if not os.path.exists(output_path):
    os.makedirs(output_path)

# COG output
write_transposed_output(output_path + '/Final_CDS_COG.tsv', 'COG', all_COGs, functions)

# GO output
write_transposed_output(output_path + '/Final_CDS_GO.tsv', 'GO', all_GOs, functions)

# EC output
write_transposed_output(output_path + '/Final_CDS_EC.tsv', 'EC', all_ECs, functions)

# KEGG_ko output
write_transposed_output(output_path + '/Final_CDS_KEGG_ko.tsv', 'KEGG_ko', all_KEGG_kos, functions)

# KEGG_Pathway output
write_transposed_output(output_path + '/Final_CDS_KEGG_Pathway.tsv', 'KEGG_Pathway', all_KEGG_Pathways, functions)

# KEGG_Module output
write_transposed_output(output_path + '/Final_CDS_KEGG_Module.tsv', 'KEGG_Module', all_KEGG_Modules, functions)

# KEGG_Reaction output
write_transposed_output(output_path + '/Final_CDS_KEGG_Reaction.tsv', 'KEGG_Reaction', all_KEGG_Reactions, functions)

# KEGG_rclass output
write_transposed_output(output_path + '/Final_CDS_KEGG_rclass.tsv', 'KEGG_rclass', all_KEGG_rclasses, functions)

# BRITE output
write_transposed_output(output_path + '/Final_CDS_BRITE.tsv', 'BRITE', all_BRITEs, functions)

# KEGG_TC output
write_transposed_output(output_path + '/Final_CDS_KEGG_TC.tsv', 'KEGG_TC', all_KEGG_TCs, functions)

# CAZy output
write_transposed_output(output_path + '/Final_CDS_CAZy.tsv', 'CAZy', all_CAZys, functions)

# BiGG_Reaction output
write_transposed_output(output_path + '/Final_CDS_BiGG_Reaction.tsv', 'BiGG_Reaction', all_BiGG_Reactions, functions)

# PFAMs output
write_transposed_output(output_path + '/Final_CDS_PFAMs.tsv', 'PFAM', all_PFAMs, functions)

