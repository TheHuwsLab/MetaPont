import os
from collections import defaultdict
import argparse
import re



################

def processFunction(input_dir, ex_taxon, minlen, prefix):

    class all_Funcs:
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
        all_CAZys = []
        all_BiGG_Reactions = []
        all_PFAMs = []


    functions = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))


    for subdir, dirs, files in os.walk(input_dir, topdown=True):
        for dir in dirs:
            try:
                current_dir = os.path.basename(dir)
                parent_of_current_dir = os.path.basename(subdir)
                if current_dir.startswith(prefix) and parent_of_current_dir == os.path.basename(input_dir):
                    print(current_dir)
                    specific_dir_path = os.path.join(subdir, dir)
                    first_file_path = os.path.join(specific_dir_path, f"{dir}_readmapped/{dir}_readmapped_cds_summary.txt")
                    second_file_path = os.path.join(specific_dir_path, f"{dir}_kraken2/{dir}_kraken2_report.txt")
                    third_file_path = os.path.join(specific_dir_path, f"{dir}_kraken2/{dir}_kraken2_report_mpa.txt")
                    forth_file_path = os.path.join(specific_dir_path, f"{dir}_kraken2/{dir}_kraken2.txt")
                    fifth_file_path = os.path.join(specific_dir_path,
                                                    f"{dir}_eggnog_mapper/{dir}_pyrodigal_eggnog_mapped.emapper.annotations")

                    if os.path.exists(first_file_path) and os.path.exists(second_file_path) and os.path.exists(third_file_path):
                        print(f"Master directory: {specific_dir_path}")
                        print(f"First file path: {first_file_path}")
                        print(f"Second file path: {second_file_path}")
                        print(f"Third file path: {third_file_path}")
                        print(f"Forth file path: {forth_file_path}")
                        print(f"Fifth file path: {fifth_file_path}")

                    ##############
                    genes = defaultdict(int)
                    with open(first_file_path, 'r') as readmap_cds_in:
                        for line in readmap_cds_in:
                            reads = line.split()[0]
                            gene = line.split()[1].replace('ID=','')
                            contig_length = int(gene.split('length_')[1].split('_cov')[0])

                            if contig_length >= int(minlen):
                                genes.update({gene:reads})

                    #############
                    lineages = {}
                    contigs = defaultdict(list)

                    ex_taxon_list = ex_taxon.split(',')

                    with open(second_file_path, 'r') as f:
                        kraken2_lines = f.readlines()

                    with open(third_file_path, 'r') as kraken_lineage_in:
                        for line in kraken_lineage_in:
                            lineages_to_check = {'d__': 'd__unknown', 'k__': 'k__unknown', 'p__': 'p__unknown',
                                                 'c__': 'c__unknown',
                                                 'o__': 'o__unknown', 'f__': 'f__unknown', 'g__': 'g__unknown',
                                                 's__': 's__unknown'}
                            if line.startswith('d__'):
                                line_data = line.split('|')
                                line_data[-1] = line_data[-1].split('\t')[0]  # check
                                # Check for missing lineages
                                for data in line_data:
                                    first_three_chars = data[:3]
                                    if first_three_chars in lineages_to_check:
                                        lineages_to_check[first_three_chars] = data

                                combined = ''
                                for key, value in lineages_to_check.items():
                                    # if value is not None and value != '' and '__unknown' not in value:
                                    combined += value + '|'
                                combined = combined[:-1]
                                # print(combined)

                                lineages[line.split('\t')[0].split('_')[-1]] = combined

                    counter = 0
                    with open(forth_file_path, 'r') as kraken_in:
                        for line in kraken_in:
                            counter += 1
                            line_data = line.strip().split('\t')
                            contig = line_data[1]
                            taxa_id = line_data[2].split(' (taxid')[0]
                            contig_length = int(line_data[3])
                            if contig_length >= minlen:
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
                                                    matching_line = next((line for line in kraken2_lines if
                                                                          re.search(escaped_pattern, line)), None)
                                                    # Loop through the lines and search for the pattern

                                                    index = kraken2_lines.index(matching_line)
                                                    line_above = kraken2_lines[index - 1]
                                                    while True:
                                                        if not any(
                                                                char.isdigit() for char in line_above.split('\t')[3]):
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
                                contigs[contig] = [line_data[2], line_data[3], full_lineage]
                                # print(counter)
                            else:
                                print("Contig length less than 2,500")
                                break

                    ###############

                    contigs_with_excluded_taxa = []
                    for contig, data in contigs.items():
                        # Check if any element in first_column_elements is in ex_taxon_list
                        if any(element in ex_taxon_list for element in data[2].split('|')):
                            contigs_with_excluded_taxa.append(contig)

                    ###############


                    with open(fifth_file_path, 'r') as emapper_in:
                        for line in emapper_in:
                            if not line.startswith('#'):
                                line_data = line.strip().split('\t')
                                current_contig = line_data[0].rsplit('_', 1)[0]
                                if current_contig in contigs_with_excluded_taxa: # This is where excluded taxa are removed
                                    continue
                                gene = line_data[0]
                                if gene in genes:
                                    gene_read_count = int(genes[gene])

                                    COGs = line_data[6]
                                    if COGs != '-':
                                        cog_list = [cog for cog in COGs]
                                        for cog in cog_list:
                                            functions[dir]['COG'][cog] += gene_read_count
                                            if cog not in all_Funcs.all_COGs:
                                                all_Funcs.all_COGs.append(cog)

                                    GOs = line_data[9]
                                    if GOs != '-':
                                        go_list = GOs.split(',')
                                        for go in go_list:
                                            functions[dir]['GO'][go] += gene_read_count
                                            if go not in all_Funcs.all_GOs:
                                                all_Funcs.all_GOs.append(go)

                                    EC = line_data[10]
                                    if EC != '-':
                                        ec_list = EC.split(',')
                                        for ec in ec_list:
                                            functions[dir]['EC'][ec] += gene_read_count
                                            if ec not in all_Funcs.all_ECs:
                                                all_Funcs.all_ECs.append(ec)

                                    KEGG_ko = line_data[11]
                                    if KEGG_ko != '-':
                                        kegg_ko_list = KEGG_ko.split(',')
                                        for kegg_ko in kegg_ko_list:
                                            functions[dir]['KEGG_ko'][kegg_ko] += gene_read_count
                                            if kegg_ko not in all_Funcs.all_KEGG_kos:
                                                all_Funcs.all_KEGG_kos.append(kegg_ko)

                                    KEGG_Pathway = line_data[12]
                                    if KEGG_Pathway != '-':
                                        kegg_pathway_list = KEGG_Pathway.split(',')
                                        for kegg_pathway in kegg_pathway_list:
                                            functions[dir]['KEGG_Pathway'][kegg_pathway] += gene_read_count
                                            if kegg_pathway not in all_Funcs.all_KEGG_Pathways:
                                                all_Funcs.all_KEGG_Pathways.append(kegg_pathway)

                                    KEGG_Module = line_data[13]
                                    if KEGG_Module != '-':
                                        kegg_module_list = KEGG_Module.split(',')
                                        for kegg_module in kegg_module_list:
                                            functions[dir]['KEGG_Module'][kegg_module] += gene_read_count
                                            if kegg_module not in all_Funcs.all_KEGG_Modules:
                                                all_Funcs.all_KEGG_Modules.append(kegg_module)

                                    KEGG_Reaction = line_data[14]
                                    if KEGG_Reaction != '-':
                                        kegg_reaction_list = KEGG_Reaction.split(',')
                                        for kegg_reaction in kegg_reaction_list:
                                            functions[dir]['KEGG_Reaction'][kegg_reaction] += gene_read_count
                                            if kegg_reaction not in all_Funcs.all_KEGG_Reactions:
                                                all_Funcs.all_KEGG_Reactions.append(kegg_reaction)

                                    KEGG_rclass = line_data[15]
                                    if KEGG_rclass != '-':
                                        kegg_rclass_list = KEGG_rclass.split(',')
                                        for kegg_rclass in kegg_rclass_list:
                                            functions[dir]['KEGG_rclass'][kegg_rclass] += gene_read_count
                                            if kegg_rclass not in all_Funcs.all_KEGG_rclasses:
                                                all_Funcs.all_KEGG_rclasses.append(kegg_rclass)

                                    BRITE = line_data[16]
                                    if BRITE != '-':
                                        brite_list = BRITE.split(',')
                                        for brite in brite_list:
                                            functions[dir]['BRITE'][brite] += gene_read_count
                                            if brite not in all_Funcs.all_BRITEs:
                                                all_Funcs.all_BRITEs.append(brite)

                                    KEGG_TC = line_data[17]
                                    if KEGG_TC != '-':
                                        kegg_tc_list = KEGG_TC.split(',')
                                        for kegg_tc in kegg_tc_list:
                                            functions[dir]['KEGG_TC'][kegg_tc] += gene_read_count
                                            if kegg_tc not in all_Funcs.all_KEGG_TCs:
                                                all_Funcs.all_KEGG_TCs.append(kegg_tc)

                                    CAZy = line_data[18]
                                    if CAZy != '-':
                                        cazy_list = CAZy.split(',')
                                        for cazy in cazy_list:
                                            functions[dir]['CAZy'][cazy] += gene_read_count
                                            if cazy not in all_Funcs.all_CAZys:
                                                all_Funcs.all_CAZys.append(cazy)

                                    BiGG_Reaction = line_data[19]
                                    if BiGG_Reaction != '-':
                                        bigg_reaction_list = BiGG_Reaction.split(',')
                                        for bigg_reaction in bigg_reaction_list:
                                            functions[dir]['BiGG_Reaction'][bigg_reaction] += gene_read_count
                                            if bigg_reaction not in all_Funcs.all_BiGG_Reactions:
                                                all_Funcs.all_BiGG_Reactions.append(bigg_reaction)

                                    PFAMs = line_data[20]
                                    if PFAMs != '-':
                                        pfams_list = PFAMs.split(',')
                                        for pfam in pfams_list:
                                            functions[dir]['PFAM'][pfam] += gene_read_count
                                            if pfam not in all_Funcs.all_PFAMs:
                                                all_Funcs.all_PFAMs.append(pfam)


                    print("Sample " + dir + " Done")



            except FileNotFoundError as e:
                print(e)

    functions = dict(sorted(functions.items()))
    return all_Funcs, functions


def write_transposed_output(output_path, category, all_items, functions):
    output = open(output_path, 'w')
    output.write('\t' + '\t'.join(functions.keys()) + '\n')

    for item in all_items:
        row_values = [item] + [str(functions[sample][category][item]) for sample in functions.keys()]
        output.write('\t'.join(row_values) + '\n')

    output.close()


def main():
    parser = argparse.ArgumentParser(description='MetaPont: Function assigned reads')
    parser.add_argument(
        "-d", "--directory", required=True,
        help="Directory containing sample directories to analyse."
    )
    parser.add_argument(
        "-p", "--prefix", required=False, default='PN',
        help="Default - 'PN': Default directory name prefix to identify sample directories to analyse.."
    )
    parser.add_argument(
        "-t", "--ex_taxon", required=False, default='k__Metazoa,k__Viridiplantae',
        help="Default - 'k__Metazoa,k__Viridiplantae': Taxon groups to exclude (provide as list e.g. 'k__Metazoa,k__Viridiplantae')."
    )
    parser.add_argument(
        "-m", "--minlen", required=False, default=2500,
        help="Default - 2500: Minimum contig length allowed for function recording (bps)."
    )

    parser.add_argument(
        "-o", "--outdir", required=False,
        help="Default - same as input dir: Location to save results directory (CDS_Final_Outputs)."
    )


    options = parser.parse_args()
    print("Running MetaPont: Function assigned reads")

    input_path = os.path.abspath(options.directory)
    if options.outdir == None:
        options.outdir = input_path
    output_path = os.path.abspath(options.outdir)
    output_path = os.path.join(output_path, 'CDS_Final_Outputs')
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    all_Funcs, functions = processFunction(input_path, options.ex_taxon, options.minlen, options.prefix)





    # COG output
    write_transposed_output(output_path + '/Final_CDS_COG.tsv', 'COG', all_Funcs.all_COGs, functions)

    # GO output
    write_transposed_output(output_path + '/Final_CDS_GO.tsv', 'GO', all_Funcs.all_GOs, functions)

    # EC output
    write_transposed_output(output_path + '/Final_CDS_EC.tsv', 'EC', all_Funcs.all_ECs, functions)

    # KEGG_ko output
    write_transposed_output(output_path + '/Final_CDS_KEGG_ko.tsv', 'KEGG_ko', all_Funcs.all_KEGG_kos, functions)

    # KEGG_Pathway output
    write_transposed_output(output_path + '/Final_CDS_KEGG_Pathway.tsv', 'KEGG_Pathway', all_Funcs.all_KEGG_Pathways, functions)

    # KEGG_Module output
    write_transposed_output(output_path + '/Final_CDS_KEGG_Module.tsv', 'KEGG_Module', all_Funcs.all_KEGG_Modules, functions)

    # KEGG_Reaction output
    write_transposed_output(output_path + '/Final_CDS_KEGG_Reaction.tsv', 'KEGG_Reaction', all_Funcs.all_KEGG_Reactions, functions)

    # KEGG_rclass output
    write_transposed_output(output_path + '/Final_CDS_KEGG_rclass.tsv', 'KEGG_rclass', all_Funcs.all_KEGG_rclasses, functions)

    # BRITE output
    write_transposed_output(output_path + '/Final_CDS_BRITE.tsv', 'BRITE', all_Funcs.all_BRITEs, functions)

    # KEGG_TC output
    write_transposed_output(output_path + '/Final_CDS_KEGG_TC.tsv', 'KEGG_TC', all_Funcs.all_KEGG_TCs, functions)

    # CAZy output
    write_transposed_output(output_path + '/Final_CDS_CAZy.tsv', 'CAZy', all_Funcs.all_CAZys, functions)

    # BiGG_Reaction output
    write_transposed_output(output_path + '/Final_CDS_BiGG_Reaction.tsv', 'BiGG_Reaction', all_Funcs.all_BiGG_Reactions, functions)

    # PFAMs output
    write_transposed_output(output_path + '/Final_CDS_PFAMs.tsv', 'PFAM', all_Funcs.all_PFAMs, functions)

if __name__ == "__main__":
    main()
    print("Complete")