
names_dmp_file  = '/Users/songweizhi/DB/taxdump/names.dmp'

tax_report_file = '/Users/songweizhi/Desktop/Adaptive_Nanopore/ref_taxon_id_uniq_lineage.txt'
output_txt      = '/Users/songweizhi/Desktop/Adaptive_Nanopore/ref_taxon_id_uniq_lineage_by_name.txt'

tax_report_file = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly_blastn_best_hits_taxon_id_uniq_tax_report.txt'
output_txt      = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly_blastn_best_hits_taxon_id_uniq_tax_report_by_name.txt'

tax_report_file = '/Users/songweizhi/Desktop/assembly_blastn_BestHit_ref_taxon_id_only_uniq_tax_report.txt'
output_txt      = '/Users/songweizhi/Desktop/assembly_blastn_BestHit_ref_taxon_id_only_uniq_tax_report_by_name.txt'


tax_report_file = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd1/000/combined_v5_best_hit_ref_id_uniq_taxon_id_only_uniq_tax_report.txt'
output_txt      = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd1/000/combined_v5_best_hit_ref_id_uniq_taxon_id_only_uniq_tax_report_by_name.txt'

tax_report_file = '/Users/songweizhi/Desktop/against_nr/assembly_blastn_best_hits_ref_id_to_taxon_id_only_uniq_tax_report.txt'
output_txt      = '/Users/songweizhi/Desktop/against_nr/assembly_blastn_best_hits_ref_id_to_taxon_id_only_uniq_tax_report_by_name.txt'


# read in names.dmp
node_id_to_name_dict = {}
for each_tax in open(names_dmp_file):
    each_tax_split = each_tax.strip().split('\t|\t')
    if 'scientific name' in each_tax:
        node_id = each_tax_split[0]
        node_name = each_tax_split[1]
        node_id_to_name_dict[node_id] = node_name

# read in tax_report_file
scientific_name_to_lineage_dict = {}
for each_tax in open(tax_report_file):
    each_tax_split = each_tax.strip().split('\t|\t')
    if len(each_tax_split) == 4:
        query_name = each_tax_split[1]
        query_lineage = each_tax_split[3]
        query_lineage_split = query_lineage.split(' ')
        query_lineage_name_list = [node_id_to_name_dict.get(i, 'NA') for i in query_lineage_split]
        query_lineage_name_str = ','.join(query_lineage_name_list[::-1])
        scientific_name_to_lineage_dict[query_name] = query_lineage_name_str

# write out
output_txt_handle = open(output_txt, 'w')
for each_tax_id in scientific_name_to_lineage_dict:
    if each_tax_id != 'taxid':
        output_txt_handle.write('%s\t%s\n' % (each_tax_id, scientific_name_to_lineage_dict[each_tax_id]))
output_txt_handle.close()


