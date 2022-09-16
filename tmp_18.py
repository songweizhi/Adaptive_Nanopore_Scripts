
def get_bin_completeness_contamination_dict(checkm_txt):
    bin_completeness_dict = {}
    bin_contamination_dict = {}
    for each_bin in open(checkm_txt):
        if not each_bin.startswith('bin\tcompleteness'):
            each_bin_split = each_bin.strip().split('\t')
            bin_completeness_dict[each_bin_split[0]] = each_bin_split[1]
            bin_contamination_dict[each_bin_split[0]] = each_bin_split[2]
    return bin_completeness_dict, bin_contamination_dict


def get_taxon_to_bin_dict(gtdb_txt):
    taxon_to_bin_dict = {}
    for each in open(gtdb_txt):
        if not each.startswith('user_genome\tclassification'):
            each_split = each.strip().split('\t')
            bin_id = each_split[0]
            taxon_str = each_split[1]
            if taxon_str not in taxon_to_bin_dict:
                taxon_to_bin_dict[taxon_str] = bin_id
            else:
                print('%s already found' % taxon_str)
    return taxon_to_bin_dict


########################################################################################################################

wd                   = '/Users/songweizhi/Desktop/CompareBinning'

bin_set_1_checkm_txt = '%s/MetaBAT2_MyCC_SemiBin_global_metawrap_50_5_bins.stats'                           % wd
bin_set_1_gtdb_txt   = '%s/MetaBAT2_MyCC_SemiBin_global_metawrap_50_5_bins.bac120.summary.tsv'              % wd

bin_set_2_checkm_txt = '%s/MetaBAT_superspecific_MyCC_SemiBin_global_metawrap_50_5_bins.stats'              % wd
bin_set_2_gtdb_txt   = '%s/MetaBAT_superspecific_MyCC_SemiBin_global_metawrap_50_5_bins.bac120.summary.tsv' % wd

bin_set_3_checkm_txt = '%s/MetaBAT_verysensitive_MyCC_SemiBin_global_metawrap_50_5_bins.stats'              % wd
bin_set_3_gtdb_txt   = '%s/MetaBAT_verysensitive_MyCC_SemiBin_global_metawrap_50_5_bins.bac120.summary.tsv' % wd

bin_set_4_checkm_txt = '%s/MetaWRAP_50_5_wd_2nd_round_refinement_metawrap_50_5_bins.stats'                  % wd
bin_set_4_gtdb_txt   = '%s/MetaWRAP_50_5_wd_2nd_round_refinement_metawrap_50_5_bins.bac120.summary.tsv'     % wd

########################################################################################################################

bin_completeness_dict_1, bin_contamination_dict_1 = get_bin_completeness_contamination_dict(bin_set_1_checkm_txt)
bin_completeness_dict_2, bin_contamination_dict_2 = get_bin_completeness_contamination_dict(bin_set_2_checkm_txt)
bin_completeness_dict_3, bin_contamination_dict_3 = get_bin_completeness_contamination_dict(bin_set_3_checkm_txt)
bin_completeness_dict_4, bin_contamination_dict_4 = get_bin_completeness_contamination_dict(bin_set_4_checkm_txt)

taxon_to_bin_dict_1 = get_taxon_to_bin_dict(bin_set_1_gtdb_txt)
taxon_to_bin_dict_2 = get_taxon_to_bin_dict(bin_set_2_gtdb_txt)
taxon_to_bin_dict_3 = get_taxon_to_bin_dict(bin_set_3_gtdb_txt)
taxon_to_bin_dict_4 = get_taxon_to_bin_dict(bin_set_4_gtdb_txt)

all_identified_taxon_set = set()
for each_taxon in taxon_to_bin_dict_1:
    all_identified_taxon_set.add(each_taxon)
for each_taxon in taxon_to_bin_dict_2:
    all_identified_taxon_set.add(each_taxon)
for each_taxon in taxon_to_bin_dict_3:
    all_identified_taxon_set.add(each_taxon)
for each_taxon in taxon_to_bin_dict_4:
    all_identified_taxon_set.add(each_taxon)

all_identified_taxon_list_sorted = sorted([i for i in all_identified_taxon_set])

print('Taxon\tbin_1\tbin_1_completeness\tbin_1_contamination\tbin_2\tbin_2_completeness\tbin_2_contamination\tbin_3\tbin_3_completeness\tbin_3_contamination\tbin_4\tbin_4_completeness\tbin_4_contamination')
for each_t in all_identified_taxon_list_sorted:
    current_taxon_bin_1 = taxon_to_bin_dict_1.get(each_t, 'NA')
    current_taxon_bin_2 = taxon_to_bin_dict_2.get(each_t, 'NA')
    current_taxon_bin_3 = taxon_to_bin_dict_3.get(each_t, 'NA')
    current_taxon_bin_4 = taxon_to_bin_dict_4.get(each_t, 'NA')
    bin_1_completeness = bin_completeness_dict_1.get(current_taxon_bin_1, 'NA')
    bin_2_completeness = bin_completeness_dict_2.get(current_taxon_bin_2, 'NA')
    bin_3_completeness = bin_completeness_dict_3.get(current_taxon_bin_3, 'NA')
    bin_4_completeness = bin_completeness_dict_4.get(current_taxon_bin_4, 'NA')
    bin_1_contamination = bin_contamination_dict_1.get(current_taxon_bin_1, 'NA')
    bin_2_contamination = bin_contamination_dict_2.get(current_taxon_bin_2, 'NA')
    bin_3_contamination = bin_contamination_dict_3.get(current_taxon_bin_3, 'NA')
    bin_4_contamination = bin_contamination_dict_4.get(current_taxon_bin_4, 'NA')

    print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_t,
                              current_taxon_bin_1, bin_1_completeness, bin_1_contamination,
                              current_taxon_bin_2, bin_2_completeness, bin_2_contamination,
                              current_taxon_bin_3, bin_3_completeness, bin_3_contamination,
                              current_taxon_bin_4, bin_4_completeness, bin_4_contamination))

