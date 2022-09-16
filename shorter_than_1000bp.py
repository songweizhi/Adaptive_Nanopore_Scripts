from Bio import SeqIO

assembly_fa             = '/Users/songweizhi/Desktop/Adaptive_Nanopore/assembly.fasta'
best_hit_txt            = '/Users/songweizhi/Desktop/Adaptive_Nanopore/combined_v5_best_hit.tab'
ref_id_txt              = '/Users/songweizhi/Desktop/Adaptive_Nanopore/ref_genbank_id.txt'
ref_to_taxon_id_txt     = '/Users/songweizhi/Desktop/Adaptive_Nanopore/ref_to_taxon_id.txt'
ref_taxon_id_uniq_txt   = '/Users/songweizhi/Desktop/Adaptive_Nanopore/ref_taxon_id_uniq.txt'

classification_aln100bp = '/Users/songweizhi/Desktop/Adaptive_Nanopore/classification_aln100bp.txt'


# get all matched refs
ref_set = set()
for each_hit in open(best_hit_txt):
    each_hit_split = each_hit.strip().split('\t')
    aln_len = int(each_hit_split[3])
    ctg_id = each_hit_split[0]
    ref_id = each_hit_split[1]
    ref_set.add(ref_id)


# write out refs
ref_id_txt_handle = open(ref_id_txt, 'w')
for each_ref in sorted([i for i in ref_set]):
    ref_id_txt_handle.write('%s\n' % each_ref)
ref_id_txt_handle.close()


# Genbank id --> taxon id
'''
cd /Users/songweizhi/Desktop/Adaptive_Nanopore
cat ref_genbank_id.txt | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId > ref_to_taxon_id.txt

esummary -db nuccore -id BA000012.4 | xtract -pattern DocumentSummary -element Caption,TaxId
CP078073.1       	1095823
XM_017150341.1      29030
'''


# get seq_len_dict
seq_len_dict = {}
for each_seq in SeqIO.parse(assembly_fa, 'fasta'):
    seq_len_dict[each_seq.id] = len(each_seq.seq)


# write out uniq taxon id
ref_taxon_id_set = set()
for each_line in open(ref_to_taxon_id_txt):
    ref_taxon_id_set.add(each_line.strip().split('\t')[1])
ref_taxon_id_uniq_txt_handle = open(ref_taxon_id_uniq_txt, 'w')
for each_ref in sorted([i for i in ref_taxon_id_set]):
    ref_taxon_id_uniq_txt_handle.write('%s\n' % each_ref)
ref_taxon_id_uniq_txt_handle.close()


# taxon id --> taxon lineage
# https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

# taxon lineage (id) --> taxon lineage (scientific name)
# python3 get_tax_lineage_by_name.py


########################################################################################################################

genbank_id_to_taxon_id_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/ref_to_taxon_id.txt'
taxon_id_to_lineage_txt    = '/Users/songweizhi/Desktop/Adaptive_Nanopore/ref_taxon_id_uniq_lineage_by_name.txt'

# read in genbank id to taxon id
genbank_id_to_taxon_id_dict = {}
for each_line in open(genbank_id_to_taxon_id_txt):
    each_line_split = each_line.strip().split('\t')
    genbank_id_to_taxon_id_dict[each_line_split[0]] = each_line_split[1]

# read in taxon id to lineage into dict
taxon_id_to_lineage_dict = {}
for each_line in open(taxon_id_to_lineage_txt):
    each_line_split = each_line.strip().split('\t')
    taxon_id_to_lineage_dict[each_line_split[0]] = each_line_split[1]

########################################################################################################################

classification_aln100bp_handle = open(classification_aln100bp, 'w')
classification_aln100bp_handle.write('contig\treference\tidentity\talignment_length\tcontig_length\taligned_percentage\ttaxonomy\n')
ctg_len_list_Bacteria_Archaea = []
ctg_len_list_Eukaryota = []

seq_num_Bacteria_Archaea = 0
seq_num_Eukaryota = 0
m = 0
for each_hit in open(best_hit_txt):
    each_hit_split = each_hit.strip().split('\t')
    aln_len = int(each_hit_split[3])
    ctg_id = each_hit_split[0]
    ref_id = each_hit_split[1]
    ctg_len = seq_len_dict[ctg_id]

    aln_pct = aln_len*100/ctg_len
    aln_pct = float("{0:.2f}".format(aln_pct))

    # get taxon lineage
    ref_no_version = ref_id
    if '.' in ref_id:
        ref_no_version = ref_id.split('.')[0]
    ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
    ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')

    if aln_len >= 100:
        if ctg_len < 100000:

            if 'Bacteria' in ref_lineage:
                ctg_len_list_Bacteria_Archaea.append(ctg_len)
            elif 'Archaea' in ref_lineage:
                ctg_len_list_Bacteria_Archaea.append(ctg_len)
            elif 'Eukaryota' in ref_lineage:
                ctg_len_list_Eukaryota.append(ctg_len)




        m += 1
        if 'Bacteria' in ref_lineage:
            seq_num_Bacteria_Archaea += 1
        elif 'Archaea' in ref_lineage:
            seq_num_Bacteria_Archaea += 1
        elif 'Eukaryota' in ref_lineage:
            seq_num_Eukaryota += 1
        else:
            pass
            #print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))

        classification_aln100bp_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], aln_pct, ref_lineage))

        # print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))

print('seq_num_Bacteria_Archaea: %s' % seq_num_Bacteria_Archaea)
print('seq_num_Eukaryota: %s' % seq_num_Eukaryota)
print(m)

classification_aln100bp_handle.close()


########################################################################################################################

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


plot_file = '/Users/songweizhi/Desktop/GTDB_vs_SILVA_confidence.png'

# similarity_list_silva = ctg_len_list_Bacteria_Archaea
# similarity_list_gtdb = ctg_len_list_Eukaryota
#
# min_value = round(min([min(similarity_list_silva), min(similarity_list_gtdb)]) - 1)
# tick_value_list = list(range(min_value, 101, 1))
#
# bin_num_1 = list(range(round(min(similarity_list_silva)), 101, 1))
# bin_num_2 = list(range(round(min(similarity_list_gtdb)), 101, 1))
#
# print(bin_num_1)
# print(bin_num_2)

plt.figure(figsize=(8,6))

#plt.hist(similarity_list_silva, bins=len(bin_num_1)-1, alpha=0.5, label="SILVA (%s)" % len(similarity_list_silva))
#plt.hist(similarity_list_gtdb, bins=len(bin_num_2)-1, alpha=0.5, label="GTDB (%s)" % len(similarity_list_gtdb))

plt.hist(ctg_len_list_Bacteria_Archaea, bins=100, alpha=0.5, label="Bacteria and Archaea")
plt.hist(ctg_len_list_Eukaryota, bins=100, alpha=0.5, label="Eukaryota")


#plt.xticks(tick_value_list, size=15)
plt.yticks(size=15)
#plt.xlabel("Alignment Length", size=16)
plt.xlabel("Contig Length", size=16)
plt.ylabel("Number of contigs", size=16)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(plot_file)
plt.close()

# print('sequence similarity >=90 (GTDB): %s' % (similarity_list_gtdb_lt_90/len(similarity_list_gtdb)))
# print('sequence similarity >=90 (SILVA): %s' % (similarity_list_silva_lt_90/len(similarity_list_silva)))

########################################################################################################################


    # if 100 <= aln_len < 1000:
    #     print(each_hit_split)
    #     m += 1
    #
    #     ref_no_version = ref_id
    #     if '.' in ref_id:
    #         ref_no_version = ref_id.split('.')[0]

        #ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
        #ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')

        #assembly_with_matched_ref_txt_handle.write('%s\n' % ctg_id)

        #if 'Eukaryota' in ref_lineage:
            #print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id,ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage ))

        #classification_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))

    #     matched_ref_list.append(ref_id)
    #     if ref_id not in ref_to_ctg_dict:
    #         ref_to_ctg_dict[ref_id] = [ctg_id]
    #     else:
    #         ref_to_ctg_dict[ref_id].append(ctg_id)
    #     n += 1
    # m += 1



#
# ref_taxon_txt       = '/Users/songweizhi/Desktop/Adaptive_Nanopore/nt_seq_taxon.txt'
# uniq_genbank_id_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/uniq_genbank_id.txt'
#
# assembly_with_matched_ref_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/assembly_with_matched_ref.txt'
#
# classification_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/classifications.txt'
#
#

# # read in genbank id to taxon id
# genbank_id_to_taxon_id_txt = '/Users/songweizhi/Desktop/55/uniq_genbank_id_taxon.txt'
# genbank_id_to_taxon_id_dict = {}
# for each_line in open(genbank_id_to_taxon_id_txt):
#     each_line_split = each_line.strip().split('\t')
#     genbank_id_to_taxon_id_dict[each_line_split[0]] = each_line_split[1]
#
# # read in taxon id to lineage into dict
# taxon_id_to_lineage_txt = '/Users/songweizhi/Desktop/55/tax_lineage_by_name.txt'
# taxon_id_to_lineage_dict = {}
# for each_line in open(taxon_id_to_lineage_txt):
#     each_line_split = each_line.strip().split('\t')
#     taxon_id_to_lineage_dict[each_line_split[0]] = each_line_split[1]
#
# classification_txt_handle = open(classification_txt, 'w')
# assembly_with_matched_ref_txt_handle = open(assembly_with_matched_ref_txt, 'w')
# m = 0
# n = 0
# matched_ref_list = []
# ref_to_ctg_dict = dict()
# for each_hit in open(best_hit_txt):
#     each_hit_split = each_hit.strip().split('\t')
#     aln_len = int(each_hit_split[3])
#     ctg_id = each_hit_split[0]
#     ref_id = each_hit_split[1]
#     if aln_len >= 1000:
#
#         ref_no_version = ref_id
#         if '.' in ref_id:
#             ref_no_version = ref_id.split('.')[0]
#
#         ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
#         ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')
#
#         assembly_with_matched_ref_txt_handle.write('%s\n' % ctg_id)
#
#         if 'Eukaryota' in ref_lineage:
#             print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id,ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage ))
#
#         classification_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))
#
#         matched_ref_list.append(ref_id)
#         if ref_id not in ref_to_ctg_dict:
#             ref_to_ctg_dict[ref_id] = [ctg_id]
#         else:
#             ref_to_ctg_dict[ref_id].append(ctg_id)
#         n += 1
#     m += 1
# assembly_with_matched_ref_txt_handle.close()
# classification_txt_handle.close()
#
#
# matched_ref_list_uniq = {i for i in matched_ref_list}
#
#
# # write out Genbank id list
# uniq_genbank_id_txt_handle = open(uniq_genbank_id_txt, 'w')
# for each_id in matched_ref_list_uniq:
#     uniq_genbank_id_txt_handle.write('%s\n' % each_id)
# uniq_genbank_id_txt_handle.close()
#

# print(m)
# print(n)
# print(len(matched_ref_list))
# print(len(matched_ref_list_uniq))
#
#
# # read in genbank id to taxon id
# genbank_id_to_taxon_id_txt = '/Users/songweizhi/Desktop/55/uniq_genbank_id_taxon.txt'
# genbank_id_to_taxon_id_dict = {}
# for each_line in open(genbank_id_to_taxon_id_txt):
#     each_line_split = each_line.strip().split('\t')
#     genbank_id_to_taxon_id_dict[each_line_split[0]] = each_line_split[1]
#
# # read in taxon id to lineage into dict
# taxon_id_to_lineage_txt = '/Users/songweizhi/Desktop/55/tax_lineage_by_name.txt'
# taxon_id_to_lineage_dict = {}
# for each_line in open(taxon_id_to_lineage_txt):
#     each_line_split = each_line.strip().split('\t')
#     taxon_id_to_lineage_dict[each_line_split[0]] = each_line_split[1]
#
# # # read into dict
# # ref_taxon_dict = dict()
# # for each_ref_taxon in open(ref_taxon_txt):
# #     each_ref_taxon_split = each_ref_taxon.strip().split('\t')
# #     ref_taxon_dict[each_ref_taxon_split[0]] = each_ref_taxon_split[1]
# #
# # for each_ref in sorted([i for i in matched_ref_list_uniq]):
# #     print('%s\t%s\t%s' % (each_ref, matched_ref_list.count(each_ref), ref_taxon_dict.get(each_ref, 'NA')))
# #
# total_len = 0
# total_len_Eukaryota = 0
# total_len_Bacteria_Archaea = 0
# total_len_Virus = 0
# total_len_unclassified = 0
# for each_ref in ref_to_ctg_dict:
#
#     matched_ctg_list = ref_to_ctg_dict[each_ref]
#
#     matched_ctg_total_len = 0
#     for each_matched_ctg in matched_ctg_list:
#         ctg_len = seq_len_dict[each_matched_ctg]
#         matched_ctg_total_len += ctg_len
#
#     total_len += matched_ctg_total_len
#
#     ref_no_version = each_ref
#     if '.' in each_ref:
#         ref_no_version = each_ref.split('.')[0]
#
#     ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
#     ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')
#
#
#     #if 'Eukaryota' in ref_lineage:
#     #    print('%s\t%s\t%s\t%s' % (each_ref, matched_ctg_total_len, ','.join(matched_ctg_list), ref_lineage))
#
#
#     if 'Bacteria' in ref_lineage:
#         total_len_Bacteria_Archaea += matched_ctg_total_len
#     elif 'Archaea' in ref_lineage:
#         total_len_Bacteria_Archaea += matched_ctg_total_len
#     elif 'Eukaryota' in ref_lineage:
#         total_len_Eukaryota += matched_ctg_total_len
#         #print('%s\t%s\t%s\t%s' % (ref_lineage, each_ref, matched_ctg_total_len, ','.join(matched_ctg_list)))
#
#     elif 'unclassified' in ref_lineage:
#         total_len_unclassified += matched_ctg_total_len
#     elif 'Virus' in ref_lineage:
#         total_len_Virus += matched_ctg_total_len
#     else:
#         print('%s\t%s\t%s\t%s' % (ref_lineage, each_ref, matched_ctg_total_len, ','.join(matched_ctg_list)))
#
#
# print('total_len :\t%sbp\t%s Mbp' % (total_len, total_len/(1024*1024)))
# print('total_len_Bacteria_Archaea :\t%sMbp\t%s'    % (total_len_Bacteria_Archaea/(1024*1024), total_len_Bacteria_Archaea/total_len))
# print('total_len_Eukaryota:\t%sMbp\t%s'            % (total_len_Eukaryota/(1024*1024), total_len_Eukaryota/total_len))
# print('total_len_unclassified:\t%sMbp\t%s'         % (total_len_unclassified/(1024*1024), total_len_unclassified/total_len))
# print('total_len_Virus:\t%sMbp\t%s'                % (total_len_Virus/(1024*1024), total_len_Virus/total_len))
#
