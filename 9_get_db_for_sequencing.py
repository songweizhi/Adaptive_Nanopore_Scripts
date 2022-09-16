from Bio import SeqIO


assembly_fa                     = '/Users/songweizhi/Desktop/assembly_max2499.fasta'
genbank_id_to_taxon_id_txt      = '/Users/songweizhi/Desktop/Adaptive_Nanopore/55/uniq_genbank_id_taxon.txt'
taxon_id_to_lineage_txt         = '/Users/songweizhi/Desktop/Adaptive_Nanopore/55/tax_lineage_by_name.txt'
best_hit_txt                    = '/Users/songweizhi/Desktop/Adaptive_Nanopore/combined_v5_best_hit.tab'


# ref_taxon_txt                   = '/Users/songweizhi/Desktop/Adaptive_Nanopore/nt_seq_taxon.txt'
# assembly_with_matched_ref_txt   = '/Users/songweizhi/Desktop/Adaptive_Nanopore/assembly_with_matched_ref.txt'
# classification_txt              = '/Users/songweizhi/Desktop/Adaptive_Nanopore/classifications.txt'



# get seq_len_dict
seq_len_dict = {}
for each_seq in SeqIO.parse(assembly_fa, 'fasta'):
    seq_len_dict[each_seq.id] = len(each_seq.seq)


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

n = 0
for each_hit in open(best_hit_txt):
    each_hit_split = each_hit.strip().split('\t')
    aln_len = int(each_hit_split[3])
    ctg_id = each_hit_split[0]
    ref_id = each_hit_split[1]

    if ctg_id in seq_len_dict:
        if aln_len >= 100:

            ref_no_version = ref_id
            if '.' in ref_id:
                ref_no_version = ref_id.split('.')[0]

            ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
            ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')
            if 'Eukaryota' in ref_lineage:
                pass
                n += 1
                print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))
            elif ('Bacteria' in ref_lineage) or ('Archaea' in ref_lineage):
                pass
                #print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))
            else:
                pass
                #print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))

            # if ('Bacteria' in ref_lineage) or ('Archaea' in ref_lineage):
            #     #print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))
            #     print(ctg_id)
            #     n += 1
            #print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id,ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage ))

print(n)


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
# matched_ref_list_uniq = {i for i in matched_ref_list}
#
#
# # read in genbank id to taxon id
# genbank_id_to_taxon_id_dict = {}
# for each_line in open(genbank_id_to_taxon_id_txt):
#     each_line_split = each_line.strip().split('\t')
#     genbank_id_to_taxon_id_dict[each_line_split[0]] = each_line_split[1]
#
#
# # read in taxon id to lineage into dict
# taxon_id_to_lineage_dict = {}
# for each_line in open(taxon_id_to_lineage_txt):
#     each_line_split = each_line.strip().split('\t')
#     taxon_id_to_lineage_dict[each_line_split[0]] = each_line_split[1]
#
