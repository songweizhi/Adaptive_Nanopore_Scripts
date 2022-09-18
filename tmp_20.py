import os
from Bio import SeqIO

########################################################################################################################

# file in
strainberry_scaffolds_info_tsv = '/Users/songweizhi/Desktop/assembly.scaffolds.info.tsv'
strainberry_scaffolds_fasta    = '/Users/songweizhi/Desktop/assembly.scaffolds.fa'
depth_txt                      = '/Users/songweizhi/Desktop/assembly.scaffolds.fa.depth.txt'

########################################################################################################################

ctg_length_dict = {}
for each_seq in SeqIO.parse(strainberry_scaffolds_fasta, 'fasta'):
    seq_id = each_seq.id
    seq_len = len(str(each_seq.seq))
    ctg_length_dict[seq_id] = seq_len

ctg_depth_dict = {}
for each_ctg in open(depth_txt):
    each_ctg_split = each_ctg.strip().split('\t')
    seq_id = each_ctg_split[0]
    seq_length = int(each_ctg_split[1])
    seq_depth = float(each_ctg_split[2])
    ctg_depth_dict[seq_id] = seq_depth

ref_to_haplotig_dict = {}
for each in open(strainberry_scaffolds_info_tsv):
    each_split = each.strip().split('\t')
    seq_id = each_split[0]
    seq_ref = each_split[4].split('=')[1]
    seq_length = ctg_length_dict[seq_id]
    seq_depth = ctg_depth_dict.get(seq_id, 0)

    hp_assign = 'unphased'
    if 'phased=false' not in each.strip():
        hp_assign = each_split[6]
    haplotig_info_to_keep = '%s__%sX__%sbp__%s' % (seq_id, seq_depth, seq_length, hp_assign)
    if seq_ref not in ref_to_haplotig_dict:
        ref_to_haplotig_dict[seq_ref] = [haplotig_info_to_keep]
    else:
        ref_to_haplotig_dict[seq_ref].append(haplotig_info_to_keep)

for each_ref in ref_to_haplotig_dict:
    print('%s\t%s' % (each_ref, '\t'.join(ref_to_haplotig_dict[each_ref])))
