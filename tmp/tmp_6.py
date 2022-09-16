from Bio import SeqIO


microbial_contigs_fa            = '/Users/songweizhi/Desktop/assembly_min2500_noneukaryota.fasta'
depth_in                        = '/Users/songweizhi/Desktop/filtered_stats_for_MetaBAT.txt'
depth_out_microbial_contigs     = '/Users/songweizhi/Desktop/filtered_stats_for_MetaBAT_assembly_min2500_noneukaryota.txt'


microbial_ctg_list = []
for each_seq in SeqIO.parse(microbial_contigs_fa, 'fasta'):
    microbial_ctg_list.append(each_seq.id)

# read in depth info
depth_dict = dict()
for each in open(depth_in):
    each_split = each.strip().split('\t')
    depth_dict[each_split[0]] = each.strip()

depth_out_microbial_contigs_handle = open(depth_out_microbial_contigs, 'w')
for each_microbial_ctg in microbial_ctg_list:
    depth_out_microbial_contigs_handle.write('%s\n' % depth_dict[each_microbial_ctg])
depth_out_microbial_contigs_handle.close()


