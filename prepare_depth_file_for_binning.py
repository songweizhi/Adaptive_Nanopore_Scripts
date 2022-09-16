from Bio import SeqIO


seq_file                = '/Users/songweizhi/Desktop/untitled_folder/assembly_min2500_nonEukaryota.fasta'
depth_file_MetaBAT_MyCC = '/Users/songweizhi/Desktop/untitled_folder/filtered_stats_MetaBAT_MyCC.txt'
depth_file_Maxbin       = '/Users/songweizhi/Desktop/untitled_folder/filtered_stats_Maxbin.txt'

depth_file_MetaBAT_MyCC_nonEukaryota = '/Users/songweizhi/Desktop/untitled_folder/filtered_stats_MetaBAT_MyCC_nonEukaryota.txt'
depth_file_Maxbin_nonEukaryota       = '/Users/songweizhi/Desktop/untitled_folder/filtered_stats_Maxbin_nonEukaryota.txt'


ctg_to_depth_dict_Maxbin = {}
for each_depth in open(depth_file_Maxbin):
    ctg_id = each_depth.strip().split('\t')[0]
    ctg_to_depth_dict_Maxbin[ctg_id] = each_depth.strip()

ctg_to_depth_dict_MetaBAT_MyCC = {}
for each_depth in open(depth_file_MetaBAT_MyCC):
    ctg_id = each_depth.strip().split('\t')[0]
    ctg_to_depth_dict_MetaBAT_MyCC[ctg_id] = each_depth.strip()

depth_file_Maxbin_nonEukaryota_handle = open(depth_file_Maxbin_nonEukaryota, 'w')
depth_file_MetaBAT_MyCC_nonEukaryota_handle = open(depth_file_MetaBAT_MyCC_nonEukaryota, 'w')
for each_seq in SeqIO.parse(seq_file, 'fasta'):
    depth_file_Maxbin_nonEukaryota_handle.write('%s\n' % ctg_to_depth_dict_Maxbin[each_seq.id])
    depth_file_MetaBAT_MyCC_nonEukaryota_handle.write('%s\n' % ctg_to_depth_dict_MetaBAT_MyCC[each_seq.id])
depth_file_Maxbin_nonEukaryota_handle.close()
depth_file_MetaBAT_MyCC_nonEukaryota_handle.close()


