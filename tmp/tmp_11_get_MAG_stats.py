import os
import glob
from Bio import SeqIO

########################################################################################################################

# cpl_ctm_str = '50_5'
#
# # input
# mag_folder      = '/Users/songweizhi/Desktop/all/metawrap_%s_bins' % cpl_ctm_str
# mag_ext         = 'fa'
# mag_gtdb_txt    = '/Users/songweizhi/Desktop/all/metawrap_%s_bins_GTDB_r207.txt' % cpl_ctm_str
# mag_checkm_txt  = '/Users/songweizhi/Desktop/all/metawrap_%s_bins.stats' % cpl_ctm_str
# mag_depth_txt   = '/Users/songweizhi/Desktop/all/metawrap_%s_bins_depth.txt' % cpl_ctm_str
#
# # output
# mag_stats_txt   = '/Users/songweizhi/Desktop/all/metawrap_%s_bins_stats.txt' % cpl_ctm_str

########################################################################################################################


# input
mag_folder      = '/Users/songweizhi/Desktop/Nanopore/rd123/get_mag_stats/metawrap_50_5_bins_manually_curated'
mag_ext         = 'fa'
mag_gtdb_txt    = '/Users/songweizhi/Desktop/Nanopore/rd123/get_mag_stats/MetaWRAP_50_5_wd_2nd_round_refinement_metawrap_50_5_bins.bac120.summary.tsv'
mag_checkm_txt  = '/Users/songweizhi/Desktop/Nanopore/rd123/get_mag_stats/MetaWRAP_50_5_wd_2nd_round_refinement_metawrap_50_5_bins_manually_curated.stats'
mag_depth_txt   = '/Users/songweizhi/Desktop/Nanopore/rd123/get_mag_stats/mag_depth.txt'

# output
mag_stats_txt   = '/Users/songweizhi/Desktop/Nanopore/rd123/get_mag_stats/metawrap_MAG_stats.txt'

########################################################################################################################

mag_taxon_dict = {}
for each_mag in open(mag_gtdb_txt):
    each_mag_split = each_mag.strip().split('\t')
    mag_taxon_dict[each_mag_split[0]]= each_mag_split[1]

mag_depth_dict = {}
for each_mag in open(mag_depth_txt):
    each_mag_split = each_mag.strip().split('\t')
    mag_id = '.'.join(each_mag_split[0].split('.')[:-1])
    mag_depth_dict[mag_id] = each_mag_split[2]

mag_cpl_dict = {}
mag_ctm_dict = {}
mag_size_dict = {}
for each_mag in open(mag_checkm_txt):
    each_mag_split = each_mag.strip().split('\t')
    mag_cpl_dict[each_mag_split[0]] = each_mag_split[1]
    mag_ctm_dict[each_mag_split[0]] = each_mag_split[2]
    mag_size_dict[each_mag_split[0]] = each_mag_split[6]


mag_file_re = '%s/*.%s' % (mag_folder, mag_ext)
mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]


mag_stats_txt_handle = open(mag_stats_txt, 'w')
mag_stats_txt_handle.write('MAG\tCompleteness\tContamination\tSize(Mbp)\tContig\tDepth(x)\tTaxonomy\n')

for each_mag in mag_file_list:
    mag_id = '.'.join(each_mag.split('.')[:-1])
    pwd_mag = '%s/%s' % (mag_folder, each_mag)
    ctg_num = 0
    for each_ctg in SeqIO.parse(pwd_mag, 'fasta'):
        ctg_num += 1

    mag_size_mbp = int(mag_size_dict[mag_id])/(1024*1024)
    mag_size_mbp = float("{0:.2f}".format(mag_size_mbp))

    mag_stats_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each_mag, mag_cpl_dict[mag_id], mag_ctm_dict[mag_id], mag_size_mbp, ctg_num, mag_depth_dict[mag_id], mag_taxon_dict.get(mag_id, 'NA')))
mag_stats_txt_handle.close()

