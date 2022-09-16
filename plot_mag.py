import os
import glob
import math
import random
import numpy as np
import seaborn as sns
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def get_color_list(color_num):

    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']

    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a', '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']

    else:
        color_num_each = math.ceil(color_num/8) + 2

        color_list_1 = sns.color_palette('Blues',  n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy',   n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green',  n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds',   n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive',  n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()

        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6, color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    return color_list_to_return_sorted


#bin_folder = '/Users/songweizhi/Desktop/metawrap_20_10_bins'
bin_folder = '/Users/songweizhi/Desktop/metawrap_50_5_bins'
bin_ext    = 'fa'
depth_file = '/Users/songweizhi/Desktop/filtered_stats_for_Maxbin.txt'
output_plot = '/Users/songweizhi/Desktop/metawrap_50_5_bins.pdf'

ctg_depth_dict = dict()
for each_ctg in open(depth_file):
    each_ctg_split = each_ctg.strip().split('\t')
    if not each_ctg.startswith('contigName	contigLen	totalAvgDepth	sample1	sample1-var'):
        ctg_depth_dict[each_ctg_split[0]] = float(each_ctg_split[1])


bin_file_re = '%s/*.%s' % (bin_folder, bin_ext)
bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]


color_list = get_color_list(len(bin_file_list))
#color_list = ['#2b7bba', '#539ecd', '#89bedc', '#000080', '#303097', '#6060ae', '#9090c5', '#ffa500', '#fdb430', '#fac35f', '#f8d28f', '#008000', '#2f972f', '#5eae5e', '#800080', '#973097', '#c590c5', '#aa1016', '#d52221', '#f44f39', '#808000', '#97972c', '#aeae59', '#5c5c5c', '#828282', '#adadad']
#color_list = ['#0b559f', '#2b7bba', '#303097', '#6060ae', '#9090c5', '#ffa500', '#fac35f', '#f8d28f', '#008000', '#2f972f', '#8dc58d', '#800080', '#973097', '#ae60ae', '#c590c5', '#aa1016', '#f44f39', '#fc8161', '#808000', '#97972c', '#aeae59', '#c4c485', '#2b2b2b', '#5c5c5c', '#828282', '#adadad']
#color_list = ['#89bedc', '#8a8acc', '#8acc8a', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c']

from matplotlib.patches import Patch

legend_elements = []
num_list_depth = []
num_list_gc = []
ctg_color_list = []
ctg_len_list = []
index = 0
for each_mag in bin_file_list:
    pwd_each_mag = '%s/%s' % (bin_folder, each_mag)
    for each_ctg in SeqIO.parse(pwd_each_mag, 'fasta'):
        ctg_depth = ctg_depth_dict[each_ctg.id]
        ctg_seq = str(each_ctg.seq).upper()
        ctg_gc = ((ctg_seq.count('G')) + (ctg_seq.count('C')))*100/len(ctg_seq)
        num_list_depth.append(ctg_depth)
        num_list_gc.append(ctg_gc)
        ctg_color_list.append(color_list[index])
        ctg_len_list.append(len(ctg_seq)/2000)

    legend_elements.append(Patch(facecolor=color_list[index], edgecolor=color_list[index], label=each_mag))

    index += 1


num_array_depth = np.array(num_list_depth)
num_array_gc = np.array(num_list_gc)
num_array_len = np.array(ctg_len_list)


plt.scatter(num_array_depth, num_array_gc, c=ctg_color_list, s=num_array_len, alpha=0.3, linewidths=0)
#plt.scatter(num_array_depth, num_array_gc, c=ctg_color_list, s=5, alpha=0.6, linewidths=0)
plt.xlabel("Depth (X)", size=12)
plt.ylabel("GC (%)", size=12)
plt.legend(handles=legend_elements, loc='center left', prop={'size': 6}, bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig(output_plot)
plt.close()

