import os
import re
import sys
import glob
import shutil
import warnings
import argparse
import itertools
import subprocess
import numpy as np
from ete3 import Tree
from time import sleep
from datetime import datetime
from string import ascii_uppercase
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqFeature as SF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram


def get_group_index_list():
    def iter_all_strings():
        size = 1
        while True:
            for s in itertools.product(ascii_uppercase, repeat=size):
                yield "".join(s)
            size += 1

    group_index_list = []
    for s in iter_all_strings():
        group_index_list.append(s)
        if s == 'ZZ':
            break
    return group_index_list


def get_node_distance(tree, node_1, node_2):
    distance = tree.get_distance(node_1, node_2)
    return distance


def get_group(n):
    leaf_name = leaf_node_list[n]
    grouping_id = bin_to_grouping_dict[leaf_name]
    return '(%s) %s' % (grouping_id, leaf_name)


def plot_clustering_dendrogram(cluster, leaf_font_size, leaf_label_func, color_threshold, dendrogram_width, dendrogram_height, pwd_png_file):

    plt.figure(figsize=(dendrogram_width, dendrogram_height))
    plt.xlabel('Distance')
    dendrogram(cluster, orientation='left', leaf_rotation=0, leaf_font_size=leaf_font_size, leaf_label_func=leaf_label_func, color_threshold=color_threshold)
    plt.axvline(x=max_d, c='k', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(pwd_png_file, dpi=300)
    plt.close()


all_distances_lol = [[], []]

# turn list of distance list into arrary
all_distances_lol_array = np.array(all_distances_lol)

# get linkage
cluster = linkage(all_distances_lol_array, method='single')

# get maximum distance for clustering
distance_list = []
for each in cluster:
    distance_list.append(each[2])

# get distance cutoff
percentile_for_distances_cutoff = 90
distance_of_percentile = np.percentile(distance_list, percentile_for_distances_cutoff)

max_d = None
if max_d is None:
    max_d = distance_of_percentile
    sleep(0.5)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Determined distance cutoff is: %s, you can change it with option "-dc"' % float("{0:.2f}".format(max_d)))
else:
    sleep(0.5)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Distance cutoff specified to %s' % max_d)

# get flat clusters
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Grouping input genomes based on above distance cutoff')
flat_clusters = fcluster(cluster, max_d, criterion='distance')

# get group number
group_index_list = []
for each_index in flat_clusters:
    if int(each_index) not in group_index_list:
        group_index_list.append(each_index)
group_number = len(group_index_list)




# define output file name with grouping included
png_file_group =            '%s_grouping_g%s.png'       % (output_prefix, group_number)
grouping_file =             '%s_grouping_g%s.txt'       % (output_prefix, group_number)
grouping_file_temp =        '%s_grouping_g%s_tmp.file'  % (output_prefix, group_number)
pwd_png_file_group =        '%s/%s'                     % (MetaCHIP_wd, png_file_group)
pwd_grouping_file =         '%s/%s'                     % (MetaCHIP_wd, grouping_file)
pwd_grouping_file_temp =    '%s/%s'                     % (MetaCHIP_wd, grouping_file_temp)

# get grouping file
group_index_list = get_group_index_list()
grouping_file_temp_handle = open(pwd_grouping_file_temp, 'w')
bin_to_grouping_dict = {}
n = 0
for each_leaf in leaf_node_list:
    leaf_cluster_index = int(flat_clusters[n])
    leaf_grouping_id = group_index_list[leaf_cluster_index - 1]
    grouping_file_temp_handle.write('%s,%s\n' % (leaf_grouping_id,each_leaf))
    bin_to_grouping_dict[each_leaf] = leaf_grouping_id
    n += 1
grouping_file_temp_handle.close()

# sort grouping file
os.system('cat %s | sort > %s; rm %s' % (pwd_grouping_file_temp, pwd_grouping_file, pwd_grouping_file_temp))

# report
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Grouping results exported to: %s' % grouping_file)

# define dendrogram leaf_font_size and plot size
dendrogram_width = 10
input_genome_num = 100
# define leaf_font_size
leaf_font_size = 0
if input_genome_num <= 50:
    leaf_font_size = 8
elif 50 < input_genome_num <= 100:
    leaf_font_size = 7
elif 100 < input_genome_num <= 150:
    leaf_font_size = 6
elif input_genome_num > 150:
    leaf_font_size = 4

# define dendrogram_height
dendrogram_height = 0
if input_genome_num <= 20:
    dendrogram_height = 6
elif 20 < input_genome_num <= 50:
    dendrogram_height = 10
elif 50 < input_genome_num <= 100:
    dendrogram_height = 15
elif 100 < input_genome_num <= 150:
    dendrogram_height = 20
elif 150 < input_genome_num <= 300:
    dendrogram_height = 25
elif 300 < input_genome_num <= 500:
    dendrogram_height = 30
elif input_genome_num > 500:
    dendrogram_height = input_genome_num * 0.06


# calculate full dendrogram
sleep(0.5)
print(datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' Get plot for visualization')
plot_clustering_dendrogram(cluster, leaf_font_size, get_group, max_d, dendrogram_width, dendrogram_height, pwd_png_file_group)


