import os
import math
from Bio import SeqIO
import os
import sys
import argparse
import kmedoids
import itertools
import numpy as np
import pandas as pd
from ete3 import Tree
from Bio import Phylo
from io import StringIO
from Bio import AlignIO
from sklearn import manifold
import multiprocessing as mp
from collections import Counter
from pyvis.network import Network
from sklearn.metrics import silhouette_score

from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def gfa_parser(gfa_file, fa_from_gfa):
    fa_from_gfa_handle = open(fa_from_gfa, 'w')
    gfa_seg_len_dict = dict()
    gfa_link_lol = []
    for each_line in open(gfa_file):
        each_line_split = each_line.strip().split('\t')

        # get gfa_seg_len_dict
        if each_line.startswith('S'):
            segment_id = each_line_split[1]
            segment_seq = each_line_split[2]
            gfa_seg_len_dict[segment_id] = len(segment_seq)

            # write out sequence
            fa_from_gfa_handle.write('>%s\n' % segment_id)
            fa_from_gfa_handle.write(segment_seq + '\n')

        # get gfa_link_lol
        elif each_line.startswith('L'):
            from_segment_id = each_line_split[1]
            from_segment_orientation = each_line_split[2]
            to_segment_id = each_line_split[3]
            to_segment_orientation = each_line_split[4]
            current_link_list = []
            current_link_list.append(from_segment_id)
            if from_segment_orientation == '+':
                current_link_list.append('r')
            elif from_segment_orientation == '-':
                current_link_list.append('l')
            current_link_list.append(to_segment_id)
            if to_segment_orientation == '+':
                current_link_list.append('l')
            elif to_segment_orientation == '-':
                current_link_list.append('r')
            gfa_link_lol.append(current_link_list)
    fa_from_gfa_handle.close()

    return gfa_seg_len_dict, gfa_link_lol


def format_mpileup_str(mpileup_str):
    mpileup_str_new = ''
    to_ignore  = 0
    for each in mpileup_str:
        if to_ignore == 0:
            if each in ['A', 'T', 'C', 'G', '*']:
                mpileup_str_new += each
            elif each in ['+', '-']:
                mpileup_str_new += ''
            elif each.isnumeric() is True:
                to_ignore = int(each)
                mpileup_str_new += ''
        else:
            to_ignore -= 1
            mpileup_str_new += ''
    return mpileup_str_new


def group_continuous_number(num_list, allowed_gap):

    # sort list/set
    num_list_sorted = sorted([i for i in num_list])

    num_list_grouped = []
    n = 0
    g = []
    while n <= (len(num_list_sorted) - 1):
        if len(g) == 0:
            g.append(num_list_sorted[n])
        else:
            if num_list_sorted[n] - g[-1] <= (allowed_gap + 1):
                g.append(num_list_sorted[n])
            else:
                num_list_grouped.append(g)
                g = [num_list_sorted[n]]
        n += 1

    num_list_grouped.append(g)

    return num_list_grouped


def get_node_distance(tree, node_1, node_2):
    distance = tree.get_distance(node_1, node_2)
    return distance


def tree_to_distance_matrix(tree_file, op_matrix):

    try:
        tree_in = Tree(tree_file, format=0)
        # get leaf node list
        leaf_node_list = []
        for leaf_node in tree_in:
            leaf_node_list.append(leaf_node.name)
        leaf_node_list = sorted(leaf_node_list)

        # get list of distance list
        all_distances_lol = []
        distance_dict = dict()
        op_matrix_handle = open(op_matrix, 'w')
        op_matrix_handle.write('\t%s\n' % '\t'.join(leaf_node_list))
        for each_node_1 in leaf_node_list:
            current_node_distance_list = []
            for each_node_2 in leaf_node_list:
                two_nodes_sorted = sorted([each_node_1, each_node_2])
                key_node_1_2 = '___'.join(two_nodes_sorted)
                if key_node_1_2 not in distance_dict:
                    distance = 0
                    if each_node_1 != each_node_2:
                        distance = get_node_distance(tree_in, each_node_1, each_node_2)
                        distance = float("{0:.3f}".format(distance))
                    distance_dict[key_node_1_2] = distance
                else:
                    distance = distance_dict[key_node_1_2]

                current_node_distance_list.append(str(distance))
            all_distances_lol.append(current_node_distance_list)
            op_matrix_handle.write('%s\t%s\n' % (each_node_1, '\t'.join([str(i) for i in current_node_distance_list])))
        op_matrix_handle.close()

    except:
        print('Problematic tree: %s' % tree_file.split('/')[-1])
        op_matrix_handle = open(op_matrix, 'w')
        op_matrix_handle.close()


# developed by iBiology (https://github.com/iBiology), based on plottree
def plottree(tree, axes, box, confidence, name, size, width, height, x, y, output):

    if os.path.isfile(tree):
        tree = Phylo.read(tree, 'newick')
    elif tree.startswith('(') and tree.endswith(';'):
        tree = Phylo.read(StringIO(tree), 'newick')
    else:
        raise TypeError('Invalid tree, tree only accepts NEWICK format file or string.')

    settings = []
    size = size or mpl.rcParams.get('font.size', 8)
    mpl.rcParams.update({'font.size': size})
    settings.append(f'fontsize (-s): {size:.1f}')

    fig, ax = plt.subplots()

    w, h = fig.get_size_inches()
    width, height = width or w, height or h
    fig.set_size_inches(width, height, forward=True)
    settings.append(f'width (-w): {width:.1f}')
    settings.append(f'height (-l): {height:.1f}')

    for clade in tree.find_clades():
        if not name and not clade.is_terminal():
            clade.name = None

    Phylo.draw(tree, do_show=False, axes=ax, show_confidence=confidence)

    x, y = x or ax.get_xlim(), y or ax.get_ylim()
    ax.set_xlim(*x), ax.set_ylim(*y)
    settings.append(f'xlim(-x): {x[0]:.2f}, {x[1]:.2f}')
    settings.append(f'ylim(-y): {y[0]:.2f}, {y[1]:.2f}')

    if axes:
        settings.append(f'Display axes ticks for further tuning: -a')
        plt.subplots_adjust(left=0.05, bottom=0.05, top=0.99, right=0.99)
    else:
        plt.subplots_adjust(left=0.01, bottom=0.01, top=0.99, right=0.99)
        ax.set_xticks([]), ax.set_yticks([])
        ax.set_xticklabels([]), ax.set_yticklabels([])
    ax.set_xlabel(''), ax.set_ylabel('')

    if box:
        settings.append('Plot a box surrounding the tree: -b')
    else:
        if axes:
            ax.spines['left'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
        else:
            for spine in ax.spines.values():
                spine.set_visible(False)
    if confidence:
        settings.append('Display confidence on the tree: -c')
    if name:
        settings.append('Display internal node name: -n')

    if output:
        settings.append(f'output figure to file (-o): {output}')
        fig.tight_layout()
        fig.savefig(output)
        plt.clf()


def executePCoA(dm_file, plot_file, drawBiplot, n_arrows, groupfile, plot_title):

    # my own way to get distance_matrix
    df = pd.read_table(dm_file, header=0, index_col=0)
    distance_matrix = df.to_numpy()
    row_name_list = list(df.index)
    n_samples = len(row_name_list)

    # execute PCoA
    mds = manifold.MDS(n_components=2, max_iter=3000, dissimilarity="precomputed", n_jobs=1)
    positions = mds.fit(distance_matrix).embedding_
    positions_with_sampleIndex = pd.DataFrame(positions, index=row_name_list)

    # General settings of the canvas
    fig = plt.figure(figsize=(12, 12))
    ax = fig.gca()
    ax.spines['right'].set_color(None)
    ax.spines['top'].set_color(None)
    #ax.spines['bottom'].set_color('gray')
    #ax.spines['left'].set_color('gray')
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    if drawBiplot:

        circle = plt.Circle((0, 0), radius=1.0, fc='none', linestyle='dashed', color='gray')
        ax.add_patch(circle)

        # compute correlations between feature vectors and data points.
        cor_pc1 = np.array([0.] * n_samples)
        cor_pc2 = np.array([0.] * n_samples)
        arrow_length = np.array([0.] * n_samples)
        # for i, current_feature in enumerate(data.index):
        #     if np.ptp(data.loc[current_feature].values) == 0:
        #         x = 0
        #         y = 0
        #     else:
        #         x = scipy.stats.pearsonr(data.loc[current_feature].values, positions[:, 0])[0]
        #         y = scipy.stats.pearsonr(data.loc[current_feature].values, positions[:, 1])[0]
        #     cor_pc1[i] = x
        #     cor_pc2[i] = y
        #     arrow_length[i] = np.sqrt(x ** 2 + y ** 2)
        arrows = pd.DataFrame(np.hstack((np.matrix(cor_pc1).T, np.matrix(cor_pc2).T, np.matrix(arrow_length).T)), index=row_name_list, columns=['x', 'y', 'len'])
        sorted_arrows = arrows.sort_values(by=['len'], ascending=False)
        # Top-{n_arrows} contributing features are drawed
        for name in sorted_arrows.index[:n_arrows]:
            ax.arrow(0.0, 0.0, arrows.loc[name, 'x'], arrows.loc[name, 'y'], ec='k', alpha=0.2)
            ax.annotate(name, xy=(arrows.loc[name, 'x'], arrows.loc[name, 'y']), xytext=(0, 0),
                        textcoords='offset points', color='k', fontsize=10)

    # draw plots using colors if samples are binned into groups
    if groupfile:
        group_names = []
        group2sample = {}
        df = pd.read_table(groupfile, header=None, index_col=0)
        for sample in df.index:
            # Use the value of second column (1) as a grouping category.
            group = df.loc[sample, 1]
            if group in group2sample:
                group2sample[group].append(sample)
            else:
                group2sample[group] = [sample]
                group_names.append(group)
        colors = itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])
        for i, current_group in enumerate(group_names):
            if len(group2sample[current_group]) == 0:
                continue
            ax.scatter(positions_with_sampleIndex.loc[group2sample[current_group], 0],
                       positions_with_sampleIndex.loc[group2sample[current_group], 1],
                       s=30, marker='o', color=next(colors), label='%s' % current_group)
        plt.legend(bbox_to_anchor=(0., 1.01, 1., 1.01), loc=3, ncol=6, mode="expand", borderaxespad=0.)
    else:
        #for i, sample_name in enumerate(data.columns):
        #    ax.annotate(sample_name, xy=(positions[i, 0], positions[i, 1]), xytext=(5, 5), textcoords='offset points', color='k', fontsize=16)
        ax.scatter(positions[:, 0], positions[:, 1], c='k', s=50)

    x_label = 'PCo1'
    y_label = 'PCo2'
    ax.annotate(x_label, xy=(0.0, -1.0), xytext=(0.0, -40.0), textcoords='offset points', ha='center', color='k', fontsize=18)
    ax.annotate(y_label, xy=(-1.0, 0.0), xytext=(-40.0, 0.0), textcoords='offset points', ha='center', color='k', fontsize=18, rotation=90)
    plt.title(plot_title)
    fig.savefig(plot_file)
    plt.close()
    fig.clf()


def get_elbow_plot(dm_txt, output_elbow_plot, bootstrap_num, plot_title):

    df = pd.read_table(dm_txt, header=0, index_col=0)
    cluster_list = list(range(1, len(df) + 1))
    cluster_list_str = [str(i) for i in cluster_list]

    loss_dict = dict()
    for each_cluster_num in cluster_list_str:
        loss_dict[each_cluster_num] = 0

    bootstrap_index = 1
    while bootstrap_index <= bootstrap_num:
        n = 1
        while n < len(df):
            current_loss = kmedoids.fasterpam(df, n).loss
            loss_dict[str(n)] += current_loss
            n += 1
        bootstrap_index += 1

    mean_loss_list = [(loss_dict[i]/bootstrap_num) for i in cluster_list_str]

    # plot
    plt.plot(cluster_list_str, mean_loss_list)
    plt.xlabel('Number of clusters')
    plt.ylabel('Loss')
    plt.title(plot_title)
    plt.tight_layout()
    plt.savefig(output_elbow_plot, bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()


def get_num_of_clusters(dm_file, silhouette_score_cutoff):

    df = pd.read_table(dm_file, header=0, index_col=0)
    df_in_arrary = df.to_numpy()

    max_score = -1
    best_cluster_num = None
    for cluster_num in range(2, df.shape[0]):

        # perform clustering
        clusters = kmedoids.fasterpam(df, cluster_num)
        label_list_1_based = [int((i + 1)) for i in clusters.labels]

        # calculate silhouette score
        score = silhouette_score(df_in_arrary, label_list_1_based)
        if score > max_score:
            max_score = score
            best_cluster_num = cluster_num

    if max_score < silhouette_score_cutoff:
        best_cluster_num = 1

    return best_cluster_num, max_score


def window_worker(arg_list):

    msa_file           = arg_list[0]
    msa_file_qc        = arg_list[1]
    dm_file            = arg_list[2]
    elbow_plot         = arg_list[3]
    tree_newick        = arg_list[4]
    tree_newick_with_c = arg_list[5]
    tree_plot_with_c   = arg_list[6]
    pcoa_plot_with_c   = arg_list[7]
    read_cluster_file  = arg_list[8]
    pos_set            = arg_list[9]
    read_id_set        = arg_list[10]
    read_to_pos_dict   = arg_list[11]
    max_N_pct_in_MSA   = arg_list[12]
    tree_format        = arg_list[13]
    FastTree_exe       = arg_list[14]
    silhouette_score_cutoff = arg_list[15]

    # sort
    pos_list_sorted     = sorted([i for i in pos_set])
    read_id_list_sorted = sorted([i for i in read_id_set])

    # write out
    current_subset_op_msa_handle = open(msa_file, 'w')
    current_subset_op_msa_qc_handle = open(msa_file_qc, 'w')
    for each_read in read_id_list_sorted:
        current_read_pos_base_dict = read_to_pos_dict.get(each_read, {})
        current_read_str = ''
        for each_pos in pos_list_sorted:
            current_pos_base = current_read_pos_base_dict.get(each_pos, 'N')
            current_read_str += current_pos_base

        current_read_str_reformatted = current_read_str.replace('*', '-')
        current_subset_op_msa_handle.write('>%s\n%s\n' % (each_read, current_read_str_reformatted))

        n_pct = 101
        if len(current_read_str_reformatted) > 0:
            n_pct = current_read_str_reformatted.count('N')*100/len(current_read_str_reformatted)

        if n_pct <= max_N_pct_in_MSA:
            current_subset_op_msa_qc_handle.write('>%s\n%s\n' % (each_read, current_read_str_reformatted))
    current_subset_op_msa_handle.close()
    current_subset_op_msa_qc_handle.close()

    # get tree
    fasttree_cmd = '%s -nt -quiet -out %s %s 2>/dev/null' % (FastTree_exe, tree_newick, msa_file_qc)
    os.system(fasttree_cmd)

    # get distance matrix
    tree_to_distance_matrix(tree_newick, dm_file)

    # get the best number of clusters
    cluster_num, cluster_score = get_num_of_clusters(dm_file, silhouette_score_cutoff)

    # perform clustering
    df = pd.read_table(dm_file, header=0, index_col=0)
    clusters = kmedoids.fasterpam(df, cluster_num)
    cluster_label_list = clusters.labels
    cluster_label_list_1_based = [int((i + 1)) for i in cluster_label_list]
    df_row_name_list = list(df.index.values)

    # write out clustering results
    read_to_cluster_dict = {}
    current_subset_read_cluster_handle = open(read_cluster_file, 'w')
    for c, r in zip(cluster_label_list_1_based, df_row_name_list):
        current_subset_read_cluster_handle.write('%s\tc%s\n' % (r, c))
        read_to_cluster_dict[r] = c
    current_subset_read_cluster_handle.close()

    # add cluster info to tree leaves and plot
    t = Tree(tree_newick, format=tree_format)
    for leaf in t:
        leaf_cluster = read_to_cluster_dict.get(leaf.name, '')
        leaf_name_new = '%s__%s' % (leaf_cluster, leaf.name)
        leaf.name = leaf_cluster
    t.write(format=tree_format, outfile=tree_newick_with_c)

    # plot tree
    plottree(tree_newick_with_c, False, False, False, False, None, None, None, None, None, tree_plot_with_c)

    # get pcoa plot
    plot_title = 'Clusters: %s, silhouette score: %s' % (cluster_num, float("{0:.2f}".format(cluster_score)))
    executePCoA(dm_file, pcoa_plot_with_c, True, 0, read_cluster_file, plot_title)

    # get elbow plot
    get_elbow_plot(dm_file, elbow_plot, 100, plot_title)


def sort_csv_by_col(file_in, file_out, col_header):
    df_in = pd.read_csv(file_in, sep='\t')
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False, sep='\t')


def filter_linkages(file_in_sorted, min_linkages, file_out):

    gnm_max_link_num_dict = {}
    file_out_handle = open(file_out, 'w')
    MarkerGene_with_assignment = set()
    GenomicSeq_best_marker_dict = {}
    with open(file_in_sorted) as fis:
        for each_match in fis:
            if each_match.startswith('W1\tW2\tLinkage'):
                file_out_handle.write(each_match)
            else:
                match_split = each_match.strip().split('\t')
                MarkerGene = match_split[0]
                GenomicSeq = match_split[1]
                linkage_num = int(match_split[2])
                if linkage_num >= min_linkages:
                    if MarkerGene not in MarkerGene_with_assignment:
                        if GenomicSeq not in GenomicSeq_best_marker_dict:
                            GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
                            gnm_max_link_num_dict[GenomicSeq] = linkage_num
                            file_out_handle.write(each_match)
                            MarkerGene_with_assignment.add(MarkerGene)
                        else:
                            MarkerGene_with_assignment.add(MarkerGene)
    file_out_handle.close()


def read_in_read_cluster_info(cluster_txt):
    window_to_clustered_reads_dict = dict()
    for each_r in open(cluster_txt):
        each_r_split = each_r.strip().split('\t')
        current_r = each_r_split[0]
        current_c = each_r_split[1]
        if current_c not in window_to_clustered_reads_dict:
            window_to_clustered_reads_dict[current_c] = {current_r}
        else:
            window_to_clustered_reads_dict[current_c].add(current_r)
    return window_to_clustered_reads_dict


def linking_windows(to_link_list, read_cluster_dir, linkage_dir, min_link_num):

    wd1_ref              = to_link_list[0]
    wd1_index            = to_link_list[1]
    wd2_ref              = to_link_list[2]
    wd2_index            = to_link_list[3]
    cluster_txt_l        = '%s/%s_wd%s_read_cluster.txt'              % (read_cluster_dir, wd1_ref, wd1_index)
    cluster_txt_r        = '%s/%s_wd%s_read_cluster.txt'              % (read_cluster_dir, wd2_ref, wd2_index)
    linkage_txt          = '%s/%s__%s---%s__%s_linkages.txt'          % (linkage_dir, wd1_ref, wd1_index, wd2_ref, wd2_index)
    linkage_txt_sorted   = '%s/%s__%s---%s__%s_linkages_sorted.txt'   % (linkage_dir, wd1_ref, wd1_index, wd2_ref, wd2_index)
    linkage_txt_filtered = '%s/%s__%s---%s__%s_linkages_filtered.txt' % (linkage_dir, wd1_ref, wd1_index, wd2_ref, wd2_index)
    read_cluster_l       = read_in_read_cluster_info(cluster_txt_l)
    read_cluster_r       = read_in_read_cluster_info(cluster_txt_r)

    # get linkages
    linkage_handle = open(linkage_txt, 'w')
    linkage_handle.write('W1\tW2\tLinkage\n')
    for lc in read_cluster_l:
        for rc in read_cluster_r:
            lc_read_set = read_cluster_l[lc]
            rc_read_set = read_cluster_r[rc]
            lc_rc_shared = set(lc_read_set).intersection(rc_read_set)
            linkage_handle.write('%s__wd%s__%s\t%s__wd%s__%s\t%s\n' % (wd1_ref, wd1_index, lc, wd2_ref, wd2_index, rc, len(lc_rc_shared)))
    linkage_handle.close()

    # sort linkages
    sort_csv_by_col(linkage_txt, linkage_txt_sorted, 'Linkage')
    os.system('rm %s' % linkage_txt)

    # fileter linkage
    filter_linkages(linkage_txt_sorted, min_link_num, linkage_txt_filtered)

    # get all linkages after filtering
    linkage_sets_after_filtering = set()
    for each_linkage in open(linkage_txt_filtered):
        each_linkage_split = each_linkage.strip().split('\t')
        linkage_str = '%s___%s' % (each_linkage_split[0], each_linkage_split[1])
        linkage_sets_after_filtering.add(linkage_str)


def visualize_linkage(linkage_txt, min_link_num, linkage_html):

    edges_list = []
    for each in open(linkage_txt):
        each_split = each.strip().split('\t')
        node_1 = each_split[0]
        node_2 = each_split[1]
        link_num = int(each_split[2])
        passed_qc = each_split[3]
        edge_tuple = tuple([node_1, node_2, link_num, passed_qc])
        edges_list.append(edge_tuple)

    net1 = Network(height='100%', width='100%', bgcolor='white', font_color='black', directed=False)
    # net1 = Network(height='750px', width='100%', bgcolor='#222222', font_color='white')

    for each_link in edges_list:
        node_1_id = each_link[0]
        node_2_id = each_link[1]
        link_weight = each_link[2]
        passed_qc = each_link[3]
        node_1_group = node_1_id.split('__')[0]
        node_2_group = node_2_id.split('__')[0]
        net1.add_node(node_1_id, shape='box', group=node_1_group)
        net1.add_node(node_2_id, shape='box', group=node_2_group)

        # hide edge if link_weight == 0
        hide_edge = False if link_weight >= min_link_num else True

        # https://visjs.github.io/vis-network/docs/network/edges.html
        #net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge, arrows='arrows:%s;%s' % (node_1_id, node_2_id))
        #net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge)
        if passed_qc == 'y':
            net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge)
        if passed_qc == 'n':
            net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge, color='#E5E7E9')

    net1.show(linkage_html)


def ReadsPhaser():

    ########################################################################################################################

    # file in
    gfa_file                = '/Users/songweizhi/Desktop/demo_bin.9/graph_subset.gfa'
    mpileup_file            = '/Users/songweizhi/Desktop/demo_bin.9/combined.mpileup'
    FastTree_exe            = 'fasttree'
    window_size             = 7500
    min_var_num             = 1
    min_var_pct             = 10
    num_threads             = 10
    max_N_pct_in_MSA        = 50
    min_link_num_filter     = 1
    min_link_num_plot       = 1
    plot_intra_ref_linkages = True
    plot_inter_ref_linkages = True
    single_cluster_silhouette_score_cutoff = 0.55

    # file out
    output_dir              = '/Users/songweizhi/Desktop/TEST'

    ########################################################################################################################

    tree_format                          = 1

    fa_from_gfa                          = '%s/fa_from_gfa.fa'               % output_dir
    aln_dir                              = '%s/01_MSAs'                      % output_dir
    tree_dir                             = '%s/02_tree'                      % output_dir
    tree_plot_dir                        = '%s/02_tree_plot'                 % output_dir
    dm_dir                               = '%s/03_distance_matrix'           % output_dir
    tree_dir_with_c                      = '%s/04_tree_with_cluster'         % output_dir
    tree_plot_dir_with_c                 = '%s/05_tree_with_cluster_plot'    % output_dir
    elbow_plot_dir                       = '%s/06_elbow_plots'               % output_dir
    read_cluster_dir                     = '%s/07_read_cluster'              % output_dir
    linkage_dir                          = '%s/08_linkage'                   % output_dir
    combined_linkage_txt                 = '%s/linkages_combined.txt'        % output_dir
    combined_linkage_txt_filtered        = '%s/linkages_filtered.txt'        % output_dir
    combined_linkage_txt_with_pass_info  = '%s/linkages.txt'                 % output_dir
    combined_linkage_html                = '%s/linkages.html'                % output_dir

    if os.path.isdir(output_dir) is True:
        os.system('rm -r %s' % output_dir)

    os.system('mkdir %s' % output_dir)
    os.system('mkdir %s' % aln_dir)
    os.system('mkdir %s' % dm_dir)
    os.system('mkdir %s' % elbow_plot_dir)
    os.system('mkdir %s' % tree_dir)
    os.system('mkdir %s' % read_cluster_dir)
    os.system('mkdir %s' % tree_dir_with_c)
    os.system('mkdir %s' % tree_plot_dir_with_c)
    os.system('mkdir %s' % linkage_dir)

    ########################################################################################################################

    # parse gfa file
    ref_seq_len_dict, gfa_link_lol = gfa_parser(gfa_file, fa_from_gfa)

    # get ref_seq_window_num_dict
    ref_seq_window_num_dict = dict()
    for each_ref in ref_seq_len_dict:
        ref_len = ref_seq_len_dict[each_ref]
        frac, whole = math.modf(ref_len/window_size)
        if frac < 0.75:
            ref_seq_window_num_dict[each_ref] = int(whole)
        else:
            ref_seq_window_num_dict[each_ref] = int(whole) + 1

    # get gfa_link_lol_by_window
    gfa_link_lol_by_window = []
    for each_gfa_link in gfa_link_lol:
        current_gfa_link_by_window = []
        current_gfa_link_by_window.append(each_gfa_link[0])
        if each_gfa_link[1] == 'l':
            current_gfa_link_by_window.append(1)
        if each_gfa_link[1] == 'r':
            current_gfa_link_by_window.append(ref_seq_window_num_dict[each_gfa_link[0]])
        current_gfa_link_by_window.append(each_gfa_link[2])
        if each_gfa_link[3] == 'l':
            current_gfa_link_by_window.append(1)
        if each_gfa_link[3] == 'r':
            current_gfa_link_by_window.append(ref_seq_window_num_dict[each_gfa_link[2]])
        gfa_link_lol_by_window.append(current_gfa_link_by_window)

    # get ref_seq_window_index_dict and ref_seq_window_pos_dict
    ref_seq_window_index_dict = dict()
    ref_seq_window_pos_dict = dict()
    for each_ref in ref_seq_window_num_dict:
        window_num = ref_seq_window_num_dict[each_ref]
        window_index_list = list(range(1, (window_num + 1)))
        ref_seq_window_index_dict[each_ref] = window_index_list
        current_ref_window_pos_dict = dict()
        for each_window in window_index_list:
            window_range_l = window_size * (each_window - 1) + 1
            window_range_r = window_size * each_window
            if each_window == window_num:
                window_range_r = ref_seq_len_dict[each_ref]
            # print('%s\t%s\t%s\t%s' % (each_ref, each_window, window_range_l, window_range_r))
            current_ref_window_pos_dict[each_window] = [window_range_l, window_range_r]
        ref_seq_window_pos_dict[each_ref] = current_ref_window_pos_dict

    ########################################################################################################################

    refseq_in_mpileup = set()
    window_to_snv_pos_dod = dict()
    window_to_reads_dict = dict()
    read_to_pos_base_dod = dict()
    for each_pos in open(mpileup_file):
        each_pos_split = each_pos.strip().split('\t')
        ref_id = each_pos_split[0]
        ref_pos = int(each_pos_split[1])
        ref_pos_dv = round(ref_pos/window_size)
        refseq_in_mpileup.add(ref_id)

        # get current pos window index
        putative_window_list = [(ref_pos_dv - 2), (ref_pos_dv - 1), ref_pos_dv, (ref_pos_dv + 1), (ref_pos_dv + 2)]
        window_index = ''
        for each_p_window in putative_window_list:
            p_window_range = ref_seq_window_pos_dict[ref_id].get(each_p_window, [])
            if len(p_window_range) == 2:
                if p_window_range[0] <= ref_pos <= p_window_range[1]:
                    window_index = each_p_window

        #################### initialize dict ####################

        # ctg level
        if ref_id not in window_to_snv_pos_dod:
            window_to_snv_pos_dod[ref_id] = dict()
        if ref_id not in window_to_reads_dict:
            window_to_reads_dict[ref_id] = dict()
        if ref_id not in read_to_pos_base_dod:
            read_to_pos_base_dod[ref_id] = dict()

        # window level
        if window_index not in window_to_snv_pos_dod[ref_id]:
            window_to_snv_pos_dod[ref_id][window_index] = set()
        if window_index not in window_to_reads_dict[ref_id]:
            window_to_reads_dict[ref_id][window_index] = set()
        if window_index not in read_to_pos_base_dod[ref_id]:
            read_to_pos_base_dod[ref_id][window_index] = dict()

        #########################################################

        cov = int(each_pos_split[3])
        mpileup_str_upper = each_pos_split[4].upper()
        mpileup_str_upper_formatted = format_mpileup_str(mpileup_str_upper)
        if cov == len(mpileup_str_upper_formatted):
            if len(set(mpileup_str_upper_formatted)) > 1:
                count_dict = Counter(mpileup_str_upper_formatted)
                count_dict_filtered = dict()
                for each_var in count_dict:
                    var_num = count_dict[each_var]
                    var_pct = var_num*100/cov
                    if (var_num >= min_var_num) and (var_pct >= min_var_pct):
                        count_dict_filtered[each_var] = var_num
                if len(count_dict_filtered) > 1:
                    read_list = each_pos_split[6].split(',')
                    #print('%s\t%s\t%s\t%sx\t%s\t%s\t%s' % (ref_id, window_index, ref_pos, cov, mpileup_str_upper_formatted, count_dict_filtered, read_list))
                    window_to_snv_pos_dod[ref_id][window_index].add(ref_pos)
                    for read_base, read_id in zip(mpileup_str_upper_formatted, read_list):
                        if read_base in count_dict_filtered:
                            if read_id not in read_to_pos_base_dod[ref_id][window_index]:
                                read_to_pos_base_dod[ref_id][window_index][read_id] = dict()
                            read_to_pos_base_dod[ref_id][window_index][read_id][ref_pos] = read_base
                            window_to_reads_dict[ref_id][window_index].add(read_id)

    # filter pos dict
    window_to_snv_pos_dod_filtered = dict()
    for each_ref in window_to_snv_pos_dod:
        current_ref_pos_dict = window_to_snv_pos_dod[each_ref]
        current_ref_pos_dict_filtered = dict()
        for each_wd in current_ref_pos_dict:
            current_ref_pos_dict_filtered[each_wd] = []
            sorted_list = sorted([i for i in current_ref_pos_dict[each_wd]])
            pos_grouped = group_continuous_number(sorted_list, 2)
            for each_pos_g in pos_grouped:
                if len(each_pos_g) == 1:
                    current_ref_pos_dict_filtered[each_wd].append(each_pos_g[0])
        window_to_snv_pos_dod_filtered[each_ref] = current_ref_pos_dict_filtered

    #################### get windows_to_link_list ####################

    windows_to_link_list = []

    # add inter ref windows to link
    if plot_inter_ref_linkages is True:
        windows_to_link_list = gfa_link_lol_by_window

    # add intra ref windows to link
    if plot_intra_ref_linkages is True:
        for each_ref in ref_seq_window_index_dict:
            if each_ref in refseq_in_mpileup:
                ref_window_list = ref_seq_window_index_dict[each_ref]
                n = 1
                while n < len(ref_window_list):
                    link_to_add = [each_ref, n, each_ref, (n + 1)]
                    windows_to_link_list.append(link_to_add)
                    n += 1

    ##################################################################

    windows_to_process_dict = dict()
    for each_to_link in windows_to_link_list:
        wd_1_ref   = each_to_link[0]
        wd_1_index = each_to_link[1]
        wd_2_ref   = each_to_link[2]
        wd_2_index = each_to_link[3]
        if wd_1_ref not in windows_to_process_dict:
            windows_to_process_dict[wd_1_ref] = set()
        if wd_2_ref not in windows_to_process_dict:
            windows_to_process_dict[wd_2_ref] = set()
        windows_to_process_dict[wd_1_ref].add(wd_1_index)
        windows_to_process_dict[wd_2_ref].add(wd_2_index)

    # get lol for window worker
    window_worker_arg_lol = []
    for each_ref in refseq_in_mpileup:
        current_ref_window_set = windows_to_process_dict[each_ref]

        for each_window in current_ref_window_set:
            msa_file           = '%s/%s_wd%s_msa.aln'          % (aln_dir, each_ref, each_window)
            msa_file_qc        = '%s/%s_wd%s_msa_qc.aln'       % (aln_dir, each_ref, each_window)
            dm_file            = '%s/%s_wd%s_dm.tab'           % (dm_dir, each_ref, each_window)
            elbow_plot         = '%s/%s_wd%s_elbow.png'        % (elbow_plot_dir, each_ref, each_window)
            tree_file          = '%s/%s_wd%s.tree'             % (tree_dir, each_ref, each_window)
            tree_file_with_c   = '%s/%s_wd%s.tree'             % (tree_dir_with_c, each_ref, each_window)
            tree_plot_with_c   = '%s/%s_wd%s_tree.png'         % (tree_plot_dir_with_c, each_ref, each_window)
            pcoa_plot_with_c   = '%s/%s_wd%s_PCoA.png'         % (tree_plot_dir_with_c, each_ref, each_window)
            read_cluster_txt   = '%s/%s_wd%s_read_cluster.txt' % (read_cluster_dir, each_ref, each_window)

            pos_set            = window_to_snv_pos_dod_filtered[each_ref][each_window]
            read_id_set        = window_to_reads_dict[each_ref][each_window]
            read_to_pos_dict   = read_to_pos_base_dod[each_ref][each_window]

            current_window_arg_list = [msa_file, msa_file_qc, dm_file, elbow_plot, tree_file, tree_file_with_c,
                                       tree_plot_with_c, pcoa_plot_with_c, read_cluster_txt, pos_set, read_id_set,
                                       read_to_pos_dict, max_N_pct_in_MSA, tree_format, FastTree_exe,
                                       single_cluster_silhouette_score_cutoff]

            window_worker_arg_lol.append(current_window_arg_list)

    # run window_worker with multiprocessing
    print('Processing %s window with %s cores' % (len(window_worker_arg_lol), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(window_worker, window_worker_arg_lol)
    pool.close()
    pool.join()

    # link windows
    for windows_to_link in windows_to_link_list:
        linking_windows(windows_to_link, read_cluster_dir, linkage_dir, min_link_num_filter)

    #################### visualize linkages ####################

    # combine linkage and filtered linkage files
    os.system('cat %s/*_linkages_sorted.txt > %s'   % (linkage_dir, combined_linkage_txt))
    os.system('cat %s/*_linkages_filtered.txt > %s' % (linkage_dir, combined_linkage_txt_filtered))

    # get all linkages after filtering
    linkage_sets_after_filtering = set()
    for each_linkage in open(combined_linkage_txt_filtered):
        if not each_linkage.startswith('W1\tW2\tLinkage'):
            each_linkage_split = each_linkage.strip().split('\t')
            linkage_str = '%s___%s' % (each_linkage_split[0], each_linkage_split[1])
            linkage_sets_after_filtering.add(linkage_str)

    # add pass info to combined linkage txt
    combined_linkage_txt_with_pass_info_handle = open(combined_linkage_txt_with_pass_info, 'w')
    for each_linkage in open(combined_linkage_txt):
        if not each_linkage.startswith('W1\tW2\tLinkage'):
            each_linkage_split = each_linkage.strip().split('\t')
            linkage_str = '%s___%s' % (each_linkage_split[0], each_linkage_split[1])
            if linkage_str in linkage_sets_after_filtering:
                combined_linkage_txt_with_pass_info_handle.write('%s\ty\n' % each_linkage.strip())
            else:
                combined_linkage_txt_with_pass_info_handle.write('%s\tn\n' % each_linkage.strip())
    combined_linkage_txt_with_pass_info_handle.close()

    # visualize linkage
    visualize_linkage(combined_linkage_txt_with_pass_info, min_link_num_plot, combined_linkage_html)

    # remove tmp files
    os.system('rm %s' % combined_linkage_txt)
    os.system('rm %s' % combined_linkage_txt_filtered)

    ############################################################

if __name__ == '__main__':

    ReadsPhaser()
