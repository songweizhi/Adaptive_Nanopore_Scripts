import os
import sys
import argparse
import kmedoids
import pandas as pd
from ete3 import Tree
from Bio import Phylo
from io import StringIO
from Bio import AlignIO
from collections import Counter
from pyvis.network import Network
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


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


def msa_to_distance_matrix(msa_file, msa_format, op_matrix):

    with open(msa_file) as msa_file_opened:
        aln = AlignIO.read(msa_file_opened, msa_format)
        calculator = DistanceCalculator('blastn')
        # Available models: identity, benner22, benner6, benner74, dayhoff, feng, genetic, gonnet1992,
        # hoxd70, johnson, jones, levin, mclachlan, mdm78, blastn, rao, risler, schneider, str, trans,
        # blosum45, blosum50, blosum62, blosum80, blosum90, pam250, pam30, pam70
        dm = calculator.get_distance(aln)

        dm_name_list = dm.names
        op_matrix_handle = open(op_matrix, 'w')
        op_matrix_handle.write('\t%s\n' % '\t'.join(dm_name_list))
        for each_name in dm_name_list:
            op_matrix_handle.write('%s\t%s\n' % (each_name, '\t'.join([str(i) for i in dm[each_name]])))
        op_matrix_handle.close()


def tree_to_distance_matrix(tree_file, op_matrix):

    # read in tree
    tree_in = Tree(tree_file, format=0)

    # get leaf node list
    leaf_node_list = []
    for leaf_node in tree_in:
        leaf_node_list.append(leaf_node.name)
    leaf_node_list = sorted(leaf_node_list)

    # get list of distance list
    all_distances_lol = []
    op_matrix_handle = open(op_matrix, 'w')
    op_matrix_handle.write('\t%s\n' % '\t'.join(leaf_node_list))
    for each_node_1 in leaf_node_list:
        current_node_distance_list = []
        for each_node_2 in leaf_node_list:
            distance = 0
            if each_node_1 != each_node_2:
                distance = get_node_distance(tree_in, each_node_1, each_node_2)
                distance = float("{0:.3f}".format(distance))
            current_node_distance_list.append(str(distance))
        all_distances_lol.append(current_node_distance_list)
        op_matrix_handle.write('%s\t%s\n' % (each_node_1, '\t'.join([str(i) for i in current_node_distance_list])))
    op_matrix_handle.close()


def get_elbow_plot(dm_txt, output_elbow_plot, bootstrap_num):

    df = pd.read_table(dm_txt, header=0, index_col=0)
    cluster_list     = list(range(1, len(df) + 1))
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
    plt.tight_layout()
    plt.savefig(output_elbow_plot, bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()


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


def visualize_linkage(linkage_txt, min_link_num, linkage_html):

    edges_list = []
    for each in open(linkage_txt):
        each_split = each.strip().split('\t')
        node_1 = each_split[0]
        node_2 = each_split[1]
        link_num = int(each_split[2])
        edge_tuple = tuple([node_1, node_2, link_num])
        edges_list.append(edge_tuple)

    net1 = Network(height='100%', width='100%', bgcolor='white', font_color='black', directed=False)
    # net1 = Network(height='750px', width='100%', bgcolor='#222222', font_color='white')

    for each_link in edges_list:
        node_1_id = each_link[0]
        node_2_id = each_link[1]
        link_weight = each_link[2]
        node_1_group = node_1_id.split('__')[0]
        node_2_group = node_2_id.split('__')[0]
        net1.add_node(node_1_id, shape='box', group=node_1_group)
        net1.add_node(node_2_id, shape='box', group=node_2_group)

        # hide edge if link_weight == 0
        hide_edge = False if link_weight >= min_link_num else True

        # https://visjs.github.io/vis-network/docs/network/edges.html
        #net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge, arrows='arrows:%s;%s' % (node_1_id, node_2_id))
        #net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge)
        net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge)

    net1.show(linkage_html)


def get_node_distance(tree, node_1, node_2):
    distance = tree.get_distance(node_1, node_2)
    return distance


def get_distance_matrix(tree_file):

    # read in tree
    tree_in = Tree(tree_file, format=0)

    # get leaf node list
    leaf_node_list = []
    for leaf_node in tree_in:
        leaf_node_list.append(leaf_node.name)
    leaf_node_list = sorted(leaf_node_list)

    # get list of distance list
    all_distances_lol = []
    for each_node_1 in leaf_node_list:
        current_node_distance_list = []
        for each_node_2 in leaf_node_list:
            distance = 0
            if each_node_1 != each_node_2:
                distance = get_node_distance(tree_in, each_node_1, each_node_2)
                distance = float("{0:.5f}".format(distance))
            current_node_distance_list.append(str(distance))
        all_distances_lol.append(current_node_distance_list)

    return all_distances_lol, leaf_node_list


def ReadsPhaser(arg):

    ########################################################################################################################

    args = {'i'           : '/Users/songweizhi/Desktop/contig_65.mpileup.subset',
            'o'           : '/Users/songweizhi/Desktop/demo',
            'ws'          : 8000,
            'min_var_num' : 15,  # 15
            'min_var_pct' : 20,  # 20
            'max_n_pct'   : 10,
            'min_link_num': 3}

    cluster_num = 2

    ########################################################################################################################

    mpileup_file =          arg['i']
    output_dir =            arg['o']
    window_size =           arg['ws']
    min_var_num =           arg['min_var_num']
    min_var_pct =           arg['min_var_pct']
    max_N_pct_in_MSA =      arg['max_n_pct']
    min_link_num_to_plot =  arg['min_link_num']
    FastTree_exe            = 'fasttree'
    msa_format              = 'fasta'
    tree_format             = 1

    aln_dir                 = '%s/01_MSAs'                      % output_dir
    aln_dir_qc              = '%s/01_MSAs_QC'                   % output_dir
    tree_dir                = '%s/02_tree'                      % output_dir
    tree_plot_dir           = '%s/02_tree_plot'                 % output_dir
    dm_dir                  = '%s/03_distance_matrix'           % output_dir
    tree_dir_with_c         = '%s/04_tree_with_cluster'         % output_dir
    tree_plot_dir_with_c    = '%s/05_tree_with_cluster_plot'    % output_dir
    elbow_plot_dir          = '%s/06_elbow_plots'               % output_dir
    read_cluster_dir        = '%s/07_read_cluster'              % output_dir
    linkage_txt             = '%s/linkages.txt'                 % output_dir
    linkage_html            = '%s/linkages.html'                % output_dir

    ########################################################################################################################

    if os.path.isdir(output_dir) is True:
        os.system('rm -r %s' % output_dir)
    os.system('mkdir %s' % output_dir)
    os.system('mkdir %s' % aln_dir)
    os.system('mkdir %s' % aln_dir_qc)
    os.system('mkdir %s' % dm_dir)
    os.system('mkdir %s' % elbow_plot_dir)
    os.system('mkdir %s' % tree_dir)
    os.system('mkdir %s' % read_cluster_dir)
    os.system('mkdir %s' % tree_dir_with_c)
    os.system('mkdir %s' % tree_plot_dir_with_c)

    ########################################################################################################################

    current_pos = 1
    pos_dict = dict()
    read_id_dict = dict()
    subset_index_set = set()
    read_to_pos_dod = dict()
    for each_pos in open(mpileup_file):

        current_subset_index = (current_pos - 1)//window_size + 1
        current_pos += 1

        if current_subset_index not in read_to_pos_dod:
            read_to_pos_dod[current_subset_index] = dict()
        if current_subset_index not in read_id_dict:
            read_id_dict[current_subset_index] = set()
        if current_subset_index not in pos_dict:
            pos_dict[current_subset_index] = set()

        subset_index_set.add(current_subset_index)

        each_pos_split = each_pos.strip().split('\t')
        ref_id = each_pos_split[0]
        ref_pos = int(each_pos_split[1])
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
                    #print('%s\t%s\t%s\t%sx\t%s\t%s\t%s' % (ref_id, current_subset_index, ref_pos, cov, mpileup_str_upper_formatted, count_dict_filtered, read_list))
                    pos_dict[current_subset_index].add(ref_pos)
                    for read_base, read_id in zip(mpileup_str_upper_formatted, read_list):
                        if read_base in count_dict_filtered:
                            if read_id not in read_to_pos_dod[current_subset_index]:
                                read_to_pos_dod[current_subset_index][read_id] = dict()
                            read_to_pos_dod[current_subset_index][read_id][ref_pos] = read_base
                            read_id_dict[current_subset_index].add(read_id)

    subset_index_list_sorted = sorted([i for i in subset_index_set])

    # needs updates (to remove this line) !!!
    subset_index_list_sorted = subset_index_list_sorted[:-1]

    for each_subset in subset_index_list_sorted:

        print('Processing window %s (%s-%sbp)' % (each_subset, ((each_subset - 1)*window_size), (each_subset*window_size)))

        current_subset_op_msa               = '%s/msa_%s.aln'               % (aln_dir, each_subset)
        current_subset_op_msa_qc            = '%s/msa_%s_qc.aln'            % (aln_dir_qc, each_subset)
        current_subset_op_dm                = '%s/dm.%s.tab'                % (dm_dir, each_subset)
        current_subset_elbow_plot           = '%s/elbow_%s.png'             % (elbow_plot_dir, each_subset)
        current_subset_tree_newick          = '%s/tree_%s.tree'             % (tree_dir, each_subset)
        current_subset_tree_newick_with_c   = '%s/tree_%s.tree'             % (tree_dir_with_c, each_subset)
        current_subset_tree_plot_with_c     = '%s/tree_plot_%s.png'         % (tree_plot_dir_with_c, each_subset)
        current_subset_read_cluster         = '%s/read_cluster_%s.txt'      % (read_cluster_dir, each_subset)
        current_subset_pos_set              = pos_dict[each_subset]
        current_subset_read_id_set          = read_id_dict[each_subset]
        current_subset_read_to_pos_dict     = read_to_pos_dod[each_subset]

        # sort
        pos_list_sorted     = sorted([i for i in current_subset_pos_set])
        read_id_list_sorted = sorted([i for i in current_subset_read_id_set])

        # write out
        current_subset_op_msa_handle = open(current_subset_op_msa, 'w')
        current_subset_op_msa_qc_handle = open(current_subset_op_msa_qc, 'w')
        for each_read in read_id_list_sorted:
            current_read_pos_base_dict = current_subset_read_to_pos_dict.get(each_read, {})
            current_read_str = ''
            for each_pos in pos_list_sorted:
                current_pos_base = current_read_pos_base_dict.get(each_pos, 'N')
                current_read_str += current_pos_base

            current_read_str_reformatted = current_read_str.replace('*', '-')
            current_subset_op_msa_handle.write('>%s\n%s\n' % (each_read, current_read_str_reformatted))
            n_pct = current_read_str_reformatted.count('N')*100/len(current_read_str_reformatted)
            if n_pct <= max_N_pct_in_MSA:
                current_subset_op_msa_qc_handle.write('>%s\n%s\n' % (each_read, current_read_str_reformatted))
        current_subset_op_msa_handle.close()
        current_subset_op_msa_qc_handle.close()

        # get tree
        fasttree_cmd = '%s -nt -quiet -out %s %s 2>/dev/null' % (FastTree_exe, current_subset_tree_newick, current_subset_op_msa_qc)
        # print(fasttree_cmd)
        os.system(fasttree_cmd)

        # get distance matrix
        #msa_to_distance_matrix(current_subset_op_msa, msa_format, current_subset_op_dm)
        tree_to_distance_matrix(current_subset_tree_newick, current_subset_op_dm)

        # get elbow plot
        get_elbow_plot(current_subset_op_dm, current_subset_elbow_plot, 50)

        # perform clustering
        df = pd.read_table(current_subset_op_dm, header=0, index_col=0)
        clusters = kmedoids.fasterpam(df, cluster_num)
        cluster_label_list = clusters.labels
        cluster_label_list_1_based = [int((i + 1)) for i in cluster_label_list]
        df_row_name_list = list(df.index.values)

        # write out clustering results
        read_to_cluster_dict = {}
        current_subset_read_cluster_handle = open(current_subset_read_cluster, 'w')
        for c, r in zip(cluster_label_list_1_based, df_row_name_list):
            current_subset_read_cluster_handle.write('c%s\t%s\n' % (c, r))
            read_to_cluster_dict[r] = c
        current_subset_read_cluster_handle.close()

        # add cluster info to tree leaves and plot
        t = Tree(current_subset_tree_newick, format=tree_format)
        for leaf in t:
            leaf_cluster = read_to_cluster_dict.get(leaf.name, '')
            leaf_name_new = '%s__%s' % (leaf_cluster, leaf.name)
            leaf.name = leaf_cluster
        t.write(format=tree_format, outfile=current_subset_tree_newick_with_c)

        # plot tree
        #plottree(current_subset_tree_newick, False, False, False, False, None, None, None, None, None, current_subset_tree_plot)
        plottree(current_subset_tree_newick_with_c, False, False, False, False, None, None, None, None, None, current_subset_tree_plot_with_c)

    # read cluster info into dict
    window_to_cluster_dict = dict()
    for each_subset in subset_index_list_sorted:
        window_to_cluster_dict[each_subset] = dict()
        pwd_cluster_txt = '%s/read_cluster_%s.txt' % (read_cluster_dir, each_subset)
        for each_r in open(pwd_cluster_txt):
            each_r_split = each_r.strip().split('\t')
            current_c = each_r_split[0]
            current_r = each_r_split[1]
            if current_c not in window_to_cluster_dict[each_subset]:
                window_to_cluster_dict[each_subset][current_c] = {current_r}
            else:
                window_to_cluster_dict[each_subset][current_c].add(current_r)

    # get linkages
    linkage_txt_handle = open(linkage_txt, 'w')
    n = 1
    while n < len(subset_index_list_sorted):
        read_cluster_l = window_to_cluster_dict[n]
        read_cluster_r = window_to_cluster_dict[n + 1]
        for lc in read_cluster_l:
            for rc in read_cluster_r:
                lc_read_set = read_cluster_l[lc]
                rc_read_set = read_cluster_r[rc]
                lc_rc_shared = set(lc_read_set).intersection(rc_read_set)
                linkage_txt_handle.write('w%s__%s\tw%s__%s\t%s\n' % (n, lc, (n + 1), rc, len(lc_rc_shared)))
        n += 1
    linkage_txt_handle.close()

    # visualize linkage
    visualize_linkage(linkage_txt, min_link_num_to_plot, linkage_html)


if __name__ == '__main__':

    # ReadsPhaser_parser = argparse.ArgumentParser()
    # ReadsPhaser_parser.add_argument('-i',               required=True,                           help='input mpileup file')
    # ReadsPhaser_parser.add_argument('-o',               required=True,                           help='output directory')
    # ReadsPhaser_parser.add_argument('-ws',              required=True, type=int,                 help='window size (in bp)')
    # ReadsPhaser_parser.add_argument('-min_var_num',     required=True, type=int,                 help='min_var_num')
    # ReadsPhaser_parser.add_argument('-min_var_pct',     required=True, type=float,               help='min_var_pct (1-100)')
    # ReadsPhaser_parser.add_argument('-max_n_pct',       required=False, type=float,  default=75, help='max_N_pct_in_MSA (1-100), default: 25')
    # ReadsPhaser_parser.add_argument('-min_link_num',    required=False, type=int,    default=3,  help='min_link_num_to_plot, default: 3')
    # args = vars(ReadsPhaser_parser.parse_args())
    args = dict
    ReadsPhaser(args)

'''

cd /Users/songweizhi/Desktop
python3 ~/PycharmProjects/Adaptive_Nanopore/ReadsPhaser.py -i contig_65.mpileup.subset -o demo -ws 8000 -min_var_num 15 -min_var_pct 20

'''
