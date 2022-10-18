import argparse
import networkx as nx
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


linkage_txt     = '/Users/songweizhi/Desktop/linkages.txt'
plot_network    = True
node_size       = 6
label_font_size = 6
output_plot     = '/Users/songweizhi/Desktop/linkages.png'


# initialize a graph
G = nx.Graph()

# add node and edge
for each in open(linkage_txt):
    each_split = each.strip().split('\t')
    print(each_split)
    node_1 = each_split[0]
    node_2 = each_split[1]
    link_num = int(each_split[2])

    G.add_node(node_1)
    G.add_node(node_2)
    G.add_edge(node_1, node_2)
    #G.add_edge(node_1, node_2, weight=link_num)





if plot_network is True:

    # specify
    graph_layout = nx.layout.kamada_kawai_layout(G)  # kamada_kawai_layout, planar_layout, fruchterman_reingold_layout

    # turn node attributes into dict
    node_attributes_dict = {}
    for node in G.nodes(data=True):
        node_attributes_dict[node[0]] = node[1]

    # plot node
    for node in G:
        nx.draw_networkx_nodes(G, graph_layout,
                               nodelist=[node],
                               node_size=node_size)

        # add customized node label
        # nx.draw_networkx_labels(g, graph_layout, nodelist=[node], font_size=8, font_color='black')

    #  all nodes label together
    nx.draw_networkx_labels(G, graph_layout, font_size=label_font_size, font_color='black')


    # plot edges
    nx.draw_networkx_edges(G, graph_layout, width=0.5, arrows=True, arrowsize=6)

    # save plot
    plt.savefig(output_plot, dpi=300)
    plt.close()



