from pyvis.network import Network


def visualize_linkage(linkage_txt, linkage_html):

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
        net1.add_node(node_1_id, shape='box', borderWidth=0, title='X reads', group=node_1_group)
        net1.add_node(node_2_id, shape='box', borderWidth=0, title='X reads', group=node_2_group)

        # hide edge if link_weight == 0
        hide_edge = True if link_weight == 0 else False

        # https://visjs.github.io/vis-network/docs/network/edges.html
        #net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge, arrows='arrows:%s;%s' % (node_1_id, node_2_id))
        net1.add_edge(node_1_id, node_2_id, value=link_weight, hidden=hide_edge)
    net1.show(linkage_html)


linkage_txt  = '/Users/songweizhi/Desktop/demo/linkages.txt'
linkage_html = '/Users/songweizhi/Desktop/demo/linkages.html'

visualize_linkage(linkage_txt, linkage_html)
