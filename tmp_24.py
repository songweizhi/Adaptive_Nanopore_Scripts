import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

dm_file = '/Users/songweizhi/Desktop/demo/op_dir/distance_matrix/dm.1.tab'

label_list = []
value_lol  = []
for each in open(dm_file):
    if not each.startswith('\t'):
        each_split = each.strip().split('\t')
        label = each_split[0]
        label_list.append(label)
        value_list = [float(i) for i in each_split[1:]]
        value_lol.append(value_list)

mat = np.array(value_lol)

print(len(mat))
print(len(label_list))

linkage_matrix = linkage(mat, "single")
print(len(linkage_matrix))
print(label_list)
dendrogram(linkage_matrix, orientation='left', labels=label_list, leaf_rotation=0, leaf_font_size=6)
plt.title("test")
plt.show()


