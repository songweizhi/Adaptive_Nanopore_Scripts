import kmedoids
import pandas as pd
from sklearn.metrics import silhouette_score


def get_num_of_clusters(dm_file):

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

    if max_score < 0.55:
        best_cluster_num = 1

    return best_cluster_num, max_score


dm_file = '/Users/songweizhi/Desktop/TEST_good/03_distance_matrix/edge_20_wd39_dm.tab'
cluster_num, score = get_num_of_clusters(dm_file)
print('Best\t%s\t%s' % (cluster_num, score))


