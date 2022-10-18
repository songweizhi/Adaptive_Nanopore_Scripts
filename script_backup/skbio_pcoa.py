import pandas as pd
# from skbio.stats.ordination import pcoa
# from skbio.stats.ordination import pcoa_biplot
# from skbio.stats.ordination import OrdinationResults

# pcoa
# http://scikit-bio.org/docs/0.5.7/generated/skbio.stats.ordination.pcoa.html#skbio.stats.ordination.pcoa
# https://github.com/biocore/scikit-bio/blob/0.5.7/skbio/stats/ordination/_principal_coordinate_analysis.py#L23

# OrdinationResults
# http://scikit-bio.org/docs/0.5.7/generated/skbio.stats.ordination.OrdinationResults.html#skbio.stats.ordination.OrdinationResults

dm_file = '/Users/songweizhi/Desktop/Irene/chromosome_op/03_distance_matrix/dm.1.tab'

df = pd.read_table(dm_file, header=0, index_col=0)


index = a_dataframe.index
a_list = list(index)



df_in_arrary = df.to_numpy()
print(df_in_arrary)
print(df.columns)
print(list(df.columns))
print(list(df))
print(len(list(df.columns)))
#print(df)
#pcoa(df)

#OrdinationResults(pcoa(df))


