import pandas as pd
from skbio.stats.ordination import pcoa

current_subset_op_dm = '/Users/songweizhi/Desktop/demo/1.tab'
output_plot          = '/Users/songweizhi/Desktop/demo/1.png'


df = pd.read_table(current_subset_op_dm, header=0, index_col=0)


pcoa(df, method='eigh')
