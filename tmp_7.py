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


bin_folder = '/Users/songweizhi/Desktop/eukaryotic_MAGs_CONCOCT'
bin_ext    = 'fa'


bin_file_re = '%s/*.%s' % (bin_folder, bin_ext)
bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]


for each_mag in bin_file_list:
    pwd_mag = '%s/%s' % (bin_folder, each_mag)
    bin_size = 0
    for each_ctg in SeqIO.parse(pwd_mag, 'fasta'):
        bin_size += len(each_ctg.seq)
    bin_size_kbp = bin_size/1024
    if bin_size_kbp >= 200:
        #print('%s\t%s' % (each_mag, bin_size_kbp))
        print('cp %s ../eukaryotic_MAGs_CONCOCT_200kbp' % each_mag)


