import matplotlib.pyplot as plt


best_hit_tab        = '/Users/songweizhi/Desktop/rd2_first1000_vs_to_reject_BestHit.tab'
best_hit_tab        = '/Users/songweizhi/Desktop/rd2_reads_first_1000bp_vs_rd1_unmapped_reads_BestHit.tab'
aln_len_cutoff_1    = 500
aln_len_cutoff_2    = 750
aln_len_cutoff_3    = 1000


iden_list_1 = []
iden_list_2 = []
iden_list_3 = []
for each in open(best_hit_tab):
    each_split = each.strip().split('\t')
    iden = float(each_split[2])
    aln_len = int(each_split[3])
    qstart = int(each_split[6])
    qend = int(each_split[7])
    sstart = int(each_split[8])
    send = int(each_split[9])
    if aln_len >= aln_len_cutoff_1:
        iden_list_1.append(iden)
    if aln_len >= aln_len_cutoff_2:
        iden_list_2.append(iden)
    if aln_len >= aln_len_cutoff_3:
        iden_list_3.append(iden)

plt.hist(iden_list_1, bins=200, alpha=0.7, color='lightgreen')
plt.hist(iden_list_2, bins=200, alpha=0.7, color='orange')
plt.hist(iden_list_3, bins=200, alpha=0.7, color='lightblue')
plt.show()


