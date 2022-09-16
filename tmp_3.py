
classification_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/classifications.txt'
filtered_stats_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/filtered_stats.txt'


ctg_cov_dict = {}
for each in open(filtered_stats_txt):
    each_split = each.strip().split('\t')
    if not each.startswith('#'):
        ctg_cov_dict[each_split[0]] = int(each_split[2])


total_cov = 0
total_cov_Bacteria_Archaea = 0
total_cov_Eukaryota = 0
total_cov_others = 0
for each_line in open(classification_txt):
    each_line_split = each_line.strip().split('\t')
    ctg_len = int(each_line_split[4])
    weighted_cov = ctg_len * ctg_cov_dict[each_line_split[0]]
    total_cov += weighted_cov
    if 'Bacteria' in each_line:
        total_cov_Bacteria_Archaea += weighted_cov
    elif 'Archaea' in each_line:
        total_cov_Bacteria_Archaea += weighted_cov
    elif 'Eukaryota' in each_line:
        total_cov_Eukaryota += weighted_cov
    else:
        total_cov_others += weighted_cov


print(total_cov_Bacteria_Archaea)
print(total_cov_Eukaryota)
print(total_cov_others)
print(total_cov_Bacteria_Archaea*100/total_cov)
print(total_cov_Eukaryota*100/total_cov)
print(total_cov_others*100/total_cov)

