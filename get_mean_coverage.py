
stats_file = '/Users/songweizhi/Desktop/filtered_stats.txt'

total_len = 0
total_cov = 0
for each in open(stats_file):
    if not each.startswith('#'):
        each_split = each.strip().split('\t')
        print(each_split)
        ctg_len = int(each_split[1])
        ctg_cov = int(each_split[2])
        total_len += ctg_len
        total_cov += (ctg_len*ctg_cov)

print(total_cov/total_len)
print(total_cov/(1024*1024*1024))