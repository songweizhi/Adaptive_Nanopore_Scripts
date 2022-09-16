
sam_file     = 'assembly_min2500_nonEukaryota.sam'
sam_file_new = 'assembly_min2500_nonEukaryota_only_mapped.sam'

sam_file_new_handle = open(sam_file_new, 'w')
for each in open(sam_file):
    if each.startswith('@'):
        sam_file_new_handle.write(each)
    else:
        each_split = each.strip().split('\t')
        ref_id = each_split[2]
        if ref_id != '*':
            sam_file_new_handle.write(each)
sam_file_new_handle.close()
