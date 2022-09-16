from Bio import AlignIO
from collections import Counter
from Bio.Phylo.TreeConstruction import DistanceCalculator


def format_mpileup_str(mpileup_str):
    mpileup_str_new = ''
    to_ignore  = 0
    for each in mpileup_str:
        if to_ignore == 0:
            if each in ['A', 'T', 'C', 'G', '*']:
                mpileup_str_new += each
            elif each in ['+', '-']:
                mpileup_str_new += ''
            elif each.isnumeric() is True:
                to_ignore = int(each)
                mpileup_str_new += ''
        else:
            to_ignore -= 1
            mpileup_str_new += ''
    return mpileup_str_new


def msa_to_distance_matrix(msa_file, msa_format):
    with open(msa_file) as msa_file_opened:
        aln = AlignIO.read(msa_file_opened, msa_format)
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        return dm


########################################################################################################################

mpileup_file    = '/Users/songweizhi/Desktop/demo.mpileup'
window_size     = 5000
min_var_num     = 5
min_var_pct     = 20
txt_out_folder  = '/Users/songweizhi/Desktop/demo'
msa_format      = 'fasta'

########################################################################################################################

subset_index_set = set()
read_to_pos_dod = dict()
read_id_dict = dict()
pos_dict = dict()
current_pos = 1
for each_pos in open(mpileup_file):

    current_subset_index = (current_pos - 1)//window_size + 1
    current_pos += 1

    if current_subset_index not in read_to_pos_dod:
        read_to_pos_dod[current_subset_index] = dict()
    if current_subset_index not in read_id_dict:
        read_id_dict[current_subset_index] = set()
    if current_subset_index not in pos_dict:
        pos_dict[current_subset_index] = set()

    subset_index_set.add(current_subset_index)

    each_pos_split = each_pos.strip().split('\t')
    ref_id = each_pos_split[0]
    ref_pos = int(each_pos_split[1])
    cov = int(each_pos_split[3])
    mpileup_str_upper = each_pos_split[4].upper()
    mpileup_str_upper_formatted = format_mpileup_str(mpileup_str_upper)
    if cov == len(mpileup_str_upper_formatted):
        if len(set(mpileup_str_upper_formatted)) > 1:
            count_dict = Counter(mpileup_str_upper_formatted)
            count_dict_filtered = dict()
            for each_var in count_dict:
                var_num = count_dict[each_var]
                var_pct = var_num*100/cov
                if (var_num >= min_var_num) and (var_pct >= min_var_pct):
                    count_dict_filtered[each_var] = var_num
            if len(count_dict_filtered) > 1:
                read_list = each_pos_split[6].split(',')
                # print('%s\t%s\t%s\t%s\t%s\t%s' % (ref_id, ref_pos, cov, mpileup_str_upper_formatted, count_dict_filtered, read_list))
                pos_dict[current_subset_index].add(ref_pos)
                for read_base, read_id in zip(mpileup_str_upper_formatted, read_list):
                    if read_base in count_dict_filtered:
                        # print('%s\t%s' % (read_base, read_id))
                        if read_id not in read_to_pos_dod[current_subset_index]:
                            read_to_pos_dod[current_subset_index][read_id] = dict()
                        read_to_pos_dod[current_subset_index][read_id][ref_pos] = read_base
                        read_id_dict[current_subset_index].add(read_id)


for each_subset in sorted([i for i in subset_index_set]):

    current_subset_op_txt           = '%s/%s.aln' % (txt_out_folder, each_subset)
    current_subset_pos_set          = pos_dict[each_subset]
    current_subset_read_id_set      = read_id_dict[each_subset]
    current_subset_read_to_pos_dict = read_to_pos_dod[each_subset]

    # sort
    pos_list_sorted     = sorted([i for i in current_subset_pos_set])
    read_id_list_sorted = sorted([i for i in current_subset_read_id_set])

    # write out
    txt_out_handle = open(current_subset_op_txt, 'w')
    for each_read in read_id_list_sorted:
        current_read_pos_base_dict = current_subset_read_to_pos_dict.get(each_read, {})
        current_read_str = ''
        for each_pos in pos_list_sorted:
            current_pos_base = current_read_pos_base_dict.get(each_pos, 'N')
            current_read_str += current_pos_base
        txt_out_handle.write('>%s\n%s\n' % (each_read, current_read_str.replace('*', '-')))
    txt_out_handle.close()



msa_file = '/Users/songweizhi/Desktop/demo/1.aln'
msa_format = 'fasta'
dm = msa_to_distance_matrix(msa_file, msa_format)
print(dm)