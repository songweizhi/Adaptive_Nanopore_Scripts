import os
from Bio import SeqIO


def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)
    return sorted(folder_list)


refined_mag             = '/Users/songweizhi/Desktop/get_unrefined_counterparts_wd/bin.20.fa'
unrefined_mag_folder    = '/Users/songweizhi/Desktop/get_unrefined_counterparts_wd/unrefined_MAGs'
output_fa               = '/Users/songweizhi/Desktop/get_unrefined_counterparts_wd/bin.20.unrefined.fa'


refined_mag_ctg_set = set()
for each_ctg in SeqIO.parse(refined_mag, 'fasta'):
    refined_mag_ctg_set.add(each_ctg.id)
print('refined_mag_ctg_set (%s):\t%s' % (len(refined_mag_ctg_set), refined_mag_ctg_set))
print()


unrefined_mag_folder_list = get_no_hidden_folder_list(unrefined_mag_folder)


ctg_seq_dict = {}
ctg_to_extract = set()
for bin_subfolder in unrefined_mag_folder_list:

    # get bin file list in each input bin subfolder
    pwd_bin_subfolder = '%s/%s' % (unrefined_mag_folder, bin_subfolder)
    bin_file_list = get_no_hidden_folder_list(pwd_bin_subfolder)

    bin_ext_list = set()
    for each_bin in bin_file_list:
        each_bin_ext = each_bin.split('.')[-1]
        bin_ext_list.add(each_bin_ext)

    if len(bin_ext_list) > 1:
        print('Program exited, please make sure all bins within %s folder have the same extension' % bin_subfolder)
        exit()

    current_unrefined_counterpart_set = set()
    for each_mag in bin_file_list:
        pwd_unrefined_mag = '%s/%s' % (pwd_bin_subfolder, each_mag)
        for each_seq in SeqIO.parse(pwd_unrefined_mag, 'fasta'):
            seq_id = each_seq.id
            seq_seq = str(each_seq.seq)
            ctg_seq_dict[seq_id] = seq_seq
            if seq_id in refined_mag_ctg_set:
                current_unrefined_counterpart_set.add(each_mag)

    for each_counterpart in current_unrefined_counterpart_set:
        pwd_counterpart = '%s/%s/%s' % (unrefined_mag_folder, bin_subfolder, each_counterpart)
        current_counterpart_ctg_set_in_refined = set()
        current_counterpart_ctg_set_not_in_refined= set()
        for each_ctg in SeqIO.parse(pwd_counterpart, 'fasta'):
            ctg_to_extract.add(each_ctg.id)
            if each_ctg.id in refined_mag_ctg_set:
                current_counterpart_ctg_set_in_refined.add(each_ctg.id)
            else:
                current_counterpart_ctg_set_not_in_refined.add(each_ctg.id)
        print('%s\t%s (in refined):     %s' % (bin_subfolder, each_counterpart, current_counterpart_ctg_set_in_refined))
        print('%s\t%s (not in refined): %s' % (bin_subfolder, each_counterpart, current_counterpart_ctg_set_not_in_refined))

print('ctg_to_extract (%s): %s' % (len(ctg_to_extract), ctg_to_extract))


output_fa_handle = open(output_fa, 'w')
for each_ctg in ctg_to_extract:
    output_fa_handle.write('>%s\n' % each_ctg)
    output_fa_handle.write('%s\n' % ctg_seq_dict[each_ctg])
output_fa_handle.close()


