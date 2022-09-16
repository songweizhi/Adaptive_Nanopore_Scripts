
########################################################################################################################

wd                                          = '/Users/songweizhi/Desktop/Adaptive_Nanopore/check_efficiency/rd1'
log_txt                                     = '%s/adaptive_sampling_FAR63609_7877d00a.csv'                                      % wd
no_decision_and_stop_receiving_reads_txt    = '%s/adaptive_sampling_FAR63609_7877d00a_no_decision_and_stop_receiving_reads.txt' % wd
unblock_reads_txt                           = '%s/adaptive_sampling_FAR63609_7877d00a_unblock_reads.txt'                        % wd
filtered_reads_txt                          = '%s/combined_pass_trimmed50bp.txt'                                                % wd
filtered_reads_txt_unlogged                 = '%s/combined_pass_trimmed50bp_unlogged.txt'                                       % wd

########################################################################################################################

# wd                                          = '/Users/songweizhi/Desktop/Adaptive_Nanopore/check_efficiency/rd2'
# log_txt                                     = '%s/adaptive_sampling_FAS36615_f4f0c7ad.csv'                                      % wd
# no_decision_and_stop_receiving_reads_txt    = '%s/adaptive_sampling_FAS36615_f4f0c7ad_no_decision_and_stop_receiving_reads.txt' % wd
# unblock_reads_txt                           = '%s/adaptive_sampling_FAS36615_f4f0c7ad_unblock_reads.txt'                        % wd
# filtered_reads_txt                          = '%s/combined_pass_rd2_trimmed50bp_noCROP.txt'                                     % wd
# filtered_reads_txt_unlogged                 = '%s/combined_pass_rd2_trimmed50bp_noCROP_unlogged.txt'                            % wd

########################################################################################################################

# wd                                          = '/Users/songweizhi/Desktop/Adaptive_Nanopore/check_efficiency/rd3'
# log_txt                                     = '%s/adaptive_sampling_FAU09901_43970497.csv'                                      % wd
# no_decision_and_stop_receiving_reads_txt    = '%s/adaptive_sampling_FAU09901_43970497_no_decision_and_stop_receiving_reads.txt' % wd
# unblock_reads_txt                           = '%s/adaptive_sampling_FAU09901_43970497_unblock_reads.txt'                        % wd
# filtered_reads_txt                          = '%s/combined_pass_rd3_trimmed50bp_noCROP.txt'                                     % wd
# filtered_reads_txt_unlogged                 = '%s/combined_pass_rd3_trimmed50bp_noCROP_unlogged.txt'                            % wd

########################################################################################################################

unblock_reads_txt_handle = open(unblock_reads_txt, 'w')
no_decision_and_stop_receiving_reads_txt_handle = open(no_decision_and_stop_receiving_reads_txt, 'w')
read_len_dict = {}
for each in open(log_txt):
    if not each.startswith('batch_time,'):
        each_split = each.strip().split(',')
        read_id = each_split[4]
        read_len = int(each_split[5])
        pore_action = each_split[6]
        read_len_dict[read_id] = read_len
        if pore_action == 'unblock':
            unblock_reads_txt_handle.write(read_id + '\n')
        elif pore_action in ['no_decision', 'stop_receiving']:
            no_decision_and_stop_receiving_reads_txt_handle.write(read_id + '\n')
        else:
            print(each_split)
unblock_reads_txt_handle.close()
no_decision_and_stop_receiving_reads_txt_handle.close()

########################################################################################################################

no_decision_and_stop_receiving_reads_set = set()
for each_read in open(no_decision_and_stop_receiving_reads_txt):
    no_decision_and_stop_receiving_reads_set.add(each_read.strip())

unblock_reads_set = set()
for each_read in open(unblock_reads_txt):
    unblock_reads_set.add(each_read.strip())

filtered_reads_set = set()
for each_read in open(filtered_reads_txt):
    filtered_reads_set.add(each_read.strip())

filtered_reads_txt_unlogged_handle = open(filtered_reads_txt_unlogged, 'w')
filtered_in_unblock_set = set()
filtered_in_no_decision_and_stop_receiving_set = set()
all_the_rest_filtered_set = set()
for each_filtered_read in filtered_reads_set:
    if each_filtered_read in unblock_reads_set:
        filtered_in_unblock_set.add(each_filtered_read)
    elif each_filtered_read in no_decision_and_stop_receiving_reads_set:
        filtered_in_no_decision_and_stop_receiving_set.add(each_filtered_read)
    else:
        all_the_rest_filtered_set.add(each_filtered_read)
        filtered_reads_txt_unlogged_handle.write(each_filtered_read + '\n')
filtered_reads_txt_unlogged_handle.close()

print('unblock_reads_set\t%s'                                   % len(unblock_reads_set))
print('no_decision_and_stop_receiving_reads_set\t%s'            % len(no_decision_and_stop_receiving_reads_set))
print('filtered_reads_set\t%s'                                  % len(filtered_reads_set))
print(len(no_decision_and_stop_receiving_reads_set) + len(unblock_reads_set))
print()
print('filtered_in_unblock_set\t%s\t%s'                         % (len(filtered_in_unblock_set),                        len(filtered_in_unblock_set)*100/len(filtered_reads_set)))
print('filtered_in_no_decision_and_stop_receiving_set\t%s\t%s'  % (len(filtered_in_no_decision_and_stop_receiving_set), len(filtered_in_no_decision_and_stop_receiving_set)*100/len(filtered_reads_set)))
print('all_the_rest_filtered_set\t%s\t%s'                       % (len(all_the_rest_filtered_set),                      len(all_the_rest_filtered_set)*100/len(filtered_reads_set)))
print('filtered_reads_set\t%s'                                  % len(filtered_reads_set))


'''
rd1 unblock                         343512
rd1 no_decision_and_stop_receiving  1229641
rd1 total                           1573153

rd2 unblock                         836325
rd2 no_decision_and_stop_receiving  694360
rd2 total                           1530685

rd3 unblock                         1656369
rd3 no_decision_and_stop_receiving  1387249
rd3 total                           3043618


rd1 filtered_in_unblock_set	                        295121	15.924768780824726
rd1 filtered_in_no_decision_and_stop_receiving_set	351828	18.984686113899052
rd1 all_the_rest_filtered_set	                    1206271	65.09054510527622   (adaptative sampling stop working)
rd1 filtered_reads_set	                            1853220

rd2 filtered_in_unblock_set	                        715905	11.496184478508997
rd2 filtered_in_no_decision_and_stop_receiving_set	305743	4.9096988162015585
rd2 all_the_rest_filtered_set	                    5205679	83.59411670528945   (adaptative sampling stop working)
rd2 filtered_reads_set	                            6227327

rd3 filtered_in_unblock_set	                        1437144	35.13910340937404
rd3 filtered_in_no_decision_and_stop_receiving_set	527792	12.90485690135459
rd3 all_the_rest_filtered_set	                    2124935	51.95603968927137   (adaptative sampling stop working)
rd3 filtered_reads_set	                            4089871
'''
