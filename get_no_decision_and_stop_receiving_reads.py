
########################################################################################################################

# file in
log_txt                                     = 'adaptive_sampling_FAU09901_43970497.csv'

# file out
no_decision_and_stop_receiving_reads_txt    = 'no_decision_and_stop_receiving_reads.txt'
unblock_reads_txt                           = 'unblock_reads.txt'


'''
# format of log file
batch_time,read_number,channel,num_samples,read_id,sequence_length,decision
1639358350.380947,117,58,4008,b1052071-672a-4e9c-af4e-764f7041925a,394,no_decision
1639358350.380947,138,59,4014,c1566468-a04d-4a3c-aa7d-17da882227d3,436,unblock

# Once you get the no_decision_and_stop_receiving_reads.txt, you can extract these reads from combined_pass.fastq with:
BioSAK select_seq -seq combined_pass_trimmed50bp.fastq -id no_decision_and_stop_receiving_reads.txt -out no_decision_and_stop_receiving_reads.fastq -option 1 -fq

'''

########################################################################################################################

# parse log file
unblock_reads_set = set()
no_decision_and_stop_receiving_reads_set = set()
for each in open(log_txt):
    if not each.startswith('batch_time,'):
        each_split = each.strip().split(',')
        read_id = each_split[4]
        pore_action = each_split[6]
        if pore_action == 'unblock':
            unblock_reads_set.add(read_id)
        elif pore_action in ['no_decision', 'stop_receiving']:
            no_decision_and_stop_receiving_reads_set.add(read_id)

# write out unblock_reads
unblock_reads_txt_handle = open(unblock_reads_txt, 'w')
for each_unblock_read in unblock_reads_set:
    unblock_reads_txt_handle.write(each_unblock_read + '\n')
unblock_reads_txt_handle.close()

# write out no_decision_and_stop_receiving_reads
no_decision_and_stop_receiving_reads_txt_handle = open(no_decision_and_stop_receiving_reads_txt, 'w')
for each_no_decision_and_stop_receiving_read in no_decision_and_stop_receiving_reads_set:
    no_decision_and_stop_receiving_reads_txt_handle.write(each_no_decision_and_stop_receiving_read + '\n')
no_decision_and_stop_receiving_reads_txt_handle.close()
