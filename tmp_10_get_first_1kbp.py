from Bio import SeqIO


fa_file             = '/srv/scratch/z5039045/Adaptive_Nanopore/rd1/adaptive_sampling_FAR63609_7877d00a_no_decision_and_stop_receiving_reads_uniq.fasta'
fa_file_first_1kbp  = '/srv/scratch/z5039045/Adaptive_Nanopore/rd1/adaptive_sampling_FAR63609_7877d00a_no_decision_and_stop_receiving_reads_uniq_first_1kbp.fasta'

fa_file_first_1kbp_handle = open(fa_file_first_1kbp, 'w')
for each_seq in SeqIO.parse(fa_file, 'fasta'):
    fa_file_first_1kbp_handle.write('>%s\n' % each_seq.id)
    if len(each_seq.seq) >= 1000:
        fa_file_first_1kbp_handle.write('%s\n'  % each_seq.seq[:1000])
    else:
        fa_file_first_1kbp_handle.write('%s\n' % each_seq.seq)
fa_file_first_1kbp_handle.close()



fa_file             = '/srv/scratch/z5039045/Adaptive_Nanopore/rd2/adaptive_sampling_FAS36615_f4f0c7ad_no_decision_and_stop_receiving_reads_uniq.fasta'
fa_file_first_1kbp  = '/srv/scratch/z5039045/Adaptive_Nanopore/rd2/adaptive_sampling_FAS36615_f4f0c7ad_no_decision_and_stop_receiving_reads_uniq_first_1kbp.fasta'

fa_file_first_1kbp_handle = open(fa_file_first_1kbp, 'w')
for each_seq in SeqIO.parse(fa_file, 'fasta'):
    fa_file_first_1kbp_handle.write('>%s\n' % each_seq.id)
    if len(each_seq.seq) >= 1000:
        fa_file_first_1kbp_handle.write('%s\n'  % each_seq.seq[:1000])
    else:
        fa_file_first_1kbp_handle.write('%s\n' % each_seq.seq)
fa_file_first_1kbp_handle.close()



fa_file             = '/srv/scratch/z5039045/Adaptive_Nanopore/rd3/adaptive_sampling_FAU09901_43970497_no_decision_and_stop_receiving_reads_uniq.fasta'
fa_file_first_1kbp  = '/srv/scratch/z5039045/Adaptive_Nanopore/rd3/adaptive_sampling_FAU09901_43970497_no_decision_and_stop_receiving_reads_uniq_first_1kbp.fasta'

fa_file_first_1kbp_handle = open(fa_file_first_1kbp, 'w')
for each_seq in SeqIO.parse(fa_file, 'fasta'):
    fa_file_first_1kbp_handle.write('>%s\n' % each_seq.id)
    if len(each_seq.seq) >= 1000:
        fa_file_first_1kbp_handle.write('%s\n'  % each_seq.seq[:1000])
    else:
        fa_file_first_1kbp_handle.write('%s\n' % each_seq.seq)
fa_file_first_1kbp_handle.close()
