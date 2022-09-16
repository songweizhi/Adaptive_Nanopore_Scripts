from Bio import SeqIO


rd1_bin24_29_fasta      = '/Users/songweizhi/Desktop/compare_giant_eukaryotic_MAGs/rd1_bin24_29.fasta'
rd2_bin28_31_fasta      = '/Users/songweizhi/Desktop/compare_giant_eukaryotic_MAGs/rd2_bin28_31.fasta'
blast_op                = '/Users/songweizhi/Desktop/compare_giant_eukaryotic_MAGs/rd2_bin28_31_vs_rd1_bin24_29.tab'
cov_cutoff_list         = [90, 50, 30, 20, 10, 5]
iden_cutoff             = 99


'''
98	90	10.7	13.98
98	50	23.27	30.27
98	30	29.9	38.68
98	20	33.11	42.5
98	10	36.26	45.87
98	5	37.62	47.2

99	90	5.7	7.5
99	50	12.05	15.82
99	30	15.1	19.8
99	20	16.59	21.64
99	10	18.08	23.27
99	5	18.72	23.96

'''

total_q_len = 0
query_seq_len_dict = {}
for each_seq in SeqIO.parse(rd2_bin28_31_fasta, 'fasta'):
    query_seq_len_dict[each_seq.id] = len(each_seq.seq)
    total_q_len += len(each_seq.seq)

total_s_len = 0
subject_seq_len_dict = {}
for each_seq in SeqIO.parse(rd1_bin24_29_fasta, 'fasta'):
    subject_seq_len_dict[each_seq.id] = len(each_seq.seq)
    total_s_len += len(each_seq.seq)

for cov_cutoff in cov_cutoff_list:
    matched_query_dict = {}
    matched_subject_dict = {}
    for each in open(blast_op):
        each_split = each.strip().split('\t')
        iden = float(each_split[2])
        aln_len = int(each_split[3])
        query_len = query_seq_len_dict[each_split[0]]
        subject_len = subject_seq_len_dict[each_split[1]]
        query_cov = aln_len*100/query_len
        subject_cov = aln_len*100/subject_len
        query_cov = float("{0:.2f}".format(query_cov))
        subject_cov = float("{0:.2f}".format(subject_cov))
        qstart = int(each_split[6])
        qend = int(each_split[7])
        sstart = int(each_split[8])
        send = int(each_split[9])

        if iden >= iden_cutoff:

            if each_split[0] not in matched_query_dict:
                matched_query_dict[each_split[0]] = set()

            if each_split[1] not in matched_subject_dict:
                matched_subject_dict[each_split[1]] = set()

            if (query_cov >= cov_cutoff) or (subject_cov >= cov_cutoff):
                #print('%s(%s-%s)(%s)\t%s(%s-%s)(%s)\t%s\t%s' % (each_split[0], qstart, qend, query_cov, each_split[1], sstart, send, subject_cov, iden, aln_len))
                for each_q_bp in list(range(min([qstart, qend]), (max([qstart, qend]) + 1))):
                    matched_query_dict[each_split[0]].add(each_q_bp)

                for each_s_bp in list(range(min([sstart, send]), (max([sstart, send]) + 1))):
                    matched_subject_dict[each_split[1]].add(each_s_bp)

    matched_query_total_len = 0
    for each_q in matched_query_dict:
        matched_query_total_len += len(matched_query_dict[each_q])

    matched_subject_total_len = 0
    for each_s in matched_subject_dict:
        matched_subject_total_len += len(matched_subject_dict[each_s])

    matched_query_pct   = float("{0:.2f}".format(matched_query_total_len*100/total_q_len))
    matched_subject_pct = float("{0:.2f}".format(matched_subject_total_len*100/total_s_len))

    print('%s\t%s\t%s\t%s' % (iden_cutoff, cov_cutoff, matched_query_pct, matched_subject_pct))
