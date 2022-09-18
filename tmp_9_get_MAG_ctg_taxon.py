from Bio import SeqIO


classification_txt  = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/flye_ctg_classifications.txt'
mag_file            = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/MyCC_20220427_0832_56mer_0.7_cov/Cluster.37.fasta'
mag_file            = '/Users/songweizhi/Desktop/MyCC_Cluster.23_56mer_0.7_cov/Cluster.9.fasta'

mag_file                    = '/Users/songweizhi/Desktop/untitled_folder/untitled_folder/Cluster.3.fasta'
mag_file_prokaryotic_seqs   = '/Users/songweizhi/Desktop/untitled_folder/untitled_folder/Cluster.3.prokaryotic.seqs.fasta'
mag_file_nonprokaryotic_seqs   = '/Users/songweizhi/Desktop/untitled_folder/untitled_folder/Cluster.3.nonprokaryotic.seqs.fasta'

ctg_id_to_ref_taxon_dict = {}
for each in open(classification_txt):
    each_split = each.strip().split('\t')
    ctg_id = each_split[0]
    ref_taxon = each_split[5]
    ctg_id_to_ref_taxon_dict[ctg_id] = ref_taxon

mag_file_prokaryotic_seqs_handle = open(mag_file_prokaryotic_seqs, 'w')
mag_file_nonprokaryotic_seqs_handle = open(mag_file_nonprokaryotic_seqs, 'w')
seq_num_total = 0
seq_num_eukaryota = 0
seq_num_bacteria = 0
seq_num_archaea = 0
seq_num_viruses = 0
seq_num_unclassified = 0
seq_num_na = 0
for each_seq in SeqIO.parse(mag_file, 'fasta'):
    seq_id = each_seq.id
    seq_ref_taxon = ctg_id_to_ref_taxon_dict.get(seq_id, 'NA')
    seq_num_total += 1
    if 'organisms,Eukaryota' in seq_ref_taxon:
        seq_num_eukaryota += 1
        mag_file_nonprokaryotic_seqs_handle.write('>%s\n' % seq_id)
        mag_file_nonprokaryotic_seqs_handle.write('%s\n' % str(each_seq.seq))
        #print(seq_ref_taxon)
    elif 'organisms,Bacteria' in seq_ref_taxon:
        seq_num_bacteria += 1
        mag_file_prokaryotic_seqs_handle.write('>%s\n' % seq_id)
        mag_file_prokaryotic_seqs_handle.write('%s\n' % str(each_seq.seq))
    elif 'organisms,Archaea,' in seq_ref_taxon:
        seq_num_archaea += 1
        mag_file_prokaryotic_seqs_handle.write('>%s\n' % seq_id)
        mag_file_prokaryotic_seqs_handle.write('%s\n' % str(each_seq.seq))
    elif 'Viruses,' in seq_ref_taxon:
        seq_num_viruses += 1
        mag_file_prokaryotic_seqs_handle.write('>%s\n' % seq_id)
        mag_file_prokaryotic_seqs_handle.write('%s\n' % str(each_seq.seq))
    elif 'unclassified entries,' in seq_ref_taxon:
        seq_num_unclassified += 1
        mag_file_prokaryotic_seqs_handle.write('>%s\n' % seq_id)
        mag_file_prokaryotic_seqs_handle.write('%s\n' % str(each_seq.seq))
    elif seq_ref_taxon == 'NA':
        seq_num_na += 1
        mag_file_nonprokaryotic_seqs_handle.write('>%s\n' % seq_id)
        mag_file_nonprokaryotic_seqs_handle.write('%s\n' % str(each_seq.seq))
    else:
        print(seq_ref_taxon)

mag_file_prokaryotic_seqs_handle.close()
mag_file_nonprokaryotic_seqs_handle.close()

print('seq_num_eukaryota\t%s\t%s'       % (seq_num_eukaryota,       float("{0:.2f}".format(seq_num_eukaryota*100/seq_num_total))))
print('seq_num_bacteria\t%s\t%s'        % (seq_num_bacteria,        float("{0:.2f}".format(seq_num_bacteria*100/seq_num_total))))
print('seq_num_archaea\t%s\t%s'         % (seq_num_archaea,         float("{0:.2f}".format(seq_num_archaea*100/seq_num_total))))
print('seq_num_viruses\t%s\t%s'         % (seq_num_viruses,         float("{0:.2f}".format(seq_num_viruses*100/seq_num_total))))
print('seq_num_unclassified\t%s\t%s'    % (seq_num_unclassified,    float("{0:.2f}".format(seq_num_unclassified*100/seq_num_total))))
print('seq_num_na\t%s\t%s'              % (seq_num_na,              float("{0:.2f}".format(seq_num_na*100/seq_num_total))))

'''
Cluster.30.fasta
seq_num_eukaryota	    378	2.51
seq_num_bacteria	    369	2.45
seq_num_archaea	        6	0.04
seq_num_viruses	        10	0.07
seq_num_unclassified	3	0.02
seq_num_na	            14323	94.92

Cluster.23.fasta
seq_num_eukaryota	408	24.21
seq_num_bacteria	39	2.31
seq_num_archaea	    0	0.0
seq_num_viruses	    0	0.0
seq_num_unclassified	2	0.12
seq_num_na	        1236	73.35
'''