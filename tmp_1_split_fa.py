from Bio import SeqIO

fa_file = '/Users/songweizhi/Desktop/assembly.fasta'

file_index = 1
wrote_seq_num = 0
file_set = set()
for each in SeqIO.parse(fa_file, 'fasta'):
    sub_file = '/Users/songweizhi/Desktop/999/sub%s.fa' % file_index

    sub_handle = open(sub_file, 'a')
    sub_handle.write('>%s\n' % each.id)
    sub_handle.write('%s\n' % each.seq)
    sub_handle.close()
    wrote_seq_num += 1

    if wrote_seq_num == 300:
        file_index += 1
        wrote_seq_num = 0

    file_name = 'sub%s' % file_index
    file_set.add(file_name)

for each_file in file_set:
    cmd = 'blastn -db /data/bio/blastv5/nt -query %s.fa -outfmt 6 -out %s.tab -num_threads 12' % (each_file, each_file)
    print(cmd)


