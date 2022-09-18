from Bio import SeqIO

fa_file   = '/Users/songweizhi/Desktop/assembly.fasta'
out_short = '/Users/songweizhi/Desktop/assembly_short.fasta'
out_long  = '/Users/songweizhi/Desktop/assembly_long.fasta'

out_short_handle = open(out_short, 'w')
out_long_handle = open(out_long, 'w')
for each_seq in SeqIO.parse(fa_file, 'fasta'):
    if len(each_seq.seq) >= 2500:
        out_long_handle.write('>%s\n' % each_seq.id)
        out_long_handle.write('%s\n' % each_seq.seq)
    else:
        out_short_handle.write('>%s\n' % each_seq.id)
        out_short_handle.write('%s\n' % each_seq.seq)
out_short_handle.close()
out_long_handle.close()
