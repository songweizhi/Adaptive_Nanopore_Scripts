import os
from Bio import SeqIO


########################################################################################################################

# file in
reads_fastq                     = '/Users/songweizhi/Desktop/bin.3_reads_min5000bp.fastq'
strainberry_scaffolds_fasta     = '/Users/songweizhi/Desktop/assembly.scaffolds.fa'
strainberry_scaffolds_info_tsv  = '/Users/songweizhi/Desktop/assembly.scaffolds.info.tsv'
modified_hp_txt                 = '/Users/songweizhi/Desktop/hp.txt'
sam_file                        = '/Users/songweizhi/Desktop/assembly.scaffolds.fa.no.secondary.sam'
sam_file_phased_ctgs            = '/Users/songweizhi/Desktop/phased_ctgs.fasta.no.secondary.sam'
modified_hp                     = True

# file out
phased_ctg_fasta                = '/Users/songweizhi/Desktop/phased_ctgs_modified_hp.fasta'
reads_subset_dir                = '/Users/songweizhi/Desktop/reads_subset_dir_modified_hp'
contig_subset_dir               = '/Users/songweizhi/Desktop/contig_subset_dir_modified_hp'

########################################################################################################################

if os.path.isdir(contig_subset_dir) is True:
    os.system('rm -r %s' % contig_subset_dir)
if os.path.isdir(reads_subset_dir) is True:
    os.system('rm -r %s' % reads_subset_dir)
os.system('mkdir %s' % reads_subset_dir)
os.system('mkdir %s' % contig_subset_dir)


hp_set = set()
seq_to_hp_dict = {}
hp_to_seq_dict = {}
if modified_hp is False:
    for each in open(strainberry_scaffolds_info_tsv):
        each_split = each.strip().split('\t')
        seq_id = each_split[0]
        hp_assign = ''
        if 'phased=false' in each.strip():
            hp_assign = 'unphased'
            if 'unphased' not in hp_to_seq_dict:
                hp_to_seq_dict['unphased'] = {seq_id}
            else:
                hp_to_seq_dict['unphased'].add(seq_id)
        else:
            hp_assign = each_split[6].split('=')[1]
            hp_set.add(hp_assign)
            if hp_assign not in hp_to_seq_dict:
                hp_to_seq_dict[hp_assign] = {seq_id}
            else:
                hp_to_seq_dict[hp_assign].add(seq_id)
        seq_to_hp_dict[seq_id] = hp_assign
else:
    for each in open(modified_hp_txt):
        each_split = each.strip().split('\t')
        seq_id = each_split[0]
        hp_assign = each_split[1]

        seq_to_hp_dict[seq_id] = hp_assign

        if hp_assign not in hp_to_seq_dict:
            hp_to_seq_dict[hp_assign] = {seq_id}
        else:
            hp_to_seq_dict[hp_assign].add(seq_id)

        if hp_assign != 'unphased':
            hp_set.add(hp_assign)


phased_ctg_set = set()
phased_ctg_fasta_handle = open(phased_ctg_fasta, 'w')
for each_seq in SeqIO.parse(strainberry_scaffolds_fasta, 'fasta'):
    seq_id = each_seq.id
    seq_hp = seq_to_hp_dict.get(seq_id)
    current_hp_ctg_fa = '%s/hp_%s.fasta' % (contig_subset_dir, seq_hp)
    current_hp_ctg_fa_handle = open(current_hp_ctg_fa, 'a')
    current_hp_ctg_fa_handle.write('>%s\n' % seq_id)
    current_hp_ctg_fa_handle.write('%s\n' % str(each_seq.seq))
    current_hp_ctg_fa_handle.close()
    if seq_hp not in ['unphased', 'NA']:
        phased_ctg_fasta_handle.write('>%s\n' % seq_id)
        phased_ctg_fasta_handle.write('%s\n' % str(each_seq.seq))
        phased_ctg_set.add(seq_id)
phased_ctg_fasta_handle.close()


reads_mapped_to_phased_ctgs = set()
for each_line in open(sam_file_phased_ctgs):
    if not each_line.startswith('@'):
        each_line_split = each_line.strip().split('\t')
        read_id = each_line_split[0]
        ref_id = each_line_split[2]
        if ref_id in phased_ctg_set:
            reads_mapped_to_phased_ctgs.add(read_id)


read_to_hp_dict = {}
for each_line in open(sam_file):
    if not each_line.startswith('@'):
        each_line_split = each_line.strip().split('\t')
        read_id = each_line_split[0]
        ref_id = each_line_split[2]
        ref_hp = seq_to_hp_dict.get(ref_id, 'NA')

        if read_id not in read_to_hp_dict:
            read_to_hp_dict[read_id] = {ref_hp}
        else:
            read_to_hp_dict[read_id].add(ref_hp)


hp_to_reads_dict_uniq = {}
for each_read in read_to_hp_dict:
    mapped_hp_set = read_to_hp_dict[each_read]
    if (len(mapped_hp_set) == 1) and (mapped_hp_set != {'unphased'}) and (mapped_hp_set != {'NA'}):
        hp_id = [i for i in mapped_hp_set][0]
        if hp_id not in hp_to_reads_dict_uniq:
            hp_to_reads_dict_uniq[hp_id] = {each_read}
        else:
            hp_to_reads_dict_uniq[hp_id].add(each_read)


for each_read in SeqIO.parse(reads_fastq, 'fastq'):
    read_id = each_read.id
    mapped_hp_set = read_to_hp_dict[read_id]
    if len(mapped_hp_set) == 1:
        mapped_hp_id = [i for i in mapped_hp_set][0]
        if mapped_hp_id == 'unphased':
            for each_hp in hp_set:
                current_hp_reads_fq = '%s/hp_%s.fastq' % (reads_subset_dir, each_hp)
                current_hp_reads_fq_handle = open(current_hp_reads_fq, 'a')
                SeqIO.write(each_read, current_hp_reads_fq_handle, 'fastq')
                current_hp_reads_fq_handle.close()
        else:
            current_hp_reads_fq = '%s/hp_%s.fastq' % (reads_subset_dir, mapped_hp_id)
            current_hp_reads_fq_handle = open(current_hp_reads_fq, 'a')
            SeqIO.write(each_read, current_hp_reads_fq_handle, 'fastq')
            current_hp_reads_fq_handle.close()
