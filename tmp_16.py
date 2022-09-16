from Bio import SeqIO

########################################################################################################################

# fa_file     = '/Users/songweizhi/Desktop/assembly_min2500_circular_contigs_no_11044_13150.fa'
# fa_folder   = '/Users/songweizhi/Desktop/assembly_min2500_circular_contigs_no_11044_13150'

# for each_seq in SeqIO.parse(fa_file, 'fasta'):
#     seq_id = each_seq.id
#     # pwd_fa_file = '%s/%s.fa' % (fa_folder, seq_id)
#     # fa_file_handle = open(pwd_fa_file, 'w')
#     # fa_file_handle.write('>%s\n' % seq_id)
#     # fa_file_handle.write('%s\n'  % each_seq.seq)
#     # fa_file_handle.close()
#     blastx_cmd = 'blastx -db /data/bio/blastv5/nr -query %s.fa -out %s_vs_nr.tab -outfmt 6 -num_threads 12' % (seq_id, seq_id)
#     print(blastx_cmd)

########################################################################################################################

# sam_file     = '/Users/songweizhi/Desktop/777/contig_13150.sam'
# reads_id_txt = '/Users/songweizhi/Desktop/777/contig_13150_read_id.txt'
# sam_file     = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/assembly_modified_to_get_archaeal_MAG_reads/ar_ctgs.sam'
# reads_id_txt = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/assembly_modified_to_get_archaeal_MAG_reads/ar_ctgs_read_id.txt'
#
# reads_id_set = set()
# for each_line in open(sam_file):
#     if not each_line.startswith('@'):
#         each_line_split = each_line.strip().split('\t')
#         reads_id_set.add(each_line_split[0])
#
# reads_id_txt_handle = open(reads_id_txt, 'w')
# for each_id in reads_id_set:
#     reads_id_txt_handle.write(each_id + '\n')
# reads_id_txt_handle.close()

########################################################################################################################

# fq_file         = '/Users/songweizhi/Desktop/777/subset.fastq'
# reads_id_txt    = '/Users/songweizhi/Desktop/777/contig_13150_read_id.txt'
# fq_file_subset  = '/Users/songweizhi/Desktop/777/subset_subset.fastq'
#
# fq_file         = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_pass_rd12_trimmed50bp.fastq'
# reads_id_txt    = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/contig_13150_read_id.txt'
# fq_file_subset  = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/contig_13150_reads.fastq'
#
# fq_file         = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_pass_rd12_trimmed50bp.fastq'
# reads_id_txt    = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/contig_11044_read_id.txt'
# fq_file_subset  = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/contig_11044_reads.fastq'

# fq_file         = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_pass_rd12_trimmed50bp.fastq'
# reads_id_txt    = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/assembly_modified_to_get_archaeal_MAG_reads/ar_ctgs_read_id.txt'
# fq_file_subset  = '/srv/scratch/z5039045/Adaptive_Nanopore/combined_rd12_flye_wd/assembly_modified_to_get_archaeal_MAG_reads/ar_ctgs_read.fastq'
#
# read_id_set = set()
# for each_id in open(reads_id_txt):
#     read_id_set.add(each_id.strip())
#
# fq_file_subset_handle = open(fq_file_subset, 'w')
# for each_read in SeqIO.parse(fq_file, 'fastq'):
#     read_id = each_read.id
#     if read_id in read_id_set:
#         SeqIO.write(each_read, fq_file_subset_handle, 'fastq')
# fq_file_subset_handle.close()

########################################################################################################################

# cmd_txt = '/Users/songweizhi/Desktop/prokka_cmds.txt'
# cmd_txt_handle = open(cmd_txt, 'w')
# for each_ctg in open('/Users/songweizhi/Desktop/assembly_circular_ctgs_2500bp_241.txt'):
#     ctg_id = each_ctg.strip()
#     prokka_cmd_1 = 'prokka --cpus 12 --kingdom Archaea --prefix %s_A --locustag %s_A --strain %s_A --outdir %s_A_prokka_wd %s.fa' % (ctg_id, ctg_id, ctg_id, ctg_id, ctg_id)
#     prokka_cmd_2 = 'prokka --cpus 12 --kingdom Bacteria --prefix %s_B --locustag %s_B --strain %s_B --outdir %s_B_prokka_wd %s.fa' % (ctg_id, ctg_id, ctg_id, ctg_id, ctg_id)
#     prokka_cmd_3 = 'prokka --cpus 12 --kingdom Mitochondria --prefix %s_M --locustag %s_M --strain %s_M --outdir %s_M_prokka_wd %s.fa' % (ctg_id, ctg_id, ctg_id, ctg_id, ctg_id)
#     prokka_cmd_4 = 'prokka --cpus 12 --kingdom Viruses --prefix %s_V --locustag %s_V --strain %s_V --outdir %s_V_prokka_wd %s.fa' % (ctg_id, ctg_id, ctg_id, ctg_id, ctg_id)
#     print(prokka_cmd_1)
#     print(prokka_cmd_2)
#     print(prokka_cmd_3)
#     print(prokka_cmd_4)
#     cmd_txt_handle.write(prokka_cmd_1 + '\n')
#     cmd_txt_handle.write(prokka_cmd_2 + '\n')
#     cmd_txt_handle.write(prokka_cmd_3 + '\n')
#     cmd_txt_handle.write(prokka_cmd_4 + '\n')
# cmd_txt_handle.close()

########################################################################################################################
