from Bio import SeqIO

####################################################### file in ########################################################

# rd1
assembly_fa                 = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd1/000/assembly.fasta'
blast_op                    = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd1/000/combined_v5_best_hit.tab'
genbank_id_to_taxon_id_txt  = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd1/000/combined_v5_best_hit_ref_id_uniq_taxon_id.txt'
taxon_id_to_lineage_txt     = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd1/000/combined_v5_best_hit_ref_id_uniq_taxon_id_only_uniq_tax_report_by_name.txt'

'''
Category	Contig_num	Contig_len	Contig_num_pct	Contig_len_pct	
Eukaryota	701	14749221	17.69	12.44
Bacteria	3187	100663825	80.44	84.88
Archaea	34	2154460	0.86	1.82
Viruses	2	42640	0.05	0.04
Unclassified	18	539727	0.45	0.46
'''

# # rd2
# assembly_fa                 = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly.fasta'
# blast_op                    = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly_blastn_BestHit.tab'
# genbank_id_to_taxon_id_txt  = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly_blastn_BestHit_ref_taxon_id.txt'
# taxon_id_to_lineage_txt     = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly_blastn_BestHit_ref_taxon_id_only_uniq_tax_report_by_name.txt'

'''
Category	Contig_num	Contig_len	Contig_num_pct	Contig_len_pct	
Eukaryota	964	33853003	19.47	22.75
Bacteria	3875	110907318	78.28	74.52
Archaea	21	2125133	0.42	1.43
Viruses	28	710724	0.57	0.48
Unclassified	45	684293	0.91	0.46
'''

########################################################################################################################

# get seq_len_dict
seq_len_dict = {}
for each_seq in SeqIO.parse(assembly_fa, 'fasta'):
    seq_len_dict[each_seq.id] = len(each_seq.seq)

# read in genbank id to taxon id
genbank_id_to_taxon_id_dict = {}
for each_line in open(genbank_id_to_taxon_id_txt):
    each_line_split = each_line.strip().split('\t')
    genbank_id_to_taxon_id_dict[each_line_split[0]] = each_line_split[1]

# read in taxon id to lineage into dict
taxon_id_to_lineage_dict = {}
for each_line in open(taxon_id_to_lineage_txt):
    each_line_split = each_line.strip().split('\t')
    taxon_id_to_lineage_dict[each_line_split[0]] = each_line_split[1]

seq_num_total = 0
seq_num_eukaryota = 0
seq_num_bacteria = 0
seq_num_archaea = 0
seq_num_viruses = 0
seq_num_unclassified = 0
seq_len_total = 0
seq_len_eukaryota = 0
seq_len_bacteria = 0
seq_len_archaea = 0
seq_len_viruses = 0
seq_len_unclassified = 0
for each_hit in open(blast_op):
    each_hit_split = each_hit.strip().split('\t')
    ctg_id = each_hit_split[0]
    ref_id = each_hit_split[1]
    aln_len = int(each_hit_split[3])
    if aln_len >= 200:
        ref_no_version = ref_id
        if '.' in ref_id:
            ref_no_version = ref_id.split('.')[0]

        ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
        ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')
        seq_len = seq_len_dict[ctg_id]
        seq_num_total += 1
        seq_len_total += seq_len
        if 'organisms,Eukaryota' in ref_lineage:
            seq_num_eukaryota += 1
            seq_len_eukaryota += seq_len

        elif 'organisms,Bacteria' in ref_lineage:
            seq_num_bacteria += 1
            seq_len_bacteria += seq_len

        elif 'organisms,Archaea,' in ref_lineage:
            seq_num_archaea += 1
            seq_len_archaea += seq_len

        elif 'Viruses,' in ref_lineage:
            seq_num_viruses += 1
            seq_len_viruses += seq_len

        elif 'unclassified entries,' in ref_lineage:
            seq_num_unclassified += 1
            seq_len_unclassified += seq_len
        else:
            print('%s\t%s' % (ref_lineage, each_hit_split))

print('Category\tContig_num\tContig_len\tContig_num_pct\tContig_len_pct\t')
print('Eukaryota\t%s\t%s\t%s\t%s'       % (seq_num_eukaryota,       seq_len_eukaryota,       float("{0:.2f}".format(seq_num_eukaryota*100/seq_num_total)),      float("{0:.2f}".format(seq_len_eukaryota*100/seq_len_total))))
print('Bacteria\t%s\t%s\t%s\t%s'        % (seq_num_bacteria,        seq_len_bacteria,        float("{0:.2f}".format(seq_num_bacteria*100/seq_num_total)),       float("{0:.2f}".format(seq_len_bacteria*100/seq_len_total))))
print('Archaea\t%s\t%s\t%s\t%s'         % (seq_num_archaea,         seq_len_archaea,         float("{0:.2f}".format(seq_num_archaea*100/seq_num_total)),        float("{0:.2f}".format(seq_len_archaea*100/seq_len_total))))
print('Viruses\t%s\t%s\t%s\t%s'         % (seq_num_viruses,         seq_len_viruses,         float("{0:.2f}".format(seq_num_viruses*100/seq_num_total)),        float("{0:.2f}".format(seq_len_viruses*100/seq_len_total))))
print('Unclassified\t%s\t%s\t%s\t%s'    % (seq_num_unclassified,    seq_len_unclassified,    float("{0:.2f}".format(seq_num_unclassified*100/seq_num_total)),   float("{0:.2f}".format(seq_len_unclassified*100/seq_len_total))))
