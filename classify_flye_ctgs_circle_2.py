from Bio import SeqIO


###################################################### file in/out #####################################################

# file in
assembly_fa                 = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd12/assembly.fasta'
blast_op                    = '/Users/songweizhi/Desktop/against_nr/combined_vs_nr_BestHit.tab'
genbank_id_to_taxon_id_txt  = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd12/assembly_blastn_best_hits_ref_id_to_taxon_id.txt'
taxon_id_to_lineage_txt     = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd12/assembly_blastn_best_hits_taxon_id_uniq_tax_report_by_name.txt'
circular_ctg_id_txt         = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd12/assembly_min2500_531_circular_contigs_no_11044_13150.txt'

# file out
classification_txt          = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd12/flye_ctg_classifications_circular.txt'



# file in
assembly_fa                 = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd12/assembly.fasta'
blast_op                    = '/Users/songweizhi/Desktop/against_nr/combined_vs_nr_BestHit.tab'
genbank_id_to_taxon_id_txt  = '/Users/songweizhi/Desktop/against_nr/assembly_blastn_best_hits_ref_id_to_taxon_id.txt'
taxon_id_to_lineage_txt     = '/Users/songweizhi/Desktop/against_nr/assembly_blastn_best_hits_ref_id_to_taxon_id_only_uniq_tax_report_by_name.txt'
circular_ctg_id_txt         = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd12/assembly_min2500_531_circular_contigs_no_11044_13150.txt'

# file out
classification_txt          = '/Users/songweizhi/Desktop/against_nr/flye_ctg_classifications_circular.txt'



########################################################################################################################

circular_ctg_set = set()
for each_line in open(circular_ctg_id_txt):
    circular_ctg_set.add(each_line.strip())


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

classification_txt_handle = open(classification_txt, 'w')
total = 0
long_aln = 0
eukaryotic_ctg_num = 0
for each_hit in open(blast_op):
    each_hit_split = each_hit.strip().split('\t')
    ctg_id = each_hit_split[0]
    ref_id = each_hit_split[1]
    aln_len = int(each_hit_split[3])

    if ctg_id in circular_ctg_set:

        if aln_len >= 33:
            long_aln += 1

            ref_no_version = ref_id
            if '.' in ref_id:
                ref_no_version = ref_id.split('.')[0]

            ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
            ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')

            if 'Viruses' in ref_lineage:
                print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id,ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage ))
                eukaryotic_ctg_num += 1
            else:
                pass
                #print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id,ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage ))

            classification_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))
        total += 1
classification_txt_handle.close()

print(total)
print(long_aln)
print(eukaryotic_ctg_num)
