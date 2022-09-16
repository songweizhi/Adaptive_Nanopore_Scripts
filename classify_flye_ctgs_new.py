from Bio import SeqIO

'''
# Genbank id --> taxon id (for one ID)
esummary -db nuccore -id BA000012.4 | xtract -pattern DocumentSummary -element Caption,TaxId
esummary -db protein -id CAB1096347.1 | xtract -pattern DocumentSummary -element Caption,TaxId

# Genbank id --> taxon id (for multiple IDs)
cat assembly_blastn_best_hits_ref_id.tab | epost -db nuccore | esummary | xtract -pattern DocumentSummary -element Caption,TaxId > assembly_blastn_best_hits_ref_id_to_taxon_id.txt
cat assembly_blastn_best_hits_ref_id_sorted_uniq.tab | epost -db protein | esummary | xtract -pattern DocumentSummary -element Caption,TaxId > assembly_blastn_best_hits_ref_id_to_taxon_id.txt

# keep only the taxon id column in combined_v5_best_hit_ref_id_uniq_taxon_id.txt

# then
cat assembly_blastn_BestHit_ref_taxon_id_only.txt| sort | uniq > assembly_blastn_BestHit_ref_taxon_id_only_uniq.txt
cat assembly_blastn_best_hits_ref_id_to_taxon_id_only.txt| sort | uniq > assembly_blastn_best_hits_ref_id_to_taxon_id_only_uniq.txt

cd /Users/songweizhi/Desktop/Adaptive_Nanopore/rd1/000
cat combined_v5_best_hit_ref_id_uniq.txt | epost -db nuccore | esummary | xtract -pattern DocumentSummary -element Caption,TaxId > combined_v5_best_hit_ref_id_uniq_taxon_id.txt
# keep only the taxon id column in combined_v5_best_hit_ref_id_uniq_taxon_id.txt
cat combined_v5_best_hit_ref_id_uniq_taxon_id_only.txt| sort | uniq > combined_v5_best_hit_ref_id_uniq_taxon_id_only_uniq.txt

# taxon id --> taxon lineage (upload the uniq ids obtained from above command, select the " full taxid lineage" and save in file)
https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

# taxon lineage (id) --> taxon lineage (scientific name)
get_tax_lineage_by_name.py

'''

'''
cd /Users/songweizhi/Desktop/Adaptive_Nanopore/rd1/000
cat combined_v5_best_hit_ref_id.txt| sort | uniq > combined_v5_best_hit_ref_id_uniq.txt
'''

####################################################### file in ########################################################

assembly_fa                 = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly.fasta'
blast_op                    = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly_blastn_best_hits.tab'

genbank_id_to_taxon_id_txt  = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly_blastn_best_hits_ref_id_to_taxon_id.txt'
# NG_015961	9606
# CP001349	460265

taxon_id_to_lineage_txt     = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/assembly_blastn_best_hits_taxon_id_uniq_tax_report_by_name.txt'
# 1001585	cellular organisms,Bacteria,Proteobacteria,Gammaproteobacteria,Pseudomonadales,Pseudomonadaceae,Pseudomonas,Pseudomonas aeruginosa group,Pseudomonas mendocina,Pseudomonas mendocina NK-01
# 100272	cellular organisms,Eukaryota,environmental samples,uncultured eukaryote

####################################################### file out #######################################################

classification_txt          = '/Users/songweizhi/Desktop/Adaptive_Nanopore/rd2/flye_ctg_classifications.txt'

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

classification_txt_handle = open(classification_txt, 'w')
total = 0
long_aln = 0
eukaryotic_ctg_num = 0
for each_hit in open(blast_op):
    each_hit_split = each_hit.strip().split('\t')
    ctg_id = each_hit_split[0]
    ref_id = each_hit_split[1]
    aln_len = int(each_hit_split[3])
    if aln_len >= 200:
        long_aln += 1

        ref_no_version = ref_id
        if '.' in ref_id:
            ref_no_version = ref_id.split('.')[0]

        ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
        ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')

        if 'Eukaryota' in ref_lineage:
            print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id,ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage ))
            eukaryotic_ctg_num += 1
        classification_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))
    total += 1
classification_txt_handle.close()

print(total)
print(long_aln)
print(eukaryotic_ctg_num)
