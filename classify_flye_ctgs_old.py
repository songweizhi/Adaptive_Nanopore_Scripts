from Bio import SeqIO


assembly_fa                     = '/Users/songweizhi/Desktop/Adaptive_Nanopore/assembly.fasta'
best_hit_txt                    = '/Users/songweizhi/Desktop/Adaptive_Nanopore/combined_v5_best_hit.tab'
ref_taxon_txt                   = '/Users/songweizhi/Desktop/Adaptive_Nanopore/nt_seq_taxon.txt'
uniq_genbank_id_txt             = '/Users/songweizhi/Desktop/Adaptive_Nanopore/uniq_genbank_id.txt'
assembly_with_matched_ref_txt   = '/Users/songweizhi/Desktop/Adaptive_Nanopore/assembly_with_matched_ref.txt'
classification_txt              = '/Users/songweizhi/Desktop/Adaptive_Nanopore/classifications.txt'


# best_hit_txt                    = '/Users/songweizhi/Desktop/Adaptive_Nanopore/combined_v5_best_hit.tab'
# ref_taxon_txt                   = '/Users/songweizhi/Desktop/Adaptive_Nanopore/nt_seq_taxon.txt'
# uniq_genbank_id_txt             = '/Users/songweizhi/Desktop/Adaptive_Nanopore/uniq_genbank_id.txt'
# assembly_with_matched_ref_txt   = '/Users/songweizhi/Desktop/Adaptive_Nanopore/assembly_with_matched_ref.txt'
# assembly_fa                     = '/Users/songweizhi/Desktop/eukaryotic_contigs.fa'
# classification_txt              = '/Users/songweizhi/Desktop/eukaryotic_contigs_classifications.txt'


# get seq_len_dict
seq_len_dict = {}
for each_seq in SeqIO.parse(assembly_fa, 'fasta'):
    seq_len_dict[each_seq.id] = len(each_seq.seq)

# read in genbank id to taxon id
genbank_id_to_taxon_id_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/55/uniq_genbank_id_taxon.txt'
genbank_id_to_taxon_id_dict = {}
for each_line in open(genbank_id_to_taxon_id_txt):
    each_line_split = each_line.strip().split('\t')
    genbank_id_to_taxon_id_dict[each_line_split[0]] = each_line_split[1]

# read in taxon id to lineage into dict
taxon_id_to_lineage_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/55/tax_lineage_by_name.txt'
taxon_id_to_lineage_dict = {}
for each_line in open(taxon_id_to_lineage_txt):
    each_line_split = each_line.strip().split('\t')
    taxon_id_to_lineage_dict[each_line_split[0]] = each_line_split[1]

classification_txt_handle = open(classification_txt, 'w')
assembly_with_matched_ref_txt_handle = open(assembly_with_matched_ref_txt, 'w')
m = 0
n = 0
matched_ref_list = []
ref_to_ctg_dict = dict()
for each_hit in open(best_hit_txt):
    each_hit_split = each_hit.strip().split('\t')
    aln_len = int(each_hit_split[3])
    ctg_id = each_hit_split[0]
    ref_id = each_hit_split[1]
    if aln_len >= 1000:

        ref_no_version = ref_id
        if '.' in ref_id:
            ref_no_version = ref_id.split('.')[0]

        ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
        ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')

        assembly_with_matched_ref_txt_handle.write('%s\n' % ctg_id)

        if 'Eukaryota' in ref_lineage:
            print('%s\t%s\t%s\t%s\t%s\t%s' % (ctg_id,ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage ))

        classification_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ctg_id, ref_id, each_hit_split[2], aln_len, seq_len_dict[ctg_id], ref_lineage))

        matched_ref_list.append(ref_id)
        if ref_id not in ref_to_ctg_dict:
            ref_to_ctg_dict[ref_id] = [ctg_id]
        else:
            ref_to_ctg_dict[ref_id].append(ctg_id)
        n += 1
    m += 1
assembly_with_matched_ref_txt_handle.close()
classification_txt_handle.close()


matched_ref_list_uniq = {i for i in matched_ref_list}


# write out Genbank id list
uniq_genbank_id_txt_handle = open(uniq_genbank_id_txt, 'w')
for each_id in matched_ref_list_uniq:
    uniq_genbank_id_txt_handle.write('%s\n' % each_id)
uniq_genbank_id_txt_handle.close()



total_len = 0
total_len_Eukaryota = 0
total_len_Bacteria_Archaea = 0
total_len_Virus = 0
total_len_unclassified = 0
for each_ref in ref_to_ctg_dict:

    matched_ctg_list = ref_to_ctg_dict[each_ref]

    matched_ctg_total_len = 0
    for each_matched_ctg in matched_ctg_list:
        ctg_len = seq_len_dict[each_matched_ctg]
        matched_ctg_total_len += ctg_len

    total_len += matched_ctg_total_len

    ref_no_version = each_ref
    if '.' in each_ref:
        ref_no_version = each_ref.split('.')[0]

    ref_tax_id = genbank_id_to_taxon_id_dict.get(ref_no_version, 'NA')
    ref_lineage = taxon_id_to_lineage_dict.get(ref_tax_id, 'NA')


    #if 'Eukaryota' in ref_lineage:
    #    print('%s\t%s\t%s\t%s' % (each_ref, matched_ctg_total_len, ','.join(matched_ctg_list), ref_lineage))


    if 'Bacteria' in ref_lineage:
        total_len_Bacteria_Archaea += matched_ctg_total_len
    elif 'Archaea' in ref_lineage:
        total_len_Bacteria_Archaea += matched_ctg_total_len
    elif 'Eukaryota' in ref_lineage:
        total_len_Eukaryota += matched_ctg_total_len
        #print('%s\t%s\t%s\t%s' % (ref_lineage, each_ref, matched_ctg_total_len, ','.join(matched_ctg_list)))

    elif 'unclassified' in ref_lineage:
        total_len_unclassified += matched_ctg_total_len
    elif 'Virus' in ref_lineage:
        total_len_Virus += matched_ctg_total_len
    else:
        print('%s\t%s\t%s\t%s' % (ref_lineage, each_ref, matched_ctg_total_len, ','.join(matched_ctg_list)))


print('total_len :\t%sbp\t%s Mbp'                  % (total_len, total_len/(1024*1024)))
print('total_len_Bacteria_Archaea :\t%sMbp\t%s'    % (total_len_Bacteria_Archaea/(1024*1024),   total_len_Bacteria_Archaea/total_len))
print('total_len_Eukaryota:\t%sMbp\t%s'            % (total_len_Eukaryota/(1024*1024),          total_len_Eukaryota/total_len))
print('total_len_unclassified:\t%sMbp\t%s'         % (total_len_unclassified/(1024*1024),       total_len_unclassified/total_len))
print('total_len_Virus:\t%sMbp\t%s'                % (total_len_Virus/(1024*1024),              total_len_Virus/total_len))

