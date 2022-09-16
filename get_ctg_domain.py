from Bio import SeqIO


ctg_taxon_file  = '/Users/songweizhi/Desktop/Adaptive_Nanopore/classification_aln100bp.txt'
bin_file        = '/Users/songweizhi/Desktop/eukaryotic_contigs.fa'
out_txt         = '/Users/songweizhi/Desktop/eukaryotic_contigs.fa.txt'


# read in contig id
bin_23_ctg_set = set()
for each_seq in SeqIO.parse(bin_file, 'fasta'):
    bin_23_ctg_set.add(each_seq.id)


# read in contig classification
ctg_to_tax_dict = dict()
for each_line in open(ctg_taxon_file):
    each_line_split = each_line.strip().split('\t')
    contig_id = each_line_split[0]
    ref_taxon = each_line_split[6]
    ctg_to_tax_dict[contig_id] = ref_taxon


out_txt_handle = open(out_txt, 'w')
num_total = 0
num_na = 0
num_Eukaryota = 0
num_Bacteria = 0
num_Archaea = 0
num_others = 0
for each_ctg in bin_23_ctg_set:
    num_total += 1
    ctg_tax = ctg_to_tax_dict.get(each_ctg, 'NA')
    if ctg_tax == 'NA':
        num_na += 1
    elif 'Eukaryota' in ctg_tax:
        num_Eukaryota += 1
    elif 'Bacteria' in ctg_tax:
        num_Bacteria += 1
    elif 'Archaea' in ctg_tax:
        num_Archaea += 1
    else:
        num_others += 1
        pass
        #print(ctg_tax)
    if ctg_tax != 'NA':
        print(ctg_tax)
        out_txt_handle.write('%s\n' % ctg_tax)

out_txt_handle.close()


print()
print('Total\t%s'     % num_total)
print('No hit\t%s'    % num_na)
print('Eukaryota\t%s' % num_Eukaryota)
print('Bacteria\t%s'  % num_Bacteria)
print('Archaea\t%s'   % num_Archaea)
print('Others\t%s'    % num_others)

