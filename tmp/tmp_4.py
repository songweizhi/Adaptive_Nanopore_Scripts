
classification_aln100bp_txt = '/Users/songweizhi/Desktop/Adaptive_Nanopore/classification_aln100bp.txt'
Eukaryota_ctg_id_txt        = '/Users/songweizhi/Desktop/Adaptive_Nanopore/Eukaryota_ctg_id.txt'
nonEukaryota_ctg_id_txt     = '/Users/songweizhi/Desktop/Adaptive_Nanopore/nonEukaryota_ctg_id.txt'


Eukaryota_ctg_id_txt_handle = open(Eukaryota_ctg_id_txt, 'w')
nonEukaryota_ctg_id_txt_handle = open(nonEukaryota_ctg_id_txt, 'w')
for each_line in open(classification_aln100bp_txt):
    each_line_split = each_line.strip().split('\t')
    print(each_line_split)
    if 'Eukaryota' in each_line:
        Eukaryota_ctg_id_txt_handle.write('%s\n' % each_line_split[0])
    else:
        nonEukaryota_ctg_id_txt_handle.write('%s\n' % each_line_split[0])


Eukaryota_ctg_id_txt_handle.close()
nonEukaryota_ctg_id_txt_handle.close()

