



print('no_overlap_set: %s' % insert_commas(len(no_overlap_set)))
print('sense_transcript_set: %s' % insert_commas(len(sense_transcript_set)))
print('antisense_transcript_set: %s' % insert_commas(len(antisense_transcript_set)))
print('intron_interior_set: %s' % insert_commas(len(intron_interior_set)))
print('intron_interior_pc_set: %s' % insert_commas(len(intron_interior_pc_set)))
print('intron_interior_lncRNA_set: %s' % insert_commas(len(intron_interior_lncRNA_set)))
print('exon_overlap_set: %s' % insert_commas(len(exon_overlap_set)))


print('primary in 3ss cluster & above 100 read count')
print('all_3ss_primary: %s' % insert_commas(len(primary_3ss_exon_id_set)))
print('no_overlap_set: %s' % insert_commas(len(no_overlap_set.intersection(primary_3ss_exon_id_set))))
print('sense_transcript_set: %s' % insert_commas(len(sense_transcript_set.intersection(primary_3ss_exon_id_set))))
print('antisense_transcript_set: %s' % insert_commas(len(antisense_transcript_set.intersection(primary_3ss_exon_id_set))))
print('intron_interior_set: %s' % insert_commas(len(intron_interior_set.intersection(primary_3ss_exon_id_set))))
print('intron_interior_pc_set: %s' % insert_commas(len(intron_interior_pc_set.intersection(primary_3ss_exon_id_set))))
print('intron_interior_lncRNA_set: %s' % insert_commas(len(intron_interior_lncRNA_set.intersection(primary_3ss_exon_id_set))))
print('exon_overlap_set: %s' % insert_commas(len(exon_overlap_set.intersection(primary_3ss_exon_id_set))))








#for threshold in [1,10,100,1000,10000]:
for threshold in [100]:
    print('threshold: %s' % (insert_commas(threshold)))
    
    print('no_overlap_set: %s' % insert_commas(above_threshold_exon_id(threshold,no_overlap_set, aggregate_exon_dict)))
    print('sense_transcript_set: %s' % insert_commas(above_threshold_exon_id(threshold,sense_transcript_set, aggregate_exon_dict)))
    print('antisense_transcript_set: %s' % insert_commas(above_threshold_exon_id(threshold,antisense_transcript_set, aggregate_exon_dict)))
    print('intron_interior_set: %s' % insert_commas(above_threshold_exon_id(threshold,intron_interior_set, aggregate_exon_dict)))
    print('intron_interior_pc_set: %s' % insert_commas(above_threshold_exon_id(threshold,intron_interior_pc_set, aggregate_exon_dict)))
    print('intron_interior_lncRNA_set: %s' % insert_commas(above_threshold_exon_id(threshold,intron_interior_lncRNA_set, aggregate_exon_dict)))
    print('exon_overlap_set: %s' % insert_commas(above_threshold_exon_id(threshold,exon_overlap_set, aggregate_exon_dict)))
    print('exon_overlap_lncRNA_coding_set: %s' % insert_commas(above_threshold_exon_id( threshold,exon_overlap_set.intersection(gencode_lncRNA_overlapping_exon_ids_set), aggregate_exon_dict) ))
    print('exon_overlap_protein_coding_set: %s' % insert_commas(above_threshold_exon_id( threshold,exon_overlap_set.intersection(gencode_pc_overlapping_exon_ids_set), aggregate_exon_dict) ))
    print('exon_overlap_set-Non_lncRNA-Non_protein_coding: %s' % insert_commas(above_threshold_exon_id( threshold,exon_overlap_set.difference(gencode_pc_overlapping_exon_ids_set).difference(gencode_lncRNA_overlapping_exon_ids_set), aggregate_exon_dict) ))
    print('Number of exon_ids that overlap both a protein coding exon and a lncRNA exon: %s' % insert_commas(len(gencode_pc_overlapping_exon_ids_set.intersection(gencode_lncRNA_overlapping_exon_ids_set))))
    
    print('\n')







exon_counts_list_dict['mRNA']
exon_counts_list_dict['lncRNA']
exon_counts_list_dict['Intergenic']
exon_counts_list_dict['Antisense']
exon_counts_list_dict['Intronic']













































