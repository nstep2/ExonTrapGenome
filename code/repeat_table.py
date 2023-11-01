


if True == False:
    
    #all exon
    ea = [[eid, aggregate_exon_dict[eid]['count'], repeat_associated_exon_ids_dict[eid]] for eid in repeat_exon_keys ]




dfam_id_to_family_dict = dict()
for ii, entry in enumerate(ea):
    for entry_2 in entry[2]:
        dfam_id_to_family_dict[entry_2['family_acc']] = entry_2['family_name']
    






reapeat_table_header = list()
reapeat_table_header.append('Dfam_ID')
reapeat_table_header.append('Repeat_genome_bases')
reapeat_table_header.append('Exon_genome_bases')
reapeat_table_header.append('exon_bases_vs_repeat_count_ratio')
reapeat_table_header.append('exon_bases_inntersect_repeat_ratio')
reapeat_table_header.append('exon_count')
reapeat_table_header.append('repeat_count')
reapeat_table_header.append('exon_count_repeat_count_ratio')
reapeat_table_header.append('average_exon_expression')
reapeat_table_header.append('ID')
reapeat_table_header.append('exon_repeat_enrichment_relative_genome_ratio')
reapeat_table_header.append('pval_sf')
reapeat_table_header.append('pval_cdf')



repeat_family_dict = dict()
for key in dfam_id_to_family_dict:
    repeat_family_dict[dfam_id_to_family_dict[key]] = list()





repeat_val_table = list()
for key in repeat_stat_rows_dict:
    row_list = list()
    
    row_list.append(key)
    
    
    if key in dfam_id_to_family_dict:
        family_id = dfam_id_to_family_dict[key]
        row_list.append(family_id)
    else:
        continue
    
    repeat = repeat_stat_rows_dict[key]
                
    repeat_genome_bases = repeat['rep_genome_bases']
    row_list.append(repeat_genome_bases)
    
    exon_bases = repeat['exon_genome_bases']
    row_list.append(exon_bases)
    
    ratio_exon_base_to_repeat_base = exon_bases/repeat_genome_bases
    row_list.append(ratio_exon_base_to_repeat_base )
    
    exon_bases_inntersect_repeat = repeat['exon_rep_intersect_bases']
    ratio_exon_base_to_repeat_base = exon_bases_inntersect_repeat/repeat_genome_bases
    row_list.append(ratio_exon_base_to_repeat_base )
    
    exon_genome_count = repeat['exon_genome_count']
    row_list.append(exon_genome_count)
    
    rep_genome_counts = repeat['rep_genome_counts']
    row_list.append(rep_genome_counts)
    
    exons_per_repeat = exon_genome_count/rep_genome_counts
    row_list.append(exons_per_repeat)    
    
    avg_exon_exp = repeat['avg_exon_exp']
    row_list.append(avg_exon_exp)
    
    
    
    ratio   = repeat_stat_ranking_dict[key][-2]
    row_list.append(ratio)
    
    sf      = repeat_stat_ranking_dict[key][-3]
    row_list.append(sf)
    
    cdf     = repeat_stat_ranking_dict[key][-1]
    row_list.append(cdf)
    
    
    
    
    repeat_val_table.append(row_list)




repeat_val_table = sorted(repeat_val_table, key=lambda x:(x[9]), reverse=True)


new_repeat_val_table = list()
for ii, entry in enumerate(repeat_val_table):
    new_repeat_val_table.append(entry)


new_repeat_val_table_strings = list()
for ii, val in enumerate(new_repeat_val_table):
    new_repeat_val_table_strings.append('\t'.join(map(str, val)) )


with open(exp_output_path.Dfam_repeats + 'repeat_genome_stats_table.txt','wt') as f:
    f.write(' \t '.join(reapeat_table_header))
    
    for ii, val in enumerate(new_repeat_val_table_strings):
        print(val)
        f.write(val + '\n')
    






'''

    repeat_stat_ranking_dict[key].append(test_hyper_3) #sf
    repeat_stat_ranking_dict[key].append(test_hyper_3_ratio)  #enrichment ratio
    repeat_stat_ranking_dict[key].append(test_hyper_4) #cdf

'''














