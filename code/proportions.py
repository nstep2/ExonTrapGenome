




pdf_stats = PdfPages(exp_output_path.proportion_bases_different_annotations+'proportion_bases_different_annotations_%d_reads.pdf' % (exon_count_build))





exon_id_cat_dict   = region_exon_id_sets_dict

genome_region_list = region_exon_id_sets_list







hg_chromosome_bases = 2 * 3088286401
ss_hg_chromosome_bases = 3088286401

exon_regions_primary_exon_id_list = list()
count_exon_regions = 0
for chrom in tree_3ss:
    for strand in tree_3ss[chrom]:
        for interval in tree_3ss[chrom][strand]:
            count_exon_regions += 1
            
            max_exon_id = el.get_max_count_exon_id_in_list(interval[2],aggregate_exon_dict)
            exon_regions_primary_exon_id_list.append(max_exon_id)
            


exon_regions_primary_exon_id_list=list(set(exon_regions_primary_exon_id_list))



exon_id_list = pc_middle_exon_id_list
exon_id_list = el.exon_id_intersection(exon_id_list, aggregate_exon_dict.keys())
counts_list = (el.get_exon_dict_counts(exon_id_list, aggregate_exon_dict))

print('Internal exons from protein-coding genes, in contrast, had an average of %s reads' % insert_commas(int(np.mean(counts_list))))


very_high_read_counts_exon_id_list = el.threshold_exon_ids(aggregate_exon_dict.keys(),20000,aggregate_exon_dict)

exons_above_20000_count = len(very_high_read_counts_exon_id_list)
mRNA_20000_count = len(set(very_high_read_counts_exon_id_list).intersection(pc_middle_exon_id_list))











region_base_count = dict()
region_dataset_dict = dict()

seleccted_regions_dict = dict()

seleccted_regions_dict['mRNA'] = pc_middle_exon_id_list
seleccted_regions_dict['lncRNA'] = lncRNA_middle_exon_ids

seleccted_regions_dict['intergenic'] = no_overlap_set
seleccted_regions_dict['antisense'] = antisense_transcript_set

seleccted_regions_dict['intronic_mRNA'] = intron_interior_pc_set.difference(intron_interior_lncRNA_set)
seleccted_regions_dict['intronic_lncRNA'] = intron_interior_lncRNA_set.difference(intron_interior_pc_set)





def get_number_chromosome_bases_exon_id_list(exon_id_list):
    bases_dict = dict()
    exon_id_list=list(set(exon_id_list))
    
    chrom_list = list()
    for chrom in range(1,23):
        chrom_list.append('chr%d' % chrom)
    chrom_list.append('chrX')
    chrom_list.append('chrY')
    chrom_list.append('chrM')
    
    sum_chromosom_bases = 0
    for chrom in chrom_list:
        chrom_base_list = list()
        for strand in ['+','-']:
            for exon_id in exon_id_list:
                ex = el.exon_id_values(exon_id)
                if ex.chrom != chrom:
                    continue
                if ex.strand != strand:
                    continue
                chrom_base_list += range(ex.start,ex.end+1)
            chrom_base_list=list(set(chrom_base_list))
            sum_chromosom_bases += len(chrom_base_list)
        
    
    return sum_chromosom_bases


for key in seleccted_regions_dict:
    exon_id_list = seleccted_regions_dict[key]
    exon_id_list = el.exon_id_intersection(exon_id_list, aggregate_exon_dict.keys())   
    total = get_number_chromosome_bases_exon_id_list(exon_id_list)
    
    region_dataset_dict[key] = exon_id_list
    region_base_count[key] = total











print('Intergenic regions (i.e. excluding any kind of mRNA or lncRNA, sense or antisense) contained an exon region every 15,403 bases on average (xx%% of the sequence): %.1f%%' % (region_base_count['intergenic']/hg_chromosome_bases*100))




non_exon_bases = hg_chromosome_bases
val_list = list()
for key in region_base_count:
    total = region_base_count[key]
    val_list.append([key, total])
    
    non_exon_bases-=total
    


val_list.append(['non exon',non_exon_bases])




non_exon_bases = hg_chromosome_bases
val_list = list()
for key in region_base_count:
    total = region_base_count[key]
    val_list.append([key, total])
    
    non_exon_bases-=total


x,y = zip(*val_list)

explode=0*np.ones(len(x))
explode[-1] = 0

fig, ax1 = plt.subplots()
patches, texts, autotexts = ax1.pie(y, explode=explode, labels=x, autopct='%.1f%%', startangle=0)

plt.legend(x)
plt.axis('equal')    
plt.tight_layout()
pdf_stats.savefig(fig)


import pandas as pd

#val_list
val_list_text = ["{:,}".format(a) for a in y]
percent_text = [a.get_text() for a in autotexts]
df = pd.DataFrame(zip(x,val_list_text, percent_text), columns=['region','count','percent'])













mRNA_exon_genome_len = 0
for chrom in gencode_pc_exon_IT:
    for strand in gencode_pc_exon_IT[chrom]:
        
        chrom_mRNA_exon_bases = list()
        for interval in gencode_pc_exon_IT[chrom][strand]:
            chrom_mRNA_exon_bases += (range(interval[0],interval[1]+1))
        
        chrom_mRNA_exon_bases = list(set(chrom_mRNA_exon_bases))
        mRNA_exon_genome_len += len(chrom_mRNA_exon_bases)

lncRNA_exon_genome_len = 0
for chrom in gencode_lncRNA_exon_IT:        
    for strand in gencode_lncRNA_exon_IT[chrom]:
        
        chrom_lncRNA_exon_bases = list()
        for interval in gencode_lncRNA_exon_IT[chrom][strand]:
            chrom_lncRNA_exon_bases += (range(interval[0],interval[1]+1))
        
        chrom_lncRNA_exon_bases = list(set(chrom_lncRNA_exon_bases))
        lncRNA_exon_genome_len += len(chrom_lncRNA_exon_bases)
        
        
        
        




unannotated_antisense_len = 0
no_overlap_genome_len = ss_hg_chromosome_bases
transcript_genome_len = 0 #includes antisense
chrom_transcript_genome_len_dict = dict()
mRNA_genome_len = 0
lncRNA_genome_len = 0
for chrom in gencode_pc_transcript_IT:
    chrom_transcript_bases = list()
    
    chrom_transcript_bases_plus = list()
    chrom_transcript_bases_minus = list()
    
    for strand in gencode_pc_transcript_IT[chrom]:
        
        chrom_mRNA_bases = list()
        chrom_lncRNA_bases = list()
        
        for interval in gencode_lncRNA_transcript_IT[chrom][strand]:
            
            chrom_transcript_bases += list(range(interval[0],interval[1]+1))
            chrom_lncRNA_bases   += list(range(interval[0],interval[1]+1))
            
            if strand == '+':
                chrom_transcript_bases_plus += list(range(interval[0],interval[1]+1))
            if strand == '-':
                chrom_transcript_bases_minus += list(range(interval[0],interval[1]+1))
        
        
        chrom_transcript_bases = list(set(chrom_transcript_bases))
        chrom_transcript_bases_minus = list(set(chrom_transcript_bases_minus))
        chrom_transcript_bases_plus = list(set(chrom_transcript_bases_plus))
        
        
        
        for interval in gencode_pc_transcript_IT[chrom][strand]:
            chrom_transcript_bases += (list(range(interval[0],interval[1]+1)))
            chrom_mRNA_bases += (list(range(interval[0],interval[1]+1)))
            if strand == '+':
                chrom_transcript_bases_plus += list(range(interval[0],interval[1]+1))
            if strand == '-':
                chrom_transcript_bases_minus += list(range(interval[0],interval[1]+1))
        
        
        chrom_transcript_bases = list(set(chrom_transcript_bases))
        chrom_transcript_bases_minus = list(set(chrom_transcript_bases_minus))
        chrom_transcript_bases_plus = list(set(chrom_transcript_bases_plus))
        
        chrom_mRNA_bases = list(set(chrom_mRNA_bases))
        chrom_lncRNA_bases = list(set(chrom_lncRNA_bases))
    
        mRNA_genome_len   +=  len(chrom_mRNA_bases)
        lncRNA_genome_len +=  len(chrom_lncRNA_bases)

        
    unannotated_antisense_len += len(set(chrom_transcript_bases_minus).difference(chrom_transcript_bases_plus))
    
    unannotated_antisense_len += len(set(chrom_transcript_bases_plus).difference(chrom_transcript_bases_minus))
    
    no_overlap_genome_len -= len(chrom_transcript_bases)
    transcript_genome_len += len(chrom_transcript_bases)
    chrom_transcript_genome_len_dict[chrom] = len(chrom_transcript_bases)
    
    










genome_regions_total_lengths = dict()

genome_regions_total_lengths['intergenic'] = no_overlap_genome_len
genome_regions_total_lengths['transcript'] = transcript_genome_len
genome_regions_total_lengths['antisense'] = unannotated_antisense_len
genome_regions_total_lengths['mRNA'] = mRNA_genome_len
genome_regions_total_lengths['lncRNA'] = lncRNA_genome_len

genome_regions_total_lengths['mRNA_exon'] = mRNA_exon_genome_len
genome_regions_total_lengths['lncRNA_exon'] = lncRNA_exon_genome_len



print('no_overlap_genome_len', insert_commas(no_overlap_genome_len))
print('len(no_overlap_set)', insert_commas(len(no_overlap_set)))
print('exon per N intergenic bases:', insert_commas(int(no_overlap_genome_len/len( region_dataset_dict['intergenic'] ))))
print('%.1f%% of the intergenic sequence that is exonic'% (region_base_count['intergenic']/no_overlap_genome_len*100))
print('Since there is a large amount of sequence, however, the absolute number of intergenic exonic regions is high (%s)' % insert_commas(len( region_dataset_dict['intergenic'] )))

#v = sense_transcript_set.intersection(exon_regions_primary_exon_id_list).difference(pc_middle_exon_id_list).difference(intron_interior_set)

print('exonic regions were also found within annotated introns, in the forward strand (Figure 3A (%s)' % insert_commas(len( set(region_dataset_dict['intronic_lncRNA']).union(region_dataset_dict['intronic_mRNA']) )))

percent = len(region_dataset_dict['intronic_lncRNA'])/len( set(region_dataset_dict['intronic_lncRNA']).union(region_dataset_dict['intronic_mRNA']) )
lncRNA_num = len(region_dataset_dict['intronic_lncRNA'])
pc_num = len(region_dataset_dict['intronic_mRNA'])
print(' %.1f%% within lncRNA introns(mRNA: %s, lncRNA: %s ).' % (percent, insert_commas(pc_num), insert_commas(lncRNA_num)))





def get_tid_list_from_IT(gencode_transcript_IT):

    tid_list = list()
    for chrom in gencode_transcript_IT:
        for strand in gencode_transcript_IT[chrom]:
            for interval in gencode_transcript_IT[chrom][strand]:
                if len(interval[2])<0: #only intron containing
                    continue
                tid = '{:}:{:}-{:}:{:}'.format(chrom, interval[0], interval[1], strand)
                tid_list.append(tid)
    
    return tid_list



pc_tid_list = get_tid_list_from_IT(gencode_pc_transcript_IT)

pc_transcript_bases_len = get_number_chromosome_bases_exon_id_list(pc_tid_list)



lncRNA_tid_list = get_tid_list_from_IT(gencode_lncRNA_transcript_IT)

lncRNA_transcript_bases_len = get_number_chromosome_bases_exon_id_list(lncRNA_tid_list)











def get_median_count_dict(exon_id_list):
    exon_id_list = el.exon_id_intersection(exon_id_list,exon_regions_primary_exon_id_list)
    exon_id_list = el.exon_id_intersection(exon_id_list, aggregate_exon_dict.keys())
    median_count = np.median(el.get_exon_dict_counts(exon_id_list, aggregate_exon_dict))
    return median_count
    







annotated_middle_bases_dict = dict()

annotated_middle_bases_dict['mRNA'] = get_number_chromosome_bases_exon_id_list(exon_id_cat_dict['mRNA'])

annotated_middle_bases_dict['lncRNA'] = get_number_chromosome_bases_exon_id_list(exon_id_cat_dict['lncRNA'])





region_base_count

genome_regions_total_lengths


genome_frac_dict = dict()
genome_frac_dict_colors = dict()

genome_frac_dict['mRNA_exon'] = region_base_count['mRNA']/genome_regions_total_lengths['mRNA_exon']

genome_frac_dict['lncRNA_exon'] = region_base_count['lncRNA']/genome_regions_total_lengths['lncRNA_exon']

genome_frac_dict['intergenic_exon'] = region_base_count['intergenic']/genome_regions_total_lengths['intergenic']

genome_frac_dict['antisense_exon'] = region_base_count['antisense']/genome_regions_total_lengths['antisense']

genome_frac_dict['lncRNA_intronic'] = region_base_count['intronic_lncRNA']/(genome_regions_total_lengths['lncRNA']-genome_regions_total_lengths['lncRNA_exon'])
 
genome_frac_dict['mRNA_intronic'] = region_base_count['intronic_mRNA']/(genome_regions_total_lengths['mRNA']-genome_regions_total_lengths['mRNA_exon'])






#output supplemental files
with open(exp_output_path.out_supplemental+'3B.txt','w') as f:
    #outstring = 'chrom\tstart\tend\tstrand\texon_id\tgenomic_category\n'
    #f.write(outstring)
    outstring = 'trapped_region\tgenomic_region\ttrapped_base_count\tgenome_base_count\n'
    f.write(outstring)
    
    
    outstring = 'mRNA\tinternal_mRNA_exons\t{:}\t{:}\n'.format(
    region_base_count['mRNA'],
    genome_regions_total_lengths['mRNA_exon'])
    
    #f.write(outstring)
    
    
    outstring = 'lncRNA\tInternal_lncRNA_exons\t{:}\t{:}\n'.format(
    region_base_count['lncRNA'],
    genome_regions_total_lengths['lncRNA_exon'])
    
    f.write(outstring)
    
    
    outstring = 'Intergenic\tIntergenic_bases\t{:}\t{:}\n'.format(
    region_base_count['intergenic'],
    genome_regions_total_lengths['intergenic'])
    
    f.write(outstring)
    
    outstring = 'Antisense_transcripts\tAntisense_bases\t{:}\t{:}\n'.format(
    region_base_count['antisense'],
    genome_regions_total_lengths['antisense'])
    
    f.write(outstring)
    
    
    outstring = 'Intronic_lncRNA\tIntronic_lncRNA_bases\t{:}\t{:}\n'.format(
    region_base_count['intronic_lncRNA'],
    (genome_regions_total_lengths['lncRNA']-genome_regions_total_lengths['lncRNA_exon']))
    
    f.write(outstring)
    
    
    outstring = 'Intronic_mRNA\tIntronic_mRNA_bases\t{:}\t{:}\n'.format(
    region_base_count['intronic_mRNA'],
    (genome_regions_total_lengths['mRNA']-genome_regions_total_lengths['mRNA_exon']))
    
    f.write(outstring)















an_middle_exon_mRNA = annotated_middle_bases_dict['mRNA']/genome_regions_total_lengths['mRNA_exon']

an_middle_exon_lncRNA = annotated_middle_bases_dict['lncRNA']/genome_regions_total_lengths['lncRNA_exon']




(region_base_count['intronic_mRNA']+ region_base_count['intronic_lncRNA'])/hg_chromosome_bases*100




genome_frac_dict_colors['lncRNA_exon'] = ['lncRNA', 1.0]
genome_frac_dict_colors['antisense_exon'] = ['Antisense', 1.0]
genome_frac_dict_colors['mRNA_exon'] = ['mRNA', 1.0]
genome_frac_dict_colors['intergenic_exon'] = ['Intergenic', 1.0]
genome_frac_dict_colors['mRNA_intronic'] = ['Intronic', 1.0]
genome_frac_dict_colors['lncRNA_intronic'] = ['Intronic', 0.6]
#genome_frac_dict_colors['mRNA_intronic'] = ['Intronic', .6]

for key in genome_frac_dict:
    print(key, genome_frac_dict[key]*100)

key_list, fracs = zip(*[(k,genome_frac_dict[k]*100) for k in genome_frac_dict])






genome_frac_dict_colors['premRNA'] = ['Alternative', 1.0]
genome_frac_dict_colors['prelncRNA'] = ['Alternative', 0.6]


def get_pre_transcript_stats(tid_list,aggregate_exon_IT):
    pre_mRNA_exon_bases_count = 0
    pre_mRNA_bases_count = 0
    for chrom in aggregate_exon_IT.keys():
        
        for strand in aggregate_exon_IT[chrom]:
            tx_exon_found_bases = list()
            tx_found_bases = list()
            
            for tid in tid_list:
                tx=el.exon_id_values(tid)
                if tx.chrom != chrom:
                    continue
                if tx.strand != strand:
                    continue
                if tx.length > 10000000000000:
                    continue
                tx_found_bases += list(range(tx.start,tx.end))

                exon_interval_list = aggregate_exon_IT[tx.chrom][tx.strand][tx.start:tx.end]
                
                for interval in exon_interval_list:
                    max_exon_id = el.get_max_count_exon_id_in_list(interval[2],aggregate_exon_dict)
                    
                    if aggregate_exon_dict[max_exon_id]['count'] < 100:
                        continue
                    ex = el.exon_id_values(max_exon_id)
                    tx_exon_found_bases += list(range(ex.start,ex.end))
                
                
                
            tx_exon_found_bases = list(set(tx_exon_found_bases))
            tx_found_bases = list(set(tx_found_bases))
            #print(chrom, len(tx_exon_found_bases)/len(tx_found_bases))
                    
            pre_mRNA_exon_bases_count += len(tx_exon_found_bases)
            pre_mRNA_bases_count += len(tx_found_bases)

    return pre_mRNA_exon_bases_count,pre_mRNA_bases_count



pre_mRNA_exon_bases_count,pre_mRNA_bases_count = get_pre_transcript_stats(pc_tid_list,aggregate_exon_IT)

genome_frac_dict['premRNA'] = pre_mRNA_exon_bases_count/pre_mRNA_bases_count

pre_lncRNA_exon_bases_count,pre_lncRNA_bases_count = get_pre_transcript_stats(lncRNA_tid_list,aggregate_exon_IT)

genome_frac_dict['prelncRNA'] = pre_lncRNA_exon_bases_count/pre_lncRNA_bases_count









pdf_stats.close()






    
    
pickle_path = exp_output_path.pickle_merged + "proportion_bases_different_annotations__%d.pickle" % (exon_count_build)

with open(pickle_path, "wb") as output_file:
    pickle.dump(genome_frac_dict, output_file)
    
    

'''

pickle_path = exp_output_path.pickle_merged + "proportion_bases_different_annotations__%d.pickle" % (exon_count_build)

with open(pickle_path, "rb") as input_file:
    genome_frac_dict = pickle.load(input_file)

'''



'''
get_number_chromosome_bases_exon_id_list(mRNA_ids)
#+
get_number_chromosome_bases_exon_id_list(intron_exon_ids)


get_number_chromosome_bases_exon_id_list(mRNA_transcript_ids)
'''

#=list(tid_dict.keys())

