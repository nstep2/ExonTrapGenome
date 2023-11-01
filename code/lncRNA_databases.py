import exon_id_library
import exon_id_library.exon_id_lib as el
import exon_id_library.bed as bl




outdir = exp_output_path.Parse_ENSEMBLE_GENCODE_exons_load_figures_2
pdf_plots = PdfPages(outdir+'Parse_ENSEMBLE_GENCODE_exons_load_figures_2.pdf')






class bed_recovry_class():
    def __init__(self, bed_file, aggregate_exon_dict, aggregate_exon_IT, genome_fasta):
        all_recovered_exon_id_list = set(aggregate_exon_dict.keys())
        
        self.bed_exon_ids                = bl.simple_bed(bed_file, aggregate_exon_IT)
        self.recovered_exon_ids          =  el.exon_id_intersection( self.bed_exon_ids.middle_exon_id_list, all_recovered_exon_id_list)
        self.fuzzy_middle_exon_ids       = (el.exon_id_intersection( self.bed_exon_ids.fuzzy_overlap_middle_exon_ids, all_recovered_exon_id_list))
        self.exact_fuzzy_middle_exon_ids = list(set(self.recovered_exon_ids).union(self.fuzzy_middle_exon_ids))
        
        self.stats = dict()
        self.stats['middle_exon_count']          = len(self.bed_exon_ids.middle_exon_id_list)
        self.stats['recovered_middle_exon_count']= len(self.recovered_exon_ids)
        self.stats['fuzzy_middle_exon_count']    = len(self.bed_exon_ids.fuzzy_overlap_middle_exon_ids)
        
    def print_stats(self):
        print_list              = ['middle_exon_count' , 'recovered_middle_exon_count', 'fuzzy_middle_exon_count', 'recovered_fuzzy_middle_exon_count']
        middle_exon_ratio       = self.stats['recovered_middle_exon_count']/self.stats['middle_exon_count']        
        fuzzy_middle_exon_ratio = len(set(self.recovered_exon_ids).union(self.fuzzy_middle_exon_ids))/self.stats['middle_exon_count']




def print_gc_ratios(query_exon_id_list, aggregate_exon_dict, aggregate_exon_IT):
    gc_results = el.get_exon_id_list_binned_GC(el.size_exon_id_list(query_exon_id_list,50,500), [20,30,40,50,60,70,80], genome_fasta)
    
    query_dict              = dict()
    query_assignment_dict   = dict()
    
    for ii, index in enumerate(gc_results):
        query_assignment_dict[index] = el.exon_assignment_class(gc_results[index], aggregate_exon_IT, aggregate_exon_dict)
    
    ratio_dict = dict()
    for key in sorted(query_assignment_dict.keys(), reverse=True):
        #print("\n*** Exons with GC content %s ***" % (key))
        query_assignment_dict[key].print_ratio()
        ratio_dict[key] = query_assignment_dict[key].recovery_ratio

    return ratio_dict









lncRNA_database_bed_dict   = dict()
all_recovered_exon_id_list = set(aggregate_exon_dict.keys())

bed_dict = dict()
bed_dict['NONCODE'] = exp_output_path.lncNRA_database_bed + 'NONCODEv6_hg38.lncAndGene.bed'
bed_dict['RNAcentral'] = exp_output_path.lncNRA_database_bed + 'homo_sapiens.GRCh38.bed'
bed_dict['lncipedia_5_2_hc_hg38'] = exp_output_path.lncNRA_database_bed + 'lncipedia_5_2_hc_hg38.bed'
bed_dict['hs.gencode+'] = exp_output_path.lncNRA_database_bed + 'hs.GENCODE+.lncRNA.bed'        



bed_dict_order_list = list()



all_other_lncRNA_database = list()
all_other_lncRNA_database_minus_gencode=list()

bed_gc_ratios = dict()
bed_recovery_dict = dict()
for key in bed_dict:
    bed_recovery_dict[key] = bed_recovry_class(bed_dict[key], aggregate_exon_dict, aggregate_exon_IT, genome_fasta)
    #print("\nlncRNA database: %s" % (key))
    #bed_recovery_dict[key].print_stats()
    minus_gencode_exon_ids = el.exon_id_difference(bed_recovery_dict[key].exact_fuzzy_middle_exon_ids, lncRNA_middle_exon_ids)
    all_other_lncRNA_database += bed_recovery_dict[key].exact_fuzzy_middle_exon_ids
    all_other_lncRNA_database_minus_gencode += minus_gencode_exon_ids
    #print('fraction of recovered exon_ids that are NOT in gencode: ',len(minus_gencode_exon_ids)/len(bed_recovery_dict[key].exact_fuzzy_middle_exon_ids), '(the rational is these exon ids might look less good (e.g. splice site, expression counts) than the gencode annotated lncRNA exons)')
    #bed_gc_ratios[key] = print_gc_ratios(bed_recovery_dict[key].bed_exon_ids.middle_exon_id_list, aggregate_exon_dict, aggregate_exon_IT)
    



all_known_exon_ids_set = set(exon_GENCODE_comp['all'])
all_other_lncRNA_annotated = set()


for key in bed_recovery_dict:
    #print(len(bed_recovery_dict[key].bed_exon_ids.middle_exon_id_list))
    all_known_exon_ids_set      = all_known_exon_ids_set.union(bed_recovery_dict[key].bed_exon_ids.middle_exon_id_list)    
    all_other_lncRNA_annotated  =all_other_lncRNA_annotated.union(bed_recovery_dict[key].bed_exon_ids.middle_exon_id_list)
    


all_other_known_lncRNA_exon_ids_set = all_known_exon_ids_set.difference(exon_GENCODE_comp['all'])



threshold_exon_it_set = set()
for exon_id in aggregate_exon_dict:
    entry = aggregate_exon_dict[exon_id]
    if entry['count'] >= 100:
        threshold_exon_it_set.add(exon_id)












bed_file            = bed_dict['hs.gencode+']
gencode_plus_hs     = bl.simple_bed(bed_file, aggregate_exon_IT)
recovered_gencode_plus_hs_exon_id   = (el.exon_id_intersection(gencode_plus_hs.middle_exon_id_list, all_recovered_exon_id_list))
lncRNA_database_bed_dict['gencode_plus_hs'] = gencode_plus_hs



bed_file    = bed_dict['NONCODE'] 
NONCODE     = bl.simple_bed(bed_file, aggregate_exon_IT)
recovered_NONCODE_exon_id           = el.exon_id_intersection(NONCODE.middle_exon_id_list, all_recovered_exon_id_list)
lncRNA_database_bed_dict['NONCODE'] = NONCODE






CHESS_gff = exp_output_path.lncNRA_database_bed + 'chess2.2.gff'

#CHESS_gff_to_bed = exp_output_path.lncNRA_database_bed + 'chess2.2_exons.bed'

chess_line_counter = 0
with open(CHESS_gff) as f:
    CHESS_genes = dict()
    for ii, line in enumerate(f):
        if line[0]=='#':
            continue
        line_split=line.split('\t')
        if line_split[1] == 'CHESS' or line_split[1] == 'Gnomon':
            chess_line_counter += 1
            ID_line = line_split[8][3:]
            chess_ID = ID_line.split(';')[0]
            CHESS_genes[chess_ID] = list()

CHESS_transcripts=dict()
chess_transcripts = dict()
with open(CHESS_gff) as f:
    for ii, line in enumerate(f):
        if line[0]=='#':
            continue
        line_split=line.strip().split('\t')
        if line_split[2] == 'transcript':
            transcript_id = line_split[0]
            ID_line = line_split[8][3:].split(';')[0]
            chess_transcripts#[transcript_id]
            CHESS_transcripts[ID_line] = list()


with open(CHESS_gff) as f:
    for ii, line in enumerate(f):
        if line[0]=='#':
            continue
        line_split=line.strip().split('\t')
        if line_split[2] == 'exon':
            if line_split[8].find('Parent=') >= 0:
                parent = line_split[8][7:].split(';')[0]
                CHESS_transcripts[parent].append('1')




bed_file        = bed_dict['RNAcentral'] 
RNAcentral      = bl.simple_bed(bed_file, aggregate_exon_IT)
#recovered_RNAcentral_exon_id            = el.exon_id_intersection(RNAcentral.middle_exon_id_list, all_recovered_exon_id_list)
recovered_RNAcentral_exon_id            = el.exon_id_intersection(RNAcentral.fuzzy_overlap_middle_exon_ids, all_recovered_exon_id_list)
lncRNA_database_bed_dict['RNAcentral']  = RNAcentral










recovered_gencode_pc_exon_id_counts      =  el.get_exon_dict_counts( el.exon_id_intersection(pc_middle_exon_id_list, all_recovered_exon_id_list), aggregate_exon_dict)
recovered_gencode_lncRNA_exon_id_counts  =  el.get_exon_dict_counts( el.exon_id_intersection(lncRNA_middle_exon_ids, all_recovered_exon_id_list), aggregate_exon_dict)
recovered_gencode_plus_hs_exon_id_counts =  el.get_exon_dict_counts(recovered_gencode_plus_hs_exon_id, aggregate_exon_dict)
recovered_NONCODE_exon_id_counts         =  el.get_exon_dict_counts(recovered_NONCODE_exon_id, aggregate_exon_dict)

d1 = el.exon_id_difference(recovered_gencode_plus_hs_exon_id, lncRNA_exon_id_list)
d1 = el.get_exon_dict_counts(d1, aggregate_exon_dict)
d2 = el.exon_id_difference(recovered_NONCODE_exon_id, lncRNA_exon_id_list)
d2 = el.get_exon_dict_counts(d2, aggregate_exon_dict)

lncRNA_key_list = list()
lncRNA_gencode_removed_counts_list = list()
for key in bed_recovery_dict:
    q = bed_recovery_dict[key]
    recovered_non_gencode_exon_ids = set(q.recovered_exon_ids).difference(lncRNA_middle_exon_ids)
    counts_list = el.get_exon_dict_counts(recovered_non_gencode_exon_ids, aggregate_exon_dict)
    lncRNA_key_list.append(key)
    lncRNA_gencode_removed_counts_list.append(counts_list)


np.mean(el.get_exon_dict_counts(bed_recovery_dict['RNAcentral'].exact_fuzzy_middle_exon_ids, aggregate_exon_dict))




counts_list_list = [recovered_gencode_pc_exon_id_counts,recovered_gencode_lncRNA_exon_id_counts]
counts_list_list += lncRNA_gencode_removed_counts_list

key_list=['mRNA','lncRNA']
key_list += lncRNA_key_list


fig = plt.figure()
plt.boxplot(counts_list_list,showfliers=False)
plt.yscale('log')
plt.ylabel('exon expression counts')
plt.xticks(range(1,len(key_list)+1), key_list, rotation=30)
plt.title('exon expression counts for different lncRNA sources')
plt.tight_layout()
pdf_plots.savefig(fig)




txt_file_path = exp_output_path.Parse_ENSEMBLE_GENCODE_exons_load_figures_2
txt_file_path = txt_file_path + 'S5A.txt'

with open(txt_file_path, 'w') as f:
    f.write('category\tcounts\n')
    
    for counter in range(len(key_list)):
        f.write('{:}\t'.format(key_list[counter]))
        for ii, val in enumerate(counts_list_list[counter]):
            f.write('{:},')
        f.write('\n')
    
    
    








exon_regions_primary_exon_id_list = list()
count_exon_regions = 0
for chrom in tree_3ss:
    for strand in tree_3ss[chrom]:
        for interval in tree_3ss[chrom][strand]:
            count_exon_regions += 1
            
            if len(interval[2]) > 1:
                max_exon_id = el.get_max_count_exon_id_in_list(interval[2],aggregate_exon_dict)
                exon_regions_primary_exon_id_list.append(max_exon_id)
            else:
                exon_regions_primary_exon_id_list.append(interval[2][0])

exon_regions_primary_exon_id_list=list(set(exon_regions_primary_exon_id_list))






##################
################## get PC exon_id recovery fraction




mRNA_recovery_percent, highly_overlapping, exact_list, has_overlapping = el.recovery_percent(pc_middle_exon_id_list, aggregate_exon_IT, aggregate_exon_dict, genome_fasta)

lncRNA_recovery_percent, highly_overlapping, exact_list, has_overlapping =  el.recovery_percent(lncRNA_middle_exon_ids, aggregate_exon_IT, aggregate_exon_dict, genome_fasta)










#calculate the number of exons from the lncRNA databases
lncRNA_database_exons = set()
lncRNA_database_fuzzy_found_exons = set()
recovery_fraction_lncRNA_database_minus_gencode_dict=dict()
label_fraction_list = list()
#label_fraction_list.append(['mRNA', mRNA_recovery_percent])
label_fraction_list.append(['', 0])
label_fraction_list.append(['GENCODE lncRNA', lncRNA_recovery_percent])

found_lists  = list()
missed_lists = list()
all_lists    = list()

label_proportion_all_list = list()
for key in bed_recovery_dict:
    
    a=bed_recovery_dict[key]
    
    print(key, len(a.bed_exon_ids.middle_exon_id_list))
    lncRNA_database_exons=lncRNA_database_exons.union(a.bed_exon_ids.middle_exon_id_list)
    lncRNA_database_fuzzy_found_exons=lncRNA_database_fuzzy_found_exons.union(a.exact_fuzzy_middle_exon_ids)
    
    
    #recovery_fraction_lncRNA_database_minus_gencode_dict[key]=len(set(el.threshold_exon_ids(a.exact_fuzzy_middle_exon_ids,100,aggregate_exon_dict)).difference(lncRNA_middle_exon_ids+pc_middle_exon_id_list))/len(set(a.bed_exon_ids.middle_exon_id_list).difference(lncRNA_middle_exon_ids+pc_middle_exon_id_list))

    

    numerator = (set(el.threshold_exon_ids(a.exact_fuzzy_middle_exon_ids,100,aggregate_exon_dict)).difference(lncRNA_middle_exon_ids+pc_middle_exon_id_list))
    denom     = (set(a.bed_exon_ids.middle_exon_id_list).difference(lncRNA_middle_exon_ids+pc_middle_exon_id_list))
    
    recovery_frac = len(numerator)/len(denom)
    recovery_fraction_lncRNA_database_minus_gencode_dict[key]=recovery_frac
    

    found_exons = set(el.threshold_exon_ids(a.exact_fuzzy_middle_exon_ids,100,aggregate_exon_dict) ).difference(lncRNA_middle_exon_ids+pc_middle_exon_id_list)
    
    found_lists.append(list(numerator))
    all_lists.append(list(denom))
    missed_lists.append(list(denom.difference(numerator)))
    #missed = set(el.threshold_exon_ids(a.exact_fuzzy_middle_exon_ids,100,aggregate_exon_dict)).difference(lncRNA_middle_exon_ids+pc_middle_exon_id_list)
    #found_lists.append()


    label_fraction_list.append([key, recovery_fraction_lncRNA_database_minus_gencode_dict[key]])
    
    label_proportion_all_list.append([key, len(set(el.threshold_exon_ids(a.exact_fuzzy_middle_exon_ids,100,aggregate_exon_dict)).difference(lncRNA_middle_exon_ids+pc_middle_exon_id_list))])




exact_include_gencode = len(set(lncRNA_database_exons).intersection( el.threshold_exon_ids(aggregate_exon_dict.keys(),100,aggregate_exon_dict) ))/len(lncRNA_database_exons)


fuzzy_include_gencode = len(set(lncRNA_database_fuzzy_found_exons).intersection(el.threshold_exon_ids(aggregate_exon_dict.keys(),100,aggregate_exon_dict) ))/len(lncRNA_database_exons)



thresholded_unannotated_exon_ids = set(el.threshold_exon_ids(aggregate_exon_dict.keys(),100,aggregate_exon_dict)).difference(lncRNA_middle_exon_ids+pc_middle_exon_id_list)




fuzzy_100_read_count = len(set(lncRNA_database_fuzzy_found_exons).intersection(thresholded_unannotated_exon_ids))/len(lncRNA_database_exons)


exact_100_read_count = len(set(lncRNA_database_exons).intersection(thresholded_unannotated_exon_ids))/len(lncRNA_database_exons)


print("fraction of lncRNA databases recovered: ")

print('(annotated) exact_include_gencode: %.1f%%' %(exact_include_gencode*100))
print('(found, including fuzzy) fuzzy_include_gencode: %.1f%%' %(fuzzy_include_gencode*100))

print('(annotated) exact_100_read_count: %.1f%%' %(exact_100_read_count*100))
print('(found, including fuzzy) fuzzy_100_read_count: %.1f%%' %(fuzzy_100_read_count*100))

print('count all found at least 1 read (including) fuzzy exons: {:,}'.format(len(lncRNA_database_fuzzy_found_exons)))

print('count all annotated  exons: {:,}'.format(len(lncRNA_database_exons)))


len(lncRNA_database_exons.union(lncRNA_database_exons))




#fuzzy_100_read_count = len(set(lncRNA_database_fuzzy_found_exons).intersection(thresholded_unannotated_exon_ids))/len(lncRNA_database_exons)

labels   =['ET recovered','not recovered']
fractions=[fuzzy_100_read_count,1-fuzzy_100_read_count]






colors = ['b','b','orange','g','r','purple']
fig = plt.figure()
plt.title('percent lncRNA database (with gencode removed) recovered by ET')
labels, fractions = zip(*label_fraction_list)
for ii, val in enumerate(fractions):
    if 'mRNA' == labels[ii]:
        plt.bar(ii,fractions[ii], color=color_dict[labels[ii]])
    elif 'GENCODE lncRNA' == labels[ii]:
        plt.bar(ii,fractions[ii], color=color_dict['lncRNA'])
    else:
        plt.bar(ii,fractions[ii], color = colors[ii])
plt.ylim(0,1)
plt.yticks(np.arange(0,1.24,.25),['%s%%' % x for x in range(0,125,25)])
plt.legend(labels)
plt.tight_layout()
pdf_plots.savefig(fig)





out_path = exp_output_path.Parse_ENSEMBLE_GENCODE_exons_load_figures_2
out_path = out_path + '3F.txt'

with open(out_path, 'w') as f:
    f.write('category\tfound_exon_ids\tmissed_exon_ids\n')
    
    for ii, val in enumerate(found_lists):
        f.write(labels[ii] + '\t')
        for jj, val in enumerate(found_lists[ii]):
            f.write('{:},'.format(val))
        
        f.write('\t')
        for jj, val in enumerate(missed_lists[ii]):
            f.write('{:},'.format(val))
        
        f.write('\n')






pdf_plots.close()




