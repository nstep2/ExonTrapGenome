





import exon_id_library.exon_id_lib as el
import exon_id_library.bed as bl



intronic_snaptron_found_exon_id_list
#Intronic_pc_other_snaptron_removed



#run /home/dude/repositories/exon_def/slim_ET (copy)/bigwig_lib.py
#gencode to get GENCODE_gene_dict


'''
gtex_file_dict = {'general':"/media/8TB_1/downloaded_data/GTEX/v8/gtex_gene_v8.txt", 'tissue':"/media/8TB_1/downloaded_data/GTEX/v8/gtex_gene_v8_all_fields.txt", 'tissue_names':"/media/8TB_1/downloaded_data/GTEX/v8/gtex_gene_v8_names.txt"}

gtex_file_dict['tissue']


def gtex_bed_line_to_coord_id(line):
    line_split = line.strip().split('\t')
    coord_id = "{:}:{:}-{:}:{:}".format(line_split[0],line_split[1],line_split[2],line_split[5])
    tissue_counts= np.fromstring(line_split[9],sep=',')
    bio_type = line_split[7]
    return coord_id, tissue_counts, bio_type, tissue_counts

all_gene_tissue_dict =dict()
all_gene_tissue_stats = list()
with open(gtex_file_dict['tissue'],'r') as f:
    next(f)
    for ii, line in enumerate(f):
        coord_id, tissue_counts, bio_type, tissue_array = gtex_bed_line_to_coord_id(line)
        line_split=line.split('\t')
        #if coord_id == 'chrY:56878350-56878390:-':
        #        print(coord_id, 'found')
            
        if bio_type == 'protein_coding':
            
            
            stat_list =[coord_id, tissue_counts, max(tissue_counts), sum(tissue_counts)]
            all_gene_tissue_stats.append(stat_list)
            all_gene_tissue_dict[coord_id] =  {'tissue_counts':tissue_counts, 'max_tissue':max(tissue_counts), 'sum_tissue':sum(tissue_counts),'start':int(line_split[1]),'end':int(line_split[2]),'chrom':line_split[0],'strand':line_split[5]}


coord_ids, tissue_array, tissue_max, tissue_sum = zip(*all_gene_tissue_stats)


gtex_IT = make_IT(all_gene_tissue_dict,gencode_exon_IT)
gtex_IT['chr1']['+'][1:100000]


plt.figure()
plt.boxplot(np.log10(tissue_max),showfliers=False)
plt.yticks([0,1,2,3,4],['1','10','100','1,000','10,000'])
plt.ylabel('Max tissue RNA-seq rpkm')


plt.figure()
plt.boxplot(np.log10(tissue_sum),showfliers=False)
plt.yticks([0,1,2,3,4],['1','10','100','1,000','10,000'])
plt.ylabel('Sum tissue RNA-seq rpkm')



for gtex_key in ['Other_recovered', 'snap_only', 'Intronic_pc_other_snaptron_removed']:
    
    
    print(gtex_key)

    sum_score_list = list()
    max_score_list = list()
    found_gtex = dict()
    for ii, exon_stuff in  enumerate(query_cons_dict[gtex_key]):
        
        exon_id = exon_stuff[0]
        if aggregate_exon_dict[exon_id]['count'] < 10000:
            continue
        
        
        ex=el.exon_id_values(exon_id)
        intervals = gtex_IT[ex.chrom][ex.strand][ex.start:ex.end]
        for interval in intervals:
            gtex_id=interval[2]
            found_gtex[gtex_id] = all_gene_tissue_dict[gtex_id]
            max_score_list.append(all_gene_tissue_dict[gtex_id]['max_tissue'])
            sum_score_list.append(all_gene_tissue_dict[gtex_id]['sum_tissue'])
        #for interval in intervals:
            #found_gtex[]

    
    max_score_list = [ found_gtex[key]['max_tissue'] for key in found_gtex ]
    plt.figure()
    plt.boxplot(np.log10(max_score_list),showfliers=False)
    plt.yticks([0,1,2,3,4],['1','10','100','1,000','10,000'])
    plt.ylabel('Sum tissue RNA-seq rpkm')
    print(np.mean(max_score_list))
    
    sum_score_list = [ found_gtex[key]['sum_tissue'] for key in found_gtex ]
    plt.figure()
    plt.boxplot(np.log10(sum_score_list),showfliers=False)
    plt.yticks([0,1,2,3,4],['1','10','100','1,000','10,000'])
    plt.ylabel('Sum tissue RNA-seq rpkm')
    print(np.mean(sum_score_list))



38.705
517.6819999999999


38.67
549.5225

33.965
393.9149999999999


'''












#key_list = region_exon_id_sets_dict.keys()
key_list = ['Intronic']
cons_exon_sets_dict = dict()
for key in key_list:
    cons_exon_sets_dict[key] = region_exon_id_sets_dict[key]










#bigwig_path = exp_output_path.bigwig_path+'phastcons_17way/hg38.phastCons17way.bw'

#bw_data_phastcons_17 = bigwig_cons(bigwig_path)



bw_phastcon_path = exp_output_path.bigwig_path+'phastcons_30way/hg38.phastCons30way.bw'
bw_30way_phastcon = bigwig_cons(bw_phastcon_path)





#protein coding genes
bbb={gene_id:GENCODE_gene_dict[gene_id] for gene_id in GENCODE_gene_dict if GENCODE_gene_dict[gene_id]['tags']['gene_type']=='lncRNA'}
#len(bbb)


#protein coding genes
aaa={gene_id:GENCODE_gene_dict[gene_id] for gene_id in GENCODE_gene_dict if GENCODE_gene_dict[gene_id]['tags']['gene_type']=='protein_coding'}
#len(aaa)

#interval tree to isolate protein coding intronic exons
def make_IT(aaa, transcript_IT): # get chroms from an existing interval tree dict
    IT = dict()
    for chrom in transcript_IT:
        IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
    
    for gene_id in aaa:
        start = aaa[gene_id]['start']
        end = aaa[gene_id]['end']
        chrom = aaa[gene_id]['chrom']
        strand = aaa[gene_id]['strand']
        
        new_inteval = intervaltree.Interval(start, end, gene_id)
        IT[chrom][strand].add(new_inteval)

    return IT

pc_gene_IT = make_IT(aaa, gencode_exon_IT)
lncRNA_gene_IT = make_IT(bbb, gencode_exon_IT)


#phylop path
#bigwig_path = '/media/8TB_1/downloaded_data/cons/phylop_7way/hg38.phyloP7way.bw'
bigwig_path = exp_output_path.bigwig_path+'phylop_7way//hg38.phyloP7way.bw'
bw_data = bigwig_cons(bigwig_path)






def generate_cons_data(exon_id_list, pc_gene_IT, bw_scorer):
    
    start_time = time.time()

    intronic_exon_conservation_list = list()
    for ii, exon_id in enumerate(exon_id_list):
        
        if ii % int(len(exon_id_list)/10) == 0 and ii != 0:
            print('%d exons in %d seconds' % (ii, time.time()-start_time))
        
        ex = el.exon_id_values(exon_id)
        
        #check that intronic exon in pc_gene
        intervals = pc_gene_IT[ex.chrom][ex.strand][ex.start:ex.end]
        if len(intervals) == 0:
            continue
        
        
        score_array = bw_scorer.get_base_scores(ex.chrom, ex.start, ex.end)
        if bw_scorer.check_if_none(score_array) == False:
            mean_score = np.mean(score_array)
            #intronic_exon_conservation_list.append([exon_id, mean_score])
            
            up_score_array = bw_scorer.get_base_scores(ex.chrom, ex.start-75, ex.start-25)
            
            
            dn_score_array = bw_scorer.get_base_scores(ex.chrom, ex.end+25, ex.end+75)
            
            if bw_scorer.check_if_none(up_score_array) == False and bw_scorer.check_if_none(dn_score_array) == False:
                
                if ex.strand == '+':
                    intronic_exon_conservation_list.append([exon_id, mean_score, np.mean(up_score_array), np.mean(dn_score_array)])
                else:
                    intronic_exon_conservation_list.append([exon_id, mean_score, np.mean(dn_score_array), np.mean(up_score_array)])

    return intronic_exon_conservation_list





from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC



def frame_stop_codon_count(exon_id, aggregate_exon_dict, genome_fasta):
    frame_counts = {0:0, 1:0, 2:0}
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    seq = str( el.get_seq_from_exon_id(exon_id, genome_fasta) )
    ex = el.exon_id_values(exon_id)
    len(seq)
    
    frame_counts[0] = str(Seq(seq).translate()).count('*')
    frame_counts[1] = str(Seq(seq[1:-2]).translate()).count('*')
    frame_counts[2] = str(Seq(seq[2:-1]).translate()).count('*')
    
    '''for ii in range(0,len(seq)-2,3):
        codon = seq[ii:ii+3]
        if codon in stop_codons:
            frame_counts[0] += 1
            #print('\t', codon, ii)
        #print(codon,ii)
    
    for ii in range(1,len(seq)-2,3):
        codon = seq[ii:ii+3]
        if codon in stop_codons:
            frame_counts[1] += 1
    
    for ii in range(2,len(seq)-2,3):
        codon = seq[ii:ii+3]
        if codon in stop_codons:
            frame_counts[2] += 1
    '''
    
    has_stop_codon = {0:0,1:0,2:0}
    for ii, val in enumerate(frame_counts):
        if frame_counts[ii] > 0:
            has_stop_codon[ii] = 1
    
    return frame_counts, has_stop_codon



exon_id = 'chr19:30493460-30493577:+'


frame_stop_codon_count(exon_id, aggregate_exon_dict, genome_fasta)














class bed_recovry_class():
    def __init__(self, bed_file, exon_id_list, aggregate_exon_IT, genome_fasta):
        all_recovered_exon_id_list = set(exon_id_list)
        
        self.bed_exon_ids = bl.simple_bed(bed_file, aggregate_exon_IT)
        
        self.recovered_exon_ids =  el.exon_id_intersection( self.bed_exon_ids.middle_exon_id_list, all_recovered_exon_id_list)
        
        self.fuzzy_middle_exon_ids = (el.exon_id_intersection( self.bed_exon_ids.fuzzy_overlap_middle_exon_ids, all_recovered_exon_id_list))
        
        self.exact_fuzzy_middle_exon_ids = list(set(self.recovered_exon_ids).union(self.fuzzy_middle_exon_ids))
        
        self.stats = dict()
        self.stats['middle_exon_count']=len(self.bed_exon_ids.middle_exon_id_list)
        self.stats['recovered_middle_exon_count']=len(self.recovered_exon_ids)
        self.stats['fuzzy_middle_exon_count']=len(self.bed_exon_ids.fuzzy_overlap_middle_exon_ids)
        
    def print_stats(self):
        
        print_list = ['middle_exon_count' , 'recovered_middle_exon_count', 'fuzzy_middle_exon_count', 'recovered_fuzzy_middle_exon_count']
        
        middle_exon_ratio = self.stats['recovered_middle_exon_count']/self.stats['middle_exon_count']
        
        fuzzy_middle_exon_ratio = len(set(self.recovered_exon_ids).union(self.fuzzy_middle_exon_ids))/self.stats['middle_exon_count']
        
        print(middle_exon_ratio)
        print(fuzzy_middle_exon_ratio)
        



exp_output_path.exon_database_bed


bed_dict = dict()
bed_dict['old_genes_bed']       = exp_output_path.exon_database_bed + 'old_UCSC_genes.txt'
bed_dict['refseq_bed']          = exp_output_path.exon_database_bed + 'ncbi_refseq.txt'
bed_dict['gencode_v40_bed']     = exp_output_path.exon_database_bed + 'gencode_v40.txt'
#bed_dict['UCSC_Alt_Events']     = exp_output_path.exon_database_bed + 'UCSC_Alt_Events.txt'
bed_dict['human_mRNAs']         = exp_output_path.exon_database_bed + 'mRNA_and_est_human_mRNAs.txt'


bed_dict['AUGUSTUS']     = exp_output_path.exon_database_bed + 'AUGUSTUS.txt'
bed_dict['est_human_mRNAs']     = exp_output_path.exon_database_bed + 'est_human_mRNAs.txt'
bed_dict['genscan']     = exp_output_path.exon_database_bed + '/genscan.txt'
bed_dict['Swiss_institute_biology']     = exp_output_path.exon_database_bed + 'Swiss_institute_biology.txt'


#bed_dict['TransMap_ESTs']     = '/media/8TB_1/downloaded_data/ucsc_genes_exons/TransMap_ESTs.txt'








#exon_list_dict = dict()


'''
key = 'old_genes_bed'
file_name_path = bed_dict[key]
exon_list_dict[key] = el.load_bed_to_exon_id_list(file_name_path)

len(exon_list_dict[key])

exon_list_dict[key][:5]

len( set(exon_list_dict[key]).intersection(intron_interior_set) )
'''




recovered_intronic_exon_list = list()
bed_recovery_dict = dict()
for key in bed_dict:
    bed_recovery_dict[key] = bed_recovry_class(bed_dict[key], intron_interior_set, aggregate_exon_IT, genome_fasta)
    print("\nDatabase: {:} ({:,})".format(key, len(set(bed_recovery_dict[key].recovered_exon_ids))))
    bed_recovery_dict[key].print_stats()
    
    recovered_intronic_exon_list += bed_recovery_dict[key].recovered_exon_ids
    recovered_intronic_exon_list = list(set(recovered_intronic_exon_list))
    
    print('total found exons {:,}'.format(len(recovered_intronic_exon_list)))



for key in bed_recovery_dict:
    entry = bed_recovery_dict[key]
    
    other_sets = list()
    for key2 in bed_recovery_dict:
        if key2 == key:
            continue
        other_sets += bed_recovery_dict[key2].recovered_exon_ids
    
    val = set(entry.recovered_exon_ids).difference(other_sets)
    print("{:}: {:} intronic found and {:,} annotation unique exons".format(key, len(entry.recovered_exon_ids), len(val)))
    
    







'''
key = 'gencode_v40_bed'
bbb = bed_recovery_dict[key]

bbb.recovered_exon_ids[:5]
len(bbb.recovered_exon_ids)
len(set(bbb.recovered_exon_ids).intersection(intron_interior_set))

len(recovered_intronic_exon_list)
len(set(bbb.recovered_exon_ids).intersection(intron_interior_set))
'''







'''

snaptron_data = dict()

snaptron_data['snaptron_found_count'] = len(intronic_snaptron_found_exon_id_list)

snaptron_data['snaptron_shared_count'] = len(set(intronic_snaptron_found_exon_id_list).intersection(recovered_intronic_exon_list))

snaptron_data['snaptron_only_count'] = snaptron_found_count - snaptron_shared_count

snaptron_data['snaptron_and_annotated_count'] = len(set(intronic_snaptron_found_exon_id_list).union(recovered_intronic_exon_list))

'''




'''
#print(snaptron_data)
for key in snaptron_data:
    print('{:}\t{:,}'.format(key, snaptron_data[key]))
'''





'''

#labels,y, up, dn = zip(*intronic_exon_conservation_list)

#conservation_dict = {entry[0]:entry[1] for entry in intronic_exon_conservation_list}

conservation_dict = {entry[0]:entry[1] for entry in intronic_exon_30way_phast_list}



all_annotation_exon_conservation = list()
for ii, exon_id in enumerate(set(conservation_dict.keys()).intersection(recovered_intronic_exon_list)):
    score = conservation_dict[exon_id]
    all_annotation_exon_conservation.append(score)


plt.figure()
plt.hist(all_annotation_exon_conservation, bins = 50)





all_intron_exon_conservation = list()
for ii, exon_id in enumerate(set(conservation_dict.keys()).intersection(intron_interior_set)):
    score = conservation_dict[exon_id]
    all_intron_exon_conservation.append(score)


plt.figure()
plt.hist(all_intron_exon_conservation, bins = 50)

'''





'''

top_phastcons_exon_ids = [key for key in conservation_dict if conservation_dict[key] > .5 and aggregate_exon_dict[key]['count']>5000]
len(top_phastcons_exon_ids)

xx_start = 0
for ii in range(xx_start,xx_start+5):
    #exon_id = list(conservation_dict.keys())[ii]
    exon_id = list(top_phastcons_exon_ids)[ii]
    plot_exon_id_conservation(exon_id, bw_30way_phastcon)

'''




#super interesting
#'chrX:148522108-148522244:+'     #counts > 5000









#Other_recovered
#snap_only





'''
other_strand = {'+':'-','-':'+'}
cons_exon_sets_dict['Antisense_pc'] = list()
for exon_id in cons_exon_sets_dict['Intronic']:
    ex = el.exon_id_values(exon_id)
    intervals = pc_gene_IT[ex.chrom][other_strand[ex.strand]][ex.start:ex.end]
    if len(intervals) == 0:
        continue
    cons_exon_sets_dict['Antisense_pc'].append(exon_id)
'''

'''
np.median(el.get_exon_dict_counts(cons_exon_sets_dict['Antisense_pc'],aggregate_exon_dict))

np.median(el.get_exon_dict_counts(cons_exon_sets_dict['Antisense_lncRNA'],aggregate_exon_dict))

np.median(el.get_exon_dict_counts(cons_exon_sets_dict['Intronic_lncRNA'],aggregate_exon_dict))


np.median(el.get_exon_dict_counts(cons_exon_sets_dict['Intronic_pc'],aggregate_exon_dict))
'''




'''
other_strand = {'+':'-','-':'+'}
cons_exon_sets_dict['Antisense_lncRNA'] = list()
for exon_id in cons_exon_sets_dict['Intronic']:
    ex = el.exon_id_values(exon_id)
    intervals = lncRNA_gene_IT[ex.chrom][other_strand[ex.strand]][ex.start:ex.end]
    if len(intervals) == 0:
        continue
    cons_exon_sets_dict['Antisense_lncRNA'].append(exon_id)
'''



cons_exon_sets_dict['Snaptron'] = intronic_snaptron_found_exon_id_list

cons_exon_sets_dict['Other_recovered'] = recovered_intronic_exon_list


cons_exon_sets_dict['Intronic_pc'] = list()
for exon_id in cons_exon_sets_dict['Intronic']:
    ex = el.exon_id_values(exon_id)
    intervals = pc_gene_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) == 0:
        continue
    cons_exon_sets_dict['Intronic_pc'].append(exon_id)



'''
cons_exon_sets_dict['Intronic_lncRNA'] = list()
for exon_id in cons_exon_sets_dict['Intronic']:
    ex = el.exon_id_values(exon_id)
    intervals = lncRNA_gene_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) == 0:
        continue
    cons_exon_sets_dict['Intronic_lncRNA'].append(exon_id)
len(cons_exon_sets_dict['Intronic_lncRNA'])
'''


#cons_exon_sets_dict['conserved_introns'] =  [entry[0] for entry in intronic_exon_conservation_list if entry[1] > .25 ]
#cons_exon_sets_dict['pc_introns'] =  [entry[0] for entry in intronic_exon_conservation_list if (2*entry[1]+entry[2]+entry[3])/4 > .5 ]

cons_exon_sets_dict['Intronic_pc_other_removed'] = list(set(cons_exon_sets_dict['Intronic_pc']).difference(recovered_intronic_exon_list))

cons_exon_sets_dict['Intronic_pc_other_snaptron_removed'] = list(set(cons_exon_sets_dict['Intronic_pc']).difference(intronic_snaptron_found_exon_id_list).difference(recovered_intronic_exon_list))




'''
high_exon_threshold = 40000

for key in list(cons_exon_sets_dict.keys()):
    cons_exon_sets_dict[key + '_high_threshold'] = el.threshold_exon_ids(set(cons_exon_sets_dict[key]).intersection(aggregate_exon_dict.keys()), high_exon_threshold, aggregate_exon_dict)

'''


#cons_exon_sets_dict['Intronic_pc_other_snaptron_removed']


#cons_exon_sets_dict['conserved_introns']

#cons_exon_sets_dict['Snaptron']

#cons_exon_sets_dict['Other_recovered']


#cons_exon_sets_dict['both_snap_other']  = list(set(cons_exon_sets_dict['Snaptron']).intersection(cons_exon_sets_dict['Other_recovered']))


cons_exon_sets_dict['Other_recovered'] = recovered_intronic_exon_list


cons_exon_sets_dict['snap_only']  = list(set(cons_exon_sets_dict['Snaptron']).difference(cons_exon_sets_dict['Other_recovered']))


cons_exon_sets_dict['Intronic_pc_other_snaptron_removed'] = list(set(cons_exon_sets_dict['Intronic_pc']).difference(intronic_snaptron_found_exon_id_list).difference(recovered_intronic_exon_list))



#cons_exon_sets_dict['other_only'] = list(set(cons_exon_sets_dict['Other_recovered']).difference(cons_exon_sets_dict['Snaptron']))

#query_sets = ['Intronic_pc_other_snaptron_removed','conserved_introns','Snaptron','Other_recovered','other_only','snap_only']

#query_sets = ['Intronic_pc_other_snaptron_removed','both_snap_other','Other_recovered','Snaptron','other_only','snap_only']

#query_sets_labels = ['category 3','category 1 & 2','category 1','category 2 \n(includes cat 1)','category 1 only','category 2 only']

#query_sets = ['Intronic_pc_other_snaptron_removed','conserved_introns','Other_recovered','Snaptron','other_only','snap_only']

#len(cons_exon_sets_dict['Intronic_pc_other_snaptron_removed'])







        





'''

fig_param = dict()
for key in cons_exon_sets_dict: 
    fig_param[key]={'ylim1':0,'ylim2':0.5,'threshold':100,'max_exons':100000}

#fig_param['Intronic_other_removed']['threshold'] = 30000

fig_param['mRNA']['ylim2'] = 1.0
fig_param['mRNA_high_threshold']['ylim2'] = 1.0

fig_param['conserved_introns']['ylim2'] = 1.0
fig_param['conserved_introns_high_threshold']['ylim2'] = 1.0
fig_param['Intronic_pc_other_removed']['ylim2'] = 1.0
fig_param['Intronic_pc_other_snaptron_removed']['ylim2'] = 1.0



for key in cons_exon_sets_dict:
    fig_param[key]['max_exons'] = 1000

for key in cons_exon_sets_dict:
    if key.find('_high_threshold') > 0:
        fig_param[key]['threshold'] = high_exon_threshold


'''









"""



#zzz = cons_exon_sets_dict['Intronic_pc_high_threshold']

zzz = list(set(cons_exon_sets_dict['Intronic_pc_other_removed']).intersection(intronic_snaptron_found_exon_id_list))

#zzz = list(set( cons_exon_sets_dict['conserved_introns_high_threshold']).difference(intronic_snaptron_found_exon_id_list).difference(recovered_intronic_exon_list))
zzz = list(set(intronic_snaptron_found_exon_id_list).difference(recovered_intronic_exon_list))
len(zzz)

#chr7:849520-849604:+

for ii in range(0,190):
    exon_id = zzz[ii]
    if aggregate_exon_dict[exon_id]['count'] < 1000:
        continue
    plot_exon_id_conservation(exon_id, bw_30way_phastcon)
    #ex = el.exon_id_values(exon_id)
    #exon_id = "%s:%d-%d:%s" % (ex.chrom,ex.start+1500, ex.end+1500, ex.strand)
    #plot_exon_id_conservation(exon_id, bw_30way_phastcon)
    print(exon_id, '\t', 'count: %d' % (aggregate_exon_dict[exon_id]['count']))
    
zzz

# snaptron only

'''
snap_only_eid = [  
'chr16:75362875-75363036:-'    ,
    'chr4:159266880-159266959:+',  #EST
    'chr1:35984470-35984550:+',    #EST
    'chr1:203150632-203150782:+',  #EST
    
    'chr5:31488982-31489082:-',    #EST
    'chr8:42536355-42536472:-',    #EST
    'chr4:151266086-151266246:-',  #EST
    ]
'''
'''
chr20:47318951-47319062:- 	 count: 7036
chr3:49748869-49748927:- 	 count: 1714  #testis
chr2:6978626-6978840:+ 	 count: 4738   #optic nerve
chrX:115019850-115019978:- 	 count: 54782  #iris adipose
chr5:31488982-31489082:- 	 count: 1118 #prostate
chr3:125192586-125192693:- 	 count: 6613  #embryotic
chr16:87389551-87389714:+ 	 count: 3583 #cerabellum 
chr8:42536355-42536472:- 	 count: 1242  #hippocampus
chr3:53046460-53046790:- 	 count: 22117 #fetal kidney
chr9:114035188-114035236:+ 	 count: 3302 #brain
chr16:75362875-75363036:- 	 count: 3239 #unkown
'''




for exon_id in intronic_only_eid:
    if exon_id in aggregate_exon_dict:
        print(exon_id)


intronic_only_eid = [   
    'chr3:28835817-28836021:+',   #maybe
    'chr10:76176175-76176302:+',  #in middle large conserved region
    'chr22:30803024-30803178:+',  #sharp 5ss, 
    'chr7:69934783-69934898:+',   #meh
    'chr9:116699593-116699681:-', #nope
    'chr5:31193947-31194142:+',  #yup
    
    'chrX:72027586-72027700:+', #nope
    
    
    ]



for exon_id in intronic_only_eid:
    plot_exon_id_conservation(exon_id, bw_30way_phastcon)
    print(exon_id, '\t', 'count: %d' % (aggregate_exon_dict[exon_id]['count']))




for exon_id in snap_only_eid:
    plot_exon_id_conservation(exon_id, bw_30way_phastcon)
    print(exon_id, '\t', 'count: %d' % (aggregate_exon_dict[exon_id]['count']))


"""




111






'''


def plot_conservation_data(top_phastcons_exon_ids, fig_param, key):
    
    import scipy
    import scipy.stats
    
    
    number_found = len(top_phastcons_exon_ids)
    
    window = fig_param[key]['max_exons']
    window = min(number_found, window)
    #window = 100
    
    pad_length = 350
    exon_divis_3 = {'3':0,'not 3':0}
    left_array = np.zeros(2*pad_length)
    right_array = np.zeros(2*pad_length)
    
    import random
    random.seed(3000)
    
    
    xx_start = 0    
    #for ii in range(xx_start,xx_start+window):
    for ii in random.sample(range(number_found),window):
        #exon_id = list(conservation_dict.keys())[ii]
        exon_id = list(top_phastcons_exon_ids)[ii]
        ex = el.exon_id_values(exon_id)
        if ex.length % 3 == 0:
            exon_divis_3['3'] += 1
        else:
            exon_divis_3['not 3'] += 1
        #exon_id = "%s:%d-%d:%s" % (ex.chrom,ex.start+1500, ex.end+1500, ex.strand)
        d = plot_exon_id_conservation(exon_id, bw_30way_phastcon, pad_length = pad_length, no_plot = True)
        #d = plot_exon_id_conservation(exon_id, bw_data, pad_length = pad_length, no_plot = True)
        #d = plot_exon_id_conservation(exon_id, bw_data)
        if len(d) != 1:
            left_array  = left_array  + np.array(d[:2*pad_length])
            right_array = right_array + np.array(d[-2*pad_length:])
    
    fig = plt.figure()
    plt.plot(range(-1*pad_length, pad_length),left_array/window )
    plt.xlabel('bp from 3SS')
    plt.ylim(fig_param[key]['ylim1'], fig_param[key]['ylim2'])
    pval = scipy.stats.binom_test(exon_divis_3['3']-1, exon_divis_3['3']+exon_divis_3['not 3'], p=0.33333333, alternative='greater')
    plt.legend(['divisible 3 = {:}; not divisible 3 = {:}\n(not 3)/(3) ratio={:.3} & pval={:.1E}'.format(exon_divis_3['3'], exon_divis_3['not 3'], exon_divis_3['not 3']/exon_divis_3['3'],pval)])
    plt.title("{:} - total found: {:,}\nThreshold:{:,}".format(key, number_found, fig_param[key]['threshold']))
    pdf_plots.savefig(fig)
    
    fig = plt.figure()
    plt.plot(range(-1*pad_length, pad_length),right_array/window )
    plt.xlabel('bp from 5SS')
    plt.ylim(fig_param[key]['ylim1'], fig_param[key]['ylim2'])
    pval = scipy.stats.binom_test(exon_divis_3['3']-1, exon_divis_3['3']+exon_divis_3['not 3'], p=0.33333333, alternative='greater')
    plt.legend(['divisible 3 = {:}; not divisible 3 = {:}\n(not 3)/(3) ratio={:.3} & pval={:.1E}'.format(exon_divis_3['3'], exon_divis_3['not 3'], exon_divis_3['not 3']/exon_divis_3['3'],pval)])
    plt.title("{:} - total found: {:,}\nThreshold:{:,}".format(key, number_found, fig_param[key]['threshold']))
    pdf_plots.savefig(fig)
    
    return range(-1*pad_length, pad_length), left_array/window, right_array/window 


'''


'''



outdir = exp_output_path.old_gene 
pdf_plots = PdfPages(outdir+'old_gene__threshold_%d.pdf' % (high_exon_threshold))

plot_5ss_dict = dict()

plot_3ss_dict = dict()

for list_key in cons_exon_sets_dict:
    print(list_key)
    exon_id_list = set(aggregate_exon_dict.keys()).intersection(cons_exon_sets_dict[list_key])
    
    top_phastcons_exon_ids = [key for key in exon_id_list if key in aggregate_exon_dict and aggregate_exon_dict[key]['count']> fig_param[list_key]['threshold'] ]
    
    len(top_phastcons_exon_ids)
    
    x, y, z = plot_conservation_data(top_phastcons_exon_ids, fig_param, list_key)
    
    #if list_key in ['conserved_introns_high_threshold','Other_recovered_high_threshold','Intronic_pc_high_threshold','lncRNA_high_threshold','mRNA_high_threshold','snap_only_high_threshold']:
    if list_key in ['conserved_introns_high_threshold','Other_recovered_high_threshold','Intronic_pc_high_threshold','lncRNA_high_threshold','snap_only_high_threshold']:
        plot_5ss_dict[list_key] = x, z
        plot_3ss_dict[list_key] = x, y


fig = plt.figure()
for key in plot_3ss_dict:
    x, y = plot_3ss_dict[key]
    plt.plot(x,y,label=key)
plt.ylim(0,.4)
plt.legend()
plt.xlabel('3ss')

fig = plt.figure()
for key in plot_5ss_dict:
    x, y = plot_5ss_dict[key]
    plt.plot(x,y,label=key)
plt.ylim(0,.4)
plt.legend()
plt.xlabel('5ss')

pdf_plots.close()

'''
















'''


#relates to recursive splicing

r_5ss_dict = dict()
r_3ss_dict = dict()

with open('/home/pdf/Downloads/mmc5_sheet_nested.csv','r') as f:
    next(f)
    for line in f:
        line_split = line.split(' ')
        
        chrom  = line_split[0]
        start  = int(line_split[1])
        end    = int(line_split[2])
        strand = line_split[3]
        count  = int(line_split[4])
        
        
        if strand == '+':
            jid = '%s:%d:%s' % (chrom, end, strand)
            r_5ss_dict[jid]=count
            jid = '%s:%d:%s' % (chrom, start, strand)
            r_3ss_dict[jid]=count
        else:
            jid = '%s:%d:%s' % (chrom, end, strand)
            r_3ss_dict[jid]=count
            jid = '%s:%d:%s' % (chrom, start, strand)
            r_5ss_dict[jid]=count





'''

'''

def search_intronic(exon_id_list, r_3ss_dict, r_5ss_dict, aggregate_exon_dict):
    count_found = 0
    seq_read_counts = dict()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        #if strand == '+':
        #for ex.strand in ['+','-']:
        for ii in range(-3,4,1):
            jid = '%s:%d:%s' % (ex.chrom, ex.start+ii, ex.strand)
            if jid in r_5ss_dict or jid in r_3ss_dict:
                count_found += 1
                seq_read_counts[exon_id] = aggregate_exon_dict[exon_id]['count']
                
            
            jid = '%s:%d:%s' % (ex.chrom, ex.end+ii, ex.strand)
            if jid in r_3ss_dict or jid in r_3ss_dict:
                count_found += 1
                seq_read_counts[exon_id] = aggregate_exon_dict[exon_id]['count']
                
    
    return count_found, seq_read_counts
                    

exon_id_list = intron_interior_set
count_found, seq_read_counts = search_intronic(exon_id_list, r_3ss_dict, r_5ss_dict, aggregate_exon_dict)
print('count_found', count_found)

np.median([seq_read_counts[x] for x in seq_read_counts])

            
'''            
            
        












def generate_cons_mean_max(exon_id_list, pc_gene_IT, bw_scorer):
    
    start_time = time.time()
    intronic_exon_conservation_list = list()
    for ii, exon_id in enumerate(exon_id_list):
        
        if ii % int(len(exon_id_list)/10) == 0 and ii != 0:
            print('%d exons in %d seconds' % (ii, time.time()-start_time))
        
        ex = el.exon_id_values(exon_id)
        
        #check that intronic exon in pc_gene
        intervals = pc_gene_IT[ex.chrom][ex.strand][ex.start:ex.end]
        if len(intervals) == 0:
            continue
        
        score_array = bw_scorer.get_base_scores(ex.chrom, ex.start, ex.end)
        if bw_scorer.check_if_none(score_array) == False:
            mean_score = np.mean(score_array)
            max_score  = max(score_array)
            intronic_exon_conservation_list.append([exon_id, mean_score, max_score])
            

    return intronic_exon_conservation_list









query_list = ['Other_recovered','snap_only', 'Intronic_pc_other_snaptron_removed']
query_cons_dict = dict()
for key in query_list:
    exon_id_list = cons_exon_sets_dict[key]
    query_cons_dict[key] =  generate_cons_mean_max(exon_id_list, pc_gene_IT, bw_30way_phastcon)











query_cons_dict = dict()
for key in query_list:
    exon_id_list = cons_exon_sets_dict[key]
    query_cons_dict[key] =  generate_cons_mean_max(exon_id_list, pc_gene_IT, bw_30way_phastcon)

#len(cons_exon_sets_dict['Other_recovered'])
#len(cons_exon_sets_dict['snap_only'])
#len(cons_exon_sets_dict['Intronic_pc_other_snaptron_removed'])

sorted_cons = sorted(query_cons_dict['Intronic_pc_other_snaptron_removed'], key=lambda x: (el.exon_id_values(x[0]).chrom, el.exon_id_values(x[0]).strand, el.exon_id_values(x[0]).start), reverse=True)

    
    
sorted_cons = sorted(query_cons_dict['Intronic_pc_other_snaptron_removed'], key=lambda x: x[1], reverse=True)




'''



from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
#str(Seq(seq).translate()

def count_AA(seq, input_dict):
    for aa in seq:
        if aa not in input_dict:
            input_dict[aa] = 0
        input_dict[aa] += 1
    return input_dict

#amino_acids_frequencies_dict = count_AA(seq, amino_acids_frequencies_dict )

amino_acids_frequencies_dict = dict()
fp = open(exp_output_path.old_gene+'amino_acid_divisible_3_intronic_exons.fa','w')
f = open(exp_output_path.old_gene+'divisible_3_intronic_exons.fa','w')
count_found = 0
#exon_ids,c1,c2 = zip(*sorted_cons[710:720])
exon_ids,c1,c2 = zip(*sorted_cons)

divisble_3_non_stop_codon_exon_ids = list()
for exon_id in exon_ids[0000:-1]:
    if aggregate_exon_dict[exon_id]['count'] < 10000 or el.exon_id_values(exon_id).length % 3 != 0:
        continue
    
    frame_counts, frame_has_stop_codon = frame_stop_codon_count(exon_id, aggregate_exon_dict, genome_fasta)
    if sum([frame_has_stop_codon[key] for key in frame_counts]) == 3:
        continue
    
    count_found += 1
    divisble_3_non_stop_codon_exon_ids.append(exon_id)
    
    #plot_exon_id_conservation(exon_id, bw_30way_phastcon)
    out_string = "{:} {:}\t{:}\t{:}".format('>',exon_id, aggregate_exon_dict[exon_id]['count'], '%.1f, %.1f' % (aggregate_exon_dict[exon_id]['3ss_score'],aggregate_exon_dict[exon_id]['5ss_score']))
    print('>',exon_id, aggregate_exon_dict[exon_id]['count'], '%.1f, %.1f' % (aggregate_exon_dict[exon_id]['3ss_score'],aggregate_exon_dict[exon_id]['5ss_score']))
    seq = str( el.get_seq_from_exon_id(exon_id, genome_fasta) )
    print(seq)
    f.write(out_string + '\n')
    f.write(seq + '\n')
    
    print(frame_has_stop_codon)
    #if frame_has_stop_codon[0] == 0:
    
    if str(Seq(seq).translate()).count('*') == 0:
        fp.write(out_string+'  :0' + '\n')
        fp.write( str(Seq(seq).translate()) + '\n')
        print(Seq(seq).translate(), 0)
        amino_acids_frequencies_dict = count_AA(str(Seq(seq).translate()), amino_acids_frequencies_dict )
    #if frame_has_stop_codon[1] == 0:
    if str(Seq(seq[1:-2]).translate()).count('*') == 0:
        fp.write(out_string+'  :1' + '\n')
        fp.write( str(Seq(seq[1:-2]).translate()) + '\n')
        print(Seq(seq[1:-2]).translate(), 1)
        amino_acids_frequencies_dict = count_AA(str(Seq(seq[1:-2]).translate()), amino_acids_frequencies_dict )
    if str(Seq(seq[2:-1]).translate()).count('*') == 0:
    #if frame_has_stop_codon[1] == 0:
        fp.write(out_string+'  :2' + '\n')
        fp.write( str(Seq(seq[2:-1]).translate()) + '\n')
        print(Seq(seq[2:-1]).translate(), 2)
        amino_acids_frequencies_dict = count_AA(str(Seq(seq[2:-1]).translate()), amino_acids_frequencies_dict )
    
    
    
    
f.close()    
fp.close()    
print('count_found', count_found)    

'''


'''


import random
random.seed(121212)
thresh_aa_pair_list = list()
tot = 0
for ii, val in enumerate(amino_acids_frequencies_dict):
    
    thresh_aa_pair_list.append([tot, val])
    tot += amino_acids_frequencies_dict[val]

def generate_AA_string(thresh_aa_pair_list,tot):
    aa_list = list()
    for ii in range(200):
        sample = random.randint(0,tot)
        #print(sample)
        
        found = False
        for jj in range(len(thresh_aa_pair_list)-1):
            if sample > thresh_aa_pair_list[jj][0] and sample < thresh_aa_pair_list[jj+1][0]:
                aa_list.append(thresh_aa_pair_list[jj][1])
                found = True
                break
        if found == False:
            aa_list.append(thresh_aa_pair_list[-1][1])
    return ''.join(aa_list)



f = open(exp_output_path.old_gene+'random_aa.fa','w')
for ii in range(2000):
    f.write('> %d\n' % (ii))
    aa_seq =     generate_AA_string(thresh_aa_pair_list,tot)    
    f.write(aa_seq + '\n')
f.close()    


'''




    


exon_ids,c1,c2 = zip(*sorted_cons)
#query_cons_phylop_dict = dict()
#for key in exon_ids:
    #exon_id_list = cons_exon_sets_dict[key]
    
#query_cons_phylop_list =  generate_cons_mean_max(exon_ids, pc_gene_IT, bw_data)





'''

amino_acid_prediction_file = '/media/8TB_1/downloaded_data/interproscan-5.56-89.0/et_data/amino_acid_divisible_3_intronic_exons.tsv'

amino_acid_prediction_dict = dict()
with open(amino_acid_prediction_file,'r') as f:
    for ii, line in enumerate(f):
        line_split = line.strip().split('\t')
        exon_id = line_split[0]
        db_name = line_split[3]
        val = line_split[8]
        name = line_split[4]
        description = line_split[5]
        if db_name not in         amino_acid_prediction_dict:
            amino_acid_prediction_dict[db_name] = list()
        amino_acid_prediction_dict[db_name].append([exon_id,name, description, val])


amino_acid_prediction_dict.keys()


for key in amino_acid_prediction_dict:
    print(key, len(amino_acid_prediction_dict[key]))


for val in amino_acid_prediction_dict['PANTHER']:
    plot_exon_id_conservation(val[0], bw_30way_phastcon)

for val in amino_acid_prediction_dict['Pfam']:
    plot_exon_id_conservation(val[0], bw_30way_phastcon)

for val in amino_acid_prediction_dict['PRINTS']:
    plot_exon_id_conservation(val[0], bw_30way_phastcon)

"""for val in amino_acid_prediction_dict['MobiDBLite']:
    plot_exon_id_conservation(val[0], bw_30way_phastcon)"""
    

for val in amino_acid_prediction_dict['ProSiteProfiles']:
    plot_exon_id_conservation(val[0], bw_30way_phastcon)    



exon_id_list,name, descript, val = zip(*amino_acid_prediction_dict['MobiDBLite'])

exon_id_list,name, descript, val = zip(*amino_acid_prediction_dict['PANTHER'])


amino_acid_prediction_dict.keys()


for key in ['PANTHER', 'ProSiteProfiles', 'Pfam', 'PRINTS']:
    print(key)
    exon_id_list,name, descript, val = zip(*amino_acid_prediction_dict[key])
    for exon_id in set(exon_id_list):
        print('exon_id:', exon_id)
        #plot_exon_id_conservation(exon_id, bw_30way_phastcon)
        plot_exon_id_conservation(exon_id, bw_data)
    print('')


'''


'''
query_cons_phylop_dict = {val[0]:val[1] for val in query_cons_phylop_list}

score_list = list()
for exon_id in exon_id_list:
    if (exon_id in query_cons_phylop_dict):
        score = query_cons_phylop_dict[exon_id]
        score_list.append(score)
        #plot_exon_id_conservation(exon_id, bw_data)
        #plot_exon_id_conservation(exon_id, bw_30way_phastcon)

plt.figure()
plt.boxplot(score_list)
plt.ylim(-.5,1)
'''




'''
score_list = list()
for exon_id in divisble_3_non_stop_codon_exon_ids:
    if (exon_id not in exon_id_list):
        score = query_cons_phylop_dict[exon_id]
        score_list.append(score)
        #plot_exon_id_conservation(exon_id, bw_data)
        #plot_exon_id_conservation(exon_id, bw_30way_phastcon)


plt.figure()
plt.boxplot(score_list)
plt.ylim(-.5,1)



for exon_id in divisble_3_non_stop_codon_exon_ids[000:100]:
    if aggregate_exon_dict[exon_id]['count'] < 10000:
        continue
    plot_exon_id_conservation(exon_id, bw_30way_phastcon)
    print(exon_id, aggregate_exon_dict[exon_id]['count'])

'''

 

'''
exon_id='chr14:41801583-41801706:+'
plot_exon_id_conservation(exon_id, bw_data)

plot_exon_id_conservation(exon_id, bw_30way_phastcon)
'''








'''
exon_id = 'chr17:79398506-79398701:-'
print(exon_id, aggregate_exon_dict[exon_id]['count'], '%.1f, %.1f' % (aggregate_exon_dict[exon_id]['3ss_score'],aggregate_exon_dict[exon_id]['5ss_score']))
dummy = plot_exon_id_conservation(exon_id, bw_data)

bw_30way_phastcon
bw_data

exon_ids = ['chr3:60615147-60615218:-',
'chr18:25174959-25175100:-',
'chr1:49328124-49328186:-']

exon_ids = ['chr1:171636193-171636388:+']

for exon_id in exon_ids:
    #exon_id = 'chr5:149749284-149749479:+'
    print(exon_id, aggregate_exon_dict[exon_id]['count'], '%.1f, %.1f' % (aggregate_exon_dict[exon_id]['3ss_score'],aggregate_exon_dict[exon_id]['5ss_score']))
    dummy = plot_exon_id_conservation(exon_id, bw_30way_phastcon)
    dummy = plot_exon_id_conservation(exon_id, bw_data)
'''







sorted_cons = sorted(query_cons_dict['Intronic_pc_other_snaptron_removed'], key=lambda x: (el.exon_id_values(x[0]).chrom, el.exon_id_values(x[0]).strand, el.exon_id_values(x[0]).start), reverse=False)

exon_ids,c1,c2 = zip(*sorted_cons)
out_bed = exp_output_path.old_gene + 'et_only_intronic_exons_counts.bed'
el.export_exon_id_list_with_counts_to_bed(exon_ids, aggregate_exon_dict, out_bed)


exon_ids,c1,c2 = zip(*sorted_cons)
out_bed = exp_output_path.old_gene + 'et_only_intronic_exons.bed'
el.export_exon_id_list_to_bed(exon_ids, out_bed)












outdir = exp_output_path.old_gene 
pdf_plots = PdfPages(outdir+'old_gene_categories.pdf')


query_list = ['Other_recovered','snap_only', 'Intronic_pc_other_snaptron_removed']


query_name_dict = {'Other_recovered':'non-GENCODE mRNA','snap_only':'Snaptron-only', 'Intronic_pc_other_snaptron_removed':'Other'}

query_name_list = ['non-GENCODE mRNA','Snaptron-only', 'Other']


#query_name_dict[]


f=open(exp_output_path.out_supplemental+'3D_left.txt','w')
outstring = 'chrom\tstart\tend\tstrand\texon_id\tgenomic_category\texon_counts\n'
f.write(outstring)

fig = plt.figure()
for ii, key in enumerate(query_list):
    exon_id_list=cons_exon_sets_dict[key]
    vals = el.get_exon_dict_counts(exon_id_list ,aggregate_exon_dict)
    
    for exon_id in exon_id_list:
        ex=el.exon_id_values(exon_id)
        
        outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id)
        
        '{:}\t{:}\n'.format(outstring, query_name_list[ii])
        f.write(outstring)
    
    step = 1
    y,x = np.histogram(np.log10(vals), bins=np.arange(0,6+step,step))
    plt.bar(x[:-1]+ii/4-.25,np.log10(y),label=query_name_dict[key], width = 0.2)

f.close()

plt.xlim(1.5,4.5)

plt.legend()
ticks = [2,3,4]
plt.xticks(ticks,['{:,} - {:,}'.format(10**x,10**(x+1)) for x in ticks])
plt.xlabel('Sequencing read counts')
plt.ylabel('Exon count')
ticks = [0,1,2,3,4,5]
plt.yticks([x for x in ticks],['{:,}'.format(10**x) for x in ticks])
plt.ylim(0,5)
plt.tight_layout()
pdf_plots.savefig(fig)





'''
with open(exp_output_path.out_supplemental+'3A.txt','w') as f:
    outstring = 'chrom\tstart\tend\tstrand\texon_id\tgenomic_category\n'
    f.write(outstring)
    
    for ii, val in enumerate(labels):
        exon_id_list = sizes[ii]
        for exon_id in exon_id_list:
            ex=el.exon_id_values(exon_id)
            
            outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id)
            
            '{:}\t{:}\n'.format(outstring, val)
            f.write(outstring)

with open(exp_output_path.out_supplemental+'3D_left.txt','w') as f:
    ex = el.exon_id_values(exon_id)
    f.write('chromosome_position\tconservation_score\n')
    
    for ii, val in enumerate(score_array):
        val
        out_string = "{:}\t{:}\n".format( ex.start+ii-200,val)
    '''









 

fig = plt.figure()
for ii, key in enumerate(query_list):
    exon_id_list=cons_exon_sets_dict[key]
    dummy, vals_1, vals_2 = zip(*query_cons_dict[key])
    step = .2
    y,x = np.histogram(vals_1, bins=np.arange(0,1+step,step))
    plt.bar(x[:-1]+ 0.2*ii/4-.0,np.log10(y),label=query_name_dict[key], width = 0.05)
#plt.yscale('log')
#plt.xlim(1.5,4.5)

plt.legend()
ticks = [0,.2,.4,.6,.8]
plt.xticks(ticks,['{:.1f} - {:.1f}'.format(x,(x+.2)) for x in ticks])
plt.xlabel('Mean phastcons')
plt.ylabel('Exon count')
ticks = [0,1,2,3,4,5]
plt.yticks([x for x in ticks],['{:,}'.format(10**x) for x in ticks])
plt.ylim(0,5)
plt.tight_layout()
pdf_plots.savefig(fig)



'''
fig = plt.figure()
for ii, key in enumerate(query_list):
    exon_id_list=cons_exon_sets_dict[key]
    dummy, vals_1, vals_2 = zip(*query_cons_dict[key])
    step = .2
    y,x = np.histogram(vals_2, bins=np.arange(0,1+step,step))
    plt.bar(x[:-1]+ 0.2*ii/4-.0,np.log10(y),label=query_name_dict[key], width = 0.05)
#plt.yscale('log')
#plt.xlim(1.5,4.5)

plt.legend()
ticks = [0,.2,.4,.6,.8]
plt.xticks(ticks,['{:.1f} - {:.1f}'.format(x,(x+.2)) for x in ticks])
plt.xlabel('Max phastcons')
plt.ylabel('Counts')
ticks = [0,1,2,3,4,5]
plt.yticks([x for x in ticks],['{:,}'.format(10**x) for x in ticks])
plt.ylim(0,5)
plt.tight_layout()
pdf_plots.savefig(fig)
'''




















def plot_exon_id_full_phastcons(exon_id, bw_data, **kwargs):
    
    pad_length = 200
    if 'pad_length' in kwargs:
        pad_length = kwargs['pad_length']
    no_plot = False
    if 'no_plot' in kwargs:
        no_plot = kwargs['no_plot']
    
    
    
    
    ex = el.exon_id_values(exon_id)
    score_array = bw_data.get_base_scores(ex.chrom, ex.start - pad_length, ex.end + pad_length)
    if ex.strand == '-':
        score_array = score_array[::-1]
    
    
    if bw_data.check_if_none(score_array) == False:
        
        if no_plot == True:
            return (score_array)
            
        fig = plt.figure()
        plt.bar(range(-1*pad_length, -1*pad_length + len(score_array)), score_array )
        plt.plot(range(-1*pad_length, -1*pad_length + len(score_array)), score_array )
        plt.title('exon_id: %s' % (exon_id))
        plt.plot([len(score_array)-2*pad_length+1,len(score_array)-2*pad_length+1],[0,1], color='r')
        plt.plot([0,0],[0,1], color='r')
        plt.xlabel("bp from 3'SS")
        plt.ylabel('Phastcons')
        plt.ylim(0,1)
        
        
        if  'output_handle'  in kwargs:
            pdf_plots = kwargs['output_handle']
            pdf_plots.savefig(fig)
        
        return (score_array)
    return [-1]





exon_id='chr21:33852138-33852267:+'
score_array = plot_exon_id_full_phastcons(exon_id, bw_30way_phastcon,output_handle=pdf_plots)

aggregate_exon_dict['chr21:33852138-33852267:+']['count']


exon_id='chr21:33852138-33852267:+'
#score_array = plot_exon_id_full_phastcons(exon_id, bw_30way_phastcon)

with open(exp_output_path.out_supplemental+'3E.txt','w') as f:
    ex = el.exon_id_values(exon_id)
    f.write('chromosome_position\tconservation_score\n')
    
    for ii, val in enumerate(score_array):
        val
        out_string = "{:}\t{:}\n".format( ex.start+ii-200,val)
    




pdf_plots.close()












'''
with open(exp_output_path.out_supplemental+'3A.txt','w') as f:
    outstring = 'chrom\tstart\tend\tstrand\texon_id\tgenomic_category\n'
    f.write(outstring)
    
    for ii, val in enumerate(labels):
        exon_id_list = sizes[ii]
        for exon_id in exon_id_list:
            ex=el.exon_id_values(exon_id)
            
            outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id)
            
            '{:}\t{:}\n'.format(outstring, val)
            f.write(outstring)

'''








'''
count_means = list()
fig = plt.figure()
for ii, key in enumerate(query_sets):
    print(ii)
    exon_id_list = cons_exon_sets_dict[key]
    counts = el.get_exon_dict_counts(exon_id_list ,aggregate_exon_dict)
    count_means.append(np.mean(counts))
    plt.boxplot([counts], positions=[ii],labels=[query_sets_labels[ii]],showfliers=False,whis=[10,90])
#plt.legend(range(len(count_means)),[ "%.1f"%x for x in count_means])
plt.legend([ "%.1f"%x for ii, x in enumerate(count_means)])
plt.yscale('log')
plt.ylim(0,100000)
plt.xticks(rotation=45)
plt.ylabel('sequencing read counts')



score_means = list()
fig = plt.figure()
for ii, key in enumerate(query_sets):
    exon_id_list = cons_exon_sets_dict[key]
    scores = el.get_exon_dict_5ss_scores(exon_id_list ,aggregate_exon_dict)
    score_means.append(np.mean(scores))
    plt.boxplot([scores], positions=[ii],labels=[query_sets_labels[ii]],showfliers=False,whis=[10,90])
plt.xlim(-.5,6.5)
plt.legend([ "%.1f"%x for ii, x in enumerate(score_means)])
plt.ylim(0,13)
plt.xticks(rotation=45)
plt.ylabel('score 5ss')



score_means = list()
fig = plt.figure()
for ii, key in enumerate(query_sets):
    exon_id_list = cons_exon_sets_dict[key]
    scores = el.get_exon_dict_3ss_scores(exon_id_list ,aggregate_exon_dict)
    score_means.append(np.mean(scores))
    print(key,np.mean(scores))
    plt.boxplot([scores], positions=[ii],labels=[query_sets_labels[ii]],showfliers=False,whis=[10,90])
plt.xlim(-.5,6.5)
plt.legend([ "%.1f"%x for ii, x in enumerate(score_means)])
plt.ylim(0,13)
plt.xticks(rotation=45)
plt.ylabel('score 3ss')



length_means = list()
fig = plt.figure()
for ii, key in enumerate(query_sets):
    exon_id_list = cons_exon_sets_dict[key]
    lengths = el.get_exon_dict_lengths(exon_id_list ,aggregate_exon_dict)
    length_means.append(np.mean(lengths))
    plt.boxplot([lengths], positions=[ii],labels=[query_sets_labels[ii]],showfliers=False,whis=[10,90])
plt.xlim(-.5,6.5)
plt.legend([ "%.1f"%x for ii, x in enumerate(length_means)])
plt.ylim(50,250)
plt.xticks(rotation=45)
plt.ylabel('length (bp)')


sample_gc_list = list()
plt.figure()
for ii, key in enumerate(query_sets):
    exon_id_list = cons_exon_sets_dict[key]
    gc_list = list()
    for exon_id in exon_id_list:
        gc = el.get_gc_for_region(exon_id, genome_fasta)
        gc_list.append(gc)
    plt.boxplot([gc_list], positions=[ii],labels=[query_sets_labels[ii]],showfliers=False,whis=[10,90])
    sample_gc_list.append(np.median(gc_list))
    
plt.xlim(-.5,6.5)
plt.legend([ "%.1f"%(100*x) for ii, x in enumerate(sample_gc_list)])
#plt.ylim(50,250)
plt.xticks(rotation=45)
plt.ylabel('GC fraction')    


'''
    
    
        
        
        







