

import sys
import intervaltree
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np


import exon_id_library.exon_id_lib as el

from exon_id_library.gff import get_annotated_exons_from_GENCODE_GFF3
from exon_id_library.gff import GENCODE_exon_class

from experiment_paths.experiment_paths import *
#exp_output_path as exp_output_path








def opposite_strand(strand):
    if strand == '+':
        return '-'
    if strand == '-':
        return '+'
    raise('"%s" is not a strand!'%d)
    

def insert_commas(value):
    return f'{value:,}'













region_sets_dict = dict()
exon_sets_dict = dict()

mRNA_exons_dict = dict()
lncRNA_exons_dict = dict()







outdir = exp_output_path.Parse_ENSEMBLE_GENCODE_exons_load
pdf_plots = PdfPages(outdir+'Parse_ENSEMBLE_GENCODE_exons_load_%d_reads.pdf' % (exon_count_build))



#This fasta file was downloaded from UCSC for the exons annotated by GENCODE
fasta_comprehensive_file = '/mnt/0862DDC375402D9E/downloaded_data/Exons/GENCODE_V33_hg38/All_GENCODE_V33_comprehensive_exon_sequence.txt'
#fasta_basic_file         = '/mnt/0862DDC375402D9E/downloaded_data/Exons/GENCODE_V33_hg38/All_GENCODE_V33_basic_exon_sequence.txt'


try:
    genome_fasta['chr1'][1000:2000]
except:
    import pyfaidx
    input_fasta = '/home/pineapple/work/indexes/downloaded/Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
    genome_fasta = pyfaidx.Fasta(input_fasta)




chromsomes_list = ['chrX',  'chrY',  'chrMT', 'chrM']
for i in range(1,23):
    chromsomes_list.append( "chr%d" % (i) )

    







GENCODE_GFF3 = exp_output_path.GENCODE_GFF3


exon_GENCODE_comp, transcript_id_count_comp, dup_exon_dict, UTR_annotation, exon_id_sets_dict, GENCODE_exon_dict_comp, exon_id_transcript_sets_dict_comp, transcript_to_gene_dict_comp, GENCODE_gene_dict_comp          = get_annotated_exons_from_GENCODE_GFF3(GENCODE_GFF3, 'comprehensive', genome_fasta)


exon_GENCODE_basic, transcript_id_count_basic, dup_exon_dict, UTR_annotation, exon_id_sets_dict_basic, GENCODE_exon_dict_basic, exon_id_transcript_sets_dict, transcript_to_gene_dict, dummy  = get_annotated_exons_from_GENCODE_GFF3(GENCODE_GFF3, 'basic', genome_fasta)

GENCODE_gene_dict=GENCODE_gene_dict_comp






mRNA_exons_dict['mRNA'] = exon_id_sets_dict['protein_coding']
mRNA_exons_dict['mRNA_basic'] = exon_id_sets_dict_basic['protein_coding']




#GENCODE_exon_dict_basic
def get_transcript_to_exon_set(GENCODE_exon_dict):
    tid_dict = dict()
    #transcript_id_to_tid = dict()
    for key in GENCODE_exon_dict:
        entry = GENCODE_exon_dict[key]
        
        tid = entry.transcript_id
        

        if tid in tid_dict:
            tid_dict[tid].append(key)
        else:
            tid_dict[tid]=[key]

    return tid_dict




def get_list_first_last_exons(tid_dict,GENCODE_exon_dict):
    first_last_exon_id_pair_list = list()
    middle_exon_id_list = list()

    for tid in tid_dict:
        t_exon_list = tid_dict[tid]
        t_exon_list = sorted(t_exon_list,key=lambda x:int(GENCODE_exon_dict[x].exon_number))

        first = t_exon_list[0]
        last  = t_exon_list[-1]

        first_exon_id = GENCODE_exon_dict[first].exon_id
        last_exon_id  = GENCODE_exon_dict[last].exon_id
        for ii in range(1,len(t_exon_list)-1):
            middle_exon_id_list.append(GENCODE_exon_dict[t_exon_list[ii]].exon_id)

        first_last_exon_id_pair_list.append([first_exon_id,last_exon_id])

    return first_last_exon_id_pair_list,middle_exon_id_list





#get list of transcript exon ids in
tid_dict = get_transcript_to_exon_set(GENCODE_exon_dict_basic)
tid_comp_dict = get_transcript_to_exon_set(GENCODE_exon_dict_comp)

#make a list of first/last exons. also have that list as paired [first,last] in case I want to get first 5'ss and last 3'ss

def get_first_last_middle_exon_ids(tid_dict,GENCODE_exon_dict_basic, mRNA_exons_dict):
    
    first_last_exon_id_pair_list,middle_exon_id_list=get_list_first_last_exons(tid_dict,GENCODE_exon_dict_basic)
    
    single_exon_genes_list = [x for x in first_last_exon_id_pair_list if x[0] == x[1]] #single exon genes
    
    first_last_exon_id_pair_list = [x for x in first_last_exon_id_pair_list if x[0] != x[1]] #remove single exon genes
        
    a,b=zip(*first_last_exon_id_pair_list)
    
    first_last_exon_ids=list(set(a).union(set(b)))
    first_last_exon_ids=set(first_last_exon_ids).difference(set(middle_exon_id_list))
    
    pc_first_last_exon_ids = list(set(first_last_exon_ids).intersection(set(mRNA_exons_dict['mRNA'])))
    pc_middle_exon_id_list = list(set(middle_exon_id_list).intersection(set(mRNA_exons_dict['mRNA'])))
    
    pc_first_last_exon_ids=list(set(pc_first_last_exon_ids).difference(set(pc_middle_exon_id_list)))  #this removes first/last exons that are also annotated as middle exons in other transcripts. 
    
    return first_last_exon_ids, middle_exon_id_list, pc_first_last_exon_ids, pc_middle_exon_id_list


first_last_exon_ids, middle_exon_id_list, pc_first_last_exon_ids, pc_middle_exon_id_list = get_first_last_middle_exon_ids(tid_dict,GENCODE_exon_dict_basic, mRNA_exons_dict)






ET_exon_ids = list(aggregate_exon_dict.keys())


recovered_annotated_middle_exon_ids = el.exon_id_intersection(ET_exon_ids, pc_middle_exon_id_list)








lncRNA_exon_id_list = el.get_lncRNA_exon_ids(GENCODE_exon_dict_basic)
lncRNA_exon_id_list = list(set(lncRNA_exon_id_list))
















lncRNA_query_IT_result_dict = el.query_IT_with_exon_id_list(lncRNA_exon_id_list, aggregate_exon_IT)

share_5ss_exon_id_list,share_3ss_exon_id_list, recovering_5ss_exon_id_list, recovering_3ss_exon_id_list, exact_match_to_alternat_dict, alternate_to_exact_mach_dict = el.get_exon_ids_that_share_5ss_3ss_in_ET_data(lncRNA_query_IT_result_dict)



recovering_3ss_exon_id_list = el.exon_id_difference(recovering_3ss_exon_id_list, lncRNA_exon_id_list)

recovering_5ss_exon_id_list = el.exon_id_difference(recovering_5ss_exon_id_list, lncRNA_exon_id_list)

counts_3ss_list=el.get_exon_dict_counts(recovering_3ss_exon_id_list, aggregate_exon_dict)

counts_5ss_list=el.get_exon_dict_counts(recovering_5ss_exon_id_list, aggregate_exon_dict)









pc_first_last_query_IT_result_dict = el.query_IT_with_exon_id_list(pc_first_last_exon_ids, aggregate_exon_IT)



#pc_first_last_query_IT_result_dict
share_5ss_exon_id_list,share_3ss_exon_id_list, recovering_5ss_exon_id_list, recovering_3ss_exon_id_list, exact_match_to_alternat_dict, alternate_to_exact_mach_dict = el.get_exon_ids_that_share_5ss_3ss_in_ET_data(pc_first_last_query_IT_result_dict)


pc_middle_query_IT_result_dict = el.query_IT_with_exon_id_list(pc_middle_exon_id_list, aggregate_exon_IT)

share_5ss_exon_id_list,share_3ss_exon_id_list, recovering_5ss_exon_id_list, recovering_3ss_exon_id_list, exact_match_to_alternat_dict, alternate_to_exact_mach_dict = el.get_exon_ids_that_share_5ss_3ss_in_ET_data(pc_middle_query_IT_result_dict)



recovering_3ss_exon_id_list=el.exon_id_intersection(recovering_3ss_exon_id_list,pc_middle_exon_id_list)
recovering_5ss_exon_id_list=el.exon_id_intersection(recovering_5ss_exon_id_list,pc_middle_exon_id_list)
both_3ss_5ss_recovering_exon_id_list=el.exon_id_intersection(recovering_3ss_exon_id_list,recovering_5ss_exon_id_list)


list_exon_id_captured_by_dual_exons, list_dual_exon_id, dual_exon_id_to_pair_individual_dict = el.get_individual_exons_from_dual(exact_match_to_alternat_dict, alternate_to_exact_mach_dict)

list_dual_individual_exon_ids_not_in_ET = el.exon_id_intersection(list_exon_id_captured_by_dual_exons,ET_exon_ids) 









exon_5ss_id_dict, exon_3ss_id_dict = el.build_single_ss_dict_from_exon_id_list(set(pc_middle_exon_id_list).union(pc_first_last_exon_ids))


'''
This section seems unused now.
'''



overlapping_found_middle_exon_id_list = list()
for exon_id in pc_middle_query_IT_result_dict['found']:
    for exon_id_2 in pc_middle_query_IT_result_dict['found'][exon_id]:
        overlapping_found_middle_exon_id_list.append(exon_id_2)
overlapping_found_middle_exon_id_list=list(set(overlapping_found_middle_exon_id_list))

count_annotated_5ss = 0
count_annotated_3ss = 0
count_unannotated_5ss = 0
count_unannotated_3ss = 0
count_both_an = 0
count_5ss_an_3ss_un_an = 0
count_3ss_an_5ss_un_an = 0
count_neither_an =0
both_an_exon_id_list = list()
an_5ss_un_an_3ss_exon_id_list = list()
an_3ss_un_an_5ss_exon_id_list = list()
neither_an_exon_id_list = list()
for exon_id in overlapping_found_middle_exon_id_list:
    ex = el.exon_id_values(exon_id)
    if ex.id_5ss in exon_5ss_id_dict:
        count_annotated_5ss += 1
    else:
        count_unannotated_5ss += 1

    if ex.id_3ss in exon_3ss_id_dict:
        count_annotated_3ss += 1
    else:
        count_unannotated_3ss += 1

    if ex.id_3ss in exon_3ss_id_dict and ex.id_5ss in exon_5ss_id_dict:
        count_both_an += 1
        both_an_exon_id_list.append(exon_id)
    if ex.id_3ss in exon_3ss_id_dict and ex.id_5ss not in exon_5ss_id_dict:
        count_3ss_an_5ss_un_an += 1
        an_3ss_un_an_5ss_exon_id_list.append(exon_id)
    if ex.id_3ss not in exon_3ss_id_dict and ex.id_5ss in exon_5ss_id_dict:
        count_5ss_an_3ss_un_an += 1
        an_5ss_un_an_3ss_exon_id_list.append(exon_id)
    if ex.id_3ss not in exon_3ss_id_dict and ex.id_5ss not in exon_5ss_id_dict:
        count_neither_an += 1
        neither_an_exon_id_list.append(exon_id)















lncRNA_first_last_exon_ids = list(set(first_last_exon_ids).intersection(set(lncRNA_exon_id_list)))
lncRNA_middle_exon_ids = list(set(middle_exon_id_list).intersection(set(lncRNA_exon_id_list)))











#exon_id_list = pc_middle_exon_id_list


print("exon_assignment_class(pc_middle_exon_id_list)")
pc_middle_exon_assignment = el.exon_assignment_class(pc_middle_exon_id_list,aggregate_exon_IT,aggregate_exon_dict)

print("exon_assignment_class(pc_middle_exon_id_list - 50,500)")
pc_middle_exon_assignment = el.exon_assignment_class(el.size_exon_id_list(pc_middle_exon_id_list,50,500),aggregate_exon_IT,aggregate_exon_dict)







#pc_middle_exon_id_list
exon_id_list = pc_middle_exon_id_list
highly_overlapping, exact_list, has_overlapping=el.get_highly_overlapping_non_exact_exon_dict(exon_id_list, aggregate_exon_IT, aggregate_exon_dict)
best_overlapping_annotated_not_exact = highly_overlapping


pc_highly_overlapping=highly_overlapping
pc_exact_list=exact_list






print('len(primary_3ss_exon_id_set)', len(primary_3ss_exon_id_set))



















unique_seq_pc_middle_exon_id_list, duplicate_seq_pc_middle_exon_id_list = el.scan_seq_to_exon_id_list_dict(pc_middle_exon_id_list,aggregate_exon_dict, genome_fasta)


deduplicate_pc_middle_exon_id_list=set(el.size_exon_id_list(pc_middle_exon_id_list,50,500)).difference(duplicate_seq_pc_middle_exon_id_list)












lncRNA_exon_tid_list = el.get_lncRNA_transcript_ids(GENCODE_exon_dict_basic)

pc_exon_tid_list = el.get_protein_coding_transcript_ids(GENCODE_exon_dict_basic)




gencode_exon_IT = dict()
for chrom in aggregate_exon_IT:
    gencode_exon_IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
    

        
gencode_transcript_IT = dict()
gencode_pc_transcript_IT = dict()
gencode_lncRNA_transcript_IT = dict()
for chrom in aggregate_exon_IT:
    gencode_transcript_IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
    gencode_pc_transcript_IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
    gencode_lncRNA_transcript_IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}


##############This except block needs fixing!!!!!

problem_t_id_list = list()

t_id = 'ENST00000359512.8'

for t_id in tid_dict:
    try:
        for ex_id in tid_dict[t_id]:
            ex = GENCODE_exon_dict_comp[ex_id]
            gencode_exon_IT[ex.chrom][ex.strand][ex.start:ex.end] = ex_id
        
        t_exons = tid_dict[t_id]
        first = GENCODE_exon_dict_comp[t_exons[0]]
        last  = GENCODE_exon_dict_comp[t_exons[-1]]
        if first.strand == '+':
            start_coord = first.start
            end_coord   = last.end
        else:
            start_coord = last.start
            end_coord   = first.end
        
        gencode_transcript_IT[first.chrom][first.strand][start_coord:end_coord] = t_exons
        
        if t_id in lncRNA_exon_tid_list:
            gencode_lncRNA_transcript_IT[first.chrom][first.strand][start_coord:end_coord] = t_exons
        
        if t_id in pc_exon_tid_list:
            gencode_pc_transcript_IT[first.chrom][first.strand][start_coord:end_coord] = t_exons
        
        
        
    except:
        1






gencode_pc_exon_IT = dict()
gencode_lncRNA_exon_IT = dict()
for chrom in aggregate_exon_IT:
    gencode_pc_exon_IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
    gencode_lncRNA_exon_IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}



for exon_id in set(pc_middle_exon_id_list).union(pc_first_last_exon_ids):
    ex=el.exon_id_values(exon_id)
    chrom  = ex.chrom
    strand = ex.strand
    start  = ex.start
    end    = ex.end
    
    gencode_pc_exon_IT[chrom][strand][start:end]=exon_id

for exon_id in set(lncRNA_middle_exon_ids).union(lncRNA_first_last_exon_ids):
    ex=el.exon_id_values(exon_id)
    
    
    gencode_lncRNA_exon_IT[ex.chrom][ex.strand][ex.start:ex.end]=exon_id

gencode_pc_overlapping_exon_ids = list()
gencode_lncRNA_overlapping_exon_ids = list()
for exon_id in aggregate_exon_dict:
    ex=el.exon_id_values(exon_id)
    intervals = gencode_lncRNA_exon_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) > 0:
        gencode_lncRNA_overlapping_exon_ids.append(exon_id)
    
    intervals = gencode_pc_exon_IT[chrom][strand][start:end]
    if len(intervals) > 0:
        gencode_pc_overlapping_exon_ids.append(exon_id)
    
    
gencode_lncRNA_overlapping_exon_ids_set=set(gencode_lncRNA_overlapping_exon_ids)

gencode_pc_overlapping_exon_ids_set=set(gencode_pc_overlapping_exon_ids)









exon_id_set = set(pc_middle_exon_id_list)
exon_id_set = exon_id_set.union(pc_first_last_exon_ids)

pc_overlapping_trapped_exons_set = list()
for exon_id in aggregate_exon_dict:
    exon = aggregate_exon_dict[exon_id]
    ex = el.exon_id_values(exon_id)
    intervals = gencode_transcript_IT[ex.chrom][ex.strand][ex.start:ex.end]
    for interval in intervals:        
        pc_overlapping_trapped_exons_set += interval[2]
        #ideally, each exon in interval[2] would be checked if it overlaps the current exon_id
pc_overlapping_trapped_exons_set=set(pc_overlapping_trapped_exons_set)



    
tid_comp_dict = get_transcript_to_exon_set(GENCODE_exon_dict_comp)

#make a list of first/last exons. also have that list as paired [first,last] in case I want to get first 5'ss and last 3'ss
first_last_exon_id_pair_list,middle_exon_id_list=get_list_first_last_exons(tid_comp_dict,GENCODE_exon_dict_comp)



#all possible exons and other transcripts
exon_id_set = exon_id_transcript_sets_dict_comp['lncRNA']
exon_id_set = exon_id_set.union(exon_id_transcript_sets_dict_comp['processed_transcript'])
exon_id_set = exon_id_set.union(exon_id_transcript_sets_dict_comp['protein_coding'])
exon_id_set = exon_id_set.union(exon_id_transcript_sets_dict_comp['transcribed_unprocessed_pseudogene'])
exon_id_set = exon_id_set.union(exon_id_transcript_sets_dict_comp['unprocessed_pseudogene'])
exon_id_set = exon_id_set.union(exon_id_transcript_sets_dict_comp['processed_pseudogene'])
exon_id_set = exon_id_set.union(exon_id_transcript_sets_dict_comp['transcribed_processed_pseudogene'])
exon_id_set = exon_id_set.union(exon_id_transcript_sets_dict_comp['pseudogene'])
exon_id_set = exon_id_set.union(exon_id_transcript_sets_dict_comp['nonsense_mediated_decay'])



stranded_gencode_overlapping_trapped_exons_set   = set()
stranded_gencode_overlapping_trapped_exons_list = list()
unstranded_gencode_overlapping_trapped_exons_set   = set()
unstranded_gencode_overlapping_trapped_exons_list = list()

for exon_id in exon_id_set:
    ex = el.exon_id_values(exon_id)
    if ex.chrom not in aggregate_exon_IT:
        continue
    intervals = aggregate_exon_IT[ex.chrom][ex.strand][ex.start:ex.end]
    for interval in intervals:
        stranded_gencode_overlapping_trapped_exons_list += interval[2]
        unstranded_gencode_overlapping_trapped_exons_list += interval[2]
    
    
    
    intervals = aggregate_exon_IT[ex.chrom][opposite_strand(ex.strand)][ex.start:ex.end]
    for interval in intervals:
        unstranded_gencode_overlapping_trapped_exons_list += interval[2]
    
    
        
        
  
stranded_gencode_overlapping_trapped_exons_set=set(stranded_gencode_overlapping_trapped_exons_list)


 
unstranded_gencode_overlapping_trapped_exons_set=set(unstranded_gencode_overlapping_trapped_exons_list)


#
#for key in exon_id_transcript_sets_dict:
#    print(key, len(exon_id_transcript_sets_dict[key]))













exon_id_list=pc_middle_exon_id_list
unique_exon_id_list, duplicate_seq_pc_middle_exon_id_list = el.scan_seq_to_exon_id_list_dict(exon_id_list,aggregate_exon_dict, genome_fasta)

pc_middle_length_array = np.zeros(501)
for exon_id in unique_exon_id_list:
    ex = el.exon_id_values(exon_id)
    if ex.length <= 500:
        pc_middle_length_array[ex.length] += 1
    

recovered_exon_plus_fuzzy_list = [x[0] for x in pc_highly_overlapping.values()] + pc_exact_list
recovered_exon_plus_fuzzy_list = set(recovered_exon_plus_fuzzy_list).difference(duplicate_seq_pc_middle_exon_id_list)

ET_pc_middle_length_array = np.zeros(501)
for exon_id in el.exon_id_intersection(unique_exon_id_list, recovered_exon_plus_fuzzy_list):
    ex = el.exon_id_values(exon_id)
    if ex.length <= 500:
        ET_pc_middle_length_array[ex.length] += 1

#https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-python-numpy-scipy
def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w



#def plot_len_recovery(pc_middle_length_list,ET_pc_middle_length_list,all_exon_lengths):
 
def plot_len_recovery(pc_middle_length_list,ET_pc_middle_length_list, **kwargs):
    
    
    fig, ax1 =plt.subplots()
    
    plt.title('recovered annotated exons vs all annotated exons by size\nmoving average window = 9')
    
    if 'gc_split' not in kwargs:
        plt.plot(range(4,len(ET_pc_middle_length_list)-4),moving_average(ET_pc_middle_length_list/pc_middle_length_list, 9),color='blue')
    
    if 'gc_split' in kwargs:
        for gc in kwargs['gc_split']:
            percent, et_lengths, an_lengths = gc
            plt.plot(range(4,len(et_lengths)-4),moving_average(et_lengths/an_lengths, 9)) #,color='blue')
    
    
    
    plt.plot([0,0],[0,0],color='red')
    plt.xlabel('length (bp)')
    plt.ylabel('proportion recovered by exon trapping')
    plt.xlim(0,500)
    plt.ylim(0,1)
    plt.legend(['Exon trapping recovered','all known interal mRNA exons'])
    ax1_twin=ax1.twinx()
    
    
    frac_an=moving_average(pc_middle_length_list,9)/max(moving_average(pc_middle_length_list,9))
    max_frac_an = max(frac_an/sum(frac_an))
    
    max_ex_count = max(moving_average(pc_middle_length_list,9)) 
    
    top = int(max_ex_count)
    mid = top/2
    
    if 'gc_split' not in kwargs:
    
        ax1_twin.plot(range(4,len(ET_pc_middle_length_list)-4),frac_an/sum(frac_an),color='red')
    
    
    if 'gc_split' in kwargs:
        for gc in kwargs['gc_split']:
            percent, et_lengths, an_lengths = gc
            sum_gc_split_an_lengths = np.zeros(len(an_lengths))
            
        for gc in kwargs['gc_split']:
            percent, et_lengths, an_lengths = gc
            sum_gc_split_an_lengths += an_lengths
            
        gc_frac_an_sum = moving_average(sum_gc_split_an_lengths,9)/(max(moving_average(sum_gc_split_an_lengths,9) ) )
        
        for gc in kwargs['gc_split']:
            percent, et_lengths, an_lengths = gc
            
            gc_frac_an=moving_average(an_lengths,9)/(max(moving_average(an_lengths,9) ) )
            
            
            
            proportion_an_gc = sum(an_lengths)/sum(sum_gc_split_an_lengths)
            ax1_twin.plot(range(4,len(an_lengths)-4),gc_frac_an/sum(gc_frac_an) *proportion_an_gc   )   #,color='red')
    
    
    
    plt.ylim(0,max_frac_an)
    #plt.yticks([0,round(.5*max_frac_an*100,1),  round(max_frac_an*100,1)], ['0%%','%d%%' % ( round(max_frac_an/2*100,1)*max_ex_count ), '%d%%' % max_ex_count*(round(max_frac_an*100,1))])
    plt.yticks([0,.5*max_frac_an, max_frac_an], ['0%%','%d' % ( mid ), '%d' % (top)])
    plt.ylabel('Exon count')
    
    
    plt.tight_layout()
    pdf_plots.savefig(fig)


'''
plot_len_recovery(pc_middle_length_array,ET_pc_middle_length_array, gc_split=gc_extra_list)

'''



#all_exon_lengths = el.get_exon_dict_lengths(primary_3ss_exon_id_set,aggregate_exon_dict)
#primary_3ss_exon_id_set


all_middle_length_array = np.zeros(501)
for exon_id in  primary_3ss_exon_id_set:
    ex = el.exon_id_values(exon_id)
    if ex.length <= 500:
        all_middle_length_array[ex.length] += 1
     
        
     
plot_len_recovery(pc_middle_length_array,ET_pc_middle_length_array)






gc_an_results = el.get_exon_id_list_binned_GC(el.size_exon_id_list(unique_exon_id_list,50,500), [40,75], genome_fasta)

gc_et_results = el.get_exon_id_list_binned_GC(el.size_exon_id_list(el.exon_id_intersection(unique_exon_id_list, recovered_exon_plus_fuzzy_list),50,500), [40,75], genome_fasta)


gc_extra_list = list()
for ii, key in enumerate(gc_an_results):
    
    
    pc_length_array = np.zeros(501)
    for exon_id in gc_an_results[key]:
        ex = el.exon_id_values(exon_id)
        if ex.length <= 500:
            pc_length_array[ex.length] += 1
    
    et_length_array = np.zeros(501)
    for exon_id in gc_et_results[key]:
        ex = el.exon_id_values(exon_id)
        if ex.length <= 500:
            et_length_array[ex.length] += 1
            
    gc_extra_list.append([key, et_length_array, pc_length_array])



#plot_len_recovery(pc_middle_length_array,ET_pc_middle_length_array, gc_split=gc_extra_list)










x = ET_pc_middle_length_array/pc_middle_length_array
max(x[~np.isnan(x)])




#<> export exon ids, exon length array

with open(exp_output_path.out_supplemental+'2C_all_mRNA.txt','w') as f:
    f.write('pc_exon_id\tlength\n')
    for exon_id in unique_exon_id_list:
        ex = el.exon_id_values(exon_id)
        if ex.length <= 500:
            f.write(exon_id + '\t')
            f.write('{:}\t{:}\n'.format(ex.length,'annotated_mRNA'))
            #pc_middle_length_array[ex.length] += 1
    
with open(exp_output_path.out_supplemental+'2C_ET_mRNA.txt','w') as f:
    f.write('pc_exon_id\tlength\n')
    for exon_id in el.exon_id_intersection(unique_exon_id_list, recovered_exon_plus_fuzzy_list):
        ex = el.exon_id_values(exon_id)
        if ex.length <= 500:
            f.write(exon_id + '\t')
            f.write('{:}\t{:}\n'.format(ex.length,'ET_mRNA'))
            #pc_middle_length_array[ex.length] += 1







pdf_plots.close()




exon_RefSeq, transcript_id_count_RefSeq = 'dummy','dummy'

protein_coding_ids = 'dummy'
protein_coding_ids_basic = 'dummy'

pc_exon_tid_list = 'dummy'

lncRNA_exon_tid_list = 'dummy'




import pickle
pickle_path = exp_output_path.pickle_merged + "Parse_ENSEMBLE_GENCODE_exons_load__%d.pickle" % (exon_count_build)

with open(pickle_path, "wb") as output_file:
    pickle.dump(tid_dict, output_file)
    pickle.dump(tid_comp_dict, output_file)
    pickle.dump(first_last_exon_id_pair_list, output_file)
    
    pickle.dump(middle_exon_id_list, output_file)
    
    pickle.dump(pc_first_last_exon_ids, output_file)
    pickle.dump(pc_middle_exon_id_list, output_file)
        
    pickle.dump(exon_GENCODE_comp, output_file)
    pickle.dump(transcript_id_count_comp, output_file)
    pickle.dump(dup_exon_dict, output_file)
    pickle.dump(UTR_annotation, output_file)
    pickle.dump(exon_id_sets_dict, output_file)
    pickle.dump(GENCODE_exon_dict_comp, output_file)
    pickle.dump(exon_id_transcript_sets_dict_comp   , output_file)       
    pickle.dump(exon_GENCODE_basic, output_file)
    pickle.dump(transcript_id_count_basic, output_file)
    pickle.dump(exon_id_sets_dict_basic, output_file)
    pickle.dump(GENCODE_exon_dict_basic, output_file)
    pickle.dump(exon_id_transcript_sets_dict , output_file) 
    pickle.dump(exon_RefSeq, output_file)
    pickle.dump(transcript_id_count_RefSeq   , output_file)        
    pickle.dump(protein_coding_ids , output_file)
    pickle.dump(protein_coding_ids_basic , output_file)
    pickle.dump(lncRNA_middle_exon_ids , output_file)
    pickle.dump(lncRNA_first_last_exon_ids , output_file)
    pickle.dump(lncRNA_exon_id_list, output_file)
    pickle.dump(pc_overlapping_trapped_exons_set, output_file)
    pickle.dump(gencode_exon_IT, output_file)
    pickle.dump(gencode_pc_exon_IT, output_file)
    pickle.dump(gencode_lncRNA_exon_IT, output_file)  
    pickle.dump(gencode_transcript_IT, output_file)  
    pickle.dump(gencode_lncRNA_overlapping_exon_ids_set, output_file) 
    pickle.dump(gencode_lncRNA_overlapping_exon_ids_set, output_file) 
    pickle.dump(gencode_pc_overlapping_exon_ids_set, output_file) 
    pickle.dump(pc_exon_tid_list, output_file) 
    pickle.dump(lncRNA_exon_tid_list, output_file) 
    
    pickle.dump(gencode_lncRNA_transcript_IT, output_file) 
    pickle.dump(gencode_pc_transcript_IT, output_file) 
    pickle.dump(stranded_gencode_overlapping_trapped_exons_set, output_file) 
    
    pickle.dump(GENCODE_gene_dict, output_file) 
    pickle.dump(pc_highly_overlapping, output_file)
    
    
#transcript_to_gene_dict


el.bed_line_to_exon_id
el.exon_id_to_bed_line









#first_last_exon_ids, middle_exon_id_list, pc_first_last_exon_ids, pc_middle_exon_id_list



exon_id_list = pc_first_last_exon_ids
exon_id_list_name = 'mRNA_first_last_exon_id_list'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)


exon_id_list = middle_exon_id_list
exon_id_list_name = 'middle_exon_id_list'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)



exon_id_list = first_last_exon_ids
exon_id_list_name = 'first_last_exon_ids'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)



exon_id_list = lncRNA_exon_id_list
exon_id_list_name = 'lncRNA_exon_id_list'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)













exon_id_list = pc_middle_exon_id_list
exon_id_list_name = 'mRNA_middle_exon_id_list'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)




exon_id_list = lncRNA_middle_exon_ids
exon_id_list_name = 'lncRNA_middle_exon_id_list'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)



exon_id_list = aggregate_exon_dict.keys()
exon_id_list_name = 'ET_exon_ids'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_with_counts_to_bed(exon_id_list, aggregate_exon_dict, out_file_path)













