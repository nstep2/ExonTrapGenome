

import numpy as np

from experiment_paths.experiment_paths import *

import exon_id_library.exon_id_lib as el

import gzip
import intervaltree


'''


collapsed_exons_list = ["1_ETF_S1_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"2_ETF_cleaned_S9_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"3_ETF_cleaned_S10_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"4_ETF_cleaned_S11_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"5_ETF_cleaned_S12_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"6_ETF_S2_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"7_ETF_cleaned_S13_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"8_ETF_cleaned_S14_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"9_ETF_S3_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"10_ETF_cleaned_S15_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"11_ETF_S4_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"12_ETF_cleaned_S16_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"13_ETF_cleaned_S17_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"14_ETF_cleaned_S18_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"15_ETF_S5_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"16_ETF_cleaned_S19_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"17_ETF_cleaned_S20_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"18_ETF_S6_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"19_ETF_cleaned_S21_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"20_ETF_cleaned_S7_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"21_ETF_cleaned_S8_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"22_ETF_cleaned_S22_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz",
"23_ETF_cleaned_S23_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz"]



#collapsed_exons_list = [exp_output_path.trimmed_fastq_input_files + x for x in collapsed_exons_list]

'''

#23_ETF_cleaned_S23_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz


input_file_name_list = list()
input_file_name_list.append( [ 1, 1, [ "SRR21497380", "SRR21497381", "SRR21497382", "SRR21497383" ] ] )
input_file_name_list.append( [ 1, 2, [ "SRR21497384", "SRR21497385", "SRR21497386", "SRR21497387" ] ] )
input_file_name_list.append( [ 1, 3, [ "SRR21497388", "SRR21497389", "SRR21497390", "SRR21497391" ] ] )
input_file_name_list.append( [ 1, 4, [ "SRR21497392", "SRR21497393", "SRR21497394", "SRR21497395" ] ] )
input_file_name_list.append( [ 1, 5, [ "SRR21497396", "SRR21497397", "SRR21497398", "SRR21497399" ] ] )
input_file_name_list.append( [ 2, 6, [ "SRR21497400", "SRR21497401", "SRR21497402", "SRR21497403" ] ] )
input_file_name_list.append( [ 2, 7, [ "SRR21497404", "SRR21497405", "SRR21497406", "SRR21497407" ] ] )
input_file_name_list.append( [ 2, 8, [ "SRR21497408", "SRR21497409", "SRR21497410", "SRR21497411" ] ] )
input_file_name_list.append( [ 2, 9, [ "SRR21497412", "SRR21497413", "SRR21497414", "SRR21497415" ] ] )
input_file_name_list.append( [ 3, 10, [ "SRR21497416", "SRR21497417", "SRR21497418", "SRR21497419" ] ] )
input_file_name_list.append( [ 3, 11, [ "SRR21497420", "SRR21497421", "SRR21497422", "SRR21497423" ] ] )
input_file_name_list.append( [ 3, 12, [ "SRR21497424", "SRR21497425", "SRR21497426", "SRR21497427" ] ] )
input_file_name_list.append( [ 3, 13, [ "SRR21497428", "SRR21497429", "SRR21497430", "SRR21497431" ] ] )
input_file_name_list.append( [ 3, 14, [ "SRR21497432", "SRR21497433", "SRR21497434", "SRR21497435" ] ] )
input_file_name_list.append( [ 4, 15, [ "SRR21497436", "SRR21497437", "SRR21497438", "SRR21497439" ] ] )
input_file_name_list.append( [ 4, 16, [ "SRR21497440", "SRR21497441", "SRR21497442", "SRR21497443" ] ] )
input_file_name_list.append( [ 4, 17, [ "SRR21497444", "SRR21497445", "SRR21497446", "SRR21497447" ] ] )
input_file_name_list.append( [ 4, 18, [ "SRR21497448", "SRR21497449", "SRR21497450", "SRR21497451" ] ] )
input_file_name_list.append( [ 5, 19, [ "SRR21497452", "SRR21497453", "SRR21497454", "SRR21497455" ] ] )
input_file_name_list.append( [ 5, 20, [ "SRR21497456", "SRR21497457", "SRR21497458", "SRR21497459" ] ] )
input_file_name_list.append( [ 5, 21, [ "SRR21497460", "SRR21497461", "SRR21497462", "SRR21497463" ] ] )
input_file_name_list.append( [ 5, 22, [ "SRR21497464", "SRR21497465", "SRR21497466", "SRR21497467" ] ] )
input_file_name_list.append( [ 5, 23, [ "SRR21497468", "SRR21497469", "SRR21497470", "SRR21497471" ] ] )


suffix = '_1.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz'
collapsed_exons_list = [x[2][0]+ suffix for ii, x in enumerate(input_file_name_list)]


collapsed_exons_list = [ [input_file_name_list[ii][1],exp_output_path.trimmed_fastq_SAM_files + x] for  ii, x in enumerate(collapsed_exons_list)]
#collapsed_exons_list = [ [input_file_name_list[ii][1],exp_output_path.trimmed_fastq_input_files + x] for  ii, x in enumerate(collapsed_exons_list)]





tmp_dict = {1:'T3', 2:'T4', 3:'T5', 4:'T5', 5:'T5'}
backbone_dict = {x[1]:tmp_dict[x[0]] for x in input_file_name_list}











#exp_output_path.trimmed_fastq_SAM_files
#test = '/mnt/hgfs/main_ssd/et_test/SAM/1_ETF_S1_L001_R1_001.fastq.gz_all_trimmed_simple_paired_selected_sorted_collapsed_reshape_sorted_collapse.txt.gz'

from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
matrix5 = load_matrix5()
matrix3 = load_matrix3()


import os,sys


a=np.zeros(24)


exon_id_dict = dict()

duplicate_count = 0
for ii, pair in enumerate(collapsed_exons_list):
    lib_num, f_path = pair
    print('load file: %s' % (os.path.basename(f_path)))
    #lib_num = int(os.path.basename(f_path).split('_')[0])
    for line in gzip.open(f_path,'rt'):
        #lib_num = os.basename(f_path)
        line_split = line.strip().split('\t')
        chrom  = line_split[1] #line_split[1] is not always == to line_split[0]
        start  = line_split[2]
        end    = line_split[3]
        strand = line_split[4]
        count  = int(line_split[5])
        exon_id = '%s:%s-%s:%s' % (chrom, start, end, strand)
        
        
        if exon_id not in exon_id_dict:
            exon_id_dict[exon_id] = np.zeros(24)
        else:
          
            duplicate_count += 1
        
        exon_id_dict[exon_id][lib_num] += count


def get_3ss_seq(exon_id, genome_fasta):
    ex=el.exon_id_values(exon_id)
    if ex.strand == '+':
        try:
            seq = str(genome_fasta[ex.chrom][ex.start-21:ex.start+2])
            score = maxent.score3(seq, matrix=matrix3)
        except:
            seq='N'
            score=-5500
    if ex.strand == '-':
        try:
            seq = str(genome_fasta[ex.chrom][ex.end-4:ex.end+19].reverse.complement)
            score=maxent.score3(seq, matrix=matrix3)
        except:
            seq='N'
            score=-5500
    return seq,score


def get_5ss_seq(exon_id, genome_fasta):
    ex=el.exon_id_values(exon_id)
    if ex.strand == '+':
        try:
            seq = str(genome_fasta[ex.chrom][ex.end-4:ex.end+5])
            score = maxent.score5(seq, matrix=matrix5)
        except:
            seq='N'
            score=-5500
    if ex.strand == '-':
        try:
            seq = str(genome_fasta[ex.chrom][ex.start-7:ex.start+2].reverse.complement)
            score=maxent.score5(seq, matrix=matrix5)
        except:
            seq='N'
            score=-5500
    return seq,score






new_exon_id_dict=dict()
for ii, exon_id in enumerate(exon_id_dict):
    if ii % 1000000==0:
        print('made {:,} new exon_ids'.format((ii)))
    count_array = exon_id_dict[exon_id]
    new_exon_id_dict[exon_id] = dict()
    
    
    ex=el.exon_id_values(exon_id)
    #exon_id_dict[exon_id][''] = 
    new_exon_id_dict[exon_id]['length'] = ex.length
    new_exon_id_dict[exon_id]['count'] = sum(count_array)
    new_exon_id_dict[exon_id]['strand'] = ex.strand
    new_exon_id_dict[exon_id]['chrom'] = ex.chrom
    new_exon_id_dict[exon_id]['start'] = ex.start
    new_exon_id_dict[exon_id]['end'] = ex.end
    #exon_id_dict[exon_id]['lib_count'] = count_array
    new_exon_id_dict[exon_id]['lib_array'] = count_array
    new_exon_id_dict[exon_id]['3ss'] = get_3ss_seq(exon_id,genome_fasta)[0]
    new_exon_id_dict[exon_id]['3ss_score'] = get_3ss_seq(exon_id,genome_fasta)[1]
    new_exon_id_dict[exon_id]['5ss'] = get_5ss_seq(exon_id,genome_fasta)[0]
    new_exon_id_dict[exon_id]['5ss_score'] = get_5ss_seq(exon_id,genome_fasta)[1]
    
    new_exon_id_dict[exon_id]['lib_count'] = {ii:val for ii, val in enumerate(count_array) if val  > 0}
    




exon_IT = dict()

for exon_id in new_exon_id_dict:
    ex = el.exon_id_values(exon_id)
    if ex.chrom not in exon_IT:
        exon_IT[ex.chrom]={'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
    exon_IT[ex.chrom][ex.strand][ex.start:ex.end]=[exon_id]
    


import pickle
pickle_path = exp_output_path.pickle_merged + "updated_merged_aggregate_exon_dict.pickle"
print('Dumped aggregate pickle path:', pickle_path)
with open(pickle_path, "wb") as output_file:
    pickle.dump(new_exon_id_dict,output_file)
    pickle.dump(exon_IT,output_file)
    pickle.dump({'empty':'empty'},output_file)
    pickle.dump({'date_join_aggregate_chromosome_data':'empty'},output_file)

    output_file.close()


    
    
    
    
    





