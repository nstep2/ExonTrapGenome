


from experiment_paths.experiment_paths import *




import pyfaidx
import time
import glob
import gzip, os, sys







from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
matrix5 = load_matrix5()
matrix3 = load_matrix3()


import Bio

import matplotlib.pyplot as plt
import intervaltree

import numpy as np
import re
import sys


print(sys.argv)

'''        Specify which library to align             '''
if len(sys.argv) > 1:
    #fasta_file_location = experiment_location()
    R1_filename = sys.argv[1]
    R2_filename = sys.argv[2]
    library_number = int(sys.argv[3])
    sample_id = sys.argv[4]
    
    if sample_id == 'T3':
        darkcycle_offset = 4
    if sample_id == 'T4':
        darkcycle_offset = 3
    if sample_id == 'T5':
        darkcycle_offset = 3
    
    R1_filename = sys.argv[1]
    library_number = int(sys.argv[3])
    sam_file_dict = dict()
    
    
    sam_file_dict[library_number]= sys.argv[1]

    sam_file_experiment_name = exp_output_path.trimmed_fastq_input_files









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




tmp_dict = {1:'T3', 2:'T4', 3:'T5', 4:'T5', 5:'T5'}
backbone_dict = {x[1]:tmp_dict[x[0]] for x in input_file_name_list}

SAM_file_list = ["{:}_1.fastq.gz_all_trimmed_simple_paired.sam.gz".format(x[2][0] ) for x in input_file_name_list ]

















def process_sam_for_exon_intervals_simple_cigar_string(sam_file_name, short_id,read_depth_complexity_check_thresholds):
    
    selected_chromosome = 'chr17'
    chr_length = 83257441
    
    paired_sam_1 = os.path.basename(sam_file_name)
    paired_sam   = sam_file_name
    
    
    synthetic_dir = exp_output_path.trimmed_fastq_SAM_files_synthetic
    
    print('BEGIN: {:}'.format(paired_sam_1))
    
    new_scrable_sam_name   =  synthetic_dir  +  paired_sam_1[:-7] + '_chr17_scrable.sam.gz'
    new_scrable_sam_name_handle = gzip.open( new_scrable_sam_name , 'wt') 
    
    new_sam_name   = synthetic_dir  +  paired_sam_1[:-7] + '_chr17.sam.gz'
    new_sam_name_handle = gzip.open( new_sam_name , 'wt') 

    
    list_5ss_scores = list()
    list_3ss_scores = list()
    
    
    unique_exons_IT = dict() # 
    stranded_unique_exons_IT = dict()
    stranded_exon_collections_IT = dict()
    
    chrom_OK_set = ['chr%d'%(x) for x in range(1,23)]
    chrom_OK_set.append('chrX')
    chrom_OK_set.append('chrY')
    chrom_OK_set.append('chrM')
    
    exon_frame_length_counts = [0,0,0]
    
    count_large_alignments = 0
    paired_fragment_length_estimates_list = list()
    unique_sites_dict = dict()
    stranded_unique_sites_dict = dict()
    count_non_overlapping_exons = 0
    unique_paired_fragment_length_estimates_list = list()
    count_non_ag_termini = 0
    count_complicated_cigar_string = 0
    
    count_leading_soft_clipped_reads = 0
    count_lagging_soft_clipped_reads = 0
    
    read_check_threshold_unique_exons_dict = dict()  #for collecting the found exon intervals at specific read depths
    
    count_annotated = 0
    
    library_stats_dict = dict()  
    
    
    
    read_position = 0
    total_processed_reads = 0
    with gzip.open(paired_sam, 'rt') as f:
        for ii, line_1 in enumerate(f):
            
            total_processed_reads +=1
            read_position += 1
            
            if ii % 100000000 == 0:
                print('Processed %d lines' % (ii))
            if ii == 500000000000000:
                break
                
            
            

            
            if line_1[0] == '@':
                new_scrable_sam_name_handle.write(line_1)
                new_sam_name_handle.write(line_1)
                continue
            
            
            split_line_1 = line_1.split('\t')
            split_line = split_line_1
            
            
            
            if int(split_line[4]) < 20:   #check read alginment is acceptable quality score
                continue
            
            
            if line_1.find("NH:i:1") >=0:
                1
            else:
                #print('skipping multi-aligning read')
                #int('b')
                continue
            
            
            
            line_2 = next(f)
            read_position += 1
            split_line_2 = line_2.split('\t')        
            
            #this field is present if the read has multiple alignments. Skip these reads.
            if line_1.find('ZS:') >= 0:
                continue
            
            
            
            #check if same chromosome. If not, don't bother with the alignment
            if split_line_1[2] != split_line_2[2]:
                continue
            
            #we are only processing a single chromosome
            if split_line_1[2] != selected_chromosome:
                continue
            
            
            #if the alignment quality score is less than 20, then skip processing it
            if int(split_line_1[4]) < 20:
                 continue
            
            chrom = split_line[2]
            if chrom not in chrom_OK_set:
                continue
            
            
            
            
            
            new_sam_name_handle.write(line_1)
            new_sam_name_handle.write(line_2)
            
            #dist between = [8]
            if int(split_line_1[3]) < int(split_line_2[3]): # (forward strand)
                scrable_strand = '+'
            else:
                scrable_strand = '-'
            
            exon_len = abs(int(split_line_1[3]) - int(split_line_2[3]) )
            
            #exon_len = np.random()
            new_left = np.random.randint(1000,chr_length-1000)
            if scrable_strand == '+':
                new_right = new_left+exon_len
            else:
                new_right = new_left-exon_len
            
            scramble_split_line_1 = line_1.split('\t')
            scramble_split_line_2 = line_2.split('\t')
            
            scramble_split_line_1[3] = str(new_left)
            scramble_split_line_2[3] = str(new_right)
            
            scramble_split_line_1[7] = str(new_right)
            scramble_split_line_2[7] = str(new_left)
            
            
            scrable_line_1 = '\t'.join(scramble_split_line_1)
            new_scrable_sam_name_handle.write(scrable_line_1)
            scrable_line_2 = '\t'.join(scramble_split_line_2)
            new_scrable_sam_name_handle.write(scrable_line_2)
            
            
            
    new_scrable_sam_name_handle.close()        
    new_sam_name_handle.close()
    
    
    print('  END: {:}'.format(paired_sam_1))

            
            
            


#####################################
#
#
#                MAIN
#
#
#####################################





read_depth_complexity_check_thresholds = []







directory_path = '/mnt/hgfs/main_ssd/SAM_files/'

directory_path = exp_output_path.trimmed_fastq_SAM_files


target_files = SAM_file_list


for file_path in target_files:
    
    sam_file_name = directory_path + file_path
    short_id = 'id'
    read_depth_complexity_check_thresholds = []
    process_sam_for_exon_intervals_simple_cigar_string(sam_file_name, short_id, read_depth_complexity_check_thresholds)





from concurrent.futures import ThreadPoolExecutor


synthetic_dir = exp_output_path.trimmed_fastq_SAM_files_synthetic


params = [( directory_path + name, ['T3'], []) for name in target_files]




#####################################
#
#
#            END  MAIN
#
#
#####################################










