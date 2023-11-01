



#from experiment_paths_synthetic.experiment_paths_synthetic import *

from experiment_paths.experiment_paths import *



import pyfaidx
import time
import glob

import gzip, os, sys



#input_fasta = 

#genome_fasta = pyfaidx.Fasta(input_fasta)


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
    
    #/home/pdf/repositories/exon_def/Process_ET/build_exon_interval_from_simple_paired_data_synthetic.py

    R1_filename = sys.argv[1]
    library_number = int(sys.argv[3])
    sam_file_dict = dict()
    #reads_dir = exp_output_path.trimmed_fastq_input_files
    #R1_name = os.path.basename(R1_filename)[:-9] + '_trimmed_simple_paired.sam.gz'
    
    sam_file_dict[library_number]= sys.argv[1]

    #sam_file_experiment_name = exp_output_path.trimmed_fastq_input_files




'''
experiment_name = '20201127'

experiment_name = "-1"

if experiment_name == '20201127':

    sam_file_dict = dict()
    sam_file_experiment_name = exp_output_path.trimmed_fastq_input_files_figures 

    fastq_files = dict()
    
    reads_dir = exp_output_path.trimmed_fastq_input_files 
    folder = exp_output_path.initial_fastq_input_files
    fastq_path = folder
    folder_files_list = glob.glob(folder+"*.fastq")
    trimmed_path = fastq_path + 'trimmed/'

    for entry in folder_files_list:
        if entry[-12:-10] == 'R1':
            lib_name = entry.split('/')[-1]
            sample = lib_name[11:13]
            number = int(lib_name[7:10])
            r1 = lib_name
            r1 = r1[:-6]+"_trimmed_simple_paired_simple_paired.sam"
            sam_file_dict[number] = reads_dir + r1
       
'''










def process_sam_for_exon_intervals_simple_cigar_string(sam_file_name, short_id,read_depth_complexity_check_thresholds):
    
    selected_chromosome = 'chr17'
    chr_length = 83257441
    
    paired_sam = sam_file_name
    
    #gzip.open(paired_sam, 'rt') 
    selected_alignments_file_name   = paired_sam[:-7] + '_selected.txt.gz'
    selected_alignments_name_handle = gzip.open( selected_alignments_file_name, 'wt',compresslevel=3) 
    
    #new_sam_name   = paired_sam[:-7] + '_chr17.sam.gz'
    #new_sam_name_handle = gzip.open( new_sam_name , 'wt') 

    
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
    
    library_stats_dict = dict()  #this is a way of saving various statistics I collect during processing the library
    
    print('ah', paired_sam)
    
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
                
            
            if ii in read_depth_complexity_check_thresholds:
                read_check_threshold_unique_exons_dict[ii] = list(unique_sites_dict.keys())
            

            
            if line_1[0] == '@':
                selected_alignments_name_handle.write(line_1)
                #new_sam_name_handle.write(line_1)
                continue
            
            
            split_line_1 = line_1.split('\t')
            split_line = split_line_1
            
            ##err.. this helped get rid of a single extra read alignment, but who knows what issue is lurking because of it.
            
            
            if int(split_line[4]) < 20:   #check read alginment is reasonable quality score
                continue
            
            
            if line_1.find("NH:i:1") >=0:
                1
            else:
                print('skipping multi-aligning read')
                int('b')
                continue
            
            
            
            line_2 = next(f)
            read_position += 1
            split_line_2 = line_2.split('\t')        
            
            #this field is present if the read has multiple alignments. Skip these reads.
            if line_1.find('ZS:') >= 0:
                continue
            
            #mate is unmapped
            if int(split_line_1[1]) & 8 != 0:
                continue
            
            if split_line_1[0] != split_line_2[0]:   #at least one time I had 3 reads with the same ID which broke the 2 read per ID assumption
                print("read_position",read_position)
                print(ii,'\n', line_1,'\n',line_2)
                print(split_line_1[0])
                print(split_line_2[0])
                int('a')
            
            #check if same chromosome. If not, don't bother with the alignment
            if split_line_1[2] != split_line_2[2]:
                continue
            
            #we are only processing a single chromosome
            #if split_line_1[2] != selected_chromosome:
            #    continue
            
            
            #if the alignment quality score is less than 20, then skip processing it
            if int(split_line_1[4]) < 20:
                 continue
            if int(split_line_2[4]) < 20:
                 continue
            
            chrom = split_line[2]
            if chrom not in chrom_OK_set:
                continue
            
            #if the alignment has a 0-length alignment (probably means mate unmapped)
            if int(split_line_1[8]) == 0:
                 continue
            
            
            
            #new_sam_name_handle.write(line_1)
            #new_sam_name_handle.write(line_2)
            
            #dist between = [8]
            if int(split_line_1[3]) < int(split_line_2[3]): # (forward strand)
                scrable_strand = '+'
            else:
                scrable_strand = '-'
            
            #exon_len = abs(int(split_line_1[3]) - int(split_line_2[3]) )
            
            #exon_len = np.random()
            #new_left = np.random.randint(1000,chr_length-1000)
            #if scrable_strand == '+':
            #    new_right = new_left+exon_len
            #else:
            #    new_right = new_left-exon_len
            
            scramble_split_line_1 = line_1.strip().split('\t')
            scramble_split_line_2 = line_2.strip().split('\t')
            
            scramble_split_line_1[9]  = '*'
            scramble_split_line_1[10] = '*'
            
            scramble_split_line_1.append(scramble_split_line_2[5])
            scramble_split_line_1.append(scramble_split_line_2[4])
            scramble_split_line_1.append('\n')
            
            
            
            #print(scramble_split_line_1)
            
            #scramble_split_line_1[3] = str(new_left)
            #scramble_split_line_2[3] = str(new_right)
            
            #scramble_split_line_1[7] = str(new_right)
            #scramble_split_line_2[7] = str(new_left)
            
            
            scrable_line_1 = '\t'.join(scramble_split_line_1)
            selected_alignments_name_handle.write(scrable_line_1)
            #print('**',scrable_line_1)
            
            #scrable_line_2 = '\t'.join(scramble_split_line_2)
            #selected_alignments_name_handle.write(scrable_line_2)
            
            
            
    selected_alignments_name_handle.close()        
    #new_sam_name_handle.close()
            
    
            


#####################################
#
#
#                MAIN
#
#
#####################################




library_unique_sites_dict = dict()
library_unique_exons_IT   = dict()

read_depth_unique_exon_recovery_dict = dict()
read_depth_complexity_check_thresholds = [800,1600,3200,6400,12800,25600, 51200,102400,204800]

#for short_id in ['3','4','6','7','8','9','10']:
for short_id in sam_file_dict.keys():
    #short_id = '3'
    sam_file_name = sam_file_dict[short_id]
    #read_depth_unique_exon_recovery_dict[short_id], library_unique_sites_dict[short_id], library_unique_exons_IT[short_id] = 

process_sam_for_exon_intervals_simple_cigar_string(sam_file_name, short_id, read_depth_complexity_check_thresholds)


#unique_non_overlapping_file_handle.close()
#reading_frame_file_handle.close()


'''
import pickle

out_filename = R1_name + '.pickle'
print('Save pickled data')
print(exp_output_path.pickle_individual + out_filename)
out_pickle = open(exp_output_path.pickle_individual + out_filename, 'wb')
pickle.dump(library_unique_sites_dict,out_pickle)
pickle.dump(library_unique_exons_IT,out_pickle)
#library_unique_exons_IT
out_pickle.close()
'''


#####################################
#
#
#            END  MAIN
#
#
#####################################

