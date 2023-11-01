#~18 hours


from experiment_paths.experiment_paths import *




try:
    len(exon_data_dict_50_to_200)
except:
    exon_data_dict_50_to_200 = dict()

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
    reads_dir = exp_output_path.trimmed_fastq_input_files
    R1_name = os.path.basename(R1_filename)[:-9] + '_trimmed_simple_paired.sam.gz'
    sam_file_dict[library_number] = exp_output_path.trimmed_fastq_SAM_files + R1_name

    sam_file_experiment_name = exp_output_path.trimmed_fastq_input_files
else:
    print('no parameters given')
    exit()









def check_seq_for_frame(seq, partial_codon_start_seq, partial_codon_end_seq):
    
    #stop codons: TAG ("amber") TAA ("ochre") TGA ("opal")
    stop_codon_list = ['TAG','TAA','TGA']
    exon_seq = partial_codon_start_seq + seq + partial_codon_end_seq
    exon_length = len(exon_seq)
    current_pos = 0
    count_stop_codons = 0
    while current_pos < exon_length-2:
        codon = seq[current_pos:current_pos+3].upper()
        current_pos += 3
        if codon in stop_codon_list:
            count_stop_codons += 1
    
    return count_stop_codons





def process_sam_for_exon_intervals_simple_cigar_string(sam_file_name, short_id,read_depth_complexity_check_thresholds):
    
    paired_sam = sam_file_name
    #library_name = sam_file_name
    
    complicated_sam_name   = paired_sam[:-4] + '_multi_exon_alingments.sam'
    complicated_sam_handle = open( complicated_sam_name , 'w') 

    soft_clipped_sam_name   = paired_sam[:-4] + '_multi_exon_alingments.sam'
    soft_clipped_sam_handle = open( soft_clipped_sam_name , 'w') 

    
    list_5ss_scores = list()
    list_3ss_scores = list()
    
    
    unique_exons_IT = dict() # intervaltree.IntervalTree()
    stranded_unique_exons_IT = dict()
    stranded_exon_collections_IT = dict()
    
    chrom_OK_set = ['chr%d'%(x) for x in range(1,23)]
    chrom_OK_set.append('chrX')
    chrom_OK_set.append('chrY')
    chrom_OK_set.append('chrM')
    
    exon_frame_length_counts = [0,0,0]
    
    count_large_alignments = 0
    paired_fragment_length_estimates_list = list()
    
    #process_unique_sites_flag = False
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
    
    print('paired_sam:', paired_sam)
    
    
    total_processed_reads = 0
    with gzip.open(paired_sam, 'rt') as f:
        for ii, line_1 in enumerate(f):
            
            total_processed_reads +=1
            
            if ii % 10000000 == 0:
                print('Processed %d lines' % (ii))
            if ii == 5000000000000:
                break
            
            
            if ii in read_depth_complexity_check_thresholds:
                read_check_threshold_unique_exons_dict[ii] = list(unique_sites_dict.keys())
            

            
            if line_1[0] == '@':
                soft_clipped_sam_handle.write(line_1)
                complicated_sam_handle.write(line_1)
                continue
            
            split_line_1 = line_1.split('\t')
            split_line = split_line_1
            
            
            
            
            if int(split_line[4]) < 20:   #check read alginment is reasonable quality score
                continue
            
            if line_1.find("NH:i:1") >=0:
                1
            else:
                print('skipping multi-aligning read')
                continue
            
            
            
            line_2 = next(f)
            split_line_2 = line_2.split('\t')        
            
            
            if line_1.find('ZS:') >= 0:
                continue

            
            #check if same chromosome. If not, don't bother with the alignment
            if split_line_1[2] != split_line_2[2]:
                continue
            
            
            #if the alignment quality score is less than 20, then skip processing it
            if int(split_line_1[4]) < 20:
                 continue
            
            chrom = split_line[2]
            if chrom not in chrom_OK_set:
                continue
            
            
            
            
            cigar_string = split_line_1[5]
            cigar_split = re.findall(r'[A-Za-z]+|\d+', cigar_string)
            
            
            if len(cigar_string) > 1 and (cigar_split[1] == 'S' or cigar_split[-1] == 'S') :
                if int(cigar_split[0]) > 0 or True:
                    if cigar_split[1] == 'S':
                        count_leading_soft_clipped_reads += 1
                    if cigar_split[-1] == 'S':
                        count_lagging_soft_clipped_reads += 1
                    
                    
                    soft_clipped_sam_handle
                    soft_clipped_sam_handle.write(line_1)
                    soft_clipped_sam_handle.write(line_2)                
                    continue
            
            
            
            if len(cigar_split) >= 6:
                
                found_N = False
                for cigar_ii in range(3,len(cigar_split), 2):
                    if cigar_split[cigar_ii] == 'N':
                        found_N = True
                
            
                if found_N == True:
                    complicated_sam_handle.write(line_1)
                    complicated_sam_handle.write(line_2)
                    count_complicated_cigar_string += 1
                    continue
            
            #### NOW REPEAT FOR READ 2 IF NEEDED
            cigar_string = split_line_2[5]
            cigar_split = re.findall(r'[A-Za-z]+|\d+', cigar_string)
            
            
            if len(cigar_string) > 1 and (cigar_split[1] == 'S' or cigar_split[-1] == 'S') :
                if int(cigar_split[0]) > 0 or True:
                    if cigar_split[1] == 'S':
                        count_leading_soft_clipped_reads += 1
                    if cigar_split[-1] == 'S':
                        count_lagging_soft_clipped_reads += 1
                        
                    soft_clipped_sam_handle
                    soft_clipped_sam_handle.write(line_1)
                    soft_clipped_sam_handle.write(line_2)                
                    continue
            
            
            if len(cigar_split) >= 6:
                
                found_N = False
                for cigar_ii in range(3,len(cigar_split), 2):
                    if cigar_split[cigar_ii] == 'N':
                        found_N = True
                
                if found_N == True:
                    complicated_sam_handle.write(line_1)
                    complicated_sam_handle.write(line_2)
                    count_complicated_cigar_string += 1
                    continue
            
            
            
            
            if int(split_line_1[1]) & 64 == 64:
                mate_read_1_flag = True
            else:
                mate_read_1_flag = False
            
            if mate_read_1_flag and abs(int(split_line[8]) ) < 1000 and int(split_line[8]) != 0 and split_line[2] != '*':
                
                exon_frame_length_counts[abs(int(split_line[8])%3)] += 1
                
                if (int(split_line[8])) > 0:
    
                    start_position_estimate = int(split_line[3])
                    end_position_estimate   = int(split_line[3]) + int(split_line[8])
                else:   
                    
                    start_position_estimate = int(split_line[7])
                    end_position_estimate   = int(split_line[7]) + abs(int(split_line[8]))
                
                
                
                if int(split_line[1]) & 16 != 16:
                    strand_flag = '-'
                    strand = '-'
                    
                    start_position_estimate = start_position_estimate - darkcycle_offset

                else:
                    strand_flag = '+'
                    strand = '+'
                    
                    end_position_estimate = end_position_estimate + darkcycle_offset

                

                stranded_unique_id = "%s:%d-%d:%s" % (split_line[2], start_position_estimate, end_position_estimate, strand)
                unique_id = stranded_unique_id
                paired_fragment_length_estimates_list.append(abs(int(split_line[8])))
                


                if unique_id not in unique_sites_dict:
                    stranded_unique_sites_dict[stranded_unique_id]  = True #this generates a list of stranded exons
                    
                    chrom = split_line[2]
                    unique_sites_dict[unique_id] = {'length':abs(int(split_line[8])), 'count':1, 'unique':True, 'read_1_seq':split_line[9], 'strand':strand_flag, 'chrom':chrom, 'start':start_position_estimate, 'end':end_position_estimate, '5ss_boundary_redetermined':False, '3ss_boundary_redetermined':False,'lib_count':{short_id:1} }
                    
                    if unique_id in exon_data_dict_50_to_200:
                        count_annotated += 1
                        
                    
                    if strand_flag == '+':
                        
                        try:
                            left_side  = start_position_estimate-1-20
                            right_side = start_position_estimate-1+3
                            left_side_N_padding=''
                            if left_side < 0:
                                left_side_N_padding = 'N'*abs(left_side)
                                left_side = 0
                            upstream_seq = left_side_N_padding + str(genome_fasta[chrom][left_side:right_side])
                        except:
                            
                            left_side  = start_position_estimate-1-20
                            right_side = start_position_estimate-1+3
                            print(chrom,left_side,right_side)
                            upstream_seq = str(genome_fasta[chrom][left_side:right_side])
                        
                        unique_sites_dict[unique_id]['3ss'] = upstream_seq
                        
                        try:
                            unique_sites_dict[unique_id]['3ss_score'] = maxent.score3(upstream_seq, matrix=matrix3)
                            list_3ss_scores.append(unique_sites_dict[unique_id]['3ss_score'])
                        except:
                            unique_sites_dict[unique_id]['3ss_score'] = -5000

                        
                        
                        downstream_seq = str(genome_fasta[chrom][end_position_estimate-5:end_position_estimate+8])
                        unique_sites_dict[unique_id]['5ss'] = downstream_seq
                        try:
                            unique_sites_dict[unique_id]['5ss_score'] = maxent.score5(downstream_seq[1:10], matrix=matrix5)
                        except:
                            unique_sites_dict[unique_id]['5ss_score'] = -100
                        
                        list_5ss_scores.append(unique_sites_dict[unique_id]['5ss_score'])

                        
                        
                        
                        exon_seq = str(genome_fasta[chrom][start_position_estimate-1:end_position_estimate-1])

                        unique_sites_dict[unique_id]['seq'] = exon_seq
                        
                        up_down_length = 400
                        fudge_len = 20

                        
                        try:
                            left_side  = start_position_estimate-4-up_down_length
                            right_side = start_position_estimate-1+up_down_length
                            left_side_N_padding = ''
                            if left_side < 0:
                                left_side_N_padding = 'N'*abs(left_side)
                                left_side = 0
                            upstream_region   = left_side_N_padding + str(genome_fasta[chrom][left_side:right_side])
                        except:
                            
                            left_side  = start_position_estimate-4-up_down_length
                            right_side = start_position_estimate-1+up_down_length
                            left_side_N_padding = ''
                            if left_side < 0:
                                left_side_N_padding = 'N'*abs(left_side)
                                left_side = 0
                            print(chrom, start_position_estimate,'\n',start_position_estimate-4-up_down_length,start_position_estimate-1+up_down_length)
                            print(left_side,right_side)
                            upstream_region   = left_side_N_padding + str(genome_fasta[chrom][left_side:right_side])
                            int('a')
                        
                        try:
                            left_side  = end_position_estimate-5-up_down_length
                            right_side = end_position_estimate+8+up_down_length
                            left_side_N_padding = ''
                            if left_side < 0:
                                left_side_N_padding = 'N'*abs(left_side)
                                left_side = 0
                            downstream_region = left_side_N_padding + str(genome_fasta[chrom][left_side:right_side])
                        except:
                            left_side  = end_position_estimate-5-up_down_length
                            right_side = end_position_estimate+8+up_down_length
                            print(chrom, start_position_estimate,'\n')
                            print(left_side,right_side)
                            downstream_region = str(genome_fasta[chrom][left_side:right_side])
                        
                        unique_sites_dict[unique_id]['upstream_region'] = upstream_region
                        unique_sites_dict[unique_id]['downstream_region'] = downstream_region
                        
                        
                        
                        
                        
                    else:
                        up_5ss_range   = start_position_estimate-10
                        down_5ss_range = start_position_estimate+3
                        if up_5ss_range < 0:
                            leadning_N = 'N'*abs(up_5ss_range)
                            up_5ss_range = 0
                            upstream_seq = leadning_N + str(genome_fasta[chrom][up_5ss_range:down_5ss_range].reverse.complement)
                        else:
                            try:
                                upstream_seq = str(genome_fasta[chrom][start_position_estimate-10:start_position_estimate+3].reverse.complement)
                            except:
                                print(start_position_estimate-10)
                                print(start_position_estimate+3)
                                len(genome_fasta[chrom])
                                upstream_seq = str(genome_fasta[chrom][start_position_estimate-10:start_position_estimate+3].reverse.complement)
                        
                        unique_sites_dict[unique_id]['5ss'] = upstream_seq
                        try:
                            unique_sites_dict[unique_id]['5ss_score'] = maxent.score5(upstream_seq[1:10], matrix=matrix5)
                        except:
                            unique_sites_dict[unique_id]['5ss_score'] = 500
                        
                        list_5ss_scores.append(unique_sites_dict[unique_id]['5ss_score'])
                        
                        downstream_seq = str(genome_fasta[chrom][end_position_estimate-1:end_position_estimate+2].reverse.complement)
                        unique_sites_dict[unique_id]['3ss'] = downstream_seq
                        if downstream_seq[1:3].upper() != 'AG' and downstream_seq[0:2].upper() != 'AG' :
                            count_non_ag_termini += 1
                        
                        downstream_seq = str(genome_fasta[chrom][end_position_estimate-1-3:end_position_estimate-1+20].reverse.complement)
                        unique_sites_dict[unique_id]['3ss'] = downstream_seq 
                        
                        try:
                            unique_sites_dict[unique_id]['3ss_score'] = maxent.score3(downstream_seq, matrix=matrix3)
                            list_3ss_scores.append(unique_sites_dict[unique_id]['3ss_score'])
                        except:
                            unique_sites_dict[unique_id]['3ss_score'] = -5000
                        
                        
                        
                        
                        
                        exon_seq = str(genome_fasta[chrom][start_position_estimate-1:end_position_estimate-1].reverse.complement)
                        unique_sites_dict[unique_id]['seq'] = exon_seq
                        
                        
                        
                        up_down_length = 400
                        fudge_len = 20
                        #try:
                        tmp_up_coord = end_position_estimate-1-up_down_length
                        tmp_dn_coord = end_position_estimate+2+up_down_length
                        if tmp_up_coord < 0:
                            tmp_up_coord=0
                        upstream_region = str(genome_fasta[chrom][tmp_up_coord:tmp_dn_coord].reverse.complement)
                        
                        
                        
                        tmp_up_coord =  start_position_estimate-10-up_down_length
                        tmp_dn_coord = start_position_estimate+3+up_down_length
                        if tmp_up_coord < 0:
                            tmp_up_coord = 0
                        downstream_region = str(genome_fasta[chrom][ tmp_up_coord: tmp_dn_coord ].reverse.complement)
                        
                        
                        unique_sites_dict[unique_id]['upstream_region'] = ''
                        unique_sites_dict[unique_id]['downstream_region'] = ''
                        
                    
                    
                    partial_codon_start_seq = ''
                    partial_codon_end_seq = ''

                    
                    
                    
                    if chrom not in unique_exons_IT:
                        unique_exons_IT[ chrom ] = intervaltree.IntervalTree()
                        stranded_unique_exons_IT[ chrom ] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
                        stranded_exon_collections_IT[ chrom ] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
                    

                    
                    if len( stranded_unique_exons_IT[ chrom ][strand].overlap( start_position_estimate, end_position_estimate) ) > 0:
                        unique_sites_dict[unique_id]['unique'] = False
                    
                    else:
                        count_non_overlapping_exons += 1

                        
                    
                    

                    unique_exons_IT[ chrom ][start_position_estimate:end_position_estimate]  = unique_id                    

                    stranded_unique_exons_IT[ chrom ][strand][start_position_estimate:end_position_estimate]  = [unique_id]
                
                else:  #the unique_id already exists so we just add it's counts
                    unique_sites_dict[unique_id]['count'] += 1
                    unique_sites_dict[unique_id]['lib_count'][short_id] += 1
                    
            
            if int(split_line[8]) > 1000 :
                count_large_alignments += 1
    
    
    
    print("reads processed for %d: %d" % (short_id, total_processed_reads))
    
    print( 'Start\n\t%s\n\tsam_file\n' % short_id )

    

    out_string = '%s\t%d\n' % (short_id, count_non_overlapping_exons)
    unique_non_overlapping_file_handle.write(out_string)
    
    
    
    
    
    
    complicated_sam_handle.close()
    soft_clipped_sam_handle.close()


    return read_check_threshold_unique_exons_dict, unique_sites_dict, stranded_unique_exons_IT






#####################################
#
#
#                MAIN
#
#
#####################################




library_unique_sites_dict = dict()
library_unique_exons_IT   = dict()

read_depth_unique_exon_recovery_dict = dict()  #this code is orphaned
read_depth_complexity_check_thresholds = [800,1600,3200,6400,12800,25600, 51200,102400,204800]

for short_id in sam_file_dict.keys():

    sam_file_name = sam_file_dict[short_id]
    read_depth_unique_exon_recovery_dict[short_id], library_unique_sites_dict[short_id], library_unique_exons_IT[short_id] = process_sam_for_exon_intervals_simple_cigar_string(sam_file_name, short_id, read_depth_complexity_check_thresholds)


unique_non_overlapping_file_handle.close()
reading_frame_file_handle.close()



import pickle

out_filename = R1_name + '.pickle'
print('Save pickled data')
print(exp_output_path.pickle_individual + out_filename)
out_pickle = open(exp_output_path.pickle_individual + out_filename, 'wb')
pickle.dump(library_unique_sites_dict,out_pickle)
pickle.dump(library_unique_exons_IT,out_pickle)

out_pickle.close()



#####################################
#
#
#            END  MAIN
#
#
#####################################












