






from experiment_paths.experiment_paths import *

import subprocess

import time
script_start_time = time.time()

import os
import glob
import sys
import gzip
import numpy as np




class experiment_location():
    fastq_gz = dict()
    fastq_location = dict()
    fastq_files = dict()
    trimmed_fastq_files = dict()
    fastq_files_pairs = dict()
    aligned_fastq_files = dict()
    fasta_pair = []
    trimmed_path = ''
    def __init__(self):
        1
                
            #print(entry)
    def command_line_experiment_file(self, R1_filename, R2_filename, library_number, sample_id):  
        self.folder = exp_output_path.initial_fastq_input_files
        #folder_files_list = glob.glob(self.folder+"*.fastq")
        folder_files_list  = [R1_filename, R2_filename]
        #print(folder_files_list)
        self.fastq_path = self.folder
        self.trimmed_path = exp_output_path.trimmed_fastq_input_files
        self.tmp_path = exp_output_path.ssd_tmp
        #self.fastq_files['3'] = ['Hughes_001_3_PCR2_e_S1_R1_001.fastq', 'Hughes_001_3_PCR2_e_S1_R2_001.fastq']
        self.fastq_files = dict()
        self.number_to_sample = dict()
        
        self.fastq_files[library_number] = [R1_filename,R2_filename]
        self.number_to_sample[library_number]=sample_id
        
    
    def select_fasta_pair(self, index):
        self.fasta_pair = self.fastq_files[index]

 




#python 3 code for reverse complement - uses C function 
#https://codereview.stackexchange.com/a/151350
rev_comp = str.maketrans('ATGC', 'TACG')
def reverse_complement(astring):
    return astring.translate(rev_comp)[::-1]


#old information for flanking exon sequence used for PCR 
molecular_barcode_length = 8
exon_1_primer_seq = 'ACCTCTGGAATGGTTCAG'  #double check this
exon_2_primer_seq = 'GAGTCGTGACGTGAGT'
exon_2_primer_seq_RC = 'ACTCACGTCACGACTC'


# new sequence infromation PCR primers. The idea here is I want to easily track 
# and maintain the set of sequences without the used of a second file.
# in the future (when mature) this code should be modified to use a second file.
class exon_primer_seq_variants():
    exon_seq_set = dict()
    primer_1 = ''
    primer_2 = ''
    primer_1_length = 0
    primer_2_length = 0
    barcode_1_len = 0
    barcode_2_len = 0
    primer_1_rc = ''
    primer_2_rc = ''
    def __init__(self):
        self.exon_seq_set['original_1'] = 'ACCTCTGGAATGGTTCAG'
        self.exon_seq_set['original_2'] = 'ACTCACGTCACGACTC'
        
        self.exon_seq_set['NZS601'] = ''#'CCCTCTTTAGAGCTCCAG'
        self.exon_seq_set['NZS602'] = 'TGGGAGCAAGCAGTT'
        self.exon_seq_set['NZS603'] = ''#'CAGCCGGATGCAG'
        self.exon_seq_set['NZS604'] = 'GCAAAAGGGCGTCT'
        self.exon_seq_set['NZS605'] = ''#'CCTCTGGAATGGTTCAG'
        self.exon_seq_set['NZS606'] = 'AGTCGTGACGTGAGT'
        
        self.exon_seq_set['NZS607'] = 'CCCTCTTTAGAGCTCCAG'
        self.exon_seq_set['NZS608'] = ''#'TGGGAGCAAGCAGTT'
        self.exon_seq_set['NZS609'] = 'CAGCCGGATGCAG'
        self.exon_seq_set['NZS610'] = ''#'GCAAAAGGGCGTCT'
        self.exon_seq_set['NZS611'] = 'CCTCTGGAATGGTTCAG'
        self.exon_seq_set['NZS612'] = ''#'AGTCGTGACGTGAGT'

    def select_first_read_second_read(self, index1, index2, barcode_1_len, barcode_2_len):
        self.primer_1 = self.exon_seq_set[index1]
        self.primer_2 = self.exon_seq_set[index2]
        self.primer_1_rc = reverse_complement( self.exon_seq_set[index1] )
        self.primer_2_rc = reverse_complement( self.exon_seq_set[index2] )
        self.primer_1_length = len(self.primer_1)
        self.primer_2_length = len(self.primer_2)
        self.barcode_1_len = barcode_1_len
        self.barcode_2_len = barcode_2_len




library_choice = 'T5'
if len(sys.argv) > 1:
    print(sys.argv)
    library_choice = sys.argv[4]


if library_choice in ['T3']:
    exon_primer_seq = exon_primer_seq_variants()
    index1 = 'NZS608'
    index2 = 'NZS607' 
    barcode_1_len = 1  #   + 5
    barcode_2_len = 8


if library_choice in ['T4']:
    exon_primer_seq = exon_primer_seq_variants()
    index1 = 'NZS610'
    index2 = 'NZS609' 
    barcode_1_len = 0  #   + 5
    barcode_2_len = 8

if library_choice in ['T5']:
    exon_primer_seq = exon_primer_seq_variants()
    index1 = 'NZS612'
    index2 = 'NZS611' 
    barcode_1_len = 1  #   + 5
    barcode_2_len = 8  


exon_primer_seq.select_first_read_second_read(index1, index2, barcode_1_len, barcode_2_len)

exon_primer_seq.primer_1
exon_primer_seq.primer_2

'''        Specify which library to align             '''
if len(sys.argv) > 1:
    fasta_file_location = experiment_location()
    R1_filename = sys.argv[1]
    R2_filename = sys.argv[2]
    library_number = int(sys.argv[3])
    sample_id = sys.argv[4]
    fasta_file_location.command_line_experiment_file(R1_filename, R2_filename, library_number, sample_id)
    print('*** command line call with arguments ***\n',R1_filename,R2_filename)
    
    
else:
    print('error, no paramaters given!')
    exit()





class sequencer_read_len():
    end_1_len = 150
    end_2_len = 150

expected_read_len = sequencer_read_len()





from difflib import SequenceMatcher #has a string similarity ratio test
from Bio import pairwise2
import matplotlib.pyplot as plt





class illumina_seq_read():
    name = ''
    qual_score = ''
    strand = ''
    seq = ''
    has_barcode = False
    barcode     = ''
    

def check_read_1_primer_seq(seq_read):
    start = exon_primer_seq.barcode_1_len
    end   = exon_primer_seq.barcode_1_len + exon_primer_seq.primer_1_length
    read = seq_read.seq[start:end]
    
    return True
    
    val = SequenceMatcher(None, read , exon_primer_seq.primer_1)
    
    
    if val.ratio() > 0.70:
        return True
    else:
        return False
    

def extract_read_mol_barcode(seq_read):
    barcode_1_len = exon_primer_seq.barcode_1_len
    mol_barcode = seq_read.seq[:barcode_1_len ]
    
    return mol_barcode


def load_fasta_pair(fasta_filename_pair, directory, *args, **kwargs):
    1
    
    
    
    if fasta_filename_pair[0].find(".gz") > 0:
        gz_flag  = True
    else:
        gz_flag = False
    
    
    if gz_flag == True:
        
        buffer_factor = 16 #*32
        
        fasta_pair_handle = dict()
        
        fasta_pair_handle['1'] = gzip.open(fasta_filename_pair[0],'rt',buffer_factor*1024) #increase the buffer read by gunzip
        fasta_pair_handle['2'] = gzip.open(fasta_filename_pair[1],'rt',buffer_factor*1024)
        
        return fasta_pair_handle
        
    else:
        
        fasta_pair_handle = dict()
        
        fasta_pair_handle['1'] = open(directory+fasta_filename_pair[0],'r', buffering=buffer_factor*1024)
        fasta_pair_handle['2'] = open(directory+fasta_filename_pair[1],'r', buffering=buffer_factor*1024)
        
        print('opened input files:\n', directory+fasta_filename_pair[0], directory+fasta_filename_pair[1])
        
        return fasta_pair_handle








def create_output_files(fasta_filename_pair, directory):
    #directory = fasta_file_location.fastq_location
    
    filename_trim_len = 6  
    
    file_extension = ".fastq"
    
    gz_flag = fasta_filename_pair[0].find(".gz") > 0
    if gz_flag  == True:
        #directory = ''
        
        fasta_filename_pair[0] = os.path.basename(fasta_filename_pair[0])
        fasta_filename_pair[1] = os.path.basename(fasta_filename_pair[1])
        filename_trim_len = 9
        file_extension = ".fastq.gz"
    
    if gz_flag == True:
        output_handle = dict()
        output_handle['1'] = dict()
        output_handle['2'] = dict()
        
        
        
        tmp_1 =  os.path.basename(fasta_filename_pair[0])
        
        print('tmp_1', tmp_1)
        
        out_name = tmp_1[:-filename_trim_len] + '_trimmed_simple_paired' + file_extension
        print('\n\n   out_name', out_name)
        output_handle['1']['simple_paired'] = gzip.open(directory+out_name,'wt',compresslevel=2)
        
        tmp_2 =  os.path.basename(fasta_filename_pair[1])
        
        out_name =   tmp_2[:-filename_trim_len] + '_trimmed_simple_paired' + file_extension
        output_handle['2']['simple_paired'] = gzip.open(directory+out_name,'wt', compresslevel=2)
        
        
        out_name = fasta_filename_pair[0][:-filename_trim_len] + '_summary.txt'
        output_handle['summary'] = open(directory+out_name,'w')
    
        tmp_3 =  os.path.basename(fasta_filename_pair[0])
        
        out_name = tmp_3[:-filename_trim_len] + '_mol_barcodes.txt'
        barcode_file_handle = open(directory+out_name,'w')
    
        return output_handle, barcode_file_handle
    
    
    else:
        output_handle = dict()
        output_handle['1'] = dict()
        output_handle['2'] = dict()
        
        
        
        
        
        
        out_name = fasta_filename_pair[0][:-filename_trim_len] + '_mol_barcodes.txt'
        barcode_file_handle = open(directory+out_name,'w')
        
        
        
        
        
        
        if gz_flag != True:
            
            out_name = fasta_filename_pair[0][:-filename_trim_len] + '_trimmed_simple_paired' + file_extension
            output_handle['1']['simple_paired'] = open(directory+out_name,'w')
            out_name = fasta_filename_pair[1][:-filename_trim_len] + '_trimmed_simple_paired' + file_extension
            output_handle['2']['simple_paired'] = open(directory+out_name,'w')
        
        
        out_name = fasta_filename_pair[0][:-filename_trim_len] + '_summary.txt'
        output_handle['summary'] = open(directory+out_name,'w')
        
        return output_handle, barcode_file_handle


def close_input_output_file_handles(fasta_pair_handle):
    for key in fasta_pair_handle:
        fasta_pair_handle[key].close()

def close_input_output_file_dict_handles(fasta_pair_handle):
    for key in fasta_pair_handle:
        if key == 'summary':
            continue
        for key2 in fasta_pair_handle[key]:
            fasta_pair_handle[key][key2].close()


def get_next_seq_read_pair(line_set_dict):
    read_1 = illumina_seq_read() #make empty read
    read_2 = illumina_seq_read()    
    
    read_1_simple = illumina_seq_read() #make empty read
    read_2_simple = illumina_seq_read()    
    
    
    #get read 1    
    val = line_set_dict['1'][0]
    read_1.name = val
    read_1_simple.name=val
    val = line_set_dict['1'][1]
    read_1.seq = val
    read_1_simple.seq=val
    val = line_set_dict['1'][2]
    read_1.strand = val
    read_1_simple.strand=val
    val = line_set_dict['1'][3]
    read_1.qual_score = val
    read_1_simple.qual_score=val
    
    
    #get read 2
    val = line_set_dict['2'][0]
    read_2.name = val
    read_2_simple.name = val
    val = line_set_dict['2'][1]
    read_2.seq = val
    read_2_simple.seq = val
    val = line_set_dict['2'][2]
    read_2.strand = val
    read_2_simple.strand = val
    val = line_set_dict['2'][3]
    read_2.qual_score = val
    read_2_simple.qual_score = val
    
    read_pair_dict = dict()
    read_pair_dict['1'] = read_1
    read_pair_dict['2'] = read_2
    
    
    read_1_start = ( exon_primer_seq.primer_1_length + exon_primer_seq.barcode_1_len )
    read_1_end = read_1_start + 45
    #read_1_end=-1    #this gets the whole read
    
    #read_1_end   = -1*( exon_primer_seq.barcode_2_len + exon_primer_seq.primer_2_length )
    read_2_start = ( exon_primer_seq.barcode_2_len + exon_primer_seq.primer_2_length )
    read_2_end = read_2_start + 45
    #read_2_end = -1 #this gets the whole read
    
    
    #read_2_end   = -1*( exon_primer_seq.primer_1_length + exon_primer_seq.barcode_1_len )
    if read_1_end == 0:
        read_1_end = -1
    if read_2_end == 0:
        read_2_end = -1
    
    read_1_simple.qual_score = read_1_simple.qual_score[read_1_start:read_1_end]
    read_1_simple.seq = read_1_simple.seq[read_1_start:read_1_end]
    read_2_simple.qual_score = read_2_simple.qual_score[read_2_start:read_2_end]
    read_2_simple.seq = read_2_simple.seq[read_2_start:read_2_end]
    
    read_pair_dict['simple_1'] = read_1_simple
    read_pair_dict['simple_2'] = read_2_simple
    
    
    return read_pair_dict

def write_seq_read_pair(read_pair_dict, fasta_out_pair_handle_dict, barcode_file_handle, seq_read_flags):
    #get read 1    
    
    read_assignment_found_flag = False
    
    
    #barcode_entry = "%s\t%s\n" % (read_pair_dict['1'].name, read_pair_dict['1'].barcode)
    #barcode_file_handle.write(barcode_entry)
    
    #if 'simple_paired' in seq_read_flags and seq_read_flags['simple_paired'] == True:
    if True == True and len(read_pair_dict['simple_1'].seq) > 15 and len(read_pair_dict['simple_2'].seq) > 15:
        test_start_truncate_len_1 = 0
        test_start_truncate_len_2 = 0
        fasta_out_pair_handle_dict['1']['simple_paired'].write(read_pair_dict['simple_1'].name + '\n')
        fasta_out_pair_handle_dict['1']['simple_paired'].write(read_pair_dict['simple_1'].seq[test_start_truncate_len_1:] + '\n')
        fasta_out_pair_handle_dict['1']['simple_paired'].write(read_pair_dict['simple_1'].strand + '\n')
        fasta_out_pair_handle_dict['1']['simple_paired'].write(read_pair_dict['simple_1'].qual_score[test_start_truncate_len_1:] + '\n')
        
        #get read 2
        fasta_out_pair_handle_dict['2']['simple_paired'].write(read_pair_dict['simple_2'].name + '\n')
        fasta_out_pair_handle_dict['2']['simple_paired'].write(read_pair_dict['simple_2'].seq[test_start_truncate_len_2:] + '\n')
        fasta_out_pair_handle_dict['2']['simple_paired'].write(read_pair_dict['simple_2'].strand + '\n')
        fasta_out_pair_handle_dict['2']['simple_paired'].write(read_pair_dict['simple_2'].qual_score[test_start_truncate_len_2:] + '\n')
    
    return read_pair_dict






#this function does a few things:
    # it finds the empty reads
    # it finds the position of the adapter
def truncate_read_1_start(seq_read):
    
    start = exon_primer_seq.barcode_1_len + exon_primer_seq.primer_1_length
    
    seq_read.seq = seq_read.seq[start:]
    seq_read.qual_score = seq_read.qual_score[start:]
    
    #
    #val1 = SequenceMatcher(None, seq_read.seq, exon_2_primer_adapter_seq_RC )
    primer_2_seq = exon_primer_seq.primer_2_rc
    #print(primer_2_seq + '\n' + seq_read.seq + '\n\n' )
    if True == False:
        val1 = SequenceMatcher(None, seq_read.seq, primer_2_seq  )
    #val3 = regex.match(seq_read.seq, exon_2_primer_seq_RC)
    
    #print(val1.get_matching_blocks())
    
    #check that the reverse adapter is in the first spot
    start =   1*(exon_primer_seq.barcode_2_len)
    end =   1*(exon_primer_seq.primer_2_length + exon_primer_seq.barcode_2_len) # since there is a barcode in front
    
    
    
    
    
    
    if True == False:
        #needs to sort by longest match (z) and have 
        max_match_block_1 = 0
        max_match_block_2 = 0
        max_match_block_3 = 0
        list_of_longest_matches = list()
        for x,y,z in val1.get_matching_blocks():
            list_of_longest_matches.append([x,y,z])
        #list_of_longest_matches = sorted(list_of_longest_matches)
        list_of_longest_matches.sort(key=lambda x: x[2])  #using the operator library version may be faster
        #the above is sorted ascending. I think the next for loop is assumed to end with the longest match
    
    
        ### ===== >>>> this code area might not do what I think. it appears that it doesn't necessarily provide the right length to use to truncate the read or it misidentifies reads to truncate?
    
        for x,y,z in list_of_longest_matches:
            if z >= max_match_block_1:
                max_match_block_3 = max_match_block_2
                max_match_block_2 = max_match_block_1
                max_match_block_1 = z
        
        pos_to_truncate_read = len(seq_read.seq)  - 0      
        if len(list_of_longest_matches) > 0:
            #print(x,y,z, x-y)
            pos_to_truncate_read = x-y - 0  # this -1 was added to truncate the read by one extra base and hopefully fix an issue with grabbing an extra nucleotide from the sequencing primers.
    
    
    
        if val2.ratio() > 0.9:
            is_empty = True
        else:
            is_empty = False
    
        
    
        if sum([max_match_block_1, max_match_block_2])/len(exon_2_primer_seq_RC) >= 0.75:
            #print('ya')
            likely_has_adapter = True
        else:
            likely_has_adapter = False   
        
        
        
        
        
        if likely_has_adapter == True:
            seq_read.seq = seq_read.seq[:pos_to_truncate_read]
            seq_read.qual_score = seq_read.qual_score[:pos_to_truncate_read]    
        
        
        seq_read, SalI_shortened_flag, linker_shortened_flag  = query_read_for_restriction_site_intron_sequence(seq_read)
        
        
    else:
        is_empty = False
        likely_has_adapter = False
        SalI_shortened_flag = False
        linker_shortened_flag = False
        
    
        
    
    seq_read_flags = dict()
    seq_read_flags['is_empty']               = is_empty
    seq_read_flags['likely_has_adapter']     = likely_has_adapter
    seq_read_flags['SalI_shortened_flag']    = SalI_shortened_flag
    seq_read_flags['linker_shortened_flag']  = linker_shortened_flag
    
    
    return seq_read, seq_read_flags
    




def truncate_read_2_start(seq_read, seq_read_flags ):
       
    
    current_read_2_primer = exon_primer_seq.primer_2
    
    seq_read_flags['read2_adapter_OK'] = True 
    
    read_2_query_start = exon_primer_seq.barcode_2_len
    read_2_query_end   = exon_primer_seq.barcode_2_len + exon_primer_seq.primer_2_length


    
    return seq_read, seq_read_flags



def query_read_for_restriction_site_intron_sequence(read):
    
    SalI_seq = 'TCGACTTCATAAGCCAT'
    linker_seq = 'GAACAGGATGCTCTATC'


    
    linker_shortened_flag = False    
    SalI_shortened_flag = False
    

    
    if len(read.seq) > 0:
        SalI_seq_score = pairwise2.align.localmd( read.seq, SalI_seq, 1, -1, -2, -2, -2, -1,score_only=False)[0]
        #print(SalI_seq_score[2])
        #print(SalI_seq_score, 1.0*SalI_seq_score/len(SalI_seq))
        if(1.0 * SalI_seq_score[2] / len(SalI_seq)) >= 0.75:
            read.seq = SalI_seq_score[0][ :SalI_seq_score[3] ] #indexing is relative to the spacing alignment in the 
            #read.seq = read.seq.translate(string.maketrans('', ''), '-')
            remove_quality_score_len = SalI_seq_score[3] - SalI_seq_score[0][ :SalI_seq_score[3] ].count('-')
            read.seq = SalI_seq_score[0][ :SalI_seq_score[3] ].replace( '-','' )
            
            read.qual_score = read.qual_score[:remove_quality_score_len]
            SalI_shortened_flag = True
            #print(read.seq)
        
        
        linker_seq_score = pairwise2.align.localmd( read.seq, linker_seq, 1, -1, -2, -2, -2, -1,score_only=False)[0]
        #print(linker_seq_score)
        #print(linker_seq_score, 1.0*SalI_seq_score/len(SalI_seq))
        if(1.0 * linker_seq_score[2] / len(linker_seq)) >= 0.75:
            read.seq = linker_seq_score[0][ :linker_seq_score[3] ].replace( '-','' )
            remove_quality_score_len = linker_seq_score[3] - linker_seq_score[0][ :linker_seq_score[3] ].count('-')
            read.qual_score = read.qual_score[:remove_quality_score_len]
            #print(read.seq)
            linker_shortened_flag = True
        
        #print('')
    
    return read, SalI_shortened_flag, linker_shortened_flag











def library_primer_details(library_choice):
    #library_choice = 'T5'
    
    #for library_choice in ['T3','T4','T5']:
    
    
    #<-I have the first 3 bases of read 1 clipped
    if library_choice in ['T3']:
        exon_primer_seq = exon_primer_seq_variants()
        index1 = 'NZS608'
        index2 = 'NZS607' 
        barcode_1_len = 1  #   + 5
        barcode_2_len = 8
    
    
    if library_choice in ['T4']:
        exon_primer_seq = exon_primer_seq_variants()
        index1 = 'NZS610'
        index2 = 'NZS609' 
        barcode_1_len = 0  #   + 5
        barcode_2_len = 8
    
    if library_choice in ['T5']:
        exon_primer_seq = exon_primer_seq_variants()
        index1 = 'NZS612'
        index2 = 'NZS611' 
        barcode_1_len = 1  #   + 5
        barcode_2_len = 8  
        
    return exon_primer_seq, index1, index2, barcode_1_len, barcode_2_len





    
    

#def main():
if 1:
    
    legend_names = list()
    
    library_set = fasta_file_location.fastq_files.keys()
    
    print(library_set)
    
    for library_number in library_set:
        
        fasta_file_location.select_fasta_pair( library_number )
               
        legend_names.append(str(library_number))
        fasta_pair_handle = load_fasta_pair( fasta_file_location.fasta_pair, fasta_file_location.fastq_path) 
        
        fasta_out_pair_handle_dict, barcode_file_handle = create_output_files(fasta_file_location.fasta_pair, exp_output_path.ssd_tmp)
        
        
        
        barcode_list = list()
        read_lengths_passed_checks_list = list()
        
        fasta_file_count_stats = dict()
        fasta_file_count_stats['count_no_exon_primer_seq_reads'] = 0
        fasta_file_count_stats['count_empty_reads'] = 0
        fasta_file_count_stats['count_non_empty_reads_containing_adapter'] = 0
        fasta_file_count_stats['count_salI_site'] = 0
        fasta_file_count_stats['count_restriction_site_linker'] = 0
        fasta_file_count_stats['reads_too_short'] = 0
        
        
        file_reads_counter = 0
        
        for line_num, line1 in enumerate(fasta_pair_handle['1']):
            file_reads_counter += 1
            if file_reads_counter % 100000000 == 0:
                print("%d reads processed" % file_reads_counter, '\t', int(time.time() - script_start_time), 'seconds')
                #break
                
            
            line1 = line1[:-1]
            line2 = next(fasta_pair_handle['1'])[:-1]
            line3 = next(fasta_pair_handle['1'])[:-1]
            line4 = next(fasta_pair_handle['1'])[:-1]
            
            line_set_dict = dict()
            line_set_dict['1'] = [line1,line2,line3,line4]
            

            
            line1 = next(fasta_pair_handle['2'])[:-1]
            line2 = next(fasta_pair_handle['2'])[:-1]
            line3 = next(fasta_pair_handle['2'])[:-1]
            line4 = next(fasta_pair_handle['2'])[:-1]
            line_set_dict['2'] = [line1,line2,line3,line4]
            

            
            read_pair_dict = get_next_seq_read_pair(line_set_dict)
            
            
            #This should be part of the called function
            mol_barcode = extract_read_mol_barcode(read_pair_dict['1'])
            read_pair_dict['1'].barcode = mol_barcode
            read_pair_dict['1'].has_barcode = True
            
            if True == False:
                barcode_list.append(mol_barcode)
            


            
            read_pair_dict['1'], seq_read_flags  = truncate_read_1_start(read_pair_dict['1']) # replace the original entry with a condensed entry
            read_pair_dict['2'], seq_read_flags  = truncate_read_2_start(read_pair_dict['2'], seq_read_flags  )
            
            
            if seq_read_flags['is_empty'] == True or seq_read_flags['likely_has_adapter'] == True:
                1
            
            
            if seq_read_flags['is_empty'] == True:
                fasta_file_count_stats['count_empty_reads'] += 1
                continue
                        
            
            if seq_read_flags['likely_has_adapter'] == True:
                fasta_file_count_stats['count_non_empty_reads_containing_adapter'] += 1
            
            if seq_read_flags['is_empty'] == False and len(read_pair_dict['1'].seq) > 18 :
                write_seq_read_pair(read_pair_dict, fasta_out_pair_handle_dict, barcode_file_handle, seq_read_flags)
                
                
            elif len(read_pair_dict['1'].seq) <= 18:
                fasta_file_count_stats['reads_too_short'] += 1
                
                #read = query_read_for_restriction_site_intron_sequence(read_pair_dict['1'])
            
            if seq_read_flags['SalI_shortened_flag'] == True:
                if seq_read_flags['linker_shortened_flag'] == False:  
                    fasta_file_count_stats['count_salI_site'] += 1

            if seq_read_flags['linker_shortened_flag'] == True:
                fasta_file_count_stats['count_restriction_site_linker'] += 1
        
            
        
            if exp_output_path.limt_reads_processed > 0 and exp_output_path.limt_reads_processed <= file_reads_counter:
                break
        
        
        
         
        ####    Finish read processing
        ####    Now process some stats
        ####    
        ####
        ####
        
        
        out_string = fasta_file_location.fasta_pair[0][:-6]
        fasta_out_pair_handle_dict['summary'].write("%s\n" % out_string)
        
        for key in fasta_file_count_stats:
            print('\t',key, fasta_file_count_stats[key], '\t', 1.0*fasta_file_count_stats[key]/file_reads_counter)
            fasta_out_pair_handle_dict['summary'].write("%s\t%d\t%f\n" % (key,fasta_file_count_stats[key], 1.0*fasta_file_count_stats[key]/file_reads_counter))
        fasta_out_pair_handle_dict['summary'].write("\n\n")
        
        close_input_output_file_handles(fasta_pair_handle)
        close_input_output_file_dict_handles(fasta_out_pair_handle_dict)
        
        barcode_file_handle.close()
        
        
        tmp_1_path = fasta_out_pair_handle_dict['1']['simple_paired'].name
        tmp_2_path = fasta_out_pair_handle_dict['2']['simple_paired'].name
        
        final_dir = fasta_file_location.trimmed_path
        
        final_1_path = final_dir + os.path.basename(tmp_1_path)
        final_2_path = final_dir + os.path.basename(tmp_2_path)
        
        command_1 = "mv %s %s" % (tmp_1_path , final_1_path)
        command_2 = "mv %s %s" % (tmp_2_path , final_2_path)
        
        p1 = subprocess.Popen(command_1, shell=True)
        p2 = subprocess.Popen(command_2, shell=True)
        
        p1.wait()
        p2.wait()
        
        print(command_1+'\n')
        print(command_2+'\n')
        
        print('closed most file handles')
        print('begin plot of read lengths (commented out. Possible issue with plotting via command line?)')
        


    fasta_out_pair_handle_dict['summary'].close()
    
    print('===>>> currently the way of counting removed reads is flawed\n')
    


print('###\ntotal elapsed time:', int(time.time() - script_start_time), "seconds")







print('###\nbbduk start time:', int(time.time() - script_start_time), "seconds")



for index in fasta_file_location.fastq_files:
    fasta_file_location.select_fasta_pair(index)
    
    filename_trim_len = 9    
    file_extension=".fastq.gz"
    
    
    directory           = fasta_file_location.trimmed_path
    fasta_filename_pair = fasta_file_location.fasta_pair
    
    
    simple_file_1 = fasta_filename_pair[0][:-filename_trim_len] + '_trimmed_simple_paired' + file_extension
    simple_file_2 = fasta_filename_pair[1][:-filename_trim_len] + '_trimmed_simple_paired' + file_extension
    
    simple_file_1 = directory + simple_file_1 
    simple_file_2 = directory + simple_file_2 
    
    
    tmp_1 = simple_file_1 + "_tmp.fastq.gz"
    tmp_2 = simple_file_2 + "_tmp.fastq.gz"
    
    
    cleanup_command_1 = "cp %s %s" % (tmp_1, simple_file_1)
    cleanup_command_2 = "cp %s %s" % (tmp_2, simple_file_2)
    cleanup_command_1_rm = "rm %s" % (tmp_1)
    cleanup_command_2_rm = "rm %s" % (tmp_2)
    
    
    #tmp_1 = simple_file_1 #+ "_tmp.fastq.gz"
    #tmp_2 = simple_file_2 #+ "_tmp.fastq.gz"


    
    command = "bbduk.sh "
    command = command + " in1=%s" % (simple_file_1)
    command = command + " in2=%s" % (simple_file_2)
    command = command + " out1=%s" % (tmp_1)
    command = command + " out2=%s" % (tmp_2)
    command = command + " ref=%s" % (exp_output_path.adapter_fastq_file_path_absolute + 'trim_fasta_15mer_both_plasmid_forward.fa')
    command = command + " ktrim=r k=13 mink=5 hdist=1 tpe tbo threads={:} overwrite=true".format(exp_output_path.bbduk_threads_count )
    
    
    
    
    
    print('\n\ncommands:\n')
    print(command+'\n')
    print(cleanup_command_1+'\n')
    print(cleanup_command_2+'\n')
    
    
    
    print('start bbduk')
    p1 = subprocess.Popen(command, shell=True)
    p1.wait()
    p2 = subprocess.Popen(cleanup_command_1, shell=True)
    p3 = subprocess.Popen(cleanup_command_2, shell=True)
    p2.wait()
    p3.wait()
    p4 = subprocess.Popen(cleanup_command_1_rm, shell=True)
    p5 = subprocess.Popen(cleanup_command_2_rm, shell=True)
    p4.wait()
    p5.wait()
    
    
    print('finish bbduk')
    
    
    
    

print('###\nbbduk end time:', int(time.time() - script_start_time), "seconds")
















