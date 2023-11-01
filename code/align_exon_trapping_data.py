from experiment_paths.experiment_paths import *

import glob, sys, os, time
import subprocess


print('BEGIN alignment python script with system args: ', sys.argv )


data_pairs = list()


import os,sys
tmp_stdout = sys.stdout
sys.stdout = open(exp_output_path.trimmed_fastq_input_files + "read_alignment_stats_%s.txt" % (os.path.basename(sys.argv[1][:-9])), "w")




'''        Specify which library to align             '''

print('#################')

if len(sys.argv) > 1:

    R1_filename = sys.argv[1]
    R2_filename = sys.argv[2]
    library_number = int(sys.argv[3])
    sample_id = sys.argv[4]
    
    reads_dir = exp_output_path.trimmed_fastq_input_files
    
    
    R1_name = os.path.basename(R1_filename)[:-9] + "_trimmed.fastq.gz"
    R2_name = os.path.basename(R2_filename)[:-9] + "_trimmed.fastq.gz"
    
    r1 = R1_name
    r2 = R2_name
    
    print('*** command line call with arguments ***\n',R1_name,R2_name)
    print('reads_dir:', reads_dir)
    
    fastq_files = dict()
    
    data_pairs = [[r1,r2]]
    fastq_files[library_number] = [r1,r2]
    
    print('generated data_pairs:', data_pairs)
    
    
    












def run_data(reads_type, reads_data_list, reads_dir, hisat2_index, intron_hisat2_index):
    
    if reads_type == 'simple_paired':
        
        print(reads_data_list)
        
        reads_data_1 = reads_data_list[0]
        reads_data_2 = reads_data_list[1]
        
        input_data_1 = reads_data_1[:-9]  + "_simple_paired.fastq.gz"
        input_data_2 = reads_data_2[:-9]  + "_simple_paired.fastq.gz"
        output_file = input_data_1[:-9] + '.sam'
        output_file_gz = input_data_1[:-9] + '.sam.gz'
        output_file_bam = 'prealign_'+input_data_1[:-9] + '.bam'
        
        exp_output_path.trimmed_fastq_SAM_files
        '''
        unaligned_output = input_data_1[:-9] + 'unaligned.fastq'
        rf_output = input_data_1[:-9] + '_rf.sam'
        ff_output = input_data_1[:-9] + '_ff.sam'
        '''
        
        
        

        hisat2_call_string = "hisat2 --max-intronlen 1500 -p {:}".format( exp_output_path.hisat2_threads_count)
        hisat2_call_string = "%s %s" % (hisat2_call_string , hisat2_index)
        hisat2_call_string = "%s -1 %s" % (hisat2_call_string , reads_dir+input_data_1)
        hisat2_call_string = "%s -2 %s" % (hisat2_call_string , reads_dir+input_data_2)


        hisat2_call_string = "%s | gzip > %s" % (hisat2_call_string , exp_output_path.trimmed_fastq_SAM_files+output_file_gz)
    
    
    print('>>> hisat2 call:\n', hisat2_call_string)

    pipes = subprocess.Popen(hisat2_call_string, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = pipes.communicate()
    
    print("***** HISAT2 output *****")
    print(std_out.decode("utf-8") )    
    print(std_err.decode("utf-8") )

    



hisat2_index=exp_output_path.hisat2_index_std_chroms

for pair in data_pairs:
    print('loop pair:', pair)
    reads_data_1 = pair[0]
    reads_data_2 = pair[1]
    
    print("processing: ", reads_data_1)


    time_start = time.time()
    #print('simple_paired')
    reads_type = 'simple_paired'
    reads_data_list = [reads_data_1, reads_data_2]
    
    run_data(reads_type, reads_data_list, reads_dir, hisat2_index, '')
    time_stop = time.time() - time_start
    print('\nAlignment took %d seconds.\n\n' % (int(time_stop)))



sys.stdout = tmp_stdout 


print('END alignment python script with system args: ', sys.argv )





















