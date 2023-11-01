from experiment_paths.experiment_paths import *


import subprocess
import os
from datetime import datetime






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





reads_path = exp_output_path.initial_fastq_input_files  





def run_process_reads_for_library(input_file_name_list, backbone_dict, reads_path, lib_index):
    
    
    
    bb, lib_number, file_string = input_file_name_list[lib_index]
    
    R1_var_1 = reads_path + "{:}_1.fastq.gz".format(file_string[0])
    R1_var_2 = reads_path + "{:}_1.fastq.gz".format(file_string[1])
    R1_var_3 = reads_path + "{:}_1.fastq.gz".format(file_string[2])
    R1_var_4 = reads_path + "{:}_1.fastq.gz".format(file_string[3])
    
    R2_var_1 = reads_path + "{:}_2.fastq.gz".format(file_string[0])
    R2_var_2 = reads_path + "{:}_2.fastq.gz".format(file_string[1])
    R2_var_3 = reads_path + "{:}_2.fastq.gz".format(file_string[2])
    R2_var_4 = reads_path + "{:}_2.fastq.gz".format(file_string[3])
    
    
    
    lib_num = lib_number
    backbone = backbone_dict[lib_number]
    
    
    
    
    R1_var = f"{R1_var_1}_all.fastq.gz"
    R2_var = f"{R2_var_1}_all.fastq.gz"
    
    start_time = datetime.now()
    
    print('begin concatenate library {:}'.format(lib_num) )
    subprocess.run(f"cat {R1_var_1} {R1_var_2} {R1_var_3} {R1_var_4} > {R1_var}", shell=True)
    subprocess.run(f"cat {R2_var_1} {R2_var_2} {R2_var_3} {R2_var_4} > {R2_var}", shell=True)
    print('end concatenate library {:}'.format(lib_num) )
    
    end_time = datetime.now()
    execution_time = end_time - start_time
    print(f"\n\nConcatenate took {execution_time.seconds} seconds to complete "
          f"({execution_time.seconds // 3600} hours {execution_time.seconds % 3600 // 60} minutes).\n\n")
    
    
    repo_path = exp_output_path.git_repo 
    print('exp_output_path.git_repo',exp_output_path.git_repo)
    subprocess.run(f"python {repo_path}prepare_reads_for_alignment_dark_cycle.py {R1_var} {R2_var} {lib_num} {backbone}", shell=True)

    subprocess.run(f"python {repo_path}align_exon_trapping_data.py {R1_var} {R2_var} {lib_num} {backbone}", shell=True)
    
    subprocess.run(f"python {repo_path}build_exon_interval_from_SAM_1.py {R1_var} {R2_var} {lib_num} {backbone}", shell=True)

    
    end_time = datetime.now()
    execution_time = end_time - start_time
    
    print(f"\n\nPython script took {execution_time.seconds} seconds to complete "
          f"({execution_time.seconds // 3600} hours {execution_time.seconds % 3600 // 60} minutes).\n\n")
    
    os.remove(R1_var)
    os.remove(R2_var)




    
lib_index=0
# single test function call
#run_process_reads_for_library(input_file_name_list, backbone_dict, reads_path, lib_index)











params = [(input_file_name_list, backbone_dict, reads_path, index) for index in range(23)]

'''

from concurrent.futures import ProcessPoolExecutor


with ProcessPoolExecutor(max_workers=4) as executor:
    results = executor.map(lambda p: run_process_reads_for_library(*p), params)


for result in results:
    1

input_file_name_list, backbone_dict, reads_path, lib_index

'''




from concurrent.futures import ThreadPoolExecutor



with ThreadPoolExecutor(max_workers=4) as executor:
    results = executor.map(lambda p: run_process_reads_for_library(*p), params)

print(list(results))










