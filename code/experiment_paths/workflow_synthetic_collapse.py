#from experiment_paths.experiment_paths import *


from experiment_paths.experiment_paths import *


import subprocess
import concurrent.futures


input_list = list()
input_list.append([1, "T3", "1_ETF_S1_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([2, "T3", "2_ETF_cleaned_S9_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([3, "T3", "3_ETF_cleaned_S10_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([4, "T3", "4_ETF_cleaned_S11_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([5, "T3", "5_ETF_cleaned_S12_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([6, "T4", "6_ETF_S2_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([7, "T4", "7_ETF_cleaned_S13_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([8, "T5", "8_ETF_cleaned_S14_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([9, "T5", "9_ETF_S3_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([10, "T5", "10_ETF_cleaned_S15_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([11, "T5", "11_ETF_S4_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([12, "T5", "12_ETF_cleaned_S16_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([13, "T5", "13_ETF_cleaned_S17_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([14, "T5", "14_ETF_cleaned_S18_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([15, "T5", "15_ETF_S5_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([16, "T5", "16_ETF_cleaned_S19_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([17, "T5", "17_ETF_cleaned_S20_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([18, "T5", "18_ETF_S6_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([19, "T5", "19_ETF_cleaned_S21_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([20, "T5", "20_ETF_cleaned_S7_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([21, "T4", "21_ETF_cleaned_S8_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([22, "T4", "22_ETF_cleaned_S22_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])
input_list.append([23, "T5", "23_ETF_cleaned_S23_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz"])


for entry in input_list:
    entry[2] = entry[2][:-7] + '_chr17_scrable.sam.gz'





#SAM_files = exp_output_path.trimmed_fastq_SAM_files

SAM_files = '/mnt/hgfs/main_ssd/et_main/SAM_synthetic/'

# Function to execute reshape_sam_read_assignments using subprocess
def execute_reshape_sam_read_assignments(entry):
    input_value = entry[1]
    file_path = SAM_files + entry[2]
    
    # Assuming reshape_sam_read_assignments is a Python script
    python_script = exp_output_path.git_repo+"reshape_sam_read_assignments.py"
    
    # Construct the subprocess command
    command = ["python", python_script, input_value, file_path]
    
    try:
        print('>>> begin processing\t%s' % entry[0] )
        subprocess.run(command, check=True)
        return "Finished #{}".format(entry[0])
    except subprocess.CalledProcessError as e:
        return "Error executing the command for entry #{}: {}".format(entry[0], e)

# Maximum number of parallel processes (adjust as needed)
max_parallel_processes = 4

# Using ThreadPoolExecutor to run the commands in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=max_parallel_processes) as executor:
    futures = [executor.submit(execute_reshape_sam_read_assignments, entry) for entry in input_list]
    
    # Wait for all the futures (commands) to complete
    for future in concurrent.futures.as_completed(futures):
        result = future.result()
        print(result)





execute_reshape_sam_read_assignments(input_list[0])



























