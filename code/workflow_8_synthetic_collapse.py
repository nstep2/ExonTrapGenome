
from experiment_paths.experiment_paths import *


import subprocess
import concurrent.futures




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

collapsed_exons_list = ['{:}_1.fastq.gz_all_trimmed_simple_paired.sam.gz.'.format(x[2][0]) for x in input_file_name_list]

tmp_dict = {1:'T3', 2:'T4', 3:'T5', 4:'T5', 5:'T5'}
backbone_dict = {x[1]:tmp_dict[x[0]] for x in input_file_name_list}

collapsed_exons_list = [ [input_file_name_list[ii][1],  backbone_dict[input_file_name_list[ii][1]] , x] for  ii, x in enumerate(collapsed_exons_list)]


input_list = collapsed_exons_list


for entry in input_list:
    entry[2] = entry[2][:-8] + '_chr17_scrable.sam.gz'





SAM_files = exp_output_path.trimmed_fastq_SAM_files_synthetic


# Function to execute reshape_sam_read_assignments using subprocess
def execute_reshape_sam_read_assignments(entry):
    input_value = entry[1]
    file_path = SAM_files + entry[2]
    
    # Assuming reshape_sam_read_assignments is a Python script
    python_script = exp_output_path.git_repo+"sam_read_assignments.py"
    
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



























