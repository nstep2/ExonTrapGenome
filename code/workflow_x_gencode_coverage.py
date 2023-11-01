
from experiment_paths.experiment_paths import *


import os
import time


exon_count_build = 100


#function to run a python script with current namespace and modify the current namespace
def run_script(script_path, current_namespace):
    with open(script_path, 'r') as file:
        exec(file.read(), current_namespace)
    return current_namespace




from datetime import timedelta

def time_hms(elapsed):  #return string of time in hh:mm:ss format, round to integer seconds
    return str(timedelta(seconds=int(elapsed)))
    

script_list = [
    'parse_gencode.py',
    'chromosome_coverage.py',
    'chromosome_coverage_prints.py',
    'proportions.py',
    
    ]



initial_time    = time.time() 
cur_time        = initial_time


for ii, key in enumerate(script_list):
    start_time = cur_time 
    script_path = script_list[ii]
    print('>>>\n>>>\t{:}\n>>>'.format(os.path.basename(script_path)))
    runfile( script_path, args='exon_count_build %d' % exon_count_build, wdir = exp_output_path.git_repo, current_namespace=True)

    cur_time=time.time()
    print('completed in {:} seconds\n'.format( time_hms(cur_time-start_time )) )



print('total time: {:}'.format( time_hms(time.time()-initial_time )) )







































