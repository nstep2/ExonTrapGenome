import os
import time


#function to run a python script with current namespace and modify the current namespace
def run_script(script_path, current_namespace):
    with open(script_path, 'r') as file:
        exec(file.read(), current_namespace)
    return current_namespace



from datetime import timedelta

def time_hms(elapsed):  #return string of time in hh:mm:ss format, round to integer seconds
    return str(timedelta(seconds=int(elapsed)))
    

##################################################################




run_script_list = [
'chr_17_1.py',
'chr_17_2.py',
'ESE_0.py',
'chr_17_3.py',
'chr_17_4.py',
'chr17_shuffle_1.py',
'chr17_shuffle_2.py'
]



####################################################################





initial_time    = time.time() 
cur_time        = initial_time


for ii, val in enumerate(run_script_list):
    start_time = cur_time 
    script_path = run_script_list[ii]
    print('>>>\t',os.path.basename(script_path))
    runfile( script_path, args='exon_count_build %d' % exon_count_build, wdir = exp_output_path.git_repo, current_namespace=True)
    
    
    cur_time=time.time()
    print('completed in {:} seconds\n'.format( time_hms(cur_time-start_time )) )



print('total time: {:}'.format( time_hms(time.time() - initial_time)) )


