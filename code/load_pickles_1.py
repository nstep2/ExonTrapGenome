import matplotlib.pyplot as plt
from experiment_paths.experiment_paths import *
import pickle

import intervaltree
import exon_id_library.exon_id_lib as el


import time
import numpy as np

def interval_tree_reducer_exon_interval(a, b):
    tmp_1 = a
    tmp_2 = b
    
    tmp = list()
    if type(tmp_1) == type(list()):
        1
    else:
        1
    
    return list(set(tmp_1).union(set(tmp_2)))
    



from sys import argv
if 'exon_count_build' in argv:
    exon_count_build = int(argv[2])
else:
    exon_count_build = 100






start_time = time.time()


pickle_path = exp_output_path.pickle_merged + "merged_2_aggregate_exon_dict__%d.pickle" % (exon_count_build)


print("start loading merged pickle")
opened_pickle = open(pickle_path,'rb')
aggregate_exon_dict             = pickle.load(opened_pickle)
print('aggregate_exon_dict loaded')
aggregate_exon_IT               = pickle.load(opened_pickle)
print('aggregate_exon_IT loaded')
output_lists_dict               = pickle.load(opened_pickle)
print('output_lists_dict loaded')
aggregate_exon_dict_meta_data   = pickle.load(opened_pickle)
print('aggregate_exon_dict_meta_data loaded')
tree_5ss                        = pickle.load(opened_pickle)
print('tree_5ss loaded')
tree_3ss                        = pickle.load(opened_pickle)
print('tree_3ss loaded')
primary_3ss_exon_id_set         = pickle.load(opened_pickle)
print('primary_3ss_exon_id_set loaded')
primary_5ss_exon_id_set         = pickle.load(opened_pickle)
print('primary_5ss_exon_id_set loaded')
min_exon_count_build            = pickle.load(opened_pickle)

opened_pickle.close()
print("finish loading merged pickle: %.1f seconds" % (time.time()-start_time))
print("this pickle was generated %s" % (aggregate_exon_dict_meta_data['date_join_aggregate_chromosome_data']))




