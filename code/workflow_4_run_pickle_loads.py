from experiment_paths.experiment_paths import *

import pickle

import time
start_time = time.time()


exon_count_build = 100
'''
exon_count_build = 1
'''



repo_path = exp_output_path.git_repo 

runfile(repo_path + 'load_pickles_1.py', args='exon_count_build %d' % exon_count_build, wdir = exp_output_path.git_repo, current_namespace=True)

runfile(repo_path + 'load_pickles_2.py', args='exon_count_build %d' % exon_count_build, wdir = exp_output_path.git_repo, current_namespace=True)

runfile(repo_path + 'load_pickles_3.py', args='exon_count_build %d' % exon_count_build, wdir = exp_output_path.git_repo, current_namespace=True)





#Accessory files

runfile('HEXEvent_1.py', args='exon_count_build %d' % exon_count_build, wdir = exp_output_path.git_repo, current_namespace=True)

pickle_path = exp_output_path.Dfam_pickle + "Dfam_pickle.pickle"
with open(pickle_path, "rb") as input_file:
    exon_ids_overlapping_repeat_list = pickle.load(input_file)
    exon_finder_unique_exon_ids_overlapping_repeat_list = pickle.load(input_file)


def delete_chrom_interval_tree(IT):
    for chrom in IT:
        for strand in IT[chrom]:
            for interval in (IT[chrom][strand].items()):
                IT[chrom][strand].remove(interval)
    return 1
    




for exon_id in aggregate_exon_dict:
    exon = aggregate_exon_dict[exon_id]
    exon['5ss']=''
    exon['3ss']=''
    







pickle_path = exp_output_path.pickle_merged + "ET_SAI_MES_exon_sets_%d.pickle" % (exon_count_build)
    
with open(pickle_path, "rb") as input_file:
    D           = pickle.load(input_file)
    ET_only     = pickle.load(input_file)
    SAI_only    = pickle.load(input_file)
    MES_only    = pickle.load(input_file)
    ALL_exon_intersect          = pickle.load(input_file)
    ALL_exon_finder_intersect   = pickle.load(input_file)




exon_finder_unique_union_list = list( set(ET_only).union(SAI_only).union(MES_only) )

def make_exon_id_IT(aaa, input_IT):
    IT = dict()
    for chrom in input_IT:
        IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
    
    for exon_id in aaa:
        ex = el.exon_id_values(exon_id)
        start = ex.start
        end = ex.end
        chrom = ex.chrom
        strand = ex.strand
        
        new_inteval = intervaltree.Interval(start, end, exon_id)
        IT[chrom][strand].add(new_inteval)

    return IT

exon_finder_unique_union_IT = make_exon_id_IT(exon_finder_unique_union_list, aggregate_exon_IT)





print("finished loading pickles: {:.2f} seconds".format(time.time()-start_time))








