


import matplotlib.pyplot as plt
from experiment_paths.experiment_paths import *
import pickle

import intervaltree
import exon_id_library.exon_id_lib as el
#import time as el

import time
import numpy as np

def interval_tree_reducer_exon_interval(a, b):
    tmp_1 = a
    tmp_2 = b
    
    
    return list(set(tmp_1).union(set(tmp_2)))
    


def insert_commas(value):
    return f'{value:,}'


#min_exon_count_build = 1
#min_exon_count_build = 100

for min_exon_count_build in [1,100]:
    
    start_time = time.time()
    
    pickle_path = exp_output_path.pickle_merged + "updated_merged_aggregate_exon_dict.pickle"
    
    
    
    print("start loading merged pickle")
    opened_pickle = open(pickle_path,'rb')
    aggregate_exon_dict             = pickle.load(opened_pickle)
    aggregate_exon_IT               = pickle.load(opened_pickle)
    output_lists_dict               = pickle.load(opened_pickle)
    aggregate_exon_dict_meta_data   = pickle.load(opened_pickle)
    opened_pickle.close()
    print("finish loading merged pickle")
    print("this pickle was generated %s" % (aggregate_exon_dict_meta_data['date_join_aggregate_chromosome_data']))
    
    
    
    below_threshold_keys = list()
    for ii, key in enumerate(aggregate_exon_dict):
        exon = aggregate_exon_dict[key]
        if exon['count'] < min_exon_count_build:
            below_threshold_keys.append(key)
    below_threshold_keys=list(set(below_threshold_keys))
        
    len(below_threshold_keys)
    #rebuild intervall tree from scratch due to rare presumed tree rebuild issues
    for ii, key in enumerate(below_threshold_keys):
        del aggregate_exon_dict[key]
        
    for chrom in aggregate_exon_IT:
        for strand in aggregate_exon_IT[chrom]:
            del aggregate_exon_IT[chrom][strand]
            aggregate_exon_IT[chrom][strand] = intervaltree.IntervalTree()
    
    for ii, key in enumerate(aggregate_exon_dict):
        ex = el.exon_id_values(key)
        aggregate_exon_IT[ex.chrom][ex.strand].addi(ex.start, ex.end, [key])
        
        
    
    
    print('removed %s exons with count below: %d' % (insert_commas(len(below_threshold_keys)),min_exon_count_build))
    
    def remove_blacklist_Tra2B_exons(aggregate_exon_dict,aggregate_exon_IT):
        keep_list = list()
        remove_list = list()
        black_list_region = "chr3:185919514-185921094:-"
       
        blacklisted_exons_intervals = aggregate_exon_IT['chr3']['-'][185919497:185921019]
        
        for entry in blacklisted_exons_intervals:
            for ii, exon_id in enumerate(entry[2]):
                overlap, len_A, len_B = el.get_overlap_exon_A_with_B(exon_id, black_list_region)
                if overlap == 0:
                    #print(ii, exon_id)
                    keep_list.append(exon_id)
                else:
                    remove_list.append(exon_id)
        
        
        rebuild_tree = intervaltree.IntervalTree()
        for exon_id in keep_list:
            ex = el.exon_id_values(exon_id)
            ex_interval = intervaltree.Interval(ex.start, ex.end, [exon_id])
            rebuild_tree.add(ex_interval)
            
        rebuild_tree.merge_overlaps(data_reducer = interval_tree_reducer_exon_interval)
        
        for interval in blacklisted_exons_intervals:
            aggregate_exon_IT['chr3']['-'].remove(interval)
        
        for interval in rebuild_tree:
            aggregate_exon_IT['chr3']['-'].add(interval)
        
        
        
        removed_exon_id_expression_counts = el.get_exon_dict_counts(remove_list,aggregate_exon_dict)
        
        
        plt.figure()
        plt.boxplot([removed_exon_id_expression_counts],showfliers=False)
        plt.title('removed Tra2B intron exon counts:\n total removed = %s\ntotal exons = %d' % ( '{:,}'.format(sum(removed_exon_id_expression_counts)), len(removed_exon_id_expression_counts)))
        plt.ylabel('expression counts')
        plt.yscale('log')
        plt.show()
        
        for exon_id in remove_list:
            del aggregate_exon_dict[exon_id]
        
    
    
    remove_blacklist_Tra2B_exons(aggregate_exon_dict,aggregate_exon_IT)
    
    
    
    
    class exon():
        def __init__(self, entry):
            self.length = entry['length']
            self.count = entry['count']
            self.unique = entry['unique']
            self.read_1_seq = entry['read_1_seq']
            self.strand = entry['strand']
            self.chrom = entry['chrom']
            self.start = entry['start']
            self.end = entry['end']
            self.redetermined_5ss_boundary = entry['5ss_boundary_redetermined']
            self.redetermined_3ss_boundary = entry['3ss_boundary_redetermined']
            self.seq_5ss = entry['5ss']
            self.score_5ss = entry['5ss_score']
            self.seq_3ss = entry['3ss']
            self.score_3ss = entry['3ss_score']
            self.seq = entry['seq']
            self.upstream_region = entry['upstream_region']
            self.downstream_region = entry['downstream_region']


    end_time = time.time() - start_time
    print("time to load: %s seconds" % ('{:,}'.format(end_time)))
    
    
    
    
    
    
    tree_5ss = dict()
    tree_3ss = dict()
    
    for chrom in aggregate_exon_IT:
        tree_5ss[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
        tree_3ss[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
        for strand in aggregate_exon_IT[chrom]:
            dict_5ss = dict()
            dict_3ss = dict()
            for interval in aggregate_exon_IT[chrom][strand]:
                exon_list = interval[2]
                
                
                for exon_id in exon_list:
                    ex = el.exon_id_values(exon_id)
                    cord_5ss = ex.pos_5ss
                    cord_3ss = ex.pos_3ss
                    
                    if cord_5ss  not in dict_5ss:
                        dict_5ss[cord_5ss]=set()
                    if cord_3ss  not in dict_3ss:
                        dict_3ss[cord_3ss]=set()
                        
                    dict_5ss[cord_5ss].add(exon_id)
                    dict_3ss[cord_3ss].add(exon_id)
                    
            for cord_5ss in dict_5ss:
                exon_id_list = list(dict_5ss[cord_5ss])
                exon_id_list = sorted(exon_id_list)
                start_list = list()
                end_list = list()
                
                for exon_id in exon_id_list:
                    ex = el.exon_id_values(exon_id)
                    start_list.append( ex.start)
                    end_list.append( ex.end)
                    
                min_start = min(start_list)
                max_end = min(end_list)
                new_inteval = intervaltree.Interval(min_start, max_end, exon_id_list)
                tree_5ss[chrom][strand].add(new_inteval)
                
            for cord_3ss in dict_3ss:
                exon_id_list = list(dict_3ss[cord_3ss])
                exon_id_list = sorted(exon_id_list)
                start_list = list()
                end_list = list()
                
                for exon_id in exon_id_list:
                    ex = el.exon_id_values(exon_id)
                    start_list.append( ex.start)
                    end_list.append( ex.end)
                    
                min_start = min(start_list)
                max_end = min(end_list)
                new_inteval = intervaltree.Interval(min_start, max_end, exon_id_list)
                tree_3ss[chrom][strand].add(new_inteval)
    
    
    
    
    
    
    above_threshold_5ss = dict()
    above_threshold_3ss = dict()
    threshold_list = [0, 10,100,1000,10000,100000]
    
    
    def threshold_interval_counts(interval, threshold_list, above_threshold_dict, aggregate_exon_dict):
        exon_id_list = interval[2]
        max_count = 0
        for exon_id in exon_id_list:
            count = aggregate_exon_dict[exon_id]['count']
            
            if count > max_count:
                max_count = count
                
        for threshold in threshold_list:
            if threshold not in above_threshold_dict:
                above_threshold_dict[threshold] = 0
                
            if max_count > threshold:
                above_threshold_dict[threshold] += 1
                
    
    
    
    
    interval_5ss_lens_list = list()
    interval_3ss_lens_list = list()
    count_5ss_intervals = 0
    count_3ss_intervals = 0
    for chrom in tree_5ss:
        for strand in tree_5ss[chrom]:
            count_5ss_intervals += len(tree_5ss[chrom][strand])
            count_3ss_intervals += len(tree_3ss[chrom][strand])
            
            for interval in tree_3ss[chrom][strand]:
                threshold_interval_counts(interval, threshold_list, above_threshold_3ss, aggregate_exon_dict)
                interval_5ss_lens_list.append(len(interval[2]))
                
                
            for interval in tree_5ss[chrom][strand]:
                threshold_interval_counts(interval, threshold_list, above_threshold_5ss, aggregate_exon_dict)
                interval_3ss_lens_list.append(len(interval[2]))
                
    
    print("number of exon regions sharing the same 3ss: %s" % (insert_commas(count_3ss_intervals)))
    print("number of exon regions sharing the same 5ss: %s" % (insert_commas(count_5ss_intervals)))
    
    
    
    for threshold in sorted(list(above_threshold_5ss.keys())):
        print("threshold %s:" % insert_commas(threshold))
        print('\texon regions with unique 5ss: %s' % insert_commas(above_threshold_5ss[threshold]))
        print('\texon regions with unique 3ss: %s\n' % insert_commas(above_threshold_3ss[threshold]))
    
    
    
    
    primary_3ss_exon_id_set = list()
    max_threshold=min_exon_count_build
    count_3ss_interval_primary_above_threshold = 0
    for chrom in tree_3ss:
        for strand in tree_3ss[chrom]:
            for interval in tree_3ss[chrom][strand]:
                max_count = 0
                for exon_id in interval[2]:
                    count = aggregate_exon_dict[exon_id]['count']
                    if count > max_count:
                        max_count = count
                        max_exon_id = exon_id
                if max_count >= max_threshold:
                    count_3ss_interval_primary_above_threshold += 1
                    primary_3ss_exon_id_set.append(max_exon_id)
    
    primary_3ss_exon_id_set = set(primary_3ss_exon_id_set)
    
    print('3ss intervals with count above %d: %d' % (max_threshold, count_3ss_interval_primary_above_threshold))
    
    
    
    
    primary_5ss_exon_id_set = list()
    max_threshold=min_exon_count_build
    count_5ss_interval_primary_above_threshold = 0
    for chrom in tree_5ss:
        for strand in tree_5ss[chrom]:
            for interval in tree_5ss[chrom][strand]:
                max_count = 0
                for exon_id in interval[2]:
                    count = aggregate_exon_dict[exon_id]['count']
                    if count > max_count:
                        max_count = count
                        max_exon_id = exon_id
                if max_count >= max_threshold:
                    count_5ss_interval_primary_above_threshold += 1
                    primary_5ss_exon_id_set.append(max_exon_id)
    
    primary_5ss_exon_id_set = set(primary_5ss_exon_id_set)
    
    print('5ss intervals with count above %d: %d' % (max_threshold, count_5ss_interval_primary_above_threshold))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    plt.figure()
    plt.hist(interval_3ss_lens_list, bins=50)
    
    plt.figure()
    plt.hist(interval_5ss_lens_list, bins=50)
    
    
    
    count_above_100 = 0
    count_above_1000 = 0
    count_above_10000 = 0
    for exon_id in aggregate_exon_dict:
        entry = aggregate_exon_dict[exon_id]
        if entry['count'] > 100:
            count_above_100 += 1
        if entry['count'] > 1000:
            count_above_1000 += 1    
        if entry['count'] > 10000:
            count_above_10000 += 1    
    
    
    
    
    
    
    
    
    
    
    
    pickle_path = exp_output_path.pickle_merged + "merged_2_aggregate_exon_dict__%d.pickle" % (min_exon_count_build)
    
    with open(pickle_path, "wb") as output_file:
        pickle.dump(aggregate_exon_dict, output_file)
        pickle.dump(aggregate_exon_IT, output_file)
        pickle.dump(output_lists_dict, output_file)
        pickle.dump(aggregate_exon_dict_meta_data, output_file)
        pickle.dump(tree_5ss, output_file)
        pickle.dump(tree_3ss, output_file)
        pickle.dump(primary_3ss_exon_id_set, output_file)
        pickle.dump(primary_5ss_exon_id_set, output_file)
        pickle.dump(min_exon_count_build, output_file)
        
        
        
    
    pickle_path = exp_output_path.pickle_merged + "merged_2_aggregate_exon_dict_only__%d.pickle" % (min_exon_count_build)
    
    with open(pickle_path, "wb") as output_file:
        pickle.dump(aggregate_exon_dict, output_file)
        pickle.dump(primary_3ss_exon_id_set, output_file)
        pickle.dump(primary_5ss_exon_id_set, output_file)
        pickle.dump(min_exon_count_build, output_file)
        






exon_count_build = min_exon_count_build


exon_id_list=aggregate_exon_dict.keys()

file_name_path = exp_output_path.out_bed_main + "ET_data_columns_%d.txt" % (min_exon_count_build)
el.export_exon_id_list_with_counts_to_bed(exon_id_list, aggregate_exon_dict, file_name_path)


file_name_path = exp_output_path.out_bed_main + "ET_exons_%d.bed" % (min_exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, file_name_path)


if exon_count_build == 100:
    import subprocess
    import os
    exon_id_list=aggregate_exon_dict.keys()
    
    
    file_name_path = exp_output_path.out_bed_main + "ET_exons.bed" 
    el.export_exon_id_list_to_bed(exon_id_list, file_name_path)
    
    out = subprocess.check_output(['md5sum',file_name_path])
    with open('{:}.md5'.format(file_name_path),'w') as f:
        f.write('{:}'.format(out.decode('UTF-8') ))
        
        
    
    file_name_path = exp_output_path.out_bed_main + "ET_data_columns.txt"
    el.export_exon_id_list_with_counts_to_bed(exon_id_list, aggregate_exon_dict, file_name_path)
    
    out = subprocess.check_output(['md5sum',file_name_path])
    with open('{:}.md5'.format(file_name_path),'w') as f:
        f.write('{:}'.format(out.decode('UTF-8') ))
    




