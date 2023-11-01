



from sys import argv
if 'exon_count_build' in argv:
    exon_count_build = int(argv[2])
else:
    exon_count_build = 100




from Bio import SeqIO

from exon_id_library.gff import get_annotated_exons_from_GENCODE_GFF3
from exon_id_library.gff import GENCODE_exon_class
from experiment_paths.experiment_paths import * #exp_output_path as exp_output_path

import exon_id_library.exon_id_lib as el

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import intervaltree



def opposite_strand(strand):
    if strand == '+':
        return '-'
    if strand == '-':
        return '+'
    raise('"%s" is not a strand!'%d)
    




def exon_id_len(exon_id):
    co = exon_id.split(':')[1]
    start = int(co.split('-')[0])
    end   = int(co.split('-')[1])
    length = abs(start-end)
    return length

def exon_id_len_from_list(exon_id_list):
    exon_id_length_list = list()
    for exon_id in exon_id_list:
        length = exon_id_len(exon_id)
        exon_id_length_list.append(length)
    return exon_id_length_list


#GENCODE_exon_dict_basic
def get_transcript_to_exon_set(GENCODE_exon_dict):
    tid_dict = dict()
    for key in GENCODE_exon_dict:
        entry = GENCODE_exon_dict[key]
        try:
            tid = entry.transcript_id
        except:
            print(key)
            print(entry)
            print(type(entry))

        if tid in tid_dict:
            tid_dict[tid].append(key)
        else:
            tid_dict[tid]=[key]

    return tid_dict


#get first,last exon_id pairs by sorting transcript exons
#return a list of pairs
def get_list_first_last_exons(tid_dict,GENCODE_exon_dict):
    first_last_exon_id_pair_list = list()
    middle_exon_id_list = list()

    for tid in tid_dict:
        t_exon_list = tid_dict[tid]
        t_exon_list = sorted(t_exon_list,key=lambda x:int(GENCODE_exon_dict[x].exon_number))

        first = t_exon_list[0]
        last  = t_exon_list[-1]

        first_exon_id = GENCODE_exon_dict[first].exon_id
        last_exon_id  = GENCODE_exon_dict[last].exon_id
        for ii in range(1,len(t_exon_list)-1):
            middle_exon_id_list.append(GENCODE_exon_dict[t_exon_list[ii]].exon_id)

        first_last_exon_id_pair_list.append([first_exon_id,last_exon_id])

    return first_last_exon_id_pair_list,middle_exon_id_list















import pickle
pickle_path = exp_output_path.pickle_merged + "Parse_ENSEMBLE_GENCODE_exons_load__%d.pickle" % (exon_count_build)

with open(pickle_path, "rb") as opened_pickle:
    tid_dict= pickle.load(opened_pickle)
    tid_comp_dict= pickle.load(opened_pickle)
    first_last_exon_id_pair_list= pickle.load(opened_pickle)
    
    middle_exon_id_list= pickle.load(opened_pickle)
    
    pc_first_last_exon_ids= pickle.load(opened_pickle)
    pc_middle_exon_id_list= pickle.load(opened_pickle)
        
    exon_GENCODE_comp= pickle.load(opened_pickle)
    transcript_id_count_comp= pickle.load(opened_pickle)
    dup_exon_dict= pickle.load(opened_pickle)
    UTR_annotation= pickle.load(opened_pickle)
    exon_id_sets_dict= pickle.load(opened_pickle)
    GENCODE_exon_dict_comp= pickle.load(opened_pickle)
    exon_id_transcript_sets_dict_comp   = pickle.load(opened_pickle)       
    exon_GENCODE_basic= pickle.load(opened_pickle)
    transcript_id_count_basic= pickle.load(opened_pickle)
    exon_id_sets_dict_basic= pickle.load(opened_pickle)
    GENCODE_exon_dict_basic= pickle.load(opened_pickle)
    exon_id_transcript_sets_dict = pickle.load(opened_pickle) 
    exon_RefSeq= pickle.load(opened_pickle)
    transcript_id_count_RefSeq   = pickle.load(opened_pickle)        
    protein_coding_ids = pickle.load(opened_pickle)
    protein_coding_ids_basic = pickle.load(opened_pickle)
    lncRNA_middle_exon_ids = pickle.load(opened_pickle)
    lncRNA_first_last_exon_ids = pickle.load(opened_pickle)
    lncRNA_exon_id_list= pickle.load(opened_pickle)
    pc_overlapping_trapped_exons_set= pickle.load(opened_pickle)
    gencode_exon_IT= pickle.load(opened_pickle)
    gencode_pc_exon_IT= pickle.load(opened_pickle)
    gencode_lncRNA_exon_IT= pickle.load(opened_pickle)  
    gencode_transcript_IT= pickle.load(opened_pickle) 
    gencode_lncRNA_overlapping_exon_ids_set=pickle.load(opened_pickle) 
    gencode_lncRNA_overlapping_exon_ids_set=pickle.load(opened_pickle) 
    gencode_pc_overlapping_exon_ids_set=pickle.load(opened_pickle) 
    
    pc_exon_tid_list=pickle.load(opened_pickle) 
    lncRNA_exon_tid_list=pickle.load(opened_pickle) 
    gencode_lncRNA_transcript_IT=pickle.load(opened_pickle) 
    gencode_pc_transcript_IT=pickle.load(opened_pickle) 
    stranded_gencode_overlapping_trapped_exons_set=pickle.load(opened_pickle) 

    GENCODE_gene_dict=pickle.load(opened_pickle)
    
    pc_highly_overlapping = pickle.load(opened_pickle)


    
    
    print('\n\nFinished loading pickles.\n\n')
    
    
    
    
    
    
    
    
    
    
    
    
pickle_path = exp_output_path.pickle_merged + "Parse_ENSEMBLE_GENCODE_exons_load_reduced.pickle"


with open(pickle_path, "wb") as output_file:
    #pickle.dump(tid_dict, output_file)
    #pickle.dump(tid_comp_dict, output_file)
    pickle.dump(first_last_exon_id_pair_list, output_file)
    
    pickle.dump(middle_exon_id_list, output_file)
    
    pickle.dump(pc_first_last_exon_ids, output_file)
    pickle.dump(pc_middle_exon_id_list, output_file)
        
    pickle.dump(exon_GENCODE_comp, output_file)
    pickle.dump(transcript_id_count_comp, output_file)
    pickle.dump(dup_exon_dict, output_file)
    pickle.dump(UTR_annotation, output_file)
    pickle.dump(exon_id_sets_dict, output_file)
    pickle.dump(GENCODE_exon_dict_comp, output_file)
    pickle.dump(exon_id_transcript_sets_dict_comp   , output_file)       
    pickle.dump(exon_GENCODE_basic, output_file)
    pickle.dump(transcript_id_count_basic, output_file)
    pickle.dump(exon_id_sets_dict_basic, output_file)
    pickle.dump(GENCODE_exon_dict_basic, output_file)
    pickle.dump(exon_id_transcript_sets_dict , output_file) 
    pickle.dump(exon_RefSeq, output_file)
    pickle.dump(transcript_id_count_RefSeq   , output_file)        
    pickle.dump(protein_coding_ids , output_file)
    pickle.dump(protein_coding_ids_basic , output_file)
    pickle.dump(lncRNA_middle_exon_ids , output_file)
    pickle.dump(lncRNA_first_last_exon_ids , output_file)
    pickle.dump(lncRNA_exon_id_list, output_file)
    pickle.dump(pc_overlapping_trapped_exons_set, output_file)
    #pickle.dump(gencode_exon_IT, output_file)
    #pickle.dump(gencode_pc_exon_IT, output_file)
    #pickle.dump(gencode_lncRNA_exon_IT, output_file)  
    #pickle.dump(gencode_transcript_IT, output_file)  
    
    pickle.dump(gencode_lncRNA_overlapping_exon_ids_set, output_file) 
    pickle.dump(gencode_pc_overlapping_exon_ids_set, output_file) 
    
    
    #pickle.dump(pc_exon_tid_list, output_file) 
    #pickle.dump(lncRNA_exon_tid_list, output_file) 
    
    #pickle.dump(gencode_lncRNA_transcript_IT, output_file) 
    #pickle.dump(gencode_pc_transcript_IT, output_file) 
    





    
    print('\n\nFinished making reduced GENCODE pickle pickles.\n\n')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    