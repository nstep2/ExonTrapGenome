

from experiment_paths.experiment_paths import *




from matplotlib.patches import Rectangle
import Bio
import Bio.SeqIO
import warnings
from Bio import BiopythonParserWarning
warnings.simplefilter('ignore', BiopythonParserWarning)

from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
matrix5 = load_matrix5()
matrix3 = load_matrix3()


embl_Dfam_consensus = exp_output_path.embl_Dfam_consensus
ouput_folder = exp_output_path.Dfam_repeats
pdf_plots = PdfPages(ouput_folder+'Dfam_repeats.pdf')


ouput_folder = exp_output_path.Dfam_repeats
pdf_repeat_consenus_plots = PdfPages(ouput_folder+'Dfam_repeats_consenus.pdf')






def insert_commas(value):
    return f'{value:,}'


def replace_slash(name):
    split = name.split('/')
    return '_'.join(split)


def flatten(t):
    return [item for sublist in t for item in sublist]


def get_5ss_scores_across_seq(seq):
    seq = str(seq).upper()
    score_5ss_pos_list = list()
    
    for ii in range(0,len(seq)-9):
        test_seq = seq[ii:ii+9]
        if test_seq.find('N') < 0:
            score_5ss = maxent.score5(test_seq, matrix=matrix5)
        else:
            score_5ss = -80
        score_5ss_pos_list.append([ii, score_5ss])
    
    return score_5ss_pos_list


def get_3ss_scores_across_seq(seq):
    seq = str(seq).upper()
    score_3ss_pos_list = list()
    
    for ii in range(0,len(seq)-23):
        test_seq = seq[ii:ii+23]
        if test_seq.find('N') < 0:
            score_3ss = maxent.score3(test_seq, matrix=matrix3)
        else:
            score_3ss = -80
        score_3ss_pos_list.append([ii, score_3ss])
        
    return score_3ss_pos_list














chromsomes_list = list()
for i in range(1,23):
    chromsomes_list.append( "chr%d" % (i) )
chromsomes_list.append('chrX')
chromsomes_list.append('chrY')
chromsomes_list.append('chrMT') #NCBI uses MT
chromsomes_list.append('chrM')

chromsomes_set = set(aggregate_exon_IT.keys())


dfam_IT = dict()
for chrom in chromsomes_list:
    dfam_IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}






print('\nstart parsing dfam genome hits\n')

dfam_hg38_hits = exp_output_path.dfam_hg38_hits



dfam_name_list = list()
with open(dfam_hg38_hits) as f:
    header = next(f)
    key_list = header[1:].strip().split('\t')
    for ii, line in enumerate(f):
        line_split = line.strip().split('\t')
        dfam_name_list.append(line_split[1])
dfam_name_list=list(set(dfam_name_list))
len(dfam_name_list)





import time
import datetime
now = datetime.datetime.now()
date_time = now.strftime("%m/%d/%Y, %H:%M:%S")
start_time = time.time()




problem_dfam_id_list    = list()
repeat_consensus_bases  = dict()
exon_consensus_bases    = dict()
repeat_descriptions     = dict()

repeat_stat_ranking_dict=dict()
repeat_exon_counts_dict = dict()
repeat_exon_bases_counts_dict = dict()
repeat_exon_bases_counts_expression_dict = dict()
repeat_genome_counts_dict = dict()
repeat_genome_seq_len_dict = dict()
repeat_consensus_name_list = list()

repeat_exon_genome_intersection_counts_dict = dict()

exon_ids_overlapping_repeat_list = list()
exon_ids_contained_repeat_list = list()


exon_finder_unique_exon_ids_overlapping_repeat_list = list()
exon_finder_unique_exon_ids_contained_repeat_list = list()




indexed_records = Bio.SeqIO.index(embl_Dfam_consensus, "embl")

repeat_associated_exon_ids_dict = dict()
all_repeat_read_count_dict = dict()

for fdam_ii, current_dfam_key in enumerate(set(dfam_name_list).intersection(indexed_records.keys())):
    #if fdam_ii == 50000000000:
    #    break
        
    total_time = str(datetime.timedelta(seconds=int(time.time()-start_time)))
    #print("\n^^^\nProcessing current_dfam_key %s (#%s\n)(time so far: %s (DD:HH:MM:SS))(average_time = %.1f second)\n\n" % (current_dfam_key,insert_commas(fdam_ii),total_time,(time.time()-start_time)/(fdam_ii+1)))
    #print('\nstart loading dfam entries for repeat consensus\n')
    
    repeat_consensus_record_dict = dict()
    repeat_subtype_dict = dict()
    species = set()
    
    record = indexed_records[current_dfam_key]
    repeat_consensus_record_dict[record.id] = record
    repeat_descriptions[record.id]          = record
    
    species=species.union(set(record.annotations['keywords']))
    line_split = record.annotations['comment'].split('\n')
    for line in line_split:
        if line[:7] == 'SubType':
            subtype = line[9:]
            if subtype not in repeat_subtype_dict:
                repeat_subtype_dict[subtype] = list()
            repeat_subtype_dict[subtype].append(record)


    
    dfam_hits_dict = dict()
    with open(dfam_hg38_hits) as f:
        header = next(f)
        key_list = header[1:].strip().split('\t')
        for ii, line in enumerate(f):
            line_split = line.strip().split('\t')
            
            #line_split[4] = 
            
            if float(line_split[4]) > 0.001:
                continue
            
            if line_split[1] != current_dfam_key:
                continue #only do one repeat at a time
                
            if line_split[0] not in chromsomes_set:
                continue
            if line_split[8] == '+':
                interval_id = "%s:%s-%s:%s" % (line_split[0],line_split[9],line_split[10],line_split[8])
            else:
                interval_id = "%s:%s-%s:%s" % (line_split[0],line_split[10],line_split[9],line_split[8])
            
            dfam_hits_dict[interval_id]=dict()
            f_list = [3,4,5,14]
            i_list = [6,7,9,10,11,12,13]
            ii_key_dict = dict([[ii, key]for ii, key in enumerate(key_list)])
            
            
            for ii in f_list:
                dfam_hits_dict[interval_id][ii_key_dict[ii]]=float(line_split[ii])
            for ii in i_list:
                dfam_hits_dict[interval_id][ii_key_dict[ii]]=int(line_split[ii])
            for ii in set(range(len(key_list))).difference(f_list).difference(i_list):
                dfam_hits_dict[interval_id][ii_key_dict[ii]]=line_split[ii]
            
            
            #pick start&end based on strand
            if line_split[8] == '+':
                dfam_hits_dict[interval_id]['start']  = dfam_hits_dict[interval_id]['ali-st']
                dfam_hits_dict[interval_id]['end']  = dfam_hits_dict[interval_id]['ali-en']
            else:
                dfam_hits_dict[interval_id]['end']  = dfam_hits_dict[interval_id]['ali-st']
                dfam_hits_dict[interval_id]['start']  = dfam_hits_dict[interval_id]['ali-en']
            
            dfam_hits_dict[interval_id]['chrom'] = dfam_hits_dict[interval_id]['seq_name']
            
    
    
    
    for key in repeat_consensus_record_dict:
        seq_len = len(repeat_consensus_record_dict[key].seq)
        
        if key != current_dfam_key:
            continue #only do one repeat at a time
        
        repeat_consensus_bases[key]=[np.zeros(seq_len+2000), np.zeros(seq_len+2000), np.zeros(seq_len+2000)]
    
        
        # index 0: repeat exon_hits
        # index 1: repeat genome bases
        # index 2: repeat genome bases with exon hit
    
    
    #get exon hits for repeat
    exon_count_threshold = 100
    
    
    for ii, repeat_region in enumerate(dfam_hits_dict):
        repeat = dfam_hits_dict[repeat_region]
        repeat_consensus_name = repeat['family_acc'] 
        repeat_consensus_name_list.append(repeat_consensus_name)
        repeat_exon_counts_dict[repeat_consensus_name] = 0
        repeat_genome_counts_dict[repeat_consensus_name] = 0
        repeat_genome_seq_len_dict[repeat_consensus_name] = 0
        repeat_exon_bases_counts_dict[repeat_consensus_name] = 0
        repeat_exon_bases_counts_expression_dict[repeat_consensus_name] = 0
        repeat_exon_genome_intersection_counts_dict[repeat_consensus_name] = 0
        
    ####################
    both_consensus_and_dfam_hit_keys = set(repeat_consensus_name_list).intersection(repeat_consensus_bases.keys())
    
    
    
    
    for ii, repeat_region in enumerate(dfam_hits_dict):
        repeat = dfam_hits_dict[repeat_region]
        rep = el.exon_id_values(repeat_region)
        
        
        intervals = aggregate_exon_IT[rep.chrom]['+'][rep.start:rep.end]
        for interval in intervals:
            for interval_exon_id in interval[2]:
                iex = el.exon_id_values(interval_exon_id)
                exon_ids_overlapping_repeat_list.append(interval_exon_id)
                if rep.start <= iex.start and rep.end >= iex.end:
                    exon_ids_contained_repeat_list.append(interval_exon_id)
                
                
        intervals = aggregate_exon_IT[rep.chrom]['-'][rep.start:rep.end]
        for interval in intervals:
            for interval_exon_id in interval[2]:
                iex = el.exon_id_values(interval_exon_id)
                exon_ids_overlapping_repeat_list.append(interval_exon_id)
                if rep.start <= iex.start and rep.end >= iex.end:
                    exon_ids_contained_repeat_list.append(interval_exon_id)
        
        
        
        
        #### Code for exon finder supplemental file
        ####
        ####
        intervals = exon_finder_unique_union_IT[rep.chrom]['+'][rep.start:rep.end]
        #for interval in intervals:
        for interval in intervals:
            #for interval_exon_id in interval[2]:
            interval_exon_id  = interval[2]
            iex = el.exon_id_values(interval_exon_id)
            exon_finder_unique_exon_ids_overlapping_repeat_list.append(interval_exon_id)
            if rep.start <= iex.start and rep.end >= iex.end:
                #exon_ids_contained_repeat_list.append(interval_exon_id)
                exon_finder_unique_exon_ids_contained_repeat_list.append(interval_exon_id)
            
        intervals = exon_finder_unique_union_IT[rep.chrom]['-'][rep.start:rep.end]
        
        for interval in intervals:
            interval_exon_id  = interval[2]
            iex = el.exon_id_values(interval_exon_id)
            exon_finder_unique_exon_ids_overlapping_repeat_list.append(interval_exon_id)
            if rep.start <= iex.start and rep.end >= iex.end:
                #exon_ids_contained_repeat_list.append(interval_exon_id)
                exon_finder_unique_exon_ids_contained_repeat_list.append(interval_exon_id)
        
        
        
        
        intervals = tree_3ss[rep.chrom][rep.strand][rep.start:rep.end] #only count primary exon of cluster
        
        repeat_consensus_name = repeat['family_acc'] 
        
        #'''
        if repeat_consensus_name in both_consensus_and_dfam_hit_keys:
            repeat_seq_bases = np.zeros(len(repeat_consensus_bases[repeat_consensus_name][1]))
            if repeat['strand'] == '+':
                repeat_seq_bases[repeat['hmm-st']+1000:repeat['hmm-en']+1000]=1
                
            else:  
                repeat_consensus_len = len(repeat_consensus_bases[ repeat['family_acc'] ][0])-2000
                repeat_seq_bases[repeat_consensus_len-repeat['hmm-st']+1000:repeat_consensus_len-repeat['hmm-en']+1000]=1


            repeat_consensus_bases[repeat_consensus_name][1] += repeat_seq_bases

            repeat_genome_counts_dict[repeat_consensus_name] += 1
            repeat_genome_seq_len_dict[repeat_consensus_name] += repeat['hmm-en']-repeat['hmm-st']

        
        for interval in intervals:
            max_count_exon_id = el.get_max_count_exon_id_in_list(interval[2], aggregate_exon_dict)

            for exon_id in [max_count_exon_id]:
                
                exon = aggregate_exon_dict[exon_id]
                ex = el.exon_id_values(exon_id)
                
                exon_bases = set(range(ex.start,ex.end))
                repeat_bases = set(range(rep.start,rep.end))
                
                #we want to know the exons repeats assocaite with
                if exon_id not in repeat_associated_exon_ids_dict:
                    repeat_associated_exon_ids_dict[exon_id] = list()
                    
                exon_repeat_intersection = exon_bases.intersection(repeat_bases)
                if len(exon_repeat_intersection) > 0:
                    repeat_associated_exon_ids_dict[exon_id].append(repeat)
                

                repeat_exon_genome_intersection_counts_dict[repeat_consensus_name] += len(exon_repeat_intersection)
                
                
                if len(exon_bases.intersection(repeat_bases)) > 0:
                    #bases of the exon along the consensus                
                    

                    if repeat['strand'] == '+':
                        bases = np.array(range(repeat['start']-ex.start+repeat['hmm-st']+1000, repeat['start']-ex.end +repeat['hmm-st']+1000))
                    else:
                        x_coord = repeat['hmm-st'] - (ex.start-repeat['end']) - ex.length
                        
                        bases = np.array(range(x_coord+1000, x_coord + ex.length +1000))
                        
                    
                    
                    if exon['count'] >= exon_count_threshold:
                        if repeat_consensus_name in both_consensus_and_dfam_hit_keys:

                            repeat_consensus_bases[repeat_consensus_name][2]+=repeat_seq_bases                            
                            repeat_exon_counts_dict[repeat_consensus_name] += 1 #count repeats with exons passing threshold for statistics
                            
                            repeat_exon_bases_counts_dict[repeat_consensus_name] += ex.length
                            repeat_exon_bases_counts_expression_dict[repeat_consensus_name] += ex.length*exon['count']     
                            
                    
                    
                    if exon['count'] < exon_count_threshold:
                        continue
                    
                    try:
                        repeat_consensus_bases[repeat_consensus_name]
                        for jj, pos in enumerate(bases):
                            repeat_consensus_bases[repeat_consensus_name][0][pos]+=1
                    except:
                        problem_dfam_id_list.append(repeat_consensus_name)
    
    
    
    for ii, repeat_region in enumerate(dfam_hits_dict):
        repeat = dfam_hits_dict[repeat_region]
        repeat_consensus_name = repeat['family_acc'] 
        if repeat_consensus_name not in repeat_consensus_bases:
            continue
        
        
    
    
    
    
    
    
    
    
    
    
    #print('\n\nCalculate 5ss and 3ss scores for consensus with exon hits')
    
    repeat_consensus_pos_score_dict = dict()
    for ii,repeat_id in enumerate(repeat_consensus_bases):
        
        if repeat_id != current_dfam_key:
            continue
        
        if np.sum(repeat_consensus_bases[repeat_id][0]) == 0:
            continue
        seq = str(repeat_consensus_record_dict[repeat_id].seq)
        
        scores = get_5ss_scores_across_seq(seq)
        if scores == -1:
            repeat_consensus_pos_score_dict[repeat_id] = -1
            continue
        a,b=zip(*scores)
        
        scores = get_3ss_scores_across_seq(seq)
        if scores == -1:
            repeat_consensus_pos_score_dict[repeat_id] = -1
            continue
        a,c=zip(*scores)
        repeat_consensus_pos_score_dict[repeat_id] = np.array([a,c,b],dtype=object) #ii,3ss,5ss
    
    
    
    
    
    
    repeat_intervals = list()
    repeat_intervals_genome = list()
    
    #print('\nstart searching for exons overlapping repeats\n')
    for jj, interval_id in enumerate(dfam_hits_dict):
        repeat = dfam_hits_dict[interval_id]
        
        repeat_intervals_genome.append([repeat['family_acc'], repeat['seq_name'], repeat['hmm-st'], repeat['hmm-en'], repeat['strand']])
        
        if repeat['strand'] == '+':
            repeat_consensus_start = repeat['hmm-st']
            repeat_consensus_end   = repeat['hmm-en'] 
        else:
            repeat_consensus_len = len(repeat_consensus_bases[ repeat['family_acc'] ][0])-2000
            
            
            repeat_consensus_start = repeat['hmm-st']
            repeat_consensus_end   = repeat['hmm-en']
        
        repeat_intervals.append([repeat_consensus_start,repeat_consensus_end])
        
    repeat_intervals=sorted(repeat_intervals, key=lambda x:(x[0],x[1]))
    
    
    
    if repeat['family_acc'] == 'DF0000317.4': #only plot for paper
        
    
        fig2,ax2 = plt.subplots()
        currentAxis2 = plt.gca()
        
        
        for jj, interval_pair in enumerate(repeat_intervals):
            
            repeat_len = interval_pair[1]-interval_pair[0]
            start_coordinate=interval_pair[0]
            
            height, width = 0.5, repeat_len
            X_coord, Y_coord = start_coordinate , jj+.1
            currentAxis2.add_patch(Rectangle((X_coord - 0.1, Y_coord - 0.1), width, height, alpha=1.0, facecolor='gray'))
        
        plt.ylim(-1, jj)
        plt.xlim(   -350, 
                    len(repeat_consensus_bases[repeat['family_acc']][0])-1650
                    +1)
        
        plt.title("all repeat fragments\n {:}: {:} ".format(repeat['family_name'], repeat['family_acc'] ))
        plt.tight_layout()
        pdf_repeat_consenus_plots.savefig(fig2)
        
    
    
    plot_counter = 0
    skip_count = 0
    for ii,key in enumerate(repeat_consensus_bases):
        
        
        if np.sum(repeat_consensus_bases[key]) == 0:
            skip_count += 1
            continue
        
        if key != current_dfam_key:
            continue
        
        if ii % 100 == 0:
            print('*** %s repeats processed.' % insert_commas(ii))
        
        repeat_stat_ranking_dict[key]=[0,0,0,0,repeat_consensus_record_dict[key].description]
        # 0 - 
        # 1 - 
        # 2 - 
        # 3 - 
        # 4 - 
        bases = repeat_consensus_bases[key][0]
        repeat_genome_bases = repeat_consensus_bases[key][1]
        repeat_bases = repeat_consensus_bases[key][2]
        
        data = bases
        n_min_max=np.nonzero(data)[0]
        #data = flatten(bases)
        #if len(data) == 0:
        #    #skip_count += 1
        #    continue
        
        #all_repeat_read_count_dict[key]
        if np.max(data) < 20:   #skip the repeats with few overlapping exons?
            skip_count += 1
            continue
        
        new_x_lim_1=n_min_max[0]
        if new_x_lim_1-1000 > 0:
            new_x_lim_1=1000
        new_x_lim_2=n_min_max[-1]
        if new_x_lim_2 < len(data)-1000:
            new_x_lim_2 = new_x_lim_2-1000
        
        
        
        if repeat['family_acc'] != 'DF0000317.4':
            continue
        
        plt.figure()
        ax1=plt.subplot()
        ax1.bar(range(len(data)), data, width=20)
        locs,labels = plt.xticks()
        plt.show()
        plt.close()
        
        
        
        
        fig, ax1=plt.subplots(1)
        fig.set_size_inches(8,5)
        ax1.bar(np.arange(len(data)), data, width=20, color='#1E88E5')
        ax1.set_title("%s:%s" % (key, repeat_consensus_record_dict[key].description))
        ax1.set_ylabel('exon sequence pileup')
        ax1.set_xlabel('repeat position, consensus starts at 0')
        new_labels = [ str(int(labels[jj].get_position()[0]-1000)) for jj,L in enumerate(labels)]
        ax1.set_xticks(locs)
        ax1.set_xticklabels(new_labels)
        
            
        x,score_3ss,score_5ss = repeat_consensus_pos_score_dict[key]
        max_data = max(data)
        factor = max_data/20
        length=1000
        
        
        ax1_twin=ax1.twinx()
        range_3ss=np.arange(length+x[0],length+x[-1]+1)
        ax1_twin.plot(range_3ss,factor*np.array(score_3ss),color='#D81B60',alpha=0.75)
        range_5ss=np.arange(length+x[0]-22+9-2,+length+x[-1])
        plt.plot(range_5ss,factor*np.array(score_5ss),color='#FFC107',alpha=1)
        plt.ylim(0,20*factor)
        plt.legend(['3ss','5ss'],loc='upper right')
        ax1.legend(['exon pileup'],loc='upper left')
        
        plt.yticks([0,.33*max_data, .66*max_data, max_data], ['0','5','10','15'])
        plt.ylabel('maxentscan score')
        plt.xlim(new_x_lim_1,new_x_lim_2)
                
        plot_counter += 1
        
        
        
        
        save_name = "repeat_pileup_%s_%s.svg"% (key, repeat_consensus_record_dict[key].description)
        save_name = replace_slash(save_name)
        save_path = exp_output_path.dfam_plots + save_name
        plt.tight_layout()
        pdf_plots.savefig(fig)
        plt.close()
        
        
        ax1.set_title("%s:%s" % (key, repeat_consensus_record_dict[key].description))
        
        save_name = "repeat_pileup_%s_%s.svg"% (key, repeat_consensus_record_dict[key].description)
        save_name = replace_slash(save_name)
        save_path = exp_output_path.dfam_plots + save_name
        plt.tight_layout()
        plt.savefig(save_path)
        
        pdf_plots.savefig(fig)
        plt.close()
        
        
        
        range_3ss=np.arange(length+x[0],length+x[-1]+1)
        
        
        
        range_3ss, factor*np.array(score_3ss)
        range_5ss,factor*np.array(score_5ss)
        np.arange(len(data)), data
        
        min_x = min(list(range_3ss)+list(range_5ss)+list(np.arange(len(data))))
        max_x = max(list(range_3ss)+list(range_5ss)+list(np.arange(len(data))))
        
        len(score_3ss)
        len(data)
        
        new_x_lim_2-new_x_lim_1
        
        with open( exp_output_path.Dfam_repeats + '6B.txt','w') as f_6B:
            f_6B.write('\t'.join(['repeat_val','3ss_val','5ss_val']) + '\n')
            for ii, val in enumerate(range(min_x, max_x)):
                
                if ii in np.arange(len(data)):
                    out_p = data[ii]
                else:
                    out_p = '.'
                    
                if ii in range_3ss:
                    #out_3 = score_3ss[ii]
                    out_3 = ii
                else:
                    out_3 = '.'
                    
                if ii in range_5ss:
                    #out_5 = score_5ss[ii]
                    out_5 = ii
                else:
                    out_5 = '.'
                
                f_6B.write('\t'.join( [str(out_p), str(out_3), str(out_5)] ) + '\n')
                    
                
        
        
        #[0] - found_repeat_bases/all_genome_repeat_bases
        #[1] - found_repeat_bases per bp of repeat
        #[2] - median of recovered non-zero found bases ()
        #[3] - max of recovered repeat bases across all positions of consensus
        
        
        
        with np.errstate(divide='ignore'):
            repeat_stat_ranking_dict[key][0] = sum(repeat_bases[1000:-1000])/(sum(repeat_genome_bases[1000:-1000]))
        
        with np.errstate(divide='ignore'):
            repeat_stat_ranking_dict[key][1] = sum(repeat_bases)/((len(repeat_genome_bases)-2000))
        
        with np.errstate(divide='ignore'):
            tmp = (repeat_bases/repeat_genome_bases)
        tmp[np.isnan(tmp)]=0
    
        repeat_stat_ranking_dict[key][2] = np.median([x for x in tmp if x >0])
    
        repeat_stat_ranking_dict[key][3] = max([x for x in tmp if x >0])
        
    
    



exon_ids_overlapping_repeat_list=list(set(exon_ids_overlapping_repeat_list))
exon_ids_contained_repeat_list=list(set(exon_ids_contained_repeat_list))





#exp_output_path.Dfam_pickle


pickle_path = exp_output_path.Dfam_pickle + "Dfam_pickle.pickle"

with open(pickle_path, "wb") as output_file:
    pickle.dump(exon_ids_overlapping_repeat_list, output_file)
    exon_ids_contained_repeat_list=list(set(exon_ids_contained_repeat_list))
    pickle.dump(exon_finder_unique_exon_ids_overlapping_repeat_list, output_file)
  
'''
with open(pickle_path, "rb") as input_file:
    exon_ids_overlapping_repeat_list=pickle.load(input_file)
'''
    
    






















r_rank = [(val,repeat_stat_ranking_dict[val]) for k, val in enumerate(repeat_stat_ranking_dict)]
ka, va = zip(*r_rank)
v1,v2,v3,v4, v5 = zip(*va)

v3=np.array(v3)
v4=np.array(v4)
ratio = v3/v4
ratio[np.isnan(ratio)]=0

key_ratio_pair=[[ka[ii],ratio[ii],v1[ii], v5[ii]] for ii, key in enumerate(ratio)]
sorted_key_ratio_pair = sorted(key_ratio_pair, key = lambda x: x[1],reverse=True)
sorted_key_avg_pair = sorted(key_ratio_pair, key = lambda x: x[2],reverse=True)
sorted_key_ratio_pair[:25]
sorted_key_avg_pair[:25]



import scipy
import scipy.stats

problem_dfam_id_list  #these don't have a consensus repeat


repeat_exon_counts_dict    # count of exons passing threshold
repeat_genome_counts_dict  # number genome instances
repeat_genome_seq_len_dict # bases in the genome




genome_instances_list = list()
spliced_repeat_ratio_list = list()
spliced_repeat_exon_list = list()
spliced_repeat_exon_exp_list = list()
spliced_repeat_exon_vs_genome_len_list = list()
spliced_repeat_exon_exp_list = list()
exon_count_pval_list = list()
for key in repeat_exon_counts_dict:
    exon_count = repeat_exon_counts_dict[key]
    genome_count = repeat_genome_counts_dict[key]
    genome_exon_len_count = repeat_exon_bases_counts_dict[key]
    spliced_repeat_exon_exp = repeat_exon_bases_counts_expression_dict[key]
    genome_len_count = repeat_genome_seq_len_dict[key]
    
    #repeat_genome_counts_dict[key]
    if genome_count == 0:
        continue
    genome_instances_list.append(genome_count)
    ratio_spliced_repeat = exon_count/genome_count
    spliced_repeat_ratio_list.append(ratio_spliced_repeat)
    spliced_repeat_exon_vs_genome_len_list.append(genome_exon_len_count/genome_len_count)
    spliced_repeat_exon_exp_list.append(spliced_repeat_exon_exp/genome_len_count)
    
    exon_rate = 1200000/6200000000
    pval = scipy.stats.binom_test(exon_count, genome_len_count, exon_rate, alternative='greater')
    exon_count_pval_list.append(pval)
    
    
    
    
repeat_exon_genome_intersection_counts_dict






'''


f_path = exp_output_path.Dfam_repeats + 'repeat_fraction_genome_fragments.txt'
with open(f_path,'w') as f:
    header = 'Dfam_id\tavg\tfraction_recovered\tDfam_description\n'
    f.write(header)
    for entry in sorted_key_avg_pair:
        #line = '\t'.join(entry)+'\n'
        line = '%s\t%f\t%f\t%s\n' % (entry[0],entry[1],entry[2],entry[3])
        f.write(line)


'''






#top repeats by fraction of repeat genome fragments recovered
sorted_key_avg_pair[:25]





pdf_repeat_consenus_plots.close()

pdf_plots.close()













    










        #[0] - found_repeat_bases/all_genome_repeat_bases
        #[1] - found_repeat_bases per bp of repeat
        #[2] - median of recovered non-zero found bases ()
        #[3] - max of recovered repeat bases across all positions of consensus
        #[4] - repeat description
        #[5] - sum(repeat exon bases) / sum(repeat genome bases)
        #[6] - sum(repeat exon bases) * sum(exon counts) / sum(repeat genome bases)
        #[7] - sum(repeat exon bases only overlapping repeat fragments (exclude area outside repeat fragment even if map to consensus region)) * sum(exon counts) / sum(repeat genome bases)
        #[7] - partial exon pval - bionomial
        #[8] - partial exon pval - heypergeometric
        
#repeat_exon_genome_intersection_counts_dict


avg_exon_exp = np.mean( el.get_exon_dict_counts(primary_3ss_exon_id_set, aggregate_exon_dict))

avg_exon_exp_bases = list()
exon_middle_exon_bases = 0
for exon_id in primary_3ss_exon_id_set:
    ex=el.exon_id_values(exon_id)
    exon_middle_exon_bases += ex.length
    avg_exon_exp_bases.append(ex.length*aggregate_exon_dict[exon_id]['count'])
avg_exon_exp_bases = np.mean(avg_exon_exp_bases)

avg_exon_exp_bases
avg_exon_exp

set(repeat_genome_seq_len_dict).difference(set(repeat_stat_ranking_dict))


import pandas as pd


['id','name','exon_rep_intersect_bases', 'exon_genome_bases', 'rep_genome_bases', 'exon_genome_count', 'rep_genome_counts','repeat_exon_bases_counts_expression']





rep_genome_stats_df=dict()
rep_genome_stats_df['exon_bases'] = exon_middle_exon_bases

rep_genome_stats_df['exon_count'] = len(primary_3ss_exon_id_set)

rep_genome_stats_df['genome_bases'] = 3088286401

rep_genome_stats_df['avg_exon_exp_bases'] = avg_exon_exp_bases
rep_genome_stats_df['avg_exon_exp']       = avg_exon_exp


repeat_stat_rows_dict = dict()

#couple of different 
#r_list = list()
for ii, key in enumerate(repeat_genome_seq_len_dict):
    
    if key not in repeat_stat_ranking_dict: #not sure why 3 keys are not in this
        continue
    
    #rep_stats_df[key]
    
    new_row = dict()
    
    new_row['rep_genome_bases'] = repeat_genome_seq_len_dict[key]
    new_row['exon_rep_intersect_bases'] = repeat_exon_genome_intersection_counts_dict[key]
    new_row['exon_genome_bases'] = repeat_exon_bases_counts_dict[key]
    new_row['exon_genome_count'] = repeat_exon_counts_dict[key]
    new_row['rep_genome_counts'] = repeat_genome_counts_dict[key]
    new_row['repeat_exon_bases_counts_expression'] = repeat_exon_bases_counts_expression_dict[key]
    
    if new_row['exon_genome_bases'] != 0:
        new_row['avg_exon_exp'] = new_row['repeat_exon_bases_counts_expression']/new_row['exon_genome_bases']
    else:
        new_row['avg_exon_exp'] = 0
    
    new_row['id'] = key
    new_row['name'] = repeat_stat_ranking_dict['DF0000402.4'][4]
    
    
    repeat_stat_rows_dict[key] = new_row

    
    
    
    #(exon_middle_exon_bases/2) is used to account for only looking at exons in the same strand as the repeat
    genome_len=3088286401
    
    #over representation
    test_hyper_3 = scipy.stats.hypergeom(M=genome_len, n=(exon_middle_exon_bases/2), N=repeat_genome_seq_len_dict[key]).logsf(repeat_exon_genome_intersection_counts_dict[key]-1) 
    
    #under representation
    test_hyper_4 = scipy.stats.hypergeom(M=genome_len, n=(exon_middle_exon_bases/2), N=repeat_genome_seq_len_dict[key]).logcdf(repeat_exon_genome_intersection_counts_dict[key]) #should not have -1?
    
    
    test_hyper_3_ratio= (repeat_exon_genome_intersection_counts_dict[key]/repeat_genome_seq_len_dict[key])/((exon_middle_exon_bases/2)/(genome_len))
    
    
    
    
    repeat_stat_ranking_dict[key].append(test_hyper_3) #sf
    repeat_stat_ranking_dict[key].append(test_hyper_3_ratio)  #enrichment ratio
    repeat_stat_ranking_dict[key].append(test_hyper_4) #cdf
    



avg_exon_exp_bases = rep_genome_stats_df['avg_exon_exp_bases']
avg_exon_exp = rep_genome_stats_df['avg_exon_exp']

e_genome_bases = rep_genome_stats_df['exon_bases'] 
r_genome_bases = rep_genome_stats_df['genome_bases']
genome_len=3088286401
for key in repeat_stat_rows_dict:
    entry = repeat_stat_rows_dict[key]
    
    r_repeat_bases = entry['rep_genome_bases']
    e_repeat_bases = entry['exon_rep_intersect_bases']
    e_exon_bases   = entry['exon_genome_bases']






with open(exp_output_path.Dfam_repeats + 'repeat_stat_ranking_dict.pickle','wb') as f:
    pickle.dump(repeat_stat_rows_dict, f)
    pickle.dump(repeat_stat_ranking_dict, f)
    pickle.dump(repeat_genome_counts_dict, f)
    pickle.dump(repeat_exon_bases_counts_dict, f)
    pickle.dump(repeat_genome_seq_len_dict, f)
    pickle.dump(repeat_exon_genome_intersection_counts_dict, f)
    pickle.dump(repeat_exon_counts_dict, f)
    pickle.dump(exon_middle_exon_bases, f)
    pickle.dump(repeat_exon_bases_counts_expression_dict, f)















#int('stop here now')










'''


from experiment_paths.experiment_paths import *exp_output_path as exp_output_path
import numpy as np
import matplotlib.pyplot as plt
import pickle
import exon_id_library.exon_id_lib as el
ouput_folder = exp_output_path.Dfam_repeats
with open(exp_output_path.Dfam_repeats + 'repeat_stat_ranking_dict.pickle','rb') as f:
    repeat_stat_rows_dict = pickle.load(f)
    repeat_stat_ranking_dict = pickle.load(f)
    repeat_genome_counts_dict = pickle.load(f)
    repeat_exon_bases_counts_dict = pickle.load(f)
    repeat_genome_seq_len_dict = pickle.load(f)
    repeat_exon_genome_intersection_counts_dict = pickle.load(f)
    repeat_exon_counts_dict = pickle.load(f)   
    exon_middle_exon_bases =  pickle.load(f)   
    repeat_exon_bases_counts_expression_dict =  pickle.load(f)   

'''









pvals_sf = [repeat_stat_ranking_dict[x][-3] for ii, x in enumerate(repeat_stat_ranking_dict)]
ratios = [repeat_stat_ranking_dict[x][-2] for ii, x in enumerate(repeat_stat_ranking_dict)]
pvals_cdf = [repeat_stat_ranking_dict[x][-1] for ii, x in enumerate(repeat_stat_ranking_dict)]

key_list = list()
pvals=list()
for ii, key in enumerate(repeat_stat_ranking_dict):
    ratio = ratios[ii]
    if ratio > 1:
        pvals.append(pvals_sf[ii])
    else:
        pvals.append(pvals_cdf[ii])
    key_list.append(key)

np.percentile(pvals,20)        







threshold_1 = int(-1*np.percentile(pvals,1))
threshold_5 = int(-1*np.percentile(pvals,15))
fold_threshold = 0 #log2 fold threshold
percent_threshold = threshold_5

blue = [[key, ratios[ii], -1*pvals[ii]] for ii, key in enumerate(repeat_stat_ranking_dict) if -1*pvals[ii] < percent_threshold or abs(np.log2(ratios[ii])) < fold_threshold]

#red corresponds to top 1% ro 5% pvals and a fold threshold
red = [[key, ratios[ii], -1*pvals[ii]] for ii, key in enumerate(repeat_stat_ranking_dict) if -1*pvals[ii] > percent_threshold and abs(np.log2(ratios[ii])) > fold_threshold]





ouput_folder = exp_output_path.Dfam_repeats
pdf_plots = PdfPages(ouput_folder+'pval_vulcanoid.pdf')





key_list, x,y = zip(*blue)
enrichment_data = list(np.log2(x))
key_list, x,y = zip(*red)
enrichment_data += list(np.log2(x))


fig, ax1=plt.subplots(1)

ax1_twin=ax1.twinx()
ax1.hist(enrichment_data, bins = np.arange(-6,7,0.25), color='gray', alpha = 0.5)

key_list, x,y = zip(*blue)
enrichment_data = list(np.log2(x))
ax1_twin.scatter(np.log2(x),y/np.log(10), marker='.')
key_list, x,y = zip(*red)
enrichment_data += list(np.log2(x))
ax1_twin.scatter(np.log2(x),y/np.log(10), marker='.', color = 'red')

plt.ylabel('-log10 pval ')
plt.xlim(-6,6)
plt.xlabel('fold enrichment (log2)')
plt.ylabel("-log10 pval")
plt.legend(['5% pval threshold: {:}\nlog2 fold-threshold: {:}'.format(threshold_5, fold_threshold)])
plt.title('hypergeometric pvals for exons bases overlapping repeat bases.')

ax1_twin.set_ylim(0,440000)
plt.tight_layout()


pdf_plots.savefig(fig)





pdf_plots.close()
























        

#DF0000317.4	 0.65	26978.0



enriched_keys = key_list












with open(exp_output_path.Dfam_repeats + 'repeat_genome_exon_pvals.pickle','wb') as f:
    pickle.dump( pvals_sf , f)
    pickle.dump( pvals_cdf , f)
    pickle.dump( ratios , f)
    pickle.dump( pvals , f)
    pickle.dump( key_list , f)
    pickle.dump( blue , f)
    pickle.dump( red , f)
    pickle.dump( threshold_1 , f)
    pickle.dump( threshold_5 , f)
    pickle.dump( fold_threshold , f)
    pickle.dump( percent_threshold , f)
    pickle.dump( enriched_keys , f)
    pickle.dump( repeat_associated_exon_ids_dict , f)


'''

with open(exp_output_path.Dfam_repeats + 'repeat_genome_exon_pvals.pickle','rb') as f:
    pvals_sf = pickle.load( f )
    pvals_cdf = pickle.load( f )
    ratios = pickle.load( f )
    pvals = pickle.load( f )
    key_list = pickle.load( f )
    blue = pickle.load( f ) #below thresholds
    red = pickle.load( f )  #above thresholds
    threshold_1 = pickle.load( f )
    threshold_5 = pickle.load( f )
    fold_threshold = pickle.load( f ) 
    percent_threshold = pickle.load( f )
    enriched_keys = pickle.load( f )
    repeat_associated_exon_ids_dict = pickle.load( f )
'''





































#repeat_associated_exon_ids_dict


repeat_exon_keys = set(repeat_associated_exon_ids_dict.keys())
pc_repeat_keys =repeat_exon_keys.intersection(pc_middle_exon_id_list)


#pc exon
ec = [[eid, aggregate_exon_dict[eid]['count'], repeat_associated_exon_ids_dict[eid]] for eid in pc_repeat_keys ]

ec = sorted(ec,key=lambda x:x[1])

#all exon
ea = [[eid, aggregate_exon_dict[eid]['count'], repeat_associated_exon_ids_dict[eid]] for eid in repeat_exon_keys ]

ea = sorted(ea,key=lambda x:x[1])






with open(exp_output_path.Dfam_repeats + 'ea_ec.pickle','wb') as f:
    pickle.dump(ea, f)
    pickle.dump(ec, f)
    

'''
with open(exp_output_path.Dfam_repeats + 'ea_ec.pickle','rb') as f:
    ea = pickle.load(f)
    ec = pickle.load(f)
'''













found_subfamily_exon_ids = dict()
for entry in ec:
    exon_id = entry[0]
    repeat_list=entry[2]
    for repeat in repeat_list:
        subfamily = repeat['family_name']
        if subfamily not in found_subfamily_exon_ids:
            found_subfamily_exon_ids[subfamily] = list()
        found_subfamily_exon_ids[subfamily].append(exon_id)





subfam_count_list = list()
for subfamily in found_subfamily_exon_ids:
    subfam_count_list.append([subfamily,len(found_subfamily_exon_ids[subfamily]), found_subfamily_exon_ids[subfamily]])


subfam_count_list = sorted(subfam_count_list,key=lambda x:x[1])





repeat_number_mRNA_exons_list = list()
repeat_number_mRNA_exons_counts_list = list()
for entry in subfam_count_list:
    repeat_number_mRNA_exons_list.append([entry[0], entry[1]])
    counts_list = list()
    for exon_id in entry[2]:
        counts_list.append(aggregate_exon_dict[exon_id]['count'])
    repeat_number_mRNA_exons_counts_list.append([entry[0],counts_list])



repeat_number_mRNA_exons_list=sorted(repeat_number_mRNA_exons_list, key=lambda x:x[1])











def plot_repeat_exon_hits(repeat_entry_list, repeat_name, repeat_consensus_bases, alpha = .85, bar_width_exon = 0.2 ):
    found_repeat_list=list()
    repeat_overlap_type_dict = {'inside':0,'overlap':0, 'exon_count':0,'repeat_all_genome_counts':0,'repeat_len':list()}
    
    for entry in repeat_entry_list:
        repeat_list=entry[2]
        
        for repeat in repeat_list:
            
            #if repeat['strand']=='-':
            #    continue
            
            if repeat['family_name'].find(repeat_name) >= 0:
                found_repeat_list.append([entry[0], entry[1],  repeat])
                
                repeat_overlap_type_dict['repeat_all_genome_counts'] = repeat_genome_counts_dict[repeat['family_acc']]
                
                repeat_overlap_type_dict['repeat_len'].append(repeat['end']-repeat['start'])
                
                exon_id = entry[0]
                repeat_bases = set(range(repeat['start'],repeat['end']))
                exon_bases = set(range(el.exon_id_values(exon_id).start, el.exon_id_values(exon_id).end))
                outside_exon_bases = exon_bases.difference(repeat_bases)
                
                
                
                if len(outside_exon_bases) > 0:
                    repeat_overlap_type_dict['overlap'] += 1
                else:
                    repeat_overlap_type_dict['inside']  += 1
                
                repeat_overlap_type_dict['exon_count'] += 1
              
    
    
    
    found_repeat_list = sorted(found_repeat_list, key=lambda x: (x[2]['hmm-st'], x[2]['hmm-en']))
    
    num_found = len(found_repeat_list)
    if num_found == 0:
        return -1
    
    
    total = repeat_overlap_type_dict['inside'] + repeat_overlap_type_dict['overlap']
    
    repeat_overlap_type_dict['inside'] = repeat_overlap_type_dict['inside']/total
    repeat_overlap_type_dict['overlap'] = repeat_overlap_type_dict['overlap']/total
    
    
    
    
    fig,ax = plt.subplots()
    currentAxis = plt.gca()
    
    for ii, entry in enumerate(found_repeat_list):
        exon_id, count, repeat = entry
         
                
        width_e, width_r, x_coor_e, x_coor_r = get_repeat_exon_box_coordinates(exon_id, repeat, repeat_consensus_bases)
        
        height, width = .4, width_r
        X_coord, Y_coord = x_coor_r, ii
        currentAxis.add_patch(Rectangle((X_coord - 0.1, Y_coord - 0.1), width, height, alpha=0.4, facecolor='gray'))
        
        height, width = bar_width_exon, width_e
        X_coord, Y_coord = x_coor_e , ii+.1
        currentAxis.add_patch(Rectangle((X_coord - 0.1, Y_coord - 0.1), width, height, alpha=alpha, facecolor='blue'))
        
    
    plt.ylim(-1, num_found)
    plt.xlim(-350, 
len(repeat_consensus_bases[repeat['family_acc']][0])-1650
+1)
    plt.title("{:}: {:}\n inside{:.2}, overlap{:.2}".format(repeat_name, repeat['family_acc'], repeat_overlap_type_dict['inside'], repeat_overlap_type_dict['overlap']))
    
    plt.tight_layout()
    
    
    
    return repeat_overlap_type_dict, fig
    














def get_repeat_exon_box_coordinates(exon_id, repeat, repeat_consensus_bases):
    repeat_consensus_len = len(repeat_consensus_bases[repeat['family_acc']][0])-2000
    ex = el.exon_id_values(exon_id)
    width_e = ex.length
    width_r = abs(repeat['end']-repeat['start'])
    if repeat['strand'] == '+':
        x_coor_e = repeat['hmm-st'] + (ex.start-repeat['start'])
        x_ccor_r = repeat['hmm-st']
        
    else:
        x_ccor_r = (-1*repeat['hmm-en'])   
        
        x_ccor_r = -1* x_ccor_r - abs(repeat['hmm-en']-repeat['hmm-st'])
        
        x_coor_e = repeat['hmm-st'] - (ex.start-repeat['end']) -ex.length#- abs(repeat['hmm-en']-repeat['hmm-st']) 
       
        
    
    return width_e, width_r, x_coor_e, x_ccor_r
















def add_repeat_true_start(repeat, repeat_consensus_bases):
    
    key = repeat['family_acc']
    
    repeat_consensus_len = len(repeat_consensus_bases[key][0])-2000
    
    if repeat['strand'] == '+':
        true_start = repeat['hmm-st']
        true_end   = repeat['hmm-en']
    else:
        true_start = repeat_consensus_len-repeat['hmm-en']
        true_end   = repeat_consensus_len-repeat['hmm-st']

    repeat['true_start'] = true_start
    repeat['true_end'] = true_end


add_repeat_true_start(repeat, repeat_consensus_bases)


for entry in ec:
    repeat_list=entry[2]
    for repeat in repeat_list:
        add_repeat_true_start(repeat, repeat_consensus_bases)
        








ouput_folder = exp_output_path.Dfam_repeats
pdf_plots = PdfPages(ouput_folder+'Dfam_repeat_individual_exon_plots.pdf')

    

mRNA_overlap = 0
mRNA_inside  = 0
all_overlap  = 0
all_inside   = 0
mRNA_repeat_subfamily_overlap_proportions  = list()
all_repeat_subfamily_overlap_proportions   = list()
mRNA_repeat_subfamily_repeat_frag_len_list = list()
all_repeat_subfamily_repeat_frag_len_list  = list()
repeat_enrichment = dict()

for key in found_subfamily_exon_ids:
    
    
    
    
    repeat_name       = key
    repeat_entry_list = ec
    result, fig = plot_repeat_exon_hits(repeat_entry_list, repeat_name, repeat_consensus_bases)
    
    
    if repeat_name == 'L1P2_5end':
        pdf_plots.savefig(fig)
    else:
        plt.close()    

    mRNA_repeat_subfamily_overlap_proportions.append(result['overlap'])
    mRNA_overlap += result['overlap'] * result['exon_count']
    mRNA_inside  += result['inside'] * result['exon_count']
    mRNA_repeat_subfamily_repeat_frag_len_list.append(result['repeat_len'])
    
    repeat_name       = key
    repeat_entry_list = ea
    result,fig = plot_repeat_exon_hits(repeat_entry_list, repeat_name, repeat_consensus_bases, 1, 0.8)
    
    if repeat_name == 'L1P2_5end':
        pdf_plots.savefig(fig)
    else:
        plt.close() 
    
    all_repeat_subfamily_overlap_proportions.append(result['overlap'])
    all_overlap += result['overlap'] * result['exon_count']
    all_inside  += result['inside'] * result['exon_count']
    all_repeat_subfamily_repeat_frag_len_list.append(result['repeat_len'])
    
    repeat_enrichment[key] = result['exon_count']/result['repeat_all_genome_counts']

    




entries_list = list()
for entry in ea:
    for entry_2 in entry[2]:
        if entry_2['family_name'] == 'L1P1_orf2':
            entries_list.append(entry_2)







#print file of repeats overlapping exons

    
found_ec_repeat_counts = list()
found_ec_repeats_list=list()
for entry in ec:
    repeat_list=entry[2]
    for repeat in repeat_list:
        if repeat['family_name'].find('L1P2_5end') >= 0:
            found_ec_repeats_list.append(repeat)
            found_ec_repeat_counts.append(aggregate_exon_dict[entry[0]]['count'])
            print("{:}:\t{:}\n{:}\t{:}\t{:}:{:}-{:}:{:}\t{:} bp\n".format(entry[0], entry[1],  repeat['family_name'], repeat['family_acc'], repeat['chrom'], repeat['start'], repeat['end'], repeat['strand'], abs(repeat['start']- repeat['end']) ))



    
found_ea_repeat_counts = list()
found_ea_repeats_list=list()
for entry in ea:
    repeat_list=entry[2]
    for repeat in repeat_list:
        if repeat['family_name'].find('L1P2_5end') >= 0:
            found_ea_repeats_list.append(repeat)
            found_ea_repeat_counts.append(aggregate_exon_dict[entry[0]]['count'])
            


























'''

for entry in ec:
    repeat_list=entry[2]
    for repeat in repeat_list:
        if repeat['family_acc'].find('DF0001008.4') >= 0:
            print("{:}:\t{:}\n{:}\t{:}\n".format(entry[0], entry[1],  repeat['family_name'], repeat['family_acc'] ))
            

'''






#requires loading the pval data and 
#ea for exon all, ec for exon coding
top_mRNA_repeats, counts = zip(*repeat_number_mRNA_exons_list[-25:])

dfam_ids_of_top_mRNAs = list()
for r_list in ec:
    if len(r_list[2]) > 0:
        repeat = r_list[2][0]
        repeat_name = repeat['family_name']        
        repeat_id = repeat['family_acc']  
        if repeat_name in top_mRNA_repeats:
            dfam_ids_of_top_mRNAs.append(repeat_id)
            
dfam_ids_of_top_mRNAs=list(set(dfam_ids_of_top_mRNAs)    )

repeat_ids_top_mRNA_highest_pvals = set(enriched_keys).intersection(dfam_ids_of_top_mRNAs)









pdf_plots.close()















