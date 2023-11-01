



#run first:
#    /home/dude/repositories/exon_def/slim_ET/ESE_count_vs_ESE.py







pdf_stats = PdfPages(exp_output_path.ESE_scan+'all_kmer_ESE_plot.pdf')






import pymannkendall as mk



def get_avg_count(exon_id_list, exon_dict):
    avg_counts_list = list()
    for exon_id in set(exon_id_list).intersection(exon_dict.keys()):
        exon = exon_dict[exon_id]
        avg_counts_list.append(exon['avg_count'])
    
    return avg_counts_list

def kmer_hash(kmer):
    kmer = kmer.upper()
    if kmer.find('N') >= 0:
        return -1
    
    base_to_num = {'A':0,'C':1,'G':2,'T':3}
    kmer_hash = 0 
    for ii in range(len(kmer)):
        kmer_hash = 4*kmer_hash
        kmer_hash += base_to_num[kmer[ii]]
    return kmer_hash

def num_hash_kmer(hash_num, kmer_len):
    kmer_str = ''
    num_to_base = {0:'A',1:'C',2:'G',3:'T'}
    ii=hash_num
    for jj in range(kmer_len):
        base = num_to_base[int(hash_num/4**(kmer_len-jj-1))%4]
        kmer_str = '%s%s' % (kmer_str, base)
    return kmer_str
    
kmer_num = kmer_hash('AGACAC')
kmer = num_hash_kmer(kmer_num, 6)
kmer_hash(kmer)









if 'avg_count' not in aggregate_exon_dict[list(aggregate_exon_dict.keys())[0]]:
    
    for exon_id in aggregate_exon_dict:
        exon = aggregate_exon_dict[exon_id]
        denom = sum(exon['lib_array']>5)
        exon['avg_count'] = exon['count']/denom











































pickle_path = exp_output_path.Dfam_pickle + "Dfam_pickle.pickle"
with open(pickle_path, "rb") as input_file:
    exon_ids_overlapping_repeat_list=pickle.load(input_file)








def get_exons_splice_score_window(exon_id_list, aggregate_exon_dict, score_range):
    
    MES_range_5ss = score_range
    MES_range_3ss = score_range
    
    exon_id_MES_range_list = list()
    
    for exon_id in exon_id_list:
        exon = aggregate_exon_dict[exon_id]
        if exon['5ss_score'] <= MES_range_5ss[0] or exon['5ss_score'] >= MES_range_5ss[1]:
            continue
        if exon['3ss_score'] <= MES_range_3ss[0] or exon['3ss_score'] >= MES_range_3ss[1]:
            continue
        exon_id_MES_range_list.append(exon_id)
    
    return exon_id_MES_range_list










#It seems MES scores correlate with abundance of some motifs.. for example, 'TTTTTT' is about 2x more common when splice site scores are above 10 and trends less common as splice site scores decrease, so this one is likely better tolerated when splice sites are good, or there is a nonlinear support offered by good splice sites.. which seems unlikely as the best MES scores have lower read counts. This appears in the other novel exons, but is less clear for mRNA. Is present in lncRNA

def get_kmer_stats(initial_exon_id_list, search_motif):
    
    
    count_MK_test={True:0,False:0}
    
    MES_window_list = [[0,4],[4,6],[6,8],[8,10],[10,15]]

    
    motif_ratio_dict_dict = dict()
    motif_ratio_dict_dict['ratio_read_count'] = dict()
    motif_ratio_dict_dict['normalized_count_fraction'] = dict()
    motif_ratio_dict_dict['read_count_ratio'] = dict()
    
    for MES_window in MES_window_list:
        key_str = "[{:}, {:}]".format(MES_window[0],MES_window[1])
        #
        motif_exon_id_dict = dict()
        motif_exon_id_dict[search_motif] = {'main':list(),'other':list()}
        exon_id_list = initial_exon_id_list
        #exon_id_list = primary_3ss_exon_id_set.intersection(pc_middle_exon_id_list)
        #exon_id_list = primary_3ss_exon_id_set.intersection(no_overlap_set)
        #exon_id_list = primary_3ss_exon_id_set.intersection(lncRNA_middle_exon_ids)
        #exon_id_list = primary_3ss_exon_id_set.intersection(intron_interior_set)
        #exon_id_list = primary_3ss_exon_id_set.intersection(antisense_transcript_set)
        
        #exon_id_list=el.threshold_exon_ids(exon_id_list, 1000, aggregate_exon_dict)
        
        exon_id_list = set(exon_id_list).difference(exon_ids_overlapping_repeat_list)
        exon_id_list = get_exons_splice_score_window(exon_id_list, aggregate_exon_dict, MES_window)
        for exon_id in exon_id_list:
            seq = el.get_seq_from_exon_id(exon_id,genome_fasta)
            seq = str(seq).upper()
            if seq.find('N') >= 0:
                continue
            
            if seq.find(search_motif) >= 0:
                motif_exon_id_dict[search_motif]['main'].append(exon_id)
            else:
                motif_exon_id_dict[search_motif]['other'].append(exon_id)
        
        
        
        
        
        scores_main_other_dict = dict()
        
        #count
        scores_main_other_dict['main'] = el.get_exon_dict_counts(motif_exon_id_dict[search_motif]['main'], aggregate_exon_dict)
        #count
        scores_main_other_dict['other'] = el.get_exon_dict_counts(motif_exon_id_dict[search_motif]['other'], aggregate_exon_dict)
        
        #avg_count
        scores_main_other_dict['main'] = get_avg_count(motif_exon_id_dict[search_motif]['main'], aggregate_exon_dict)
        #avg_count
        scores_main_other_dict['other'] = get_avg_count(motif_exon_id_dict[search_motif]['other'], aggregate_exon_dict)
    
        
    
        ratio_read_count = np.median(scores_main_other_dict['main'])*( np.median(scores_main_other_dict['main'])/np.median(scores_main_other_dict['other']))
    
        normalized_count_fraction = ( np.median(scores_main_other_dict['main'])/np.median(scores_main_other_dict['other']))

        read_count_ratio = len(scores_main_other_dict['main'])/(len(scores_main_other_dict['main'])+len(scores_main_other_dict['other']))
        
        motif_ratio_dict_dict['normalized_count_fraction'][key_str] = normalized_count_fraction
        motif_ratio_dict_dict['ratio_read_count'][key_str] = ratio_read_count
        motif_ratio_dict_dict['read_count_ratio'][key_str] = read_count_ratio
    
    data_items = motif_ratio_dict_dict['read_count_ratio'].items()
    key_list, val_list = zip(*data_items)

    
    trend_true = 'not tested'
    trend = 'N/A'

    if val_list[2] >-1:
        
        data_items = motif_ratio_dict_dict['normalized_count_fraction'].items()
        key_list, val_list = zip(*data_items)
        medians_read_count=val_list
        mk_result = mk.original_test(val_list)
        mk_read_count = mk_result

            
        data_items = motif_ratio_dict_dict['ratio_read_count'].items()
        key_list, val_list = zip(*data_items)
        
        

        
        data_items = motif_ratio_dict_dict['read_count_ratio'].items()
        key_list, val_list = zip(*data_items)
        medians_fraction=val_list
        mk_result = mk.original_test(val_list)
        mk_fraction = mk_result
        count_MK_test[mk_result.h] += 1
        trend = mk_result.trend
        trend_true = mk_result.h

    return trend_true, trend, mk_fraction, mk_read_count, medians_fraction, medians_read_count, key_list






def pval_pair(kmer_mk_results, kmer):
    return [kmer_mk_results[kmer]['mk_fraction'].p, kmer_mk_results[kmer]['mk_read_count'].p]


def thresh_pvals(kmer_mk_results, kmer, thresh):
    return kmer_mk_results[kmer]['mk_fraction'].p <= thresh # and kmer_mk_results[kmer]['mk_read_count'].p <= thresh



def plot_kmer_medians(kmer, kmer_medians_list, kmer_mk_results, *args, **kwargs):
    
    write_to_pdf_flag = 'pdf_handle' in kwargs
    'pdf_handle'
    
    if 'exon_category' in kwargs:
        exon_category = '\n'+kwargs['exon_category']
    else:
        exon_category=''
    
    medians_list = kmer_medians_list[kmer]
    
    mk_result=kmer_mk_results[kmer]['mk_fraction']
    val_list = medians_list['fractions']
    fig = plt.figure()
    plt.bar(range(len(val_list)), val_list, color='red')
    plt.xticks(range(len(key_list)), key_list)
    plt.title('fraction of exons with kmer: {:}\nMann_kendall test \ntrend={:}, test={:},pval={:.2},z={:.2}'.format(kmer, mk_result.trend, mk_result.h, mk_result.p, mk_result.z+0.0 ))
    plt.ylabel('fraction of exon that contain kmer')
    plt.tight_layout()
    if write_to_pdf_flag == True:
        kwargs['pdf_handle'].savefig(fig)
        
    return val_list, key_list
        
    








search_motif = 'GAAGAA' 


exon_id_list = primary_3ss_exon_id_set.intersection(no_overlap_set)

count_MK_test = {True:0,False:0,'not tested':0}
trend_MK_test = {'increasing':0,'decreasing':0,'no trend':0,'N/A':0,}
kmer_with_trend_dict=dict()
kmer_with_trend_dict['increasing'] = list()
kmer_with_trend_dict['decreasing'] = list()
kmer_with_trend_dict['any'] = list()
kmer_len=6
kmer_mk_results = dict()
kmer_medians_list = dict()

for ii in range(0,4**kmer_len):

    if  kmer_hash('GAAGAA' )  != ii:
        continue

    if ii % int(4**kmer_len/5) == 0:
        print('processed {:,}/{:,}'.format(ii, 4**kmer_len))
    search_kmer = num_hash_kmer(ii,kmer_len)
    trend_true, trend, mk_fraction, mk_read_count,medians_fraction, medians_read_count, key_list = get_kmer_stats(exon_id_list, search_kmer)
    count_MK_test[trend_true] += 1
    trend_MK_test[trend] += 1
    if trend_true == True:
        kmer_with_trend_dict['any'].append(search_kmer)
        kmer_with_trend_dict[trend].append(search_kmer)
    kmer_mk_results[search_kmer] = {'mk_fraction':mk_fraction, 'mk_read_count':mk_read_count}
    kmer_medians_list[search_kmer]={'read_counts':medians_read_count,'fractions':medians_fraction,'key_list':key_list}

    



pdf_freq_stats = PdfPages(exp_output_path.ESE_scan+'count_freq_examples_paper.pdf')

search_motif = 'GAAGAA'
val_list, key_list = plot_kmer_medians(search_motif, kmer_medians_list, kmer_mk_results, pdf_handle = pdf_freq_stats)

pdf_freq_stats.close()





with open(exp_output_path.ESE_scan+'4F.txt', 'w') as f:
    f.write('MaxEntScan_range\tfraction\n')
    for ii, val in enumerate(val_list):
        f.write('{:}\t{:}\n'.format(key_list[ii], val_list[ii]) )
    
    
    
    
    
                     









