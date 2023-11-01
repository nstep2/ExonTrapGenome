"""

This code requires the full dataset depth (1 read thresholdq)





'''
00B9F1 	lncRNA
00A875	Antisense
0072BC	mRNA
F15A22	Intergenic
DA6FAB	Intronic
'''


"""


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


pdf_stats = PdfPages(exp_output_path.rewrite_exon_number_expectation+'rewrite_exon_number_expectation_%d_reads.pdf' % (exon_count_build))

#https://davidmathlogic.com/colorblind/#%23D81B60-%231E88E5-%23FFC107-%23004D40









def exon_ids_in_exp_bin(exon_id_list,aggregate_exon_dict, **kwargs):
    exon_id_list = el.exon_id_intersection(exon_id_list, aggregate_exon_dict.keys())
    
    #exp_bin_list = [0]
    #exp_bin_list = [int(10**x) for x in np.arange(0,6,0.5)]
    exp_bin_list = [int(10**x) for x in np.arange(0,6,0.1)]
    if 'exp_bin_list' in kwargs:
        exp_bin_list = kwargs['exp_bin_list']
    #exp_bin_list =  [int(2**x) for x in np.arange(0,17+0.1,.5)]
    #exp_bin_list=np.arange(0,100000,10)
    
    exon_id_exp_bin_dict = dict()
    
    for nn, val in enumerate(exp_bin_list):
        if nn +1 == len(exp_bin_list): 
            continue
        #key = '%d-%d' % (exp_bin_list[nn],exp_bin_list[nn+1])
        key = exp_bin_list[nn]
        exon_id_exp_bin_dict[key] = list()
        for ii, exon_id in enumerate(exon_id_list):
            count = aggregate_exon_dict[exon_id]['count']
            if count >= exp_bin_list[nn] and count < exp_bin_list[nn+1]:
                exon_id_exp_bin_dict[key].append(exon_id)
        
    return exon_id_exp_bin_dict
                
        
    
    
    
    
    














pair_list = list()

pair_list.append([pc_middle_exon_id_list, 'mRNA'])
pair_list.append([lncRNA_middle_exon_ids, 'lncRNA'])
pair_list.append([antisense_transcript_set, 'Antisense'])
pair_list.append([no_overlap_set, 'Intergenic'])
pair_list.append([intron_interior_set, 'Intronic'])


exon_annotation_bins_dict = dict()
exon_annotation_bins_low_res_dict = dict()
exon_annotation_bins_linear_dict = dict()
for pair in pair_list:
    exon_id_list = pair[0]
    exon_annotation_bins_dict[pair[1]] = exon_ids_in_exp_bin(pair[0],aggregate_exon_dict)
    exp_bin_list = [int(10**x) for x in np.arange(0, 6.0, 0.5)]
    exon_annotation_bins_low_res_dict[pair[1]] = exon_ids_in_exp_bin(pair[0],aggregate_exon_dict, exp_bin_list= exp_bin_list)
    exp_bin_list = np.arange(100,100000,1000)
    exon_annotation_bins_linear_dict[pair[1]] = exon_ids_in_exp_bin(pair[0],aggregate_exon_dict, exp_bin_list= exp_bin_list)





def make_stats_3ss_dict(exon_annotation_bins_dict, aggregate_exon_dict):
        
    stats_3ss_dict = dict()
    for exon_id_list_name in exon_annotation_bins_dict:
        exon_sets_dict=exon_annotation_bins_dict[exon_id_list_name]
        
        score_list = list()
        score_list_stats = list()
        counts_list = list()
        for key in sorted(list(exon_sets_dict.keys())):
            
            exon_id_list = exon_sets_dict[key]
            if len(exon_id_list) == 0:
                #score_list.append([key,[]])
                #counts_list.append([key, []])
                score_list.append([key,0])
                counts_list.append([key, 0])
                continue
            score_list_3ss = el.get_exon_dict_3ss_scores(exon_id_list,aggregate_exon_dict)
            score_list.append([key,score_list_3ss])
            score_list_stats.append([key,np.median(score_list_3ss),np.percentile(score_list_3ss,25),np.percentile(score_list_3ss,75)])
            counts_list.append([key, len(score_list_3ss)])
            
        
        
        stats_3ss_dict[exon_id_list_name] = [score_list_stats,counts_list]
    
    
    return stats_3ss_dict

stats_3ss_dict = make_stats_3ss_dict(exon_annotation_bins_dict, aggregate_exon_dict)

stats_3ss_low_res_dict = make_stats_3ss_dict(exon_annotation_bins_low_res_dict, aggregate_exon_dict)

stats_3ss_linear_dict = make_stats_3ss_dict(exon_annotation_bins_linear_dict, aggregate_exon_dict)





stats_5ss_dict = dict()
for exon_id_list_name in exon_annotation_bins_dict:
    exon_sets_dict=exon_annotation_bins_dict[exon_id_list_name]
    
    score_list = list()
    score_list_stats = list()
    counts_list = list()
    for key in sorted(list(exon_sets_dict.keys())):
        exon_id_list = exon_sets_dict[key]
        if len(exon_id_list) == 0:
            continue
        score_list_5ss = el.get_exon_dict_5ss_scores(exon_id_list,aggregate_exon_dict)
        score_list.append([key,score_list_5ss])
        score_list_stats.append([key,np.median(score_list_5ss),np.percentile(score_list_5ss,25),np.percentile(score_list_5ss,75)])
        counts_list.append([key, len(score_list_5ss)])
    
    
    stats_5ss_dict[exon_id_list_name] = [score_list_stats,counts_list]
    
    
    
    



#Stacked plots
fig = plt.figure()
key_list = list()
y_list = list()
color_order_list = list()
max_x = [0]
for ii, exon_id_list_name in enumerate(stats_3ss_dict):
    key_list.append(exon_id_list_name)
    score_list_stats,counts_list = stats_3ss_dict[exon_id_list_name]
    
    x,y = zip(*counts_list)
    if len(x) > len(max_x):
        max_x = x
    y_list.append(y[:47])  #length depends on whether there are exons at the higher counts
    color_order_list.append(color_dict[exon_id_list_name])


# Pad the y lists to make them of equal length
max_length = max(len(yi) for yi in y_list)
padded_y_list = [np.pad(yi, (0, max_length - len(yi)), mode='constant') for yi in y_list]

# Convert the padded y list into a NumPy array
y_stack = np.array(padded_y_list, dtype=float)


x_stack = np.array(np.log10(x), dtype=float)
y_stack = np.array([np.array(yi) for ii, yi in enumerate(y_list)], dtype=float)


import numpy as np
import matplotlib.pyplot as plt





fig = plt.figure()
plt.stackplot(x_stack[:len(y_stack[0])], y_stack,colors=color_order_list, labels=key_list)
plt.ylim(0,20000)
plt.xticks(np.log10([10,100,1000,10000,100000]),['{:}'.format(k) for ii, k in enumerate([10,100,1000,10000,100000]) ])
plt.xlim(4,5)
plt.legend(reversed(plt.legend().legendHandles), reversed(key_list))







fig, ax = plt.subplots()
key_list = ['mRNA', 'lncRNA','Intergenic']


score_list_stats,counts_list = stats_3ss_dict['mRNA']
x,y,q1,q2 = zip(*score_list_stats)
plt.plot(np.log10(x),y, color = '#FFC20A', marker ='D')
ax.fill_between(np.log10(x), (q1), (q2), color = '#FFC20A', alpha=.2)

score_list_stats,counts_list = stats_3ss_dict['lncRNA']
x,y,q1,q2 = zip(*score_list_stats)
plt.plot(np.log10(x),y, color = '#0C7BDC', linestyle='--')
ax.fill_between(np.log10(x), (q1), (q2), color = '#0C7BDC', alpha=.1)

score_list_stats,counts_list = stats_3ss_dict['Intergenic']
x,y,q1,q2 = zip(*score_list_stats)
plt.plot(np.log10(x),y, color = 'r', linestyle='solid')
ax.fill_between(np.log10(x), (q1), (q2), color = 'r', alpha=.1)


plt.xticks(np.log10([10,100,1000,10000,100000]),["{:,}".format(k) for ii, k in enumerate([10,100,1000,10000,100000]) ])

plt.xlim(0,5)
plt.legend(key_list)
plt.ylabel('3ss scores (maxentscan)')
plt.ylim(0,12)
plt.tight_layout()
pdf_stats.savefig(fig)

plt.xlim(2,5)

plt.tight_layout()
pdf_stats.savefig(fig)





fig, ax = plt.subplots()
key_list = ['mRNA', 'lncRNA','Intergenic']

score_list_stats,counts_list = stats_5ss_dict['mRNA']
x,y,q1,q2 = zip(*score_list_stats)
plt.plot(np.log10(x),y, color = '#FFC20A', linestyle='solid', marker ='D')
ax.fill_between(np.log10(x), (q1), (q2), color = '#FFC20A', alpha=.1)

score_list_stats,counts_list = stats_5ss_dict['lncRNA']
x,y,q1,q2 = zip(*score_list_stats)
plt.plot(np.log10(x),y, color = '#0C7BDC', linestyle='--')
ax.fill_between(np.log10(x), (q1), (q2), color = '#0C7BDC', alpha=.1)

score_list_stats,counts_list = stats_5ss_dict['Intergenic']
x,y,q1,q2 = zip(*score_list_stats)
plt.plot(np.log10(x),y, color = 'r', linestyle='solid')
ax.fill_between(np.log10(x), (q1), (q2), color = 'r', alpha=.1)


plt.xticks(np.log10([10,100,1000,10000,100000]),['{:}'.format(k) for ii, k in enumerate([10,100,1000,10000,100000]) ])
plt.xlim(0,5)
plt.legend(key_list)
plt.ylim(0,12)
plt.ylabel('5ss scores (maxentscan)')
plt.tight_layout()
pdf_stats.savefig(fig)













fig = plt.figure()
key_list = list()
for ii, exon_id_list_name in enumerate(stats_3ss_low_res_dict):
    key_list.append(exon_id_list_name)
    score_list_stats,counts_list = stats_3ss_low_res_dict[exon_id_list_name]
    
    x,y = zip(*counts_list)
    plt.plot(np.log10(x),y, color = color_dict[exon_id_list_name], marker='o',markersize=1.5)






plt.xticks(np.log10([1,10,100,1000,10000,100000]),['{:}'.format(k) for ii, k in enumerate([1, 10,100,1000,10000,100000]) ])
plt.tight_layout()
plt.ylabel('count')
plt.xlabel('log axis')
plt.legend(key_list)
plt.xlim(np.log10(1),5)
plt.tight_layout()
pdf_stats.savefig(fig)









pdf_stats.close()















