

stop_codons = ['TAG', 'TAA','TGA']






outdir = exp_output_path.library_specifics 
pdf_plots = PdfPages(outdir+'library_specifics.pdf')
pdf_scatter_plots = PdfPages(outdir+'library_specifics_scatter.pdf')


def insert_commas(value):
    return f'{value:,}'

import matplotlib.pyplot as plt







library_codes = dict()
code_to_libraries = dict()

library_codes[1]='1'
library_codes[2]='1'
library_codes[3]='1'
library_codes[4]='1'
library_codes[5]='1'
code_to_libraries['1'] = [1,2,3,4,5]

library_codes[6]='2'
library_codes[7]='2'
library_codes[21]='2'
library_codes[22]='2'
code_to_libraries['2'] = [6,7,21,22]

library_codes[8]='3'
library_codes[9]='3'
library_codes[10]='3'
library_codes[11]='3'
library_codes[23]='3'
code_to_libraries['3'] = [8,9,10,11,23]

library_codes[12]='4'
library_codes[13]='4'
library_codes[14]='4'
library_codes[15]='4'
code_to_libraries['4'] = [12,13,14,15]

library_codes[16]='5'
library_codes[17]='5'
library_codes[18]='5'
library_codes[19]='5'
library_codes[20]='5'
code_to_libraries['5'] = [16,17,18,19,20]

backbone_counts = dict()
backbone_counts['1']=0
backbone_counts['2']=0
backbone_counts['3']=0
backbone_counts['4']=0
backbone_counts['5']=0

backbone_specific_counts = dict()
backbone_specific_counts['1']=0
backbone_specific_counts['2']=0
backbone_specific_counts['3']=0
backbone_specific_counts['4']=0
backbone_specific_counts['5']=0

backbone_unique_exon_ids = dict()
backbone_unique_exon_ids['1']=list()
backbone_unique_exon_ids['2']=list()
backbone_unique_exon_ids['3']=list()
backbone_unique_exon_ids['4']=list()
backbone_unique_exon_ids['5']=list()

for exon_id in aggregate_exon_dict:
    exon=aggregate_exon_dict[exon_id]
    for lib_id in exon['lib_count']:
        count = exon['lib_count'][lib_id]
        backbone_counts[library_codes[lib_id]] += count
        
    for backbone in backbone_counts:
        backbone_count = 0
        non_backbone_count = 0
        for lib_id in exon['lib_count']:
            count = exon['lib_count'][lib_id]
            if backbone == library_codes[lib_id]:
                backbone_count += count
            else:
                non_backbone_count += count
        
        if backbone_count > 0 and non_backbone_count == 0:
            backbone_unique_exon_ids[backbone].append(exon_id)
            backbone_specific_counts[backbone] += backbone_count




def above_threshold_exon_id_list(threshold, exon_id_set, aggregate_exon_dict):
    count_above_threshold = 0
    exon_ids_passing_threshold=list()
    for exon_id in exon_id_set:
        if aggregate_exon_dict[exon_id]['count'] >= threshold:
            count_above_threshold += 1
            exon_ids_passing_threshold.append(exon_id)
            #print(exon_id)
    return exon_ids_passing_threshold


backbone_exon_ids = dict()
backbone_exon_ids['1']=0
backbone_exon_ids['2']=0
backbone_exon_ids['3']=0
backbone_exon_ids['4']=0
backbone_exon_ids['5']=0

backbone_counts_list_dict = dict()
backbone_counts_list_list = list()

threshold = 100
backbone_index_list = list()
for backbone_id in backbone_unique_exon_ids:
    
    exon_id_set = backbone_unique_exon_ids[backbone_id]
    threshold_exon_ids = above_threshold_exon_id_list(threshold, exon_id_set, aggregate_exon_dict)
    count = len(threshold_exon_ids)
    
    backbone_exon_ids[backbone_id] = threshold_exon_ids
    
    print("backbone unique exon_ids for %s:" % backbone_id, insert_commas(count))
    
    #count_list = el.get_exon_dict_counts(threshold_exon_ids, aggregate_exon_dict)
    
    count_list = el.get_exon_dict_counts(exon_id_set, aggregate_exon_dict)
    backbone_counts_list_dict[backbone_id]=count_list
    backbone_counts_list_list.append(count_list)
    backbone_index_list.append(backbone_id)

for backbone_id in backbone_unique_exon_ids:
    count = backbone_specific_counts[backbone_id]
    print("backbone specific counts for %s:" % backbone_id, insert_commas(count))






for backbone_id in backbone_unique_exon_ids:
    scores = el.get_exon_dict_3ss_scores(el.threshold_exon_ids(backbone_unique_exon_ids[backbone_id],100,aggregate_exon_dict),aggregate_exon_dict)
    print("exons above 100 and backbone specific 3ss score median for %s: %.1f" % (backbone_id, np.median(scores)))

for backbone_id in backbone_unique_exon_ids:
    scores = el.get_exon_dict_5ss_scores(el.threshold_exon_ids(backbone_unique_exon_ids[backbone_id],100,aggregate_exon_dict),aggregate_exon_dict)
    print("exons above 100 and backbone specific 5ss score median for %s: %.1f" % (backbone_id, np.median(scores)))






#backbone_counts
for backbone_id in backbone_counts:
    count = backbone_counts[backbone_id]
    print("backbone total read counts for %s:" % backbone_id, insert_commas(count))



exon_library_scatter_dict = dict()
for backbone in {'1':0,'2':0,'3':0,'4':0,'5':0}:
    exon_library_scatter_dict[backbone] = list()

exon_id_list = el.exon_id_intersection(pc_middle_exon_id_list, aggregate_exon_dict.keys())
for ii, exon_id in enumerate(aggregate_exon_dict):

    exon = aggregate_exon_dict[exon_id]
    exon['backbone_counts'] = {'1':0,'2':0,'3':0,'4':0,'5':0}
    
    count_dict = dict()
    for backbone in exon['backbone_counts']:

        count = 0
        for code in code_to_libraries[backbone]:
            if code in exon['lib_count']:
                count += exon['lib_count'][code]
        count = count/len(code_to_libraries[backbone])
        exon['backbone_counts'][backbone] = count
        count_dict[backbone] = count
    
    threshold_all_backbones = 100
    above_threshold_all_backbones = True

    for backbone in count_dict:
        if count_dict[backbone] >= threshold_all_backbones:
            1
            #above_threshold_all_backbones = True
        else:
            above_threshold_all_backbones = False
            1
    if above_threshold_all_backbones == True:
        for backbone in count_dict:
            exon_library_scatter_dict[backbone].append(count_dict[backbone] )





min_reads_threshold = 9
exon_id_all_backbones = list()
for exon_id in aggregate_exon_dict:
    exon = aggregate_exon_dict[exon_id]
    backbone_count_number = 0
    for backbone in exon['backbone_counts']:
        if exon['backbone_counts'][backbone] > min_reads_threshold:
            backbone_count_number += 1
    if backbone_count_number == 5:
        exon_id_all_backbones.append(exon_id)




a = len(set(exon_id_all_backbones).intersection(pc_middle_exon_id_list))/len(exon_id_all_backbones)

b = len(set(exon_id_all_backbones).intersection(lncRNA_middle_exon_ids))/len(exon_id_all_backbones)

print('count exons in all backbones: %s' % (insert_commas(len(exon_id_all_backbones))))

print('exons in all libraries - percent that are annotated: %.1f%% -- %s' % (100*a,insert_commas(int(a*len(exon_id_all_backbones)))))

print('exons in all libraries - percent that are lncRNA: %.1f%% -- %s' % (100*b,insert_commas(int(b*len(exon_id_all_backbones)))))





exon_id_list = exon_id_all_backbones 
a = el.get_exon_dict_counts(exon_id_list, aggregate_exon_dict)


exon_id_list = set(exon_id_all_backbones).intersection(pc_middle_exon_id_list)
b = el.get_exon_dict_counts(exon_id_list, aggregate_exon_dict)


exon_id_list = set(exon_id_all_backbones).intersection(lncRNA_middle_exon_ids)
c = el.get_exon_dict_counts(exon_id_list, aggregate_exon_dict)






exon_id_list = el.exon_id_intersection(exon_id_all_backbones, pc_middle_exon_id_list)

mod_3_exon_id_codon_counts=dict()
mod_3_exon_id_codon_counts_seq=dict()
for exon_id in exon_id_list:
    exon = aggregate_exon_dict[exon_id]
    #if exon['count'] < 1000:
    #    continue
    
    seq = el.get_seq_from_exon_id(exon_id,genome_fasta)
    mod_3_exon_id_codon_counts_seq[exon_id]=seq
    if len(seq) % 3 != 0:
        continue
    mod_3_exon_id_codon_counts[exon_id]={0:0,1:0,2:0,3:0}
    
    
    seq = seq.upper()
    for ii in range(0,len(seq),3):
        codon = seq[ii:ii+3]
        #print(codon)
        if codon in stop_codons:
            mod_3_exon_id_codon_counts[exon_id][0] += 1
            mod_3_exon_id_codon_counts[exon_id][3] += 1

    for ii in range(1,len(seq),3): # T4 is -1 frame
        codon = seq[ii:ii+3]
        #print(ii,codon)
        if codon in stop_codons:
            mod_3_exon_id_codon_counts[exon_id][2] += 1
    
    for ii in range(2,len(seq),3): # T3 is -2 frame
        codon = seq[ii:ii+3]
        #print(ii,codon)
        if codon in stop_codons:
            mod_3_exon_id_codon_counts[exon_id][1] += 1



with open('2B','w') as f:
    f.write('chrom\tstart\tend\strand\exon_id\tframe 0 count\tframe 1 count\tframe 2 count\n')
    for exon_id in mod_3_exon_id_codon_counts:
        ex = el.exon_id_values(exon_id)
        f.write(ex.chrom+'\t')
        f.write(str(ex.start)+'\t')
        f.write(str(ex.end)+'\t')
        f.write(ex.strand+'\t')
        f.write(exon_id+'\t')
        f.write(str(mod_3_exon_id_codon_counts[exon_id][0])+'\t')
        f.write(str(mod_3_exon_id_codon_counts[exon_id][1]) +'\t')
        f.write(str(mod_3_exon_id_codon_counts[exon_id][2]) +'\t')
        
 
    
 
    
 
    
 
    
 
    
 
    
#1,2,3
frame_to_backbone_num = {0:2,1:0,2:1}
in_frame = {'3':list(),'1':list(),'2':list()}
out_frame = {'3':list(),'1':list(),'2':list()}
out_frame_above_4 = {'3':list(),'1':list(),'2':list()}

scatter_no_stop_vs_stop_list_0_1 = list()
scatter_no_stop_vs_stop_list_0_2 = list()


scatter_no_stop_vs_stop_list_1_0 = list()
scatter_no_stop_vs_stop_list_1_2 = list()


scatter_no_stop_vs_stop_list_2_0 = list()
scatter_no_stop_vs_stop_list_2_1 = list()


for ii, exon_id in enumerate(mod_3_exon_id_codon_counts):
    exon=aggregate_exon_dict[exon_id]
    bb_count = exon['backbone_counts']
    
    #check T5 backbone
    if mod_3_exon_id_codon_counts[exon_id][3] == 0:
        in_frame['3'].append(bb_count['3']) #3 is T5
        #for jj in ['1','2']:
    if mod_3_exon_id_codon_counts[exon_id][3] > 0:
        out_frame['3'].append(bb_count['3']) #3 is T5
    if mod_3_exon_id_codon_counts[exon_id][3] > 4:
        out_frame_above_4['3'].append(bb_count['3']) #3 is T5
    
    if mod_3_exon_id_codon_counts[exon_id][3] == 0:
        if mod_3_exon_id_codon_counts[exon_id][1] > 0:
            pair = [bb_count['3'], bb_count['1']]
            scatter_no_stop_vs_stop_list_0_1.append(pair)
    if mod_3_exon_id_codon_counts[exon_id][3] == 0:
        if mod_3_exon_id_codon_counts[exon_id][2] > 0:
            pair = [bb_count['3'], bb_count['2']]
            scatter_no_stop_vs_stop_list_0_2.append(pair)
    
        
    #check T3 backbone
    if mod_3_exon_id_codon_counts[exon_id][1] == 0:
        in_frame['1'].append(bb_count['1'])
    if mod_3_exon_id_codon_counts[exon_id][1] > 0:
        out_frame['1'].append(bb_count['1'])
    if mod_3_exon_id_codon_counts[exon_id][1] > 4:
        out_frame_above_4['1'].append(bb_count['1'])
    
    if mod_3_exon_id_codon_counts[exon_id][1] == 0:
        if mod_3_exon_id_codon_counts[exon_id][2] > 0:
            pair = [bb_count['1'], bb_count['2']]
            scatter_no_stop_vs_stop_list_1_2.append(pair)
            
    if mod_3_exon_id_codon_counts[exon_id][1] == 0:
        if mod_3_exon_id_codon_counts[exon_id][3] > 0:
            pair = [bb_count['1'], bb_count['3']]
            scatter_no_stop_vs_stop_list_1_0.append(pair)

        
    
    #check T4 backbone
    if mod_3_exon_id_codon_counts[exon_id][2] == 0:
        in_frame['2'].append(bb_count['2'])
    if mod_3_exon_id_codon_counts[exon_id][2] > 0:
        out_frame['2'].append(bb_count['2'])
    if mod_3_exon_id_codon_counts[exon_id][2] > 4:
        out_frame_above_4['2'].append(bb_count['2'])

    if mod_3_exon_id_codon_counts[exon_id][2] == 0:
        if mod_3_exon_id_codon_counts[exon_id][1] > 0:
            pair = [bb_count['2'], bb_count['1']]
            scatter_no_stop_vs_stop_list_2_1.append(pair)
            
    if mod_3_exon_id_codon_counts[exon_id][2] == 0:
        if mod_3_exon_id_codon_counts[exon_id][3] > 0:
            #print(exon_id)
            pair = [bb_count['2'], bb_count['3']]
            scatter_no_stop_vs_stop_list_2_0.append(pair)        
       

import scipy
import scipy.stats

color_val = ['tab:blue','tab:orange','tab:green']
fig = plt.figure()
for ii, index_val in enumerate(['3','1','2']):
    a = in_frame[index_val]
    b = out_frame[index_val]
    c = out_frame_above_4[index_val]
    

    statstic_b,pval_b = scipy.stats.ranksums(a,b)
    statstic_c,pval_c = scipy.stats.ranksums(a,c)

    
    print('Wilcoxon rank-sum backcbone %s ' % (index_val))
    print('pval a-b: ', (pval_b))

    width  = 0.3
    
    plt.bar([ii],[np.median(a)], yerr=[scipy.stats.sem(a)], capsize=3, width=width, color=color_val[ii], alpha=1)
    plt.bar([ii+width*1.25],[np.median(b)], yerr=[scipy.stats.sem(b)], capsize=3, width=width, color=color_val[ii], alpha=.8)
    


plt.title('median backbone counts')
plt.legend(['reading frame 0','reading frame 1','reading frame 2'])
ax = plt.gca()
leg = ax.get_legend()
leg.legendHandles[0].set_color(color_val[0])
leg.legendHandles[1].set_color(color_val[1])
leg.legendHandles[2].set_color(color_val[2])

plt.xticks([0,width*1.25, 1+0,1+width*1.25, 2+0,2+width*1.25], ['0','1+','0','1+','0','1+'])

plt.xlabel('# in-frame stop codons')
plt.ylabel('median splicing inclusion values')
plt.tight_layout()
pdf_plots.savefig(fig)







fig = plt.figure()
x,y=zip(*scatter_no_stop_vs_stop_list_0_1)
data = np.log2(np.array(x)/np.array(y))
plt.hist(data, bins=np.arange(-8,8,.25))
plt.legend(['median %.1f' % np.median(data)])
plt.xlabel('T5:T4')
plt.ylabel('count')
plt.title('ratio: in frame no stop codon vs out of frame with stop codon')
plt.tight_layout()
pdf_plots.savefig(fig)


fig = plt.figure()
x,y=zip(*scatter_no_stop_vs_stop_list_0_2)
data = np.log2(np.array(x)/np.array(y))
plt.hist(data, bins=np.arange(-9,9,.5))
plt.legend(['median %.1f' % np.median(data)])
plt.xlabel('T5:T3')
plt.ylabel('count')
plt.title('ratio: in frame no stop codon vs out of frame with stop codon')
plt.tight_layout()
pdf_plots.savefig(fig)


fig = plt.figure()
x,y=zip(*scatter_no_stop_vs_stop_list_1_2)
data = np.log2(np.array(x)/np.array(y))
plt.hist(data, bins=np.arange(-9,9,.5))
plt.legend(['median %.1f' % np.median(data)])
plt.xlabel('T3:T4')
plt.ylabel('count')
plt.title('ratio: in frame no stop codon vs out of frame with stop codon')
plt.tight_layout()
pdf_plots.savefig(fig)


fig = plt.figure()
x,y=zip(*scatter_no_stop_vs_stop_list_1_0)
data = np.log2(np.array(x)/np.array(y))
plt.hist(data, bins=np.arange(-9,9,.5))
plt.legend(['median %.1f' % np.median(data)])
plt.xlabel('T3:T5 - Note that T3 tends to have more counts than the other backbones. Needs normalizaton?')
plt.ylabel('count')
plt.title('ratio: in frame no stop codon vs out of frame with stop codon')
plt.tight_layout()
pdf_plots.savefig(fig)


fig = plt.figure()
x,y=zip(*scatter_no_stop_vs_stop_list_2_1)
data = np.log2(np.array(x)/np.array(y))
plt.hist(data, bins=np.arange(-9,9,.5))
plt.legend(['median %.1f' % np.median(data)])
plt.xlabel('T4:T3')
plt.ylabel('count')
plt.title('ratio: in frame no stop codon vs out of frame with stop codon')
plt.tight_layout()
pdf_plots.savefig(fig)


fig = plt.figure()
x,y=zip(*scatter_no_stop_vs_stop_list_2_0)
data = np.log2(np.array(x)/np.array(y))
plt.hist(data, bins=np.arange(-9,9,.5))
plt.legend(['median %.1f' % np.median(data)])
plt.xlabel('T4:T5')
plt.ylabel('count')
plt.title('ratio: in frame no stop codon vs out of frame with stop codon')
plt.show()
plt.tight_layout()
pdf_plots.savefig(fig)



pdf_plots.close()
pdf_scatter_plots.close()





scatter_no_stop_vs_stop_list_2_0

scatter_no_stop_vs_stop_list_0_1






















