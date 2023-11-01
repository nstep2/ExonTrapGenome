from experiment_paths.experiment_paths import *

import numpy as np
import pyfaidx
import matplotlib.pyplot as plt
import pandas as pd



input_fasta = exp_output_path.chr17_shuffle_input_folder + 'chr17_shuffle.fa'
genome_fasta = pyfaidx.Fasta(input_fasta)

len_chr_17 = len(genome_fasta['chr17_shuffle'])

chr17_spliceAI_scores_text_file = exp_output_path.chr17_shuffle_folder+'chr17_shuffle_+_spliceAI_computations.txt'

83258000

spliceAI_score_list_3ss = np.zeros(len_chr_17,dtype=np.float32)
spliceAI_score_list_5ss = np.zeros(len_chr_17,dtype=np.float32)

with open(chr17_spliceAI_scores_text_file,'rt') as f:
    for ii, line in enumerate(f):
        split = line.strip().split('\t')
        pos = int(split[0])
        score_3ss = float(split[1])
        score_5ss = float(split[2])
        
        if ii >= len_chr_17:  #predictions made for 'N' bases
            continue

        spliceAI_score_list_3ss[ii] = score_3ss
        spliceAI_score_list_5ss[ii] = score_5ss

print("loaded spliceAI chr17 shuffle splice site scores")




upper_search_window = 222  #90th
lower_search_window = 63   #10th


splice_AI_exon_id_list = list()
splice_AI_score_pairs = list()
for ii, score_5 in enumerate(spliceAI_score_list_5ss):
    if ii % 10000000 == 0 and ii != 0:
        print('{:,} bp processed\t with {:,} exons found'.format(ii, len(splice_AI_exon_id_list)))
    if score_5 >= 0.2:
        subset_3ss = spliceAI_score_list_3ss[ii-upper_search_window:ii-lower_search_window]
        max_3ss = max(subset_3ss)
        if max_3ss >= 0.2:
            #pair_5ss_3ss.append([list_scores_5ss[ii], max_3ss])
            
            index_3ss = np.argmax(subset_3ss) 
            
            
            
            SAI_exon_id = "%s:%d-%d:%s" % ('chr17_shuffle', ii-upper_search_window+index_3ss-(500-1) , ii-(500-2) ,'+')
            splice_AI_exon_id_list.append(SAI_exon_id)
            splice_AI_score_pairs.append([min(max_3ss,score_5),[max_3ss,score_5],SAI_exon_id])

            
print('{:,} bp processed\t with {:,} exons found'.format(ii, len(splice_AI_exon_id_list)))




df = pd.read_csv(exp_output_path.spliceAI_chr_17_analysis+'chr17_SAI_exon_ids.txt', sep='\t', lineterminator='\n', header=None)

SAI_17_unannotated = list()

exon_ids       = list(df[0])
score_3ss_list = list(df[1])
score_5ss_list = list(df[2])


scores = list((zip(score_3ss_list,score_5ss_list)))
mean_scores=np.mean(scores,axis=1)



import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pdf_stats = PdfPages(exp_output_path.chr17_shuffle_folder+'SAI_chr17_and_shuffle.pdf')

plt.figure()
r=plt.hist(mean_scores,bins=np.arange(0,1.01,0.02),density=True)

a,b,c = zip(*splice_AI_score_pairs)
fig = plt.figure()
plt.plot(r[1][1:],r[0],color='orange')
plt.hist([np.mean(x) for x in b], bins=np.arange(0,1.01,0.02),density=True,color='tab:blue')
plt.xlabel("chr 17 scramble SAI exon scores MEAN(3ss,5ss)")
plt.title('thresholded SAI exons based on MEAN score of pair')
plt.ylabel('count')
plt.legend(['chr17','chr17 shuffle{:,} exons found'.format(len(a))])

pdf_stats.savefig(fig)

pdf_stats.close()
















