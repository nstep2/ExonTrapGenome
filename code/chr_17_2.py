#Needs:
# pc_middle_exon_id_list    
# /mnt/hgfs/main_ssd/slim_test/code/exon_def/et_main/submit_aug_2023/ch17_shuffle_2.py for chr17_shuffle_+_spliceAI_computations.txt

# > 24 hrs


from matplotlib.backends.backend_pdf import PdfPages


import numpy as np


len_chr_17 = len(genome_fasta['chr17'])


chr_17_1
ch17_shuffle_1
ch17_shuffle_2



def count_gc(seq):
    seq=str(seq).upper()
    gc_count = 0
    seq_len = len(seq)
    for ii, b in enumerate(seq):
        if b == 'C' or b == 'G':
            gc_count += 1
    
    return gc_count, gc_count/seq_len
            







#upper_search_window = 222  90th
#lower_search_window = 63   10th

pc_exon_lengths = [el.exon_id_values(exon_id).length for exon_id in pc_middle_exon_id_list]
upper_search_window = int(np.quantile(pc_exon_lengths,.9))
lower_search_window = int(np.quantile(pc_exon_lengths,.1))


exp_output_path.chr17_folder
chr17_spliceAI_scores_text_file = exp_output_path.chr17_folder+'chr17_+_spliceAI_computations.txt'



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



def get_chr17_plus_exon_ids(exon_id_list):
    
    
    chr17_pc_middle_exons_10_200 = list()
    chr17_pc_middle_exons = list()
    
    for ii, exon_id in enumerate(exon_id_list):
        ex=el.exon_id_values(exon_id)
    
        if ex.strand == '+' and ex.chrom=='chr17':
            chr17_pc_middle_exons.append(exon_id)
            
            if ex.length >= lower_search_window and ex.length <= upper_search_window:
                chr17_pc_middle_exons_10_200.append(exon_id)
    
                
    return chr17_pc_middle_exons_10_200,chr17_pc_middle_exons



chr17_pc_middle_exons_10_200,chr17_pc_middle_exons = get_chr17_plus_exon_ids(pc_middle_exon_id_list)

chr17_lncRNA_middle_exons_10_200,chr17_lcnRNA_middle_exons = get_chr17_plus_exon_ids(lncRNA_middle_exon_ids)





SAI_power_exponent = 1
SAI_3ss_AG_count = 0
SAI_3ss_other_count = 0
SAI_GT_count = 0
SAI_AT_count = 0
SAI_other_count = 0
splice_AI_exon_id_list = list()


splice_AI_score_pairs = list()

for ii, score_5 in enumerate(spliceAI_score_list_5ss):
    #if ii == 10000000:
    #    break
    if ii % 10000000 == 0 and ii != 0:
        print('{:,} bp processed\t with {:,} exons found'.format(ii, len(splice_AI_exon_id_list)))
    if score_5 >= 0.2:#**(SAI_power_exponent):
        subset_3ss = spliceAI_score_list_3ss[ii-upper_search_window:ii-lower_search_window]
        max_3ss = max(subset_3ss)
        if max_3ss >= 0.2:#**(SAI_power_exponent):
            #pair_5ss_3ss.append([list_scores_5ss[ii], max_3ss])
            
            index_3ss = np.argmax(subset_3ss) 
            
            index_3ss_list = subset_3ss >= 0.2#**(SAI_power_exponent)
            for jj, val in enumerate(index_3ss_list):
                if val == 0:
                    continue
                index_3ss = jj
                pos_3ss = ii-upper_search_window + index_3ss
                pos_5ss = ii-(500-2)
                
                seq_5ss_dint = (str(genome_fasta['chr17'][pos_5ss-1:pos_5ss+1]).upper())
                seq_3ss_dint = (str(genome_fasta['chr17'][pos_3ss-500+1-3:pos_3ss-500+1-1]).upper())
                if seq_5ss_dint == 'GT':
                    SAI_GT_count += 1
                elif seq_5ss_dint == 'AT':
                    SAI_AT_count += 1
                    continue
                else:
                    SAI_other_count += 1
                    continue
                
                #print(seq_3ss_dint)
                if seq_3ss_dint == 'AG':
                    SAI_3ss_AG_count += 1
                else:
                    SAI_3ss_other_count += 1
                    continue
            
            
            SAI_exon_id = "%s:%d-%d:%s" % ('chr17', ii-upper_search_window+index_3ss-(500-1) , ii-(500-2) ,'+')
            splice_AI_exon_id_list.append(SAI_exon_id)
            splice_AI_score_pairs.append([min(max_3ss,score_5),[max_3ss,score_5],SAI_exon_id])                

len(splice_AI_exon_id_list)





gencode_overlap_exon_ids_list = list()
for ii, exon_id in enumerate(splice_AI_exon_id_list):
    ex = el.exon_id_values(exon_id)

    intervals=gencode_exon_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) > 0:
        gencode_overlap_exon_ids_list.append(exon_id)



    
unannotated_exon_ids = list(set(splice_AI_exon_id_list).difference(gencode_overlap_exon_ids_list))
len(unannotated_exon_ids)



with open(exp_output_path.spliceAI_chr_17_analysis+'chr17_SAI_exon_ids.txt','w') as f:
    SAI_score_pairs_dict={entry[2]:entry[1] for entry in splice_AI_score_pairs}
    
    for exon_id in unannotated_exon_ids:
        entry=SAI_score_pairs_dict[exon_id]
        out_string = '{:}\t{:}\t{:}\n'.format(exon_id,entry[0],entry[1])
        f.write(out_string)











        
        





import matplotlib.pyplot as plt
import numpy as np






maxentscan_5ss_file = exp_output_path.chr17_folder +'chr17_+_maxentscan_5ss.txt'

maxentscan_3ss_file = exp_output_path.chr17_folder +'chr17_+_maxentscan_3ss.txt'


list_scores_5ss = np.zeros(len_chr_17)
list_ratio_5ss  = np.zeros(len_chr_17)
with open(maxentscan_5ss_file, 'rt') as f:
    for ii, line in enumerate(f):
        split = line.strip().split('\t')
        list_scores_5ss[ii] = float(split[1])
        list_ratio_5ss[ii]  = float(split[2])
        

list_scores_3ss = np.zeros(len_chr_17)
list_ratio_3ss  = np.zeros(len_chr_17)
with open(maxentscan_3ss_file, 'rt') as f:
    for ii, line in enumerate(f):
        split = line.strip().split('\t')
        list_scores_3ss[ii] = float(split[1])
        list_ratio_3ss[ii] = float(split[2])



MES_3ss_AG_count = 0
MES_3ss_other_count = 0
MES_GT_count=0
MES_AT_count=0
MES_other_count=0
MES_threshold_score = 6
pair_5ss_3ss = list()
pair_exon_id_list = list()


splice_AI_score_pairs = list()

for ii, val in enumerate(list_scores_5ss):
    #if ii == 1000000:
    #    break
    ii=ii-4
    if list_scores_5ss[ii] >= MES_threshold_score:
        subset_3ss = list_scores_3ss[ii-upper_search_window-21-4:ii-lower_search_window-21-4]
        max_3ss = max(subset_3ss)
        if max_3ss >= MES_threshold_score:
            pair_5ss_3ss.append([list_scores_5ss[ii], max_3ss])
            #index_3ss = np.argmax(subset_3ss) 
            
            index_3ss_list = subset_3ss >= MES_threshold_score
            for jj, val in enumerate(index_3ss_list):
                if val == 0:
                    continue
                index_3ss = jj
                
                pos_3ss = ii-upper_search_window -21-4 + index_3ss
                genome_pos_5ss = ii+4
                seq_5ss_dint = (str(genome_fasta['chr17'][genome_pos_5ss-1:genome_pos_5ss+1]).upper())
                #print(seq_5ss_dint)
                seq_3ss_dint = (str(genome_fasta['chr17'][pos_3ss-3+21:pos_3ss+21-1]).upper())
                if seq_5ss_dint == 'GT':
                    MES_GT_count += 1
                elif seq_5ss_dint == 'AT':
                    MES_AT_count += 1
                    continue
                else:
                    MES_other_count += 1
                    continue
                if seq_3ss_dint == 'AG':
                    MES_3ss_AG_count += 1
                else:
                    MES_3ss_other_count += 1
                    continue
                max_ent_exon_id = "%s:%d-%d:%s" % ('chr17', pos_3ss+21, ii+4,'+')
                pair_exon_id_list.append(max_ent_exon_id)

            
        
        




pair_5ss_3ss=np.array(pair_5ss_3ss)


maxent_pc_middle_list = (el.exon_id_intersection(pair_exon_id_list, chr17_pc_middle_exons_10_200))

SAI_pc_middle_list = (el.exon_id_intersection(splice_AI_exon_id_list, chr17_pc_middle_exons_10_200))

SAI_maxent_list = el.exon_id_intersection(pair_exon_id_list, splice_AI_exon_id_list)




primary_3ss_exon_id_set  #already thresholded by 100
ET_3ss_AG_count = 0
ET_3ss_other_count = 0
ET_GT_count=0
ET_AT_count=0
ET_other_count=0
ET_chr_17_keys = list()
#for exon_id in aggregate_exon_dict:
for exon_id in el.threshold_exon_ids(aggregate_exon_dict.keys(),  100, aggregate_exon_dict):
    #el.threshold_exon_ids()
#for ii, exon_id in enumerate(primary_3ss_exon_id_set):
    #if ii == 1000:
    #    break
    ex=el.exon_id_values(exon_id)
    if ex.strand == '+' and ex.chrom == 'chr17' :
        if aggregate_exon_dict[exon_id]['count'] < 100:
            continue
        if ex.length >= lower_search_window and ex.length <= upper_search_window:
            
            seq_5ss_dint = (str(genome_fasta['chr17'][ex.pos_5ss-1:ex.pos_5ss+1]).upper())
            seq_3ss_dint = (str(genome_fasta['chr17'][ex.pos_3ss-3:ex.pos_3ss+-1]).upper())
            #print(seq_5ss_dint)
            if seq_5ss_dint == 'GT':
                ET_GT_count += 1
            elif seq_5ss_dint == 'AT':
                ET_AT_count += 1
                continue
            else:
                ET_other_count += 1
                continue
            if seq_3ss_dint == 'AG':
                ET_3ss_AG_count += 1
            else:
                ET_3ss_other_count += 1
                continue
            
            ET_chr_17_keys.append(exon_id)
    
len(ET_chr_17_keys)

ET_pc_middle_list =(el.exon_id_intersection(ET_chr_17_keys,chr17_pc_middle_exons_10_200))
len(ET_pc_middle_list)


















def get_ESE_count_rate(exon_id_list, genome_fasta, MES_range_3ss, MES_range_5ss, PESE):
    
    
    gc_list=list()
    ESE_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        
        score_5ss = ex.score_5ss(genome_fasta)[0]
        score_3ss = ex.score_3ss(genome_fasta)[0]
        if score_5ss <= MES_range_5ss[0] or score_5ss >= MES_range_5ss[1]:
            continue
        if score_3ss <= MES_range_3ss[0] or score_3ss >= MES_range_3ss[1]:
            continue
        
        
        seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        
        found_PESE = PESE.score_seq(seq)
    
        if found_PESE == -1:
            continue
        
        ESE_list.append([len(found_PESE), len(found_PESE)/ex.length, score_3ss, score_5ss, exon_id])
        

    return ESE_list












ouput_folder = exp_output_path.spliceAI_chr_17_analysis
pdf_venn = PdfPages(ouput_folder+'pdf_venn.pdf')




import venn




from pylab import  cm
cmap = cm.get_cmap('Reds', 56405)


import matplotlib
cmap = matplotlib.cm.get_cmap('Reds')
rgba = cmap(0.5)
print(rgba)

color_nums = [44,35,6789,293,773,2492,4687,635,2356,3057,3081,21,31,3166,56405]

log_color_nums = [[x,np.log10(x)] for ii, x in enumerate(color_nums)]

for val in log_color_nums:
    print("{:>6} \t\t{:.1f}".format(val[0],val[1]))



log_color_nums = [[x,np.log10(x)] for ii, x in enumerate(color_nums)]
for val in log_color_nums:
    a = val[1]
    c = a/8
    hex_rgb = matplotlib.colors.rgb2hex( (1,1-c,1-c) )
    print("{:>6} \t\t{:.1f}\t\t{:}".format(val[0],val[1], hex_rgb))







color_nums = [11078,1215,787,87,801,918,920,144,61223,1756,475,3429,584,19,16]

log_color_nums = [[x,np.log10(x)] for ii, x in enumerate(color_nums)]
for val in log_color_nums:
    a = val[1]
    c = a/8
    hex_rgb = matplotlib.colors.rgb2hex( (1,1-c,1-c) )
    print("{:>6} \t\t{:.1f}\t\t{:}".format(val[0],val[1], hex_rgb))











import venn


import venn
labels = venn.get_labels([set(chr17_pc_middle_exons_10_200), set(ET_chr_17_keys), set(splice_AI_exon_id_list), set(pair_exon_id_list)], fill=['number', 'logic'])
#fig, ax = venn.venn(labels, names=['annotated', 'ET', 'SpliceAI', 'MaxEntScan'])
#fig.show()



from venn import generate_petal_labels, draw_venn, generate_colors








dataset = {'annotated':set(chr17_pc_middle_exons_10_200+chr17_lncRNA_middle_exons_10_200), 'ET':set(ET_chr_17_keys), 'SpliceAI':set(splice_AI_exon_id_list), 'MaxEntScan':set(pair_exon_id_list)}

'''
for key in dataset:
    exon_id_list = dataset[key]
    unique = el.check_exon_ids_by_blat(exon_id_list, genome_fasta)
    print('{:} unique/all: {:}/{:} '.format(key, len(unique) ,len(exon_id_list)))
    dataset[key] = set(unique)
'''





import venn
fig,ax = plt.subplots(figsize=(12,8))
dataset = {'annotated':set(chr17_pc_middle_exons_10_200+chr17_lncRNA_middle_exons_10_200), 'ET':set(ET_chr_17_keys), 'SpliceAI':set(splice_AI_exon_id_list), 'MaxEntScan':set(pair_exon_id_list)}
datasets=dict()
for key in dataset:
    new_key = '%s:  %d' % (key, len(dataset[key]))
    datasets[new_key] = dataset[key]
#, fill=['number', 'logic'}
#venn.venn(datasets,ax=ax,cmap=['#D6D6D6','#EBE837','#31BBED','#FBAA28'])
venn.venn(datasets,ax=ax,cmap=['#D6D6D6','#31BBED','#FBAA28','#EBE837'])
plt.title('chr17+ pc + lncRNA exons 4 group venn')
plt.tight_layout()
pdf_venn.savefig(fig)



#supplemental 5A
with open(exp_output_path.out_supplemental+'5A.txt','w') as f:
    outstring = 'chrom\tstart\tend\tstrand\texon_id\texon_finder\n'
    f.write(outstring)
    
    exon_id_list = set()
    for exon_id_list_name in dataset:
        exon_id_list = exon_id_list.union(dataset[exon_id_list_name])
        
    for exon_id in exon_id_list:
        ex=el.exon_id_values(exon_id)
        
        outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id,exon_id_list_name)
        
        val = ''
        for exon_id_list_name in dataset:
            if exon_id in dataset[exon_id_list_name]:
                val = '{:},{:}'.format(val,exon_id_list_name)
        val = val[1:] #remove starting ',' if present
        
        outstring = '{:}\t{:}\n'.format(outstring, val)
        f.write(outstring)

    






dataset_3ss = dict()
for key in dataset:
    dataset_3ss[key] = list()

for key in dataset:
    for exon_id in dataset[key]:
        ex = el.exon_id_values(exon_id)
        dataset_3ss[key].append(ex.start)

for key in dataset:
    dataset_3ss[key] = set(dataset_3ss[key])

import venn
fig,ax = plt.subplots(figsize=(12,8))

datasets=dict()
for key in dataset_3ss:
    new_key = '%s:  %d' % (key, len(dataset_3ss[key]))
    datasets[new_key] = dataset_3ss[key]

venn.venn(datasets,ax=ax,cmap=['#D6D6D6','#31BBED','#FBAA28','#EBE837'])
plt.title('chr17+ pc + lncRNA exons 4 group venn, share 3ss')
plt.tight_layout()
pdf_venn.savefig(fig)



#supplemental S5A
with open(exp_output_path.out_supplemental+'S7A.txt','w') as f:
    outstring = 'chrom\tstart\tend\tstrand\texon_id\texon_finder\n'
    f.write(outstring)
    
    exon_id_list = set()
    for exon_id_list_name in dataset:
        exon_id_list = exon_id_list.union(dataset[exon_id_list_name])
        
    for exon_id in exon_id_list:
        ex=el.exon_id_values(exon_id)
        
        outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id,exon_id_list_name)
        
        val = ''
        for exon_id_list_name in dataset_3ss:
            exon_end = ex.start
            if exon_end in dataset_3ss[exon_id_list_name]:
                
                val = '{:},{:}'.format(val,exon_id_list_name)
        val = val[1:] #remove starting ',' if present
        
        outstring = '{:}\t{:}\n'.format(outstring, val)
        f.write(outstring)






#dataset = {'annotated':set(chr17_pc_middle_exons_10_200+chr17_lncRNA_middle_exons_10_200), 'ET':set(ET_chr_17_keys), 'SpliceAI':set(splice_AI_exon_id_list), 'MaxEntScan':set(pair_exon_id_list)}

dataset_5ss = dict()
for key in dataset:
    dataset_5ss[key] = list()

for key in dataset:
    for exon_id in dataset[key]:
        ex = el.exon_id_values(exon_id)
        dataset_5ss[key].append(ex.end)

for key in dataset:
    dataset_5ss[key] = set(dataset_5ss[key])

import venn
fig,ax = plt.subplots(figsize=(12,8))

datasets=dict()
for key in dataset_5ss:
    new_key = '%s:  %d' % (key, len(dataset_5ss[key]))
    datasets[new_key] = dataset_5ss[key]

venn.venn(datasets,ax=ax,cmap=['#D6D6D6','#31BBED','#FBAA28','#EBE837'])
plt.title('chr17+ pc + lncRNA exons 4 group venn, share 5ss')
plt.tight_layout()
pdf_venn.savefig(fig)

















pdf_venn.close()























#write maxentscan exon_ids

m_exon_id_file_path = exp_output_path.chr17_folder +'chr17_+_maxentscan_exon_ids_score_%d.txt' % (MES_threshold_score)

with open(m_exon_id_file_path,'wt') as f:
    for exon_id in pair_exon_id_list:
        f.write(exon_id)
    


























'''



int('stop')








x,y = zip(*pair_5ss_3ss)

counts, bins = np.histogram(y, bins = np.arange(-50,15,.1))
plt.figure()
plt.hist(bins[:-1],bins,weights=counts)
plt.plot([threshold_score,threshold_score],[0,max(counts)])
plt.title('max 3ss score when 5ss above 0')



pair_threshold = 12

for pair_threshold in [6,7,8,9,10,12,14,16,18,20]:
    pair_above_thresh_list = list()
    both_above_thresh_list = list()
    for ii, pair in enumerate(pair_5ss_3ss):
        if sum(pair) >= pair_threshold:
            pair_above_thresh_list.append(pair)
        if pair[0] >= pair_threshold and pair[1] >= pair_threshold:
            both_above_thresh_list.append(pair)
            
    print("pair threshold above %d has %d exons" % (pair_threshold,len(pair_above_thresh_list)))
    
    print("both threshold above %d has %d exons" % (pair_threshold,len(both_above_thresh_list)))
    
    
    x,y = zip(*pair_above_thresh_list)
    
    counts, bins = np.histogram(y, bins = np.arange(-50,15,.1))
    plt.figure()
    plt.hist(bins[:-1],bins,weights=counts)
    #plt.plot([threshold_score,threshold_score],[0,max(counts)])
    plt.title('max 3ss score when 3ss+5ss above %d' % (pair_threshold))
    plt.xlabel('3ss score (maxentscan)')
    
    
    counts, bins = np.histogram(x, bins = np.arange(-50,15,.1))
    plt.figure()
    plt.hist(bins[:-1],bins,weights=counts)
    #plt.plot([threshold_score,threshold_score],[0,max(counts)])
    plt.title('max 5ss score when 3ss+5ss above %d' % (pair_threshold))
    plt.xlabel('5ss score (maxentscan)')
>>>>>>> f25229cec6ecd847ff0f5cab9be40ebea3d7e829

















'''











