











class PESE_class():
    def __init__(self, chasin_PESS,PESE_threshold, PESE_name):
        self.motif_len = len(list(chasin_PESS.keys())[0])
        all_kmer = dict()
        num_to_base = {0:'A',1:'C',2:'G',3:'T'}
        for ii in range(4**self.motif_len):
            kmer_str = ''
            for jj in range(self.motif_len):
                base = num_to_base[int(ii/4**(jj))%4]
                kmer_str = '%s%s' % (kmer_str, base)
            all_kmer[kmer_str] = -1000000
        for key in chasin_PESS:
            all_kmer[key] = chasin_PESS[key]
        self.PESE = all_kmer# chasin_PESS
        self.threshold = PESE_threshold
        self.name = PESE_name
    def score_seq(self, query_seq, **args):
        seq = query_seq.upper()
        if seq.find('N') != -1:
            return -1
        result_list = list()
        motif_window = 50
        
        #full range
        base_positions = range(0,len(seq)-self.motif_len)
        #left end
        base_positions_left = range(0,len(seq[:motif_window])-self.motif_len)
        # right end
        base_positions_right = range(len(seq)-len(seq[-motif_window:])-self.motif_len,len(seq)-self.motif_len)
        
        for ii in base_positions:
            oct_seq = seq[ii:ii+self.motif_len]
            #if oct_seq in self.PESE:
            
            score = self.PESE[oct_seq]
            if 'silencer' not in args:
                if score >= self.threshold:
                    result_list.append([ii, score, oct_seq])
            else:
                if score <= self.threshold and score > -1000000:
                    result_list.append([ii, score, oct_seq])
            
                
        
        return result_list
    
    def score(self,octomer):
        if octomer in self.PESE:
            return self.PESE[octomer]
        return 0
    def generate_ESE_empty_counts_dict(self):
        ESE_empty_counts_dict = dict()
        for ese_nmer in self.PESE:
            ESE_empty_counts_dict[ese_nmer] = 0
        return ESE_empty_counts_dict
    def add_results_to_ESE_counts_dict(self, result_list, ESE_counts_dict):
        for result in result_list:
            ESE_nmer = result[2]
            ESE_counts_dict[ESE_nmer] += 1
    
    


















#https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-python-numpy-scipy
def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w




#repeat_option=exon_ids_overlapping_repeat_list
#repeat_option=[]
repeat_option = exon_ids_overlapping_repeat_list + exon_finder_unique_exon_ids_overlapping_repeat_list








def flanking_PESE_exon_id_list(exon_id_list, PESE, PESE_len, genome_fasta, **args):
    flank_len = 800
    count_PESE_array_3ss = np.zeros(2*flank_len+1)
    count_PESE_array_5ss = np.zeros(2*flank_len+1)
    
    total_seq_len = 0
    exon_len_vs_count_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        if ex.chrom == 'chrM':
            continue
        
        if ex.strand == '+': 
            query_id = '%s:%d-%d:%s' % (ex.chrom, ex.start-flank_len, ex.start+flank_len+1 + PESE_len, ex.strand)
        else:
            query_id = '%s:%d-%d:%s' % (ex.chrom, ex.end-flank_len-1-PESE_len, ex.end+flank_len , ex.strand)
            
        
        #if ex.strand == '-':
        #    continue
        
        seq = el.get_seq_from_exon_id(query_id, genome_fasta)
        
        #print(query_id, seq)
        for ii in range(0, len(seq) - PESE_len):
            query_seq = seq[ii:ii+PESE_len+1].upper()
            if 'silencer' not in args:
                found_PESE = PESE.score_seq(query_seq)
            else:
                found_PESE = PESE.score_seq(query_seq,silencer='')
                #print(ii, found_PESE)
                
            #print(ii, query_seq)
            if found_PESE != -1 and len(found_PESE) != 0:
                #print(ii, found_PESE)
                count_PESE_array_3ss[ii] += 1
        

                
        if ex.strand == '+': 
            query_id = '%s:%d-%d:%s' % (ex.chrom, ex.end-flank_len, ex.end+flank_len+1 + PESE_len, ex.strand)
        else:
            query_id = '%s:%d-%d:%s' % (ex.chrom, ex.start-flank_len-1-PESE_len, ex.start+flank_len , ex.strand)
        
        
        seq = el.get_seq_from_exon_id(query_id, genome_fasta)
    
        #print(query_id, seq)
        for ii in range(0, len(seq) - PESE_len):
            query_seq = seq[ii:ii+PESE_len+1].upper()
            if 'silencer' not in args:
                found_PESE = PESE.score_seq(query_seq)
            else:
                found_PESE = PESE.score_seq(query_seq,silencer='')
                #print(ii, found_PESE)
                
            #print(ii, query_seq)
            if found_PESE != -1 and len(found_PESE) != 0:
                #print(ii, found_PESE)
                count_PESE_array_5ss[ii] += 1
                        
        
        
        
    return count_PESE_array_3ss, count_PESE_array_5ss





exon_id_list=pc_middle_exon_id_list[:1000]
PESE_len = 6
count_PESE_array_3ss, count_PESE_array_5ss = flanking_PESE_exon_id_list(exon_id_list, PESE, PESE_len, genome_fasta)







exon_id_list=list(no_overlap_set)[:1000]
PESE_len = 6
count_PESE_array_3ss, count_PESE_array_5ss = flanking_PESE_exon_id_list(exon_id_list, PESE, PESE_len, genome_fasta)






pickle_path = exp_output_path.pickle_merged + "ET_SAI_MES_exon_sets_%d.pickle" % (exon_count_build)

finder_unique=dict()
with open(pickle_path, "rb") as input_file:
    D = pickle.load(input_file)
    finder_unique['ET_only'] = pickle.load(input_file)
    finder_unique['SAI_only'] = pickle.load(input_file)
    finder_unique['MES_only'] = pickle.load(input_file)
    finder_unique['mRNA'] = pc_middle_exon_id_list
    finder_unique['intergenic'] = no_overlap_set
    


import matplotlib.pyplot as plt



















ISS_dict_3 = dict()
ISS_dict_5 = dict()
ISS_dict = dict()
ISS_list = list()

ISS_path = '/mnt/hgfs/main_ssd/Linux_2TB/ESE/ISS_PMID_20685814.txt'

with open( ISS_path, 'r') as f:
    next(f)
    for line in f:
        line_split = line.split()
        #print(line.strip())
        if line_split[1] == '5ISS':
            ISS_list.append(line_split[0])
            ISS_dict_5[line_split[0]] = -1

        if line_split[1] == '3ISS':
            ISS_list.append(line_split[0])
            ISS_dict_3[line_split[0]] = -1

        if line_split[1] == '5ISS' or line_split[1] == '3ISS':
            ISS_list.append(line_split[0])
            ISS_dict[line_split[0]] = -1




chasin_ESE_filepath=exp_output_path.ESE_input_folder+'chasin_8mer_scores_wayback_machine.txt'
chasin_PESE=dict()
chasin_PESS=dict()
chasin_PISS=dict()
chasin_PISE=dict()
PESS_threshold = -2.62
PESE_threshold = 2.62
with open(chasin_ESE_filepath,'r') as f:
    for ii in range(23):
        #next(f)
        print(next(f).strip())
    
    
    for jj, line in enumerate(f):
        line = line.strip()
        line_split = line.split('\t')
        if line_split[0]=='':   #empty lines occur 1-2 times
            continue
        #exonic
        P_index = float(line_split[1])
        octomer=line_split[0].upper()
        if P_index >= PESE_threshold:
            chasin_PESE[octomer]=P_index
        if P_index <= PESS_threshold:
            chasin_PESS[octomer]=P_index
        
        #intronic
        P_index = float(line_split[2])
        octomer=line_split[0].upper()
        if P_index >= PESE_threshold:
            chasin_PISE[octomer]=P_index
        if P_index <= PESS_threshold:
            chasin_PISS[octomer]=P_index
        
        

ISS_scorer = PESE_class(chasin_PISS, 0, 'PISS')


ISE_scorer = PESE_class(chasin_PISE, 0, 'PISE')



ISS_scorer_5 = PESE_class(ISS_dict_5, 0, 'PISS')

ISS_scorer_3 = PESE_class(ISS_dict_3, 0, 'PISS')


ISS_scorer = PESE_class(ISS_dict, 0, 'PISS')




xc = [ {a:ISS_scorer.PESE[a]} for a in ISS_scorer.PESE if ISS_scorer.PESE[a] > -10000]






if True == False:
    
    splice_enhance_left_dict_3 = dict()
    splice_enhance_right_dict_3 = dict()
    
    
    for ii, key in enumerate(finder_unique):
        
        exon_id_list=set(finder_unique[key]).difference(repeat_option)
        PESE_len = 8
        count_PESE_array_3ss, count_PESE_array_5ss = flanking_PESE_exon_id_list(exon_id_list, ISE_scorer, PESE_len, genome_fasta)
        splice_enhance_left_dict_3[key] = count_PESE_array_3ss, count_PESE_array_5ss
        
    
    for ii, key in enumerate(finder_unique):
        
        exon_id_list=set(finder_unique[key]).difference(repeat_option)
        PESE_len = 8
        count_PESE_array_3ss, count_PESE_array_5ss = flanking_PESE_exon_id_list(exon_id_list, ISS_scorer, PESE_len, genome_fasta,silencer=True)
        splice_enhance_right_dict_3[key] = count_PESE_array_3ss, count_PESE_array_5ss
    
    


splice_enhance_left_dict_3 = dict()
splice_enhance_right_dict_3 = dict()


for ii, key in enumerate(finder_unique):
    
    exon_id_list=set(finder_unique[key]).difference(repeat_option)

    PESE_len = 6

    splice_enhance_left_dict_3[key] = [],[]


for ii, key in enumerate(finder_unique):
    
    exon_id_list=set(finder_unique[key]).difference(repeat_option)

    PESE_len = 6
    count_PESE_array_3ss, count_PESE_array_5ss = flanking_PESE_exon_id_list(exon_id_list, ISS_scorer, PESE_len, genome_fasta,silencer=True)
    splice_enhance_right_dict_3[key] = count_PESE_array_3ss, count_PESE_array_5ss















if True == False:
    
    splice_enhance_left_dict_4 = dict()
    splice_enhance_right_dict_4 = dict()
    
    
    
    
    for ii, key in enumerate(finder_unique):
        
        exon_id_list=set(finder_unique[key]).difference(repeat_option)
        PESE_len = 8
        count_PESE_array_3ss, count_PESE_array_5ss = flanking_PESE_exon_id_list(exon_id_list, ISE_scorer, PESE_len, genome_fasta)
        
        splice_enhance_left_dict_4[key] = count_PESE_array_3ss, count_PESE_array_5ss
    
    for ii, key in enumerate(finder_unique):
        
        exon_id_list=set(finder_unique[key]).difference(repeat_option)
        PESE_len = 8
        count_PESE_array_3ss, count_PESE_array_5ss = flanking_PESE_exon_id_list(exon_id_list, ISS_scorer, PESE_len, genome_fasta,silencer=True)
        
        splice_enhance_right_dict_4[key] = count_PESE_array_3ss, count_PESE_array_5ss







splice_enhance_left_dict_4 = dict()
splice_enhance_right_dict_4 = dict()




for ii, key in enumerate(finder_unique):
    
    exon_id_list=set(finder_unique[key]).difference(repeat_option)
    PESE_len = 6

    splice_enhance_left_dict_4[key] = [],[]


for ii, key in enumerate(finder_unique):
    
    exon_id_list=set(finder_unique[key]).difference(repeat_option)

    PESE_len = 6
    count_PESE_array_3ss, count_PESE_array_5ss = flanking_PESE_exon_id_list(exon_id_list, ISS_scorer, PESE_len, genome_fasta,silencer=True)
    
    splice_enhance_right_dict_4[key] = count_PESE_array_3ss, count_PESE_array_5ss
    



















outdir = exp_output_path.review_response_ESE 
pdf_plots = PdfPages(outdir+'review_response_ESE.pdf')
















color_set_dict = {'annotated':'#D6D6D6','ET_only':'#31BBED','SAI_only':'#FBAA28','MES_only':'#EBE837'}


color_set_dict['mRNA'] = color_dict['mRNA']
color_set_dict['intergenic'] = color_dict['Intergenic']




fig, axs = plt.subplots(1, 1, figsize=(8,3) )
fig.suptitle('3ss')


for ii, key in enumerate(['mRNA', 'intergenic']):   
    exon_id_list=set(finder_unique[key]).difference(repeat_option)
    PESE_len = 8
    count_PESE_array_3ss, count_PESE_array_5ss = splice_enhance_right_dict_3[key]  
    
    axs.plot(range(len(count_PESE_array_3ss)-18), moving_average(count_PESE_array_3ss/len(exon_id_list),19), label = key, c=color_set_dict[key])    
    axs.set_ylabel('Density')

    xax_dist=int( (len(count_PESE_array_5ss)-1)/2)
    xax_range = range(-xax_dist,xax_dist+1,400)
    axs.set_xticks([x+xax_dist for ii,x in enumerate(xax_range)])
    axs.set_xticklabels(["{:,}".format(x) for ii, x in enumerate(xax_range)])
    axs.set_ylim(0.005,0.04)
    
axs.set_title('ISS')    
axs.legend()    

plt.tight_layout()
plt.show()


pdf_plots.savefig(fig)



























fig, axs = plt.subplots(1, 1, figsize=(8,3) )
fig.suptitle('3ss')



for ii, key in enumerate(['mRNA', 'ET_only', 'SAI_only', 'MES_only']):   
    exon_id_list=set(finder_unique[key]).difference(repeat_option)
    PESE_len = 8
    count_PESE_array_3ss, count_PESE_array_5ss = splice_enhance_right_dict_3[key]  
    
    axs.plot(range(len(count_PESE_array_3ss)-18), moving_average(count_PESE_array_3ss/len(exon_id_list),19), label = key, c=color_set_dict[key])    
    axs.set_ylabel('Density')

    xax_dist=int( (len(count_PESE_array_5ss)-1)/2)
    xax_range = range(-xax_dist,xax_dist+1,400)
    axs.set_xticks([x+xax_dist for ii,x in enumerate(xax_range)])
    axs.set_xticklabels(["{:,}".format(x) for ii, x in enumerate(xax_range)])
    axs.set_ylim(0.005,0.04)
    


axs.set_title('ISS')    

axs.legend()    

plt.tight_layout()
plt.show()



pdf_plots.savefig(fig)












fig, axs = plt.subplots(1, 1, figsize=(8,3) )
fig.suptitle('5ss')


for ii, key in enumerate(['mRNA', 'intergenic']):   
    exon_id_list=set(finder_unique[key]).difference(repeat_option)
    PESE_len = 8
    count_PESE_array_3ss, count_PESE_array_5ss = splice_enhance_right_dict_4[key]  
    
    axs.plot(range(len(count_PESE_array_5ss)-18), moving_average(count_PESE_array_5ss/len(exon_id_list),19), label = key, c=color_set_dict[key])    
    axs.set_ylabel('Density')

    xax_dist=int( (len(count_PESE_array_5ss)-1)/2)
    xax_range = range(-xax_dist,xax_dist+1,400)
    axs.set_xticks([x+xax_dist for ii,x in enumerate(xax_range)])
    axs.set_xticklabels(["{:,}".format(x) for ii, x in enumerate(xax_range)])
    
axs.set_ylim(0.005,0.04)
axs.set_title('ISS')    
axs.legend()    

plt.tight_layout()
plt.show()



pdf_plots.savefig(fig)





fig, axs = plt.subplots(1, 1, figsize=(8,3) )
fig.suptitle('5ss')



for ii, key in enumerate(['mRNA', 'ET_only', 'SAI_only', 'MES_only']):   
    exon_id_list=set(finder_unique[key]).difference(repeat_option)
    PESE_len = 8
    count_PESE_array_3ss, count_PESE_array_5ss = splice_enhance_right_dict_4[key]  
    
    axs.plot(range(len(count_PESE_array_5ss)-18), moving_average(count_PESE_array_5ss/len(exon_id_list),19), label = key, c=color_set_dict[key])    
    axs.set_ylabel('Density')

    xax_dist=int( (len(count_PESE_array_5ss)-1)/2)
    xax_range = range(-xax_dist,xax_dist+1,400)
    axs.set_xticks([x+xax_dist for ii,x in enumerate(xax_range)])
    axs.set_xticklabels(["{:,}".format(x) for ii, x in enumerate(xax_range)])
    axs.set_ylim(0.005,0.04)
    


axs.set_title('ISS')    
axs.legend()    

plt.tight_layout()
plt.show()


pdf_plots.savefig(fig)










pdf_plots.close()









