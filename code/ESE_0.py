
"""


AWS recommended color blind safe colors
(Information taken from cookbook-r.com/Graphs/Colors_(ggplot2)/)
00000
F7931D
00B9F1 	lncRNA
00A875	Antisense
ECDE38	
0072BC	mRNA
F15A22	Intergenic
DA6FAB	Intronic


"""


import exon_id_library.exon_id_lib as el


exon_count_threshold = 100


PESE_name = 'chasin_ESEseq'







if 'avg_count' not in aggregate_exon_dict[list(aggregate_exon_dict.keys())[0]]:
    
    for exon_id in aggregate_exon_dict:
        exon = aggregate_exon_dict[exon_id]
        denom = sum(exon['lib_array']>5)
        exon['avg_count'] = exon['count']/denom
    





# https://genome.cshlp.org/content/21/8/1360.long

exp_output_path.ESE_input_folder
burge_ESE_file = exp_output_path.ESE_input_folder +'burge/ESE_burge.txt'
ESESeq_filepath = exp_output_path.ESE_input_folder +'ESEseq/Supplemental_Table_1_main.txt'
rescue_ESE_filepath = exp_output_path.ESE_input_folder +'rescue-ese_wayback_machine.txt'




ESE_seq_threshold = 0
ESE_seq_score=dict()
with open(ESESeq_filepath,'r') as f:
    for ii in range(10):
        line = next(f)
    
        
    ESE_score_column = 23 
    
    for jj, line in enumerate(f):
        line = line.strip()
        line_split = line.split('\t')
        if line_split[24] == 'E':
           hexamer = line_split[1] 
           score = float(line_split[ESE_score_column])
           if score >= ESE_seq_threshold:
               ESE_seq_score[hexamer] = score
           
        if line_split[0] == '4096':
            break
        
        




ESE_seq_threshold = 0
ESE_seq_silencer_score=dict()
with open(ESESeq_filepath,'r') as f:
    for ii in range(10):
        line = next(f)
        
        
    
    ESE_score_column = 23 
    
    
    for jj, line in enumerate(f):
        line = line.strip()
        line_split = line.split('\t')
        if line_split[24] == 'S':
           hexamer = line_split[1] 
           score = float(line_split[ESE_score_column])
           if score < ESE_seq_threshold:
               ESE_seq_silencer_score[hexamer] = score
           
        if line_split[0] == '4096':
            break







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
    
    







try:
    type(ESE_seq)
except:
    ESE_seq          = PESE_class(ESE_seq_score,ESE_seq_threshold, 'ESEseq')    
    ESE_silencer_seq = PESE_class(ESE_seq_silencer_score,ESE_seq_threshold, 'ESEseq')
    
   
    
#Alternate PESE scorers - 6mers
if PESE_name == 'chasin_ESEseq':
    PESE=ESE_seq








pickle_path = exp_output_path.pickle_main
#with open(exp_output_path.pickle_main+'chasin_PESE.pickle', 'wb') as outfile:
#    pickle.dump(PESE_class(chasin_PESE,0,'chasin_PESE'),outfile)

with open(exp_output_path.pickle_main+'chasin_ESEseq.pickle', 'wb') as outfile:
    pickle.dump(PESE_class(ESE_seq_score,0, 'ESEseq') ,outfile)
    pickle.dump(PESE_class(ESE_seq_silencer_score,0, 'ESEseq_silencers') ,outfile)
    

#with open(exp_output_path.pickle_main+'burge_ESE.pickle', 'wb') as outfile:
#    pickle.dump(burge_ESE,outfile)













def fraction_covered(found_PESE,exon_len):
    exon_range = np.zeros(exon_len)
    
    for hit in found_PESE:
        exon_range[hit[0]:hit[0]+len(hit[2])] += 1
    
    fraction_covered = sum(exon_range>0)/exon_len
    return fraction_covered


def get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta):
    count_PESE = 0
    total_seq_len = 0
    exon_len_vs_count_list = list()
    exon_fraction_covered_list = list()
    for exon_id in exon_id_list:
        if el.exon_id_values(exon_id).chrom == 'chrM':
            continue
        
            
        seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        found_PESE = PESE.score_seq(seq)
        if found_PESE != -1:
            total_seq_len+=len(seq)
            count_PESE += len(found_PESE)
            exon_len_vs_count_list.append([len(seq),len(found_PESE)])
            frac_covered = fraction_covered(found_PESE,len(seq))
            exon_fraction_covered_list.append(frac_covered)
        
        if total_seq_len != 0:
            found_rate_100_bp = 100*count_PESE/total_seq_len
        else:
            found_rate_100_bp =0
    
    return found_rate_100_bp, exon_len_vs_count_list, exon_fraction_covered_list











def get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_id_list, PESE, genome_fasta, aggregate_exon_dict):
    exon_id_list = el.exon_id_intersection(exon_id_list, aggregate_exon_dict.keys())
    count_PESE = 0
    
    total_seq_len = 0
    exon_len_vs_count_list = list()
    exon_inclusion_vs_count_list = list()
    exon_fraction_covered_list = list()
    for exon_id in exon_id_list:
        seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        
        found_PESE = PESE.score_seq(seq)
        if found_PESE != -1:
            total_seq_len+=len(seq)
            count_PESE += len(found_PESE)
            
            frac_covered = fraction_covered(found_PESE,len(seq))
            exon_fraction_covered_list.append(frac_covered)
            exon_len_vs_count_list.append([len(seq),len(found_PESE)])
            
            exon_inclusion_vs_count_list.append([aggregate_exon_dict[exon_id]['count'],len(found_PESE)*100/len(seq)])

        if total_seq_len != 0:
            found_rate_100_bp = 100*count_PESE/total_seq_len
        else:
            found_rate_100_bp =0
    
    return found_rate_100_bp, exon_len_vs_count_list, exon_inclusion_vs_count_list, exon_fraction_covered_list






def get_PESE_ESE_count_dict(exon_id_list, genome_fasta, PESE):
    ESE_counts_dict = PESE.generate_ESE_empty_counts_dict()
    
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        if ex.chrom == 'chrM':
            continue
        try: 
            seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        except:
            print(exon_id)
            seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        found_PESE = PESE.score_seq(seq)
    
        if found_PESE == -1:
            continue
        for result in found_PESE:
            PESE.add_results_to_ESE_counts_dict([result], ESE_counts_dict)
            
    
    return  ESE_counts_dict






def get_exon_id_ESE_scores(exon_id_list, ESE_scorer, genome_fasta,aggregate_exon_dict, exon_count_threshold):
    
    exon_id_list = el.exon_id_intersection(exon_id_list,primary_3ss_exon_id_set)
    
    exon_ids_1000_plus = el.threshold_exon_ids(exon_id_list,exon_count_threshold,aggregate_exon_dict)
    
    ESE_100_bp, exon_len_vs_count_list, exon_inclusion_vs_count_list, exon_fraction_covered_list   =  get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_ids_1000_plus, ESE_scorer, genome_fasta,aggregate_exon_dict)
    
    return exon_ids_1000_plus, ESE_100_bp, exon_len_vs_count_list, exon_inclusion_vs_count_list, exon_fraction_covered_list










def plot_2d_hist_ESE(x,y,xlabel,ylabel,title,x_lower_lim,y_lower_lim,num_bins, pdf_stats, PESE_name):
        if type(y_lower_lim) == type(list) and  len(y_lower_lim)==2:
            y_upper_lim = y_lower_lim[1]
            y_lower_lim = y_lower_lim[0]
        else:
            y_upper_lim=70
        x_c = np.log10(x)
        y_c = (np.array(y) ) #*np.array(z))
        nx, ny = (num_bins, num_bins) #31
        xb = np.linspace(x_lower_lim, 5, nx)
        yb = np.linspace(y_lower_lim, 100, ny)
        
        
        heatmap, xedges, yedges = np.histogram2d(x_c, y_c, bins=(xb,yb))
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        
        
        
        bins_step = 0.166
        min_data_point=10**20
        for ii, val in enumerate(x_c):
            if x_c[ii] < min_data_point:
                min_data_point = x_c[ii]
        
        
        x_lower_lim = min_data_point
        
        x_bins = np.arange(min_data_point,5,bins_step)
        y_vals = list()
        for bb in x_bins:
            y_vals.append(list())
            
        for ii, val in enumerate(x_c):
            for jj, x_val in enumerate(x_bins):
            
                if  x_c[ii] >= x_val and x_c[ii] < x_val + bins_step:
                    y_vals[jj].append(y_c[ii])
        
        
        data_list = list()
        
        for ii, y_list in enumerate(y_vals):
            if len(y_list) != 0:
                median = np.median(y_list)
        
                quantiles = np.quantile(y_list,[0.25,0.75])
                spread = np.quantile(y_list,[0.1,0.9])
                data_list.append([x_bins[ii],median,quantiles[0],quantiles[1],spread[0],spread[1]])

        data_x,data_y,quantile_1,quantile_2,spread_1,spread_2 = zip(*data_list)

        return data_x,data_y
        
    





















