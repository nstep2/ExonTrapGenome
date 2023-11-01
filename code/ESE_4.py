

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





This plot ignores exons on chromosome 'ChrM'

"""
import copy

import exon_id_library.exon_id_lib as el


from matplotlib.backends.backend_pdf import PdfPages

exon_count_threshold = 100

PESE_name = 'chasin_ESEseq'


path_name = exp_output_path.ESE_scan+'ESE_scan_ESE_database_%s_exon_count_thresh_%d_%d_reads.pdf' % (PESE_name,exon_count_threshold,exon_count_build)
pdf_stats = PdfPages(path_name)


if 'avg_count' not in aggregate_exon_dict[list(aggregate_exon_dict.keys())[0]]:
    
    for exon_id in aggregate_exon_dict:
        exon = aggregate_exon_dict[exon_id]
        denom = sum(exon['lib_array']>5)
        exon['avg_count'] = exon['count']/denom
    




ESESeq_filepath = exp_output_path.ESE_folder + 'ESEseq/Supplemental_Table_1_main.txt'



rescue_ESE_filepath = exp_output_path.ESE_folder + 'rescue-ese_wayback_machine.txt'
rescue_ESE=dict()
with open(rescue_ESE_filepath,'r') as f:
    for line in f:
        line = line.strip()
        rescue_ESE[line.upper()] = True



chasin_ESE_filepath=exp_output_path.ESE_folder + 'chasin_8mer_scores_wayback_machine.txt'
chasin_PESE=dict()
chasin_PESS=dict()
PESS_threshold = -2.62
PESE_threshold = 2.62
with open(chasin_ESE_filepath,'r') as f:
    for ii in range(23):
        next(f)

    
    
    for jj, line in enumerate(f):
        line = line.strip()
        line_split = line.split('\t')
        if line_split[0]=='':   #empty lines occur 1-2 times
            continue
        P_index = float(line_split[1])
        octomer=line_split[0].upper()
        if P_index >= PESE_threshold:
            chasin_PESE[octomer]=P_index
        if P_index <= PESS_threshold:
            chasin_PESS[octomer]=P_index
      


ESE_seq_threshold = 0
ESE_seq_score=dict()
with open(ESESeq_filepath,'r') as f:
    for ii in range(10):
        line = next(f)
        #print(line)
        
    
    ESE_score_column = 23 
    
    #ESE_seq_score
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
        fig = plt.figure()
        heatmap, xedges, yedges = np.histogram2d(x_c, y_c, bins=(xb,yb))
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]




        
        #bins_step = 0.1
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
                #print(y_list)
                quantiles = np.quantile(y_list,[0.25,0.75])
                spread = np.quantile(y_list,[0.1,0.9])
                data_list.append([x_bins[ii],median,quantiles[0],quantiles[1],spread[0],spread[1]])
            

        fig, ax = plt.subplots()
        data_x,data_y,quantile_1,quantile_2,spread_1,spread_2 = zip(*data_list)
        ax.plot(data_x,data_y)
        locs,labels = plt.xticks()
        ax.fill_between(data_x, spread_1, spread_2, color='r', alpha=0.1)
        ax.fill_between(data_x, quantile_1, quantile_2, color='b', alpha=.1)
        
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title+'\n25-75%% and 10-90%% quantiles\n%s'%(PESE_name))
        plt.legend(['median'])
        
        
        xt = ax.get_xticks()
        new_labels = [ insert_commas(int(10**float(locs[jj]) )) for jj,L in enumerate(locs)]
        plt.xticks(locs,new_labels)
        plt.ylim(y_lower_lim,y_upper_lim)
        plt.xlim(x_lower_lim,5)
        plt.tight_layout()
        pdf_stats.savefig(fig)
        
        
        return data_x,data_y




class burge_ESE_class():
    def __init__(self, hex_list,PESE_name):
        self.PESE = set(hex_list)
        self.threshold = 1
        self.motif_len = len(hex_list[0])
        self.name = PESE_name
    def score_seq(self, query_seq):
        seq = query_seq.upper()
        if seq.find('N') != -1:
            return -1
        result_list = list()
        motif_window = 40

        base_positions = range(0,len(seq)-self.motif_len)

        for ii in base_positions:
            hex_seq = seq[ii:ii+self.motif_len]
            if hex_seq in self.PESE:

                result_list.append([ii, 1, hex_seq])
        return result_list
    def generate_ESE_empty_counts_dict(self):
        ESE_empty_counts_dict = dict()
        for ese_nmer in self.PESE:
            ESE_empty_counts_dict[ese_nmer] = 0
        return ESE_empty_counts_dict
    def add_results_to_ESE_counts_dict(self, result_list, ESE_counts_dict):
        for result in result_list:
            ESE_nmer = result[2]
            ESE_counts_dict[ESE_nmer] += 1





test_seq ='acacagagaatatatgtgtgcgcgatcagttacgagtcgtacgctagctagtcgtatcgagTTTTTGACtcatgcgtagTCGAGCTGtcgacgtagtcgtatgcgtatgcacgtgtacacagatagatatatgacacacacagatatgtgtctactatcatcggcttcgggtagcTCGAGCTGgtatgcgctcgatcgcggcgcgatatctggcattgcgatcgagcgatcgatcgagctgtacTTTTTGACTTTTTGAC'
test_octomer='TCGAGCTG'



class PESE_class():
    def __init__(self, chasin_PESS,PESE_threshold, PESE_name):
        self.PESE = chasin_PESS
        self.threshold = PESE_threshold
        self.motif_len = len(list(chasin_PESS.keys())[0])
        self.name = PESE_name
    def score_seq(self, query_seq):
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
        #outer
        #base_positions= set(list(base_positions_left) + list(base_positions_right)) 
        # middle
        #base_positions = range(len(seq[:motif_window]),                               len(seq)-len(seq[-motif_window:])-self.motif_len)
        for ii in base_positions:
            oct_seq = seq[ii:ii+self.motif_len]
            if oct_seq in self.PESE:
                score = self.PESE[oct_seq]
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
    
    




#initial PESE scorer
PESE = PESE_class(chasin_PESE,PESE_threshold,PESE_name)
intron_interior_exon_ids_10000_plus = el.threshold_exon_ids(region_exon_id_sets_dict['Intronic'],10000,aggregate_exon_dict)




try:
    type(ESE_seq)
except:
    ESE_seq = PESE_class(ESE_seq_score,ESE_seq_threshold, 'ESEseq')    
    
    
##initial PESE scorer
#if PESE_name == 'chasin_PESE':
#    PESE = PESE_class(chasin_PESE,PESE_threshold,PESE_name)

#Alternate PESE scorers - 6mers
if PESE_name == 'chasin_ESEseq':
    PESE=ESE_seq

#burge fairbrother
#if PESE_name == 'burge':
#    PESE=burge_ESE  


pickle_path = exp_output_path.pickle_main
with open(exp_output_path.pickle_main+'chasin_PESE.pickle', 'wb') as outfile:
    pickle.dump(PESE_class(chasin_PESE,0,'chasin_PESE'),outfile)

with open(exp_output_path.pickle_main+'chasin_ESEseq.pickle', 'wb') as outfile:
    pickle.dump(PESE_class(ESE_seq_score,0, 'ESEseq') ,outfile)

#with open(exp_output_path.pickle_main+'burge_ESE.pickle', 'wb') as outfile:
#    pickle.dump(burge_ESE,outfile)






def get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta):
    count_PESE = 0
    total_seq_len = 0
    exon_len_vs_count_list = list()
    for exon_id in exon_id_list:
        if el.exon_id_values(exon_id).chrom == 'chrM':
            continue
        
            
        seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        found_PESE = PESE.score_seq(seq)
        if found_PESE != -1:
            total_seq_len+=len(seq)
            count_PESE += len(found_PESE)
            exon_len_vs_count_list.append([len(seq),len(found_PESE)])
        
        if total_seq_len != 0:
            found_rate_100_bp = 100*count_PESE/total_seq_len
        else:
            found_rate_100_bp =0
    
    return found_rate_100_bp, exon_len_vs_count_list











def get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_id_list, PESE, genome_fasta, aggregate_exon_dict):
    exon_id_list = el.exon_id_intersection(exon_id_list, aggregate_exon_dict.keys())
    count_PESE = 0
    #count_PESE_list = list()
    total_seq_len = 0
    exon_len_vs_count_list = list()
    exon_inclusion_vs_count_list = list()
    for exon_id in exon_id_list:
        seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        
        found_PESE = PESE.score_seq(seq)
        if found_PESE != -1:
            total_seq_len+=len(seq)
            count_PESE += len(found_PESE)
            
            exon_len_vs_count_list.append([len(seq),len(found_PESE)])

            exon_inclusion_vs_count_list.append([aggregate_exon_dict[exon_id]['count'],len(found_PESE)*100/len(seq)])
        
        
        
        if total_seq_len != 0:
            found_rate_100_bp = 100*count_PESE/total_seq_len
        else:
            found_rate_100_bp =0
    
    return found_rate_100_bp, exon_len_vs_count_list, exon_inclusion_vs_count_list






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







'''

if True == False:
    
    ESE_counts_dict = get_PESE_ESE_count_dict(region_exon_id_sets_dict['mRNA'], genome_fasta, PESE)
    
    
    offset_ESE_counts_dict = get_PESE_ESE_count_dict(pc_offset_250_bp, genome_fasta, PESE)
    
    ESE_counts_dict = get_PESE_ESE_count_dict(region_exon_id_sets_dict['Intergenic'], genome_fasta, PESE)
    
    
    offset_ESE_counts_dict = get_PESE_ESE_count_dict(no_overlap_offset_250_bp, genome_fasta, PESE)
    
    
    
    
    ESE_list_results = [(key, ESE_counts_dict[key]) for key in ESE_counts_dict]
    
    offset_ESE_list_results = [(key, offset_ESE_counts_dict[key]) for key in offset_ESE_counts_dict]
    
    ESE_list_results=sorted(ESE_list_results, key=lambda x:x[0]) #sort by kmer
    offset_ESE_list_results=sorted(offset_ESE_list_results, key=lambda x:x[0]) #sort by kmer
    
    x,y = zip(*ESE_list_results)
    print(sum(y))
    
    plt.figure()
    plt.hist(y,bins=(1000))
    plt.xlim(0,7000)
    plt.ylim(0,175)
    plt.legend(['median %.1f' % (np.median(y))])
    plt.title('25%%-75%% = [%.1f, %.1f]' % (np.percentile(y,25), np.percentile(y,75)))
    
    x,y = zip(*offset_ESE_list_results)
    print(sum(y))
    
    plt.figure()
    plt.hist(y,bins=(1000))
    plt.xlim(0,7000)
    plt.ylim(0,175)
    plt.legend(['median %.1f' % (np.median(y))])
    plt.title('25%%-75%% = [%.1f, %.1f]' % (np.percentile(y,25), np.percentile(y,75)))
    
    
    ESE_list_results[:10]
    offset_ESE_list_results[:10]
    
    x,y_main = zip(*ESE_list_results)
    x,y_offset = zip(*offset_ESE_list_results)
    ratio_main_offset = list(zip(x,(np.array(y_main)/(np.array(y_offset) +1))))
    difference_main_offset = list(zip(x,(np.array(y_main) - np.array(y_offset))))
    
    x,y = zip(*ratio_main_offset)
    
    plt.figure()
    plt.hist(np.log2(y),bins=(1000))
    #plt.xscale('log')
    #plt.xlim(0,7000)
    #plt.ylim(0,175)
    plt.legend(['median %.1f' % (np.median(y))])
    plt.title('25%%-75%% = [%.1f, %.1f]' % (np.percentile(y,25), np.percentile(y,75)))
    
    
    
    x,y = zip(*difference_main_offset)
    
    plt.figure()
    plt.hist(y,bins=(1000))
    #plt.xscale('log')
    #plt.xlim(0,7000)
    #plt.xlim(-1000,1000)
    #plt.ylim(0,175)
    plt.legend(['median %.1f' % (np.median(y))])
    plt.title('25%%-75%% = [%.1f, %.1f]' % (np.percentile(y,25), np.percentile(y,75)))
    
    
    #ratio_main_offset
    #difference_main_offset
    x,y = zip(*difference_main_offset)
    thresh_75 = np.percentile(y,50)
    top_percent = [val[0] for val in difference_main_offset if val[1] >= thresh_75]
    
    #ratio_main_offset
    
    
    PESE = burge_ESE_class(top_percent,'top 50% PESE')

'''





exon_id_list = region_exon_id_sets_dict['lncRNA']
lncRNA_PESE_100_bp, lncRNA_exon_len_vs_count_list, lncRNA_exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_id_list, PESE, genome_fasta,aggregate_exon_dict)



exon_id_list = region_exon_id_sets_dict['mRNA']
pc_PESE_100_bp, pc_exon_len_vs_count_list, pc_exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_id_list, PESE, genome_fasta,aggregate_exon_dict)





exon_id_list = region_exon_id_sets_dict['mRNA']
exon_id_list = el.exon_id_intersection(exon_id_list,aggregate_exon_dict.keys())
pc_exon_ids_1000_plus = el.threshold_exon_ids(exon_id_list,exon_count_threshold,aggregate_exon_dict)
pc_1000_PESE_100_bp, pc_1000_exon_len_vs_count_list, exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(pc_exon_ids_1000_plus, PESE, genome_fasta,aggregate_exon_dict)



exon_id_list = region_exon_id_sets_dict['lncRNA']
exon_id_list = el.exon_id_intersection(exon_id_list,aggregate_exon_dict.keys())
lncRNA_exon_ids_1000_plus = el.threshold_exon_ids(exon_id_list,exon_count_threshold,aggregate_exon_dict)
lncRNA_1000_PESE_100_bp, lncRNA_1000_exon_len_vs_count_list, exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(lncRNA_exon_ids_1000_plus, PESE, genome_fasta,aggregate_exon_dict)






exon_id_list = region_exon_id_sets_dict['Intronic']
dummy, dummy, intron_interior_exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_id_list, PESE, genome_fasta,aggregate_exon_dict)


intron_interior_exon_ids_1000_plus = el.threshold_exon_ids(region_exon_id_sets_dict['Intronic'],exon_count_threshold,aggregate_exon_dict)
exon_id_list = intron_interior_exon_ids_10000_plus
intron_interior_PESE_100_bp, intron_interior_1000_exon_len_vs_count_list, exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_id_list, PESE, genome_fasta,aggregate_exon_dict)




exon_id_list = region_exon_id_sets_dict['Intergenic']
dummy, dummy, no_overlap_exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_id_list, PESE, genome_fasta,aggregate_exon_dict)


no_overlap_exon_ids_1000_plus = el.threshold_exon_ids(region_exon_id_sets_dict['Intergenic'],exon_count_threshold,aggregate_exon_dict)
exon_id_list = no_overlap_exon_ids_1000_plus
no_overlap_PESE_100_bp, no_overlap_1000_exon_len_vs_count_list, exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(exon_id_list, PESE, genome_fasta,aggregate_exon_dict)



def generate_250_bp_offset_exon_ids(exon_id_list, offset):
    offset_exon_id_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        new_ex = "%s:%d-%d:%s" % (ex.chrom,ex.start+offset,ex.end+offset,ex.strand)
        offset_exon_id_list.append(new_ex)
    return offset_exon_id_list


no_overlap_offset_250_bp = generate_250_bp_offset_exon_ids(no_overlap_exon_ids_1000_plus,250)


intron_interior_offset_250_bp = generate_250_bp_offset_exon_ids(intron_interior_exon_ids_1000_plus,250)


pc_offset_250_bp = generate_250_bp_offset_exon_ids(region_exon_id_sets_dict['mRNA'],250)


lncRNA_offset_250_bp = generate_250_bp_offset_exon_ids(region_exon_id_sets_dict['lncRNA'],250)



pc_1000_offset_250_bp = generate_250_bp_offset_exon_ids(pc_exon_ids_1000_plus,250)

pc_offset_250_bp = generate_250_bp_offset_exon_ids(region_exon_id_sets_dict['mRNA'],250)


lncRNA_1000_offset_250_bp = generate_250_bp_offset_exon_ids(lncRNA_exon_ids_1000_plus,250)



exon_id_list = no_overlap_offset_250_bp
no_overlap_offset_exon_ids_1000_plus, no_overlap_offset_exon_len_vs_count_list  = get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta)



exon_id_list = intron_interior_offset_250_bp
intron_interior_offset_exon_ids_1000_plus, intron_interior_offset_exon_len_vs_count_list  = get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta)





exon_id_list = pc_offset_250_bp
pc_all_offset_exon_ids, pc_all_offset_exon_len_vs_count_list  = get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta)


exon_id_list = lncRNA_offset_250_bp
lncRNA_offset_exon_ids, lncRNA_offset_exon_len_vs_count_list  = get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta)




exon_id_list = pc_1000_offset_250_bp
pc_offset_exon_ids_1000_plus, pc_1000_offset_exon_len_vs_count_list  = get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta)


exon_id_list = lncRNA_1000_offset_250_bp
lncRNA_offset_exon_ids_1000_plus, lncRNA_1000_offsetexon_len_vs_count_list  = get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta)
















#try:
if 'a'=='a':
    HEXEvent_exon_id_list_constitutive
    HEXEvent_exon_id_list_cassette
    
    
    
    
    exon_id_list = HEXEvent_exon_id_list_cassette
    cassette_exon_ids = el.exon_id_intersection(exon_id_list,aggregate_exon_dict.keys())
    #cassette_exon_ids_1000_plus = el.threshold_exon_ids(exon_id_list,100,aggregate_exon_dict)
    cassette_PESE_100_bp, cassette_PESE_exon_len_vs_count_list, cassette_exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(cassette_exon_ids, PESE, genome_fasta,aggregate_exon_dict)
    
    
    exon_id_list = HEXEvent_exon_id_list_constitutive
    constitutive_exon_ids = el.exon_id_intersection(exon_id_list,aggregate_exon_dict.keys())
    #cassette_exon_ids_1000_plus = el.threshold_exon_ids(exon_id_list,100,aggregate_exon_dict)
    constitutive_PESE_100_bp, constitutive_exon_len_vs_count_list, constitutive_exon_inclusion_vs_count_list  = get_PESE_per_100_bp_vs_inclusion_exon_id_list(constitutive_exon_ids, PESE, genome_fasta,aggregate_exon_dict)
        
        

    HexEvent_bar_list = list()
    HexEvent_bar_list.append(pc_PESE_100_bp)
    HexEvent_bar_list.append(constitutive_PESE_100_bp)
    HexEvent_bar_list.append(cassette_PESE_100_bp)
    HexEvent_bar_list.append(lncRNA_PESE_100_bp)
    
    



'''
fig = plt.figure()
plt.title('median ESE count vs exon expression\n%s'%(PESE_name))
plt.ylabel('ESE count')
plt.xlabel('exon expression counts')
plt.plot(pce_plot[0],pce_plot[1])

plt.plot(case_plot[0],case_plot[1])
plt.plot(lnc_plot[0],lnc_plot[1])
plt.plot(no_plot[0],no_plot[1])
plt.plot(ii_plot[0],ii_plot[1])

plt.legend(['pc','cassette','lncRNA','intergenic','intronic'])
plt.xticks(range(0,6),[insert_commas(10**x)for x in range(0,6)])
plt.xlim(2,5)
plt.ylim(0,40)
plt.tight_layout()
pdf_stats.savefig(fig)
'''
    




bar_25_75_list = list()
dummy,m = zip(*pc_exon_len_vs_count_list)
bar_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75)])
dummy,m = zip(*lncRNA_exon_len_vs_count_list)
bar_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75)])
dummy,m = zip(*intron_interior_1000_exon_len_vs_count_list)
bar_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75)])
dummy,m = zip(*no_overlap_1000_exon_len_vs_count_list)
bar_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75)])


bar_offset_25_75_list = list()


dummy,m = zip(*pc_1000_offset_exon_len_vs_count_list)
bar_offset_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75)])


dummy,m = zip(*lncRNA_offset_exon_len_vs_count_list)
bar_offset_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75)])

dummy,m = zip(*intron_interior_offset_exon_len_vs_count_list)
bar_offset_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75)])


dummy,m = zip(*no_overlap_offset_exon_len_vs_count_list)
bar_offset_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75)])



annotated_bar_list = list()
annotated_bar_list.append(pc_PESE_100_bp)
annotated_bar_list.append(lncRNA_PESE_100_bp)

count_1000_bar_list = list()
count_1000_bar_list.append(pc_1000_PESE_100_bp)
count_1000_bar_list.append(lncRNA_1000_PESE_100_bp)
count_1000_bar_list.append(intron_interior_PESE_100_bp)
count_1000_bar_list.append(no_overlap_PESE_100_bp)


offset_1000_bar_list = list()
offset_1000_bar_list.append(pc_offset_exon_ids_1000_plus)
offset_1000_bar_list.append(lncRNA_offset_exon_ids_1000_plus)
offset_1000_bar_list.append(intron_interior_offset_exon_ids_1000_plus)
offset_1000_bar_list.append(no_overlap_offset_exon_ids_1000_plus)








m1,m2,m3 = zip(*bar_25_75_list)
n1,n2,n3 = zip(*bar_offset_25_75_list)

m1,m2,m3 = np.array(m1),np.array(m2),np.array(m3)
n1,n2,n3 = np.array(n1),np.array(n2),np.array(n3)


x = np.arange(len(bar_offset_25_75_list))  # the label locations
width = 0.35  # the width of the bars



fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, m2, width, label='exon',yerr = [m2-m1,m3-m2])
rects2 = ax.bar(x + width/2, n2, width, label='offset',yerr = [n2-n1,n3-n2])
plt.xticks(range(4),['pc','lncRNA','intron_interior', 'intergenic'])
plt.ylabel('PESE counts per 100 bp\n%s octomers with threshold %.2f' % (insert_commas(len(PESE.PESE)),PESE.threshold))
plt.title('offset right 250 bases & threshold %s PESE counts\n%s, exon_thresh: %s' % (insert_commas(exon_count_threshold),PESE_name,insert_commas(exon_count_threshold)))
plt.ylim(0,max(m3)+5)
ax.legend()
plt.tight_layout()
pdf_stats.savefig(fig)




























pdf_stats.close()
























pc_exon_lengths = [el.exon_id_values(exon_id).length for exon_id in region_exon_id_sets_dict['mRNA']]
upper_search_window = int(np.quantile(pc_exon_lengths,.9))
lower_search_window = int(np.quantile(pc_exon_lengths,.1))


def get_inclusion_vs_ESE_count_for_MES_range(exon_id_list, aggregate_exon_dict, MES_range_3ss, MES_range_5ss, PESE):
    
    exon_id_list = el.size_exon_id_list(exon_id_list, 63,222) # get 10-90% mRNA exon lengths range 
    
    exon_id_list = el.exon_id_intersection(exon_id_list, aggregate_exon_dict.keys())
    
    gc_list=list()
    exon_count_vs_PESE_list = list()
    exon_countper_len_vs_PESE_list = list()
    exon_avg_count_vs_PESE_list = list()
    exon_id_MES_range_list=list()
    for exon_id in exon_id_list:
        exon = aggregate_exon_dict[exon_id]
        if exon['5ss_score'] <= MES_range_5ss[0] or exon['5ss_score'] >= MES_range_5ss[1]:
            continue
        if exon['3ss_score'] <= MES_range_3ss[0] or exon['3ss_score'] >= MES_range_3ss[1]:
            continue
        exon_id_MES_range_list.append(exon_id)
    
        #get_PESE_ESE_count_dict([exon_id],aggregate_exon_dict, PESE)
        seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        
        found_PESE = PESE.score_seq(seq)
    
        if found_PESE == -1:
            continue
        
        exon_count_vs_PESE_list.append([len(found_PESE), exon['count']])
        exon_avg_count_vs_PESE_list.append([len(found_PESE), exon['avg_count']])
        exon_countper_len_vs_PESE_list.append([len(found_PESE)/exon['length']*100, exon['count']])
        
        gc_list.append(el.get_gc_for_region(exon_id,genome_fasta))

    return exon_count_vs_PESE_list,  exon_avg_count_vs_PESE_list, exon_countper_len_vs_PESE_list, gc_list
















def plot_ESE_linear_y(x,y, title, xlabel, ylabel, x_bin_lims, y_bin_lims, *args, **kwargs):
    
    write_to_pdf_flag = 'pdf_handle' in kwargs
    'pdf_handle'
    
    x_c = np.log10(x)
    y_c = np.array(y)  

    bins_step = .25
    
    min_data_point = 2
    x_lower_lim = min_data_point
    
    x_bins_1=np.arange(min_data_point,x_bin_lims[1],x_bin_lims[3])
    x_bins_2=np.arange(x_bin_lims[1],x_bin_lims[2],x_bin_lims[4])
    
    x_bins = np.concatenate([x_bins_1,x_bins_2])

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
        
    
    
    fig, ax = plt.subplots()
    data_x,data_y,quantile_1,quantile_2,spread_1,spread_2 = zip(*data_list)
    ax.plot(data_x,data_y)
    locs,labels = plt.xticks()
    ax.fill_between(data_x, spread_1, spread_2, color='r', alpha=0.1)
    ax.fill_between(data_x, quantile_1, quantile_2, color='b', alpha=.1)
    
    plt.title('%s' % (title))
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.tight_layout()
    
    if write_to_pdf_flag == True:
        kwargs['pdf_handle'].savefig(fig)









def plot_ESE_linear_x(x,y, title, xlabel, ylabel, x_bin_lims, y_bin_lims, *args, **kwargs):
    
    write_to_pdf_flag = 'pdf_handle' in kwargs
    'pdf_handle'
    
    x_c = np.array(x)
    y_c = np.array(y)  

    bins_step = 4
    
    
    
    min_data_point = 0
    x_lower_lim = min_data_point
    
    x_bins_1=np.arange(min_data_point,x_bin_lims[1],x_bin_lims[3])
    x_bins_2=np.arange(x_bin_lims[1],x_bin_lims[2],x_bin_lims[4])
    
    x_bins = np.concatenate([x_bins_1,x_bins_2])
    
    
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
        

    fig, ax = plt.subplots()
    data_x,data_y,quantile_1,quantile_2,spread_1,spread_2 = zip(*data_list)
    ax.plot(data_x,data_y)
    locs,labels = plt.xticks()
    ax.fill_between(data_x, spread_1, spread_2, color='r', alpha=0.1)
    ax.fill_between(data_x, quantile_1, quantile_2, color='b', alpha=.1)
    plt.yscale('log', base=2)
    plt.title('%s' % (title))
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.tight_layout()
    
    if write_to_pdf_flag == True:
        kwargs['pdf_handle'].savefig(fig)











MES_range_3ss = [0,6]
MES_range_5ss = [0,6]

exon_id_list = region_exon_id_sets_dict['Intergenic']





def plot_MES_range_ESE_vs_inclusion(exon_id_list, aggregate_exon_dict, PESE, PESE_name, exon_id_list_name,  *args, **kwargs):
    
    write_to_pdf_flag = 'pdf_handle' in kwargs
    #print('pdf_handle',write_to_pdf_flag)
    #args['pdf_handle']
    
    median_GC_pair_MES_list = list()
    median_ESE_vs_inclusion_list = list()
    
    for ii, pair in enumerate([[0,4],[4,6],[6,8],[8,10],[10,15]]):
        MES_range_3ss = pair
        MES_range_5ss = pair
        
        
        exon_count_vs_PESE_list,  exon_avg_count_vs_PESE_list, exon_countper_len_vs_PESE_list, gc_list = get_inclusion_vs_ESE_count_for_MES_range(exon_id_list, aggregate_exon_dict, MES_range_3ss, MES_range_5ss, PESE)
        
        x,y = zip(*exon_countper_len_vs_PESE_list)
        
        
        
        x_data = x
        y_data = y
        
        x_bin_lims = [0,50,120,2,10] # lower, mid, upper, step lower, step upper
        y_bin_lims = [0, max(x_data)]
        
        
        
        x_data = y
        y_data = x
        
        y_bin_lims = [0,60] # lower, mid, upper, step lower, step upper
        x_bin_lims = [2, 4, 5.1, .2, .2]
        
        
        
        
        median_ESE_vs_inclusion_list.append([ii, pair,np.median(x),np.median(y)])
        median_GC_pair_MES_list.append([ii, pair,np.median(gc_list), np.mean(gc_list), np.percentile(gc_list,25), np.percentile(gc_list,75)])
    
    index, pair_list, med_x, med_y = zip(*median_ESE_vs_inclusion_list)
    
    
    title = exon_id_list_name
    fig, ax = plt.subplots()
    plt.scatter(med_y, med_x)
    plt.title(title)
    plt.xlabel('read counts')
    plt.xscale('log', base=10)
    plt.ylabel('median ESE')
    for ii, pair in enumerate(pair_list):
        ax.annotate('  %s' % str(pair), (med_y[ii], med_x[ii]))
    plt.xlim(100,100000)
    plt.tight_layout()
    if write_to_pdf_flag == True:
        kwargs['pdf_handle'].savefig(fig)
    
    
    
    
    return median_ESE_vs_inclusion_list, median_GC_pair_MES_list








exon_id_list_list = list()
exon_id_list_list.append(['pc_middle_exon_id_list', region_exon_id_sets_dict['mRNA']])
exon_id_list_list.append(['lncRNA_middle_exon_ids', region_exon_id_sets_dict['lncRNA']])
exon_id_list_list.append(['no_overlap_set', region_exon_id_sets_dict['Intergenic']])
exon_id_list_list.append(['intron_interior_set', region_exon_id_sets_dict['Intronic']])
exon_id_list_list.append(['antisense_transcript_set', region_exon_id_sets_dict['Antisense']])






#load exon_id list of exons associated with repeats
pickle_path = exp_output_path.Dfam_pickle + "Dfam_pickle.pickle"
with open(pickle_path, "rb") as input_file:
    exon_ids_overlapping_repeat_list=pickle.load(input_file)







#PESE_name=PESE.name
ESE_count_vs_inclusion_group = list()
ESE_gc_group = list()
pdf_ESE = PdfPages(exp_output_path.ESE_scan+'ESE_scan_MES_scores_%s_%d_reads.pdf' % (PESE.name,exon_count_build))
for ii, lib_pair in enumerate(exon_id_list_list):
    exon_id_list = set(lib_pair[1]).difference(exon_ids_overlapping_repeat_list)
    exon_id_list_name = lib_pair[0]

    median_ESE_vs_inclusion_list, median_GC_pair_MES_list = plot_MES_range_ESE_vs_inclusion(exon_id_list, aggregate_exon_dict, PESE, PESE.name, exon_id_list_name, pdf_handle = pdf_ESE)
    ESE_count_vs_inclusion_group.append([ii, exon_id_list_name, median_ESE_vs_inclusion_list])
    ESE_gc_group.append([ii, exon_id_list_name, median_GC_pair_MES_list])
    
    
    





index, names,ESE_vs_inclusion = zip(*ESE_count_vs_inclusion_group)
dummy,MES_range, dummy, dummy = zip(*ESE_vs_inclusion[0])
fig = plt.figure(figsize=(7,5))
for ii, val in enumerate(ESE_vs_inclusion):
    val_list = list()
    for jj, dummy in enumerate(ESE_vs_inclusion):
        val_list.append(ESE_vs_inclusion[jj][ii][2])
    plt.bar(np.arange(0,len(index))+ii*.15, val_list,width = .15)
plt.bar([len(index)],[0])
plt.legend(MES_range)
plt.xticks(index, names, rotation=30)
plt.title('ESE counts in different maxentscan ranges (legend)\n%s' % (PESE.name))
plt.ylabel('Median ESE count')
plt.tight_layout()
pdf_ESE.savefig(fig)





#1



out_path = exp_output_path.ESE_scan + 'S6B.txt'

with open(out_path,'w') as f:
    
    f.write( '{:}\t{:}\t{:}\t{:}\t{:}\t{:}\n'.format('category', '[0,4]', '[4,6]', '[6,8]', '[8,10]', '[10,15]') )
    
    for ii, val in enumerate(ESE_count_vs_inclusion_group):
        #val_list = list()
        print(val[1])
        f.write(val[1])
        for jj, val2 in enumerate(val[2]):
            print('\t{:}'.format(val2[2]))
            f.write('\t{:}'.format(val2[2]))
        f.write('\n')
    
    

















pdf_ESE.close()
































