

spliceAI_path = exp_output_path.spliceAI_model_path




outdir = exp_output_path.Compute_spliceAI_exon_id_scratch 
pdf_plots = PdfPages(outdir+'Compute_spliceAI_exon_id_scratch%d_reads.pdf' % (exon_count_build))



import time
import math


int('caution')




import tensorflow.compat.v1 as tf

tf.config.threading.set_intra_op_parallelism_threads(6)
tf.config.threading.set_inter_op_parallelism_threads(6)



from tensorflow.compat.v1.keras.models import load_model
tf.disable_v2_behavior()

from tensorflow.keras.models import load_model





models = list()

for ii in range(5):
    models.append(load_model('%sspliceai%d.h5' % (spliceAI_path, ii+1)) )




import scipy
import numpy
a = ""
num_to_ch = {1:'A', 2:'C',3:'G',4:'T'}
for i in range(11001):
    #a[i, int(scipy.rand(1)*100)%4 ] = 1
    a = a + num_to_ch[int(numpy.random.rand(1)*100)%4+1 ]

seq = a

import numpy as np
def one_hot_encode(seq):

    map = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')

    return map[np.fromstring(seq, np.int8) % 5]




try:
    genome_fasta['chr1'][1000:2000]
except:
    import pyfaidx
    input_fasta = '/home/pineapple/work/indexes/downloaded/Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
    genome_fasta = pyfaidx.Fasta(input_fasta)



#models[0].compile(optimizer="Adam", loss="mse", metrics=["mae"])



'''
########### T4 reporter seq
'''
up_seq = "CGTTTAGTGAACCGTCAGATCGCCTGGAGACGCCATCCACGCTGTGCTAGNTTTCCGGTAAGCTTACAAGTTTGTACAAAAAAGCAGGCTCAATGGCATCATATCCCTGTCATCAACACGCGAGTGCTTTTGACCAAGCAGCCCGGTCTCGAGGACATAATAACCGACGCACAGCTCTTCGCCCTCGACGGCAACAGGAGGCGACAGAGGTCCGACCGGAGCAAAAGATGCCTACGCTGCTTCGGGTATATATTGACGGCCCGCATGGTATGGGCAAGACGACAACTACGCAACTGCTGGTCGCGTTGGGTAGCCGAGATGACATCGTCTATGTCCCAGAGCCGATGACATACTGGCGGGTGCTCGGAGCAAGCGAGACTATAGCGAACATATACACAACTCAACATCGACTTGATCAAGGCGAAATATCCGCCGGCGATGCGGCGGTGGTGATGACGTCAGCGCAAATTACGATGGGTATGCCCTACGCCGTTACAGATGCGGTTCTTGCTCCTCACATTGGCGGCGAAGCGGGGAGTTCACATGCGCCACCGCCGGCTTTGACGCTTATTTTCGACCGACACCCCATTGCGGCGCTTCTGTGTTACCCGGCGGCAAGGTACCTCATGGGCTCTATGACGCCCCAAGCCGTACTGGCATTTGTCGCCCTTATTCCCCCTACTCTGCCCGGTACGAACATAGTACTTGGAGCGTTGCCTGAAGATAGGCACATTGATCGCCTGGCAAAGAGACAACGGCCCGGCGAACGACTCGATTTGGCTATGCTCGCTGCGATTAGGCGGGTGTACGGGCTGCTGGCGAACACGGTGCGATACCTCCAGTGCGGTGGTAGTTGGAGGGAGGATTGGGGTCAACTCAGCGGCACTGCAGTACCCCCCCAAGGCGCCGAGCCGCAGTCTAATGCCGGCCCACGACCCCACATCGGCGACACACTCTTTACCCTCTTTAGAGCTCCAGAACTGCTTGCTCCCAACGGTGACCTTTACAACGTCTTCGCGTGGGCACTGGACGTCCTGGCGAAGCGACTTCGCTCCATGCACGTGTTCATTCTTGATTACGACCAGTCACCAGCCGGATGCAGgtaagttatgaggtagaaaggtcaacgtctgaatctcagtacagtgtacagagcagcactatgcttaagttttagcctttgttaaaagctttgctttagtacttttactcgagaaatgcagagtagaaagttaaatcagattctaaggttctgtgtcttaatcaatttaatgtgtaactggatttctgttcaaatattctagaatataaatgggatttaaatgataacttggaatattaataataaaatgtgttttatttttaaacattgaatgaggccaggcgcagtagctcacgcctgtaatcccaacacttcaggaggccaaggcaggcagatcacctgagatccagcctggccaacatggtgaaccccgtctctaccagaaataaaaatttagccggttgtggtggtgcacgcctgtagtcccagctactggcaggtgaattgcttgaacccggcaggcagaggttgcagtgaatttagattgcaccactgcactccactgtaggcaacagagcgagagactcaaaaaaaaaaataaataaatgaatacatagtgaaagcaaaattataggctgaggcaggggaattgcttgagcccgggaatttgagaccagcctgagcagcatagggaacctccacaaaaaatgaggcaggaagccggggaggttgaagctgcagtgagccaagtacgctcccacctgagtctcctcccccagagcaagaccccaccttaaaagaaaaaattatagcctctgtcttagccttttaaaagcaaatttataatggcttataaaagcaaatAAGATacgcgTGACCCTGAATCACGT"

down_seq = "ACAGGATGCTCTATCGCTAGCtcgacttcataagccattttggaaagccaacacttttcaggcagtttcttgtaatctcttaatgcatgcatcgtacactgagggcaggtaaaatggatgccttgaaagagtccagagaataattttaaagatgttaatttagtatgacaagggtatgtggctttttatttttgccaaaaatttggtatttcggggtttgggtagggaaatgtcaatcatagtttaaaaataaagatatctttgagaatacttccaaagataaaaatgaagtaaaattgtgacattttggtcctcccttgtcaaaatgtgtaatttgacaagtcaaataatgtagttagtaagtattgaccacagttgtcattgataatggtattgaggtatatggatgtctattttgtgcattatttggcaatgtttaagaagttaacagttgtcagtcattgtctgaaagtgagtctgtaagtgtattgttcagcattttatcacccatgaaaaccctttttttatttaaagcatcactctattgaatggtactaaaggtacaaatttgttggggctgagatttttttcttagtttgatttttttttttaagagcatttacttaggaatgttctgtgttctaatactacttatttcagccaccagagagaaccctatgccagtttctaaaaagggtggcaaatcacttgacctgcaagaaacatctcattcatgctctatgaaagtttttattttactaaatgaatgacataattttttaaaacttacaaactgaatgcgtgttgtgttgcctttctgtcttcacagAGACGCCCTTTTGCAGTTGACCTCTGGAATGGTTCAGACTCACGTCACGACTCCAGGTTCTATACCAACGATCTGTGATCTCGCTCGAACATTCGCGAGAGAAATGGGTGAGGCCAACACCCAGCTTTCTTGTACAAAGTGGTGTGATCTAGAGGATCCCGGGTGGCTTCGAATTCTCGACCTCGAGACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAGATCCACTTTGGCCGCGGCTCGAGGGGGTTGGGGTTGCGCCTTTTCCAAGGCAGCCCTGGGTTTGCGCAGGGACGCGGCTGCTCTGGGCGTGGTTCCGGGAAACGCAGCGGCGCCGACCCTGGGTCTCGCACATTCTTCACGTCCGTTCGCAGCGTCACCCGGATCTTCGCCGCTACCCTTGTGGGCCCCCCGGCGACGCTTCCTGCTCCGCCCCTAAGTCGGGAAGGTTCCTTGCGGTTCGCGGCGTGCCGGACGTGACAAACGGAAGCCGCACGTCTC"



def one_hot_encode_seq(seq):
    x_ref = one_hot_encode(seq)[None, :] 
    return x_ref



insert = "tgccggtatattccataataaggagaaaatatgggtggttttattttctcttgtttcaagcacattaaccataaaacataaaatgtgtatccttggtattaaaagagaaaaaaatggacgtgggaaacaattactactgtacatagtaatgtctactgcatttctcttcttgcttataagatatttctatttctttgcagCAATTTGGTAAAATCTTAGATGTTGAAATTATTTTTAATGAGCGAGGCTCAAAGgtaagcaacttaattcactttagaattgtttagtttatttttaaatgtttgttatgaataaggggaaacaactgcactgtcttaatcagctcattgtgaaaatacacttcagtgggaattcctttagagaccaccttccgttttcaaaatggtgggaaagacatcgggttcagtttaaaatgtataaggggataattttaaattgctgttttatatatatatcaagtttccagaggaaaattatacatgc"






def get_exon_flank_seq(exon_id, genome_fasta, flank_len):
    
    ex = el.exon_id_values(exon_id)
    
    new_region_id = "%s:%d-%d:%s" % (ex.chrom, ex.start-flank_len, ex.end+flank_len, ex.strand)
    
    seq = genome_fasta[ex.chrom][ex.start-flank_len:ex.end+flank_len]
    
    if ex.strand == '-':
        seq = seq.reverse.complement
    
    pos_3ss = flank_len-1
    pos_5ss = flank_len + ex.length-1
    
    region_seq = str(seq)
    exon_seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
    
    return exon_id, region_seq, exon_seq, pos_3ss, pos_5ss, flank_len
    

flank_len=500






exon_id_list_pairs = list()
exon_id_list = pc_middle_exon_id_list
exon_id_list_name = 'mRNA'
exon_id_list_pairs.append([exon_id_list,exon_id_list_name])
exon_id_list = lncRNA_middle_exon_ids
exon_id_list_name = 'lncRNA'
exon_id_list_pairs.append([exon_id_list,exon_id_list_name])
exon_id_list = no_overlap_set
exon_id_list_name = 'Intergenic'
exon_id_list_pairs.append([exon_id_list,exon_id_list_name])
exon_id_list = intron_interior_set
exon_id_list_name = 'Intronic'
exon_id_list_pairs.append([exon_id_list,exon_id_list_name])
exon_id_list = antisense_transcript_set
exon_id_list_name = 'Antisense'
exon_id_list_pairs.append([exon_id_list,exon_id_list_name])






pickle_path = exp_output_path.Dfam_pickle + "Dfam_pickle.pickle"
with open(pickle_path, "rb") as input_file:
    exon_ids_overlapping_repeat_list=pickle.load(input_file)


for ii, pair in enumerate(exon_id_list_pairs):
    exon_id_list_pairs[ii][0] = list(set(exon_id_list_pairs[ii][0]).difference(exon_ids_overlapping_repeat_list))









exon_id_list_spliceAI_scores_dict = dict()


for exon_id_pair in exon_id_list_pairs:
    exon_id_list = exon_id_pair[0]
    exon_id_list_name = exon_id_pair[1]
    
    
    exon_count_vs_spliceAI_scores_list = list()

    
    exon_id_list = el.exon_id_intersection(exon_id_list, aggregate_exon_dict.keys())
    
    test_set = list(exon_id_list)[:10000]
    for jj, exon_id in enumerate(test_set):
        
        ex = el.exon_id_values(exon_id)
        
        exon_id, region_seq, exon_seq, pos_3ss, pos_5ss, flank_len_used  = get_exon_flank_seq(exon_id, genome_fasta, flank_len)
        fasta_line = '>%s_%d_%d_%d\t%s\n'  % (exon_id, pos_3ss, pos_5ss, flank_len_used, region_seq)
        
        fasta_entry = fasta_line
        
        test_seq = fasta_entry.strip().split('\t')[1]
        pos_3ss = fasta_entry.strip().split('\t')[0].split('_')[1]
        pos_3ss = int(pos_3ss)
        pos_5ss = fasta_entry.strip().split('\t')[0].split('_')[2]
        pos_5ss = int(pos_5ss)
        
        test_seq_pad_1 = up_seq + test_seq + down_seq
        #test_seq = str_seq
        
        left_pad_len  = ( 5500 - len(up_seq)   - len(test_seq)//2 )
        right_pad_len = ( 5501 - len(down_seq) - math.ceil( len(test_seq)/2.0 )  )
        
        test_seq_pad_2 = "N" * left_pad_len + test_seq_pad_1 + "N"*right_pad_len
        
        
        x_ref = one_hot_encode_seq(test_seq_pad_2)
        
        output_list = list()
        for ii, model in enumerate(models):
        
            start_time = time.time()
            output = model.predict(x_ref)
            end_time = time.time() - start_time
            #print("Model #%d:\t%.2f seconds" % (ii+1, end_time))
            output_list.append( np.array(output[0]) )
        
        avg_output =  (output_list[0] + output_list[1] + output_list[2] + output_list[3] + output_list[4])/5
        #predict_threshold = 0.01
        
        if ex.strand == '+':
            adjust_pos_5ss = -1
            adjust_pos_3ss = 0
        else:
            adjust_pos_5ss = 1
            adjust_pos_3ss = 2
            
        pos_5ss_score = int((1001-len(test_seq))/2)+pos_5ss + adjust_pos_5ss
        pos_3ss_score = int((1001-len(test_seq))/2)+pos_3ss+adjust_pos_3ss
        score_5ss = avg_output[pos_5ss_score][2]
        score_3ss = avg_output[pos_3ss_score][1]
        
        max_score_3ss = max(avg_output[pos_3ss_score-2:pos_3ss_score+2,1])
        max_score_5ss  = max(avg_output[pos_5ss_score-2:pos_5ss_score+2,2])
        
        counts = aggregate_exon_dict[exon_id]['count']
        exon_count_vs_spliceAI_scores_list.append([counts,max_score_3ss,max_score_5ss, exon_id])
        
    
    exon_id_list_spliceAI_scores_dict[exon_id_list_name] = exon_count_vs_spliceAI_scores_list
    
    avg_output[pos_3ss_score-2:pos_3ss_score+2,1]
    avg_output[pos_5ss_score-5:pos_5ss_score+5,2]
    
    x,y,z,d = zip(*exon_count_vs_spliceAI_scores_list)

    
    

def plot_2d_hist(x,y,xlabel,ylabel,title,x_lower_lim,y_lower_lim,num_bins,pdf_plots):
    x_c = np.log10(x)
    y_c = np.log10(np.array(y) ) #*np.array(z))
    nx, ny = (num_bins, num_bins) #31
    xb = np.linspace(2, 5, nx)
    yb = np.linspace(y_lower_lim, 0, ny)
    fig = plt.figure()
    heatmap, xedges, yedges = np.histogram2d(x_c, y_c, bins=(xb,yb))
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.clf()
    plt.imshow(heatmap.T, cmap='cividis', extent=extent, origin='lower', aspect='auto') #,norm=LogNorm(vmin=0.01, vmax=1))
    plt.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    pdf_plots.savefig(fig)
    
    bins_step = 0.25
    x_bins = np.arange(2,5,bins_step)
    y_vals = list()
    for bb in x_bins:
        y_vals.append(list())
        
    for ii, val in enumerate(x_c):
        for jj, x_val in enumerate(x_bins):
        
            if  x_c[ii] >= x_val and x_c[ii] < x_val + bins_step:
                y_vals[jj].append(y_c[ii])
    
    
    data_list = list()
    fig = plt.figure()
    for ii, y_list in enumerate(y_vals):
        median = np.median(y_list)
        quantiles = np.quantile(y_list,[0.25,0.75])
        spread = np.quantile(y_list,[0.1,0.9])
        data_list.append([x_bins[ii],median,quantiles[0],quantiles[1],spread[0],spread[1]])
      

    fig, ax = plt.subplots()
    data_x,data_y,quantile_1,quantile_2,spread_1,spread_2 = zip(*data_list)
    ax.plot(data_x,data_y)
    
    ax.fill_between(data_x, spread_1, spread_2, color='r', alpha=0.1)
    ax.fill_between(data_x, quantile_1, quantile_2, color='b', alpha=.1)
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title+'\n25-75%% and 10-90%% quantiles')
    plt.ylim(y_lower_lim,0)
    plt.xlim(x_lower_lim,5)
    ax.plot([0,5],[np.log10(0.2),np.log10(0.2)], linestyle='dashed', alpha =0.5, color='orange')

    plt.legend(['median','log10 of spliceAI 0.2 threshold'])
    pdf_plots.savefig(fig)



pickle_path = exp_output_path.SpliceAI_ET_exons + "SpliceAI_ET_exons.pickle"

with open(pickle_path, "wb") as output_file:
    pickle.dump(exon_id_list_spliceAI_scores_dict, output_file)
    








if "load" == "no":
    from experiment_paths.experiment_paths import *
    import pickle
    
    pickle_path = exp_output_path.SpliceAI_ET_exons + "SpliceAI_ET_exons.pickle"
    
    with open(pickle_path, "rb") as input_file:
        exon_id_list_spliceAI_scores_dict = pickle.load( input_file)
        










pickle_path = exp_output_path.SpliceAI_ET_exons + "SpliceAI_ET_exons.pickle"

with open(pickle_path, "rb") as input_file:
    exon_id_list_spliceAI_scores_dict = pickle.load(input_file)
    





for exon_id_list_name in exon_id_list_spliceAI_scores_dict:
    exon_count_vs_spliceAI_scores_list = exon_id_list_spliceAI_scores_dict[exon_id_list_name]  
    
    
    x,y,z = zip(*exon_count_vs_spliceAI_scores_list)
        
    
    
    
    
    xlabel='exon trapping counts (log10)'
    ylabel='spliceAI score 3ss score (log10)'
    title='%s scatter exon trapping counts vs SpliceAI score' % (exon_id_list_name)
    x_lower_lim = 2
    y_lower_lim = -10
    num_bins=29
    plot_2d_hist(x, np.array(y), xlabel, ylabel, title, x_lower_lim, y_lower_lim, num_bins, pdf_plots)
    
    
    
    xlabel='exon trapping counts (log10)'
    ylabel='spliceAI score 5ss score (log10)'
    title='%s scatter exon trapping counts vs SpliceAI score' % (exon_id_list_name)
    x_lower_lim = 2
    y_lower_lim = -10
    num_bins=29
    plot_2d_hist(x, np.array(z), xlabel, ylabel, title, x_lower_lim, y_lower_lim, num_bins, pdf_plots)
    
    
    
    
    xlabel='exon trapping counts (log10)'
    ylabel='spliceAI score 3ss score (log10)'
    title='%s scatter exon trapping counts vs SpliceAI score\ncounts above 100' % (exon_id_list_name)
    x_lower_lim = 2
    y_lower_lim = -2.5
    num_bins=15
    plot_2d_hist(x, np.array(y), xlabel, ylabel, title, x_lower_lim, y_lower_lim, num_bins, pdf_plots)
    
    xlabel='exon trapping counts (log10)'
    ylabel='spliceAI score 5ss score (log10)'
    title='%s scatter exon trapping counts vs SpliceAI score\ncounts above 100' % (exon_id_list_name)
    x_lower_lim = 2
    y_lower_lim = -2.5
    num_bins=15
    plot_2d_hist(x, np.array(z), xlabel, ylabel, title, x_lower_lim, y_lower_lim, num_bins, pdf_plots)
    
    














pdf_plots.close()








def get_bin_median(x,y,xlabel,ylabel,title,x_lower_lim,y_lower_lim,num_bins,pdf_plots):
    x_c = np.log10(x)
    y_c = np.log10(np.array(y) ) #*np.array(z))
    nx, ny = (num_bins, num_bins) #31
    xb = np.linspace(2, 5, nx)
    yb = np.linspace(y_lower_lim, 0, ny)
    
    bins_step = 0.25
    x_bins = np.arange(2,4.75,bins_step)
    y_vals = list()
    for bb in x_bins:
        y_vals.append(list())
        
    for ii, val in enumerate(x_c):
        for jj, x_val in enumerate(x_bins):
        
            if  x_c[ii] >= x_val and x_c[ii] < x_val + bins_step:
                y_vals[jj].append(y_c[ii])

    
    import scipy
    import scipy.stats
    
    res = scipy.stats.linregress(np.log10(x), np.array(y))
    
    
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
    
    X = sm.add_constant(np.log10(x), prepend=False)
    Y=np.array(y)
    
    df = pd.DataFrame({'x':np.log10(x), 'y':np.array(y)})
    model =smf.ols(formula='y ~ x',data = df).fit()
    
    
    data_list = list()
    fig = plt.figure()
    for ii, y_list in enumerate(y_vals):
        median = np.median(y_list)
        quantiles = np.quantile(y_list,[0.25,0.75])
        spread = np.quantile(y_list,[0.1,0.9])
        data_list.append([x_bins[ii],median,quantiles[0],quantiles[1],spread[0],spread[1]])
      
    data_x,data_y,quantile_1,quantile_2,spread_1,spread_2 = zip(*data_list)
    
    
    return data_x,data_y, [res, model]







outdir = exp_output_path.Compute_spliceAI_exon_id_scratch 
pdf_plots = PdfPages(outdir+'combined_Compute_spliceAI_exon_id_scratch%d_reads.pdf' % (exon_count_build))



x,y,z = zip(*exon_id_list_spliceAI_scores_dict['Intronic'])




median_lines=dict()
median_lines['Intergenic']= dict()
median_lines['Intronic'] = dict()
median_lines['mRNA']     = dict()
median_lines['lncRNA']   = dict()
median_lines['Antisense']   = dict()


for exon_id_list_name in exon_id_list_spliceAI_scores_dict:
    exon_count_vs_spliceAI_scores_list = exon_id_list_spliceAI_scores_dict[exon_id_list_name]  
    x,y,z = zip(*exon_count_vs_spliceAI_scores_list)
            
    
    
    xlabel='exon trapping counts (log10)'
    ylabel='spliceAI score 3ss score (log10)'
    title='%s scatter exon trapping counts vs SpliceAI score\ncounts above 100' % (exon_id_list_name)
    x_lower_lim = 2
    y_lower_lim = -2.5
    num_bins=25
    data_x,data_y,r_sq = get_bin_median(x, np.array(y), xlabel, ylabel, title, x_lower_lim, y_lower_lim, num_bins, pdf_plots)
    median_lines[exon_id_list_name]['3ss'] = data_x,data_y
    median_lines[exon_id_list_name]['3ss_rsq'] = r_sq
    
    xlabel='exon trapping counts (log10)'
    ylabel='spliceAI score 5ss score (log10)'
    title='%s scatter exon trapping counts vs SpliceAI score\ncounts above 100' % (exon_id_list_name)
    x_lower_lim = 2
    y_lower_lim = -2.5
    num_bins=25
    data_x,data_y,r_sq = get_bin_median(x, np.array(z), xlabel, ylabel, title, x_lower_lim, y_lower_lim, num_bins, pdf_plots)
    
    median_lines[exon_id_list_name]['5ss'] = data_x,data_y
    median_lines[exon_id_list_name]['5ss_rsq'] = r_sq





print('\nStandard Error of the Regression for 3ss (unsure correctly implemented)')
for exon_id_list_name in median_lines:
    r_sq, model = median_lines[exon_id_list_name]['3ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, np.sqrt(model.scale)) )


print('\nStandard Error of the Regression for 3ss (unsure correctly implemented)')
for exon_id_list_name in median_lines:
    r_sq, model = median_lines[exon_id_list_name]['5ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, np.sqrt(model.scale)) )



print('\nrsq for OLS model for 3ss (unsure correctly implemented)')
for exon_id_list_name in median_lines:
    r_sq, model = median_lines[exon_id_list_name]['3ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, np.sqrt(model.rsquared)) )

print('\nrsq for OLS model for 5ss (unsure correctly implemented)')
for exon_id_list_name in median_lines:
    r_sq, model = median_lines[exon_id_list_name]['5ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, np.sqrt(model.rsquared)) )


print('\nrsq for scipy linear fit for 3ss')
for exon_id_list_name in median_lines:
    res, model = median_lines[exon_id_list_name]['3ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, res.rvalue**2) )


print('\nrsq for OLS model for 5ss')
for exon_id_list_name in median_lines:
    res, model = median_lines[exon_id_list_name]['5ss_rsq']
    print( "{:}: {:.2}".format(exon_id_list_name, res.rvalue**2) )








fig,ax = plt.subplots()
for key in median_lines:
    data_x,data_y = median_lines[key]['3ss']
    plt.plot(data_x,10**np.array(data_y),label=key, color=color_dict[key])
#plt.plot([2,5],[(0.2),(0.2)], linestyle='dashed', alpha =0.5, color='orange',label='spliceAI 0.2 threshold')
plt.title('median 5ss score SpliceAI')
ax.legend()
plt.xticks([2,3,4,5],['100','1,000','10,000','100,000',])
plt.ylabel('spliceAI 3\'ss score')
plt.xlabel('exon read count')
plt.ylim(0,1)
pdf_plots.savefig(fig)


fig,ax = plt.subplots()
for key in median_lines:
    data_x,data_y = median_lines[key]['5ss']
    plt.plot(data_x,10**np.array(data_y),label=key, color=color_dict[key])
#plt.plot([2,5],[(0.2),(0.2)], linestyle='dashed', alpha =0.5, color='orange',label='spliceAI 0.2 threshold')
plt.title('median 5ss score SpliceAI')
ax.legend()
plt.xticks([2,3,4,5],['100','1,000','10,000','100,000',])
plt.ylabel('spliceAI 5\'ss score')
plt.xlabel('exon read count')
plt.ylim(0,1)
pdf_plots.savefig(fig)



fig,ax = plt.subplots()
for key in ['mRNA','lncRNA','Intergenic']:
    res, model = median_lines[key]['3ss_rsq']
    r_squared = res.rvalue**2
    data_x,data_y = median_lines[key]['3ss']
    plt.plot(data_x,10**np.array(data_y),label=key+" {:.2}".format(r_squared), color=color_dict[key])
#plt.plot([2,5],[(0.2),(0.2)], linestyle='dashed', alpha =0.5, color='orange',label='spliceAI 0.2 threshold')
         
plt.title('median 5ss score SpliceAI')
ax.legend()
plt.xticks([2,3,4,5],['100','1,000','10,000','100,000',])
plt.ylabel('spliceAI 3\'ss score')
plt.xlabel('exon read count')
plt.ylim(0,1)
pdf_plots.savefig(fig)






fig,ax = plt.subplots()
for key in ['mRNA','lncRNA','Intergenic']:
    res, model = median_lines[key]['5ss_rsq']
    r_squared = res.rvalue**2
    data_x,data_y = median_lines[key]['5ss']
    plt.plot(data_x,10**np.array(data_y),label=key+" {:.2}".format(r_squared), color=color_dict[key])
#plt.plot([2,5],[(0.2),(0.2)], linestyle='dashed', alpha =0.5, color='orange',label='spliceAI 0.2 threshold')
plt.title('median 5ss score SpliceAI')
ax.legend()
plt.xticks([2,3,4,5],['100','1,000','10,000','100,000',])
plt.ylabel('spliceAI 5\'ss score')
plt.xlabel('exon read count')
plt.ylim(0,1)
pdf_plots.savefig(fig)



pdf_plots.close()


