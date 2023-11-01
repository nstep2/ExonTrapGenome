









PESE_name = 'chasin_PESE'
PESE_name = 'chasin_ESEseq'
#PESE_name = 'burge'



#initial PESE scorer
if PESE_name == 'chasin_PESE':
    PESE = PESE_class(chasin_PESE,PESE_threshold,'chasin_PESE')    

#Alternate PESE scorers - 6mers
if PESE_name == 'chasin_ESEseq':
    PESE=ESE_seq

#burge fairbrother
if PESE_name == 'burge':
    PESE=burge_ESE  








pc_exon_lengths = [el.exon_id_values(exon_id).length for exon_id in pc_middle_exon_id_list]
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

    data_x,data_y,quantile_1,quantile_2,spread_1,spread_2 = zip(*data_list)

    
    return data_x,data_y






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



MES_range_3ss = [0,6]
MES_range_5ss = [0,6]

exon_id_list = no_overlap_set





def plot_MES_range_ESE_vs_inclusion(exon_id_list, aggregate_exon_dict, PESE, PESE_name, exon_id_list_name,  *args, **kwargs):
    
    write_to_pdf_flag = 'pdf_handle' in kwargs

    median_GC_pair_MES_list = list()
    median_ESE_vs_inclusion_list = list()
    
    main_figure_line = [0,0,'']
    main_figure_line = dict()
    
    for ii, pair in enumerate([[0,4],[4,6],[6,8],[8,10],[10,15]]):
        MES_range_3ss = pair
        MES_range_5ss = pair

        
        exon_count_vs_PESE_list,  exon_avg_count_vs_PESE_list, exon_countper_len_vs_PESE_list, gc_list = get_inclusion_vs_ESE_count_for_MES_range(exon_id_list, aggregate_exon_dict, MES_range_3ss, MES_range_5ss, PESE)

        x,y = zip(*exon_countper_len_vs_PESE_list)

        
        num_bins=50
        x_c = np.array(x)
        y_c = (np.log10(y) ) #*np.array(z))
        nx, ny = (num_bins, num_bins) #31
        xb = np.linspace(0, 50, nx)
        yb = np.linspace(min(y_c), 5, ny)

        heatmap, xedges, yedges = np.histogram2d(x_c, y_c, bins=(xb,yb))
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        line_x = np.percentile(x,5)
        
        xlabel = 'ESE count'
        ylabel = 'exon read count'
        title = '%s - ESE:%s total_data: %s\nx_median %.1f (25-75%%: %.1f-%.1f), y_median %.1f\n 3ss range [%.1f-%.1f] 5ss range [%.1f-%.1f]' % (exon_id_list_name, PESE_name, insert_commas(len(x)), np.median(x),np.percentile(x,25), np.percentile(x,75), np.median(y), MES_range_3ss[0],MES_range_3ss[1], MES_range_5ss[0],MES_range_5ss[1])

        
        
        x_data = x
        y_data = y
        
        x_bin_lims = [0,50,52,2,2] # lower, mid, upper, step lower, step upper
        y_bin_lims = [0, max(x_data)]

        
        num_bins=50
        y_c = np.array(x)
        x_c = (np.log10(y) ) #*np.array(z))
        nx, ny = (num_bins, num_bins) #31
        yb = np.linspace(0, 50, nx)
        xb = np.linspace(min(x_c), 5, ny)
        

        heatmap, xedges, yedges = np.histogram2d(x_c, y_c, bins=(xb,yb))
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        line_y = np.percentile(y_c,5)

                
        ylabel = 'ESE count'
        xlabel = 'exon read count'
        title = '%s - ESE:%s total_data: %s\nx_median %.1f (25-75%%: %.1f-%.1f), y_median %.1f\n 3ss range [%.1f-%.1f] 5ss range [%.1f-%.1f]' % (exon_id_list_name, PESE_name, insert_commas(len(y)), np.median(y),np.percentile(y,25), np.percentile(y,75), np.median(x), MES_range_3ss[0],MES_range_3ss[1], MES_range_5ss[0],MES_range_5ss[1])
        
        x_data = y
        y_data = x
        
        y_bin_lims = [0,50] # lower, mid, upper, step lower, step upper
        x_bin_lims = [2, 4, 5.1, .2, .2]
        
        if write_to_pdf_flag == True:
            data_x,data_y = plot_ESE_linear_y(x_data,y_data,title,xlabel, ylabel, x_bin_lims, y_bin_lims, pdf_handle =  kwargs['pdf_handle']) 
            
        else:
            data_x,data_y = plot_ESE_linear_y(x_data,y_data,title,xlabel, ylabel, x_bin_lims, y_bin_lims) 

        if exon_id_list_name not in main_figure_line:
            main_figure_line[exon_id_list_name] = dict()
        
        splice_site_bin = '{:d}-{:d}'.format(pair[0],pair[1])
        main_figure_line[exon_id_list_name][splice_site_bin] = [data_x,data_y]
        
            
        
            
        
        median_ESE_vs_inclusion_list.append([ii, pair,np.median(x),np.median(y)])
        median_GC_pair_MES_list.append([ii, pair,np.median(gc_list), np.mean(gc_list), np.percentile(gc_list,25), np.percentile(gc_list,75)])
    
    index, pair_list, med_x, med_y = zip(*median_ESE_vs_inclusion_list)
    
    
    return median_ESE_vs_inclusion_list, median_GC_pair_MES_list, main_figure_line








exon_id_list=pc_middle_exon_id_list
exon_id_list_name='mRNA'
median_ESE_vs_inclusion_list, median_GC_pair_MES_list, main_figure_line_mRNA = plot_MES_range_ESE_vs_inclusion(exon_id_list, aggregate_exon_dict, PESE, PESE.name, exon_id_list_name)        
median_ESE_vs_inclusion_list, median_GC_pair_MES_list, main_figure_line      = plot_MES_range_ESE_vs_inclusion(exon_id_list, aggregate_exon_dict, PESE, PESE.name, exon_id_list_name)        


exon_id_list=no_overlap_set
exon_id_list_name='Intergenic'
median_ESE_vs_inclusion_list, median_GC_pair_MES_list, main_figure_line_intergenic = plot_MES_range_ESE_vs_inclusion(exon_id_list, aggregate_exon_dict, PESE, PESE.name, exon_id_list_name)        


exon_id_list=antisense_transcript_set
exon_id_list_name='Antisense'
median_ESE_vs_inclusion_list, median_GC_pair_MES_list, main_figure_line_antisense = plot_MES_range_ESE_vs_inclusion(exon_id_list, aggregate_exon_dict, PESE, PESE.name, exon_id_list_name)        






exon_id_list_list = list()
exon_id_list_list.append(['mRNA', pc_middle_exon_id_list])
exon_id_list_list.append(['lncRNA', lncRNA_middle_exon_ids])
exon_id_list_list.append(['Intergenic', no_overlap_set])
exon_id_list_list.append(['Intronic', intron_interior_set])
exon_id_list_list.append(['Antisense', antisense_transcript_set])





pickle_path = exp_output_path.Dfam_pickle + "Dfam_pickle.pickle"
with open(pickle_path, "rb") as input_file:
    exon_ids_overlapping_repeat_list=pickle.load(input_file)


ESE_count_vs_inclusion_group = list()
ESE_gc_group = list()
pdf_ESE = PdfPages(exp_output_path.ESE_scan+'ESE_scan_MES_scores_%s_%d_reads.pdf' % (PESE.name,exon_count_build))
for ii, lib_pair in enumerate(exon_id_list_list):
    exon_id_list = set(lib_pair[1]).difference(exon_ids_overlapping_repeat_list)
    exon_id_list_name = lib_pair[0]

    median_ESE_vs_inclusion_list, median_GC_pair_MES_list, main_figure_line = plot_MES_range_ESE_vs_inclusion(exon_id_list, aggregate_exon_dict, PESE, PESE.name, exon_id_list_name, pdf_handle = pdf_ESE)    
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
    



index, names,ESE_gc = zip(*ESE_gc_group)
dummy,MES_range, dummy, dummy, dummy, dummy = zip(*ESE_gc[0])
fig = plt.figure(figsize=(8,5))
for ii, val in enumerate(ESE_gc):
    val_list = list()
    yerror = list()
    for jj, dummy in enumerate(ESE_gc):
        med = ESE_gc[jj][ii][2]
        val_list.append(med)
        yerror.append(np.array([med-ESE_gc[jj][ii][4], ESE_gc[jj][ii][5] -med ]))
    
    yerror=np.array(yerror).T
    plt.bar(np.arange(0,len(index))+ii*.15, val_list, yerr=yerror,width = .15)
plt.bar([len(index)],[0])
plt.legend(MES_range)
plt.xticks(index, names, rotation=30)
plt.title('GC percent in different maxentscan ranges (legend)\n%s' % (PESE.name))
plt.ylabel('Median GC fraction\nyerror =25%%-75%% range')
plt.tight_layout()
pdf_ESE.savefig(fig)




for region_key in ['Intergenic']:
    if region_key == 'mRNA':
        x,y=main_figure_line_mRNA[region_key]['8-10']
    if region_key == 'Intergenic':
        x,y=main_figure_line_intergenic['Intergenic']['8-10']
    
    fig,ax = plt.subplots()
    plt.title('%s - median ESE count vs exon expression\n%s'%(region_key,PESE_name))
    plt.ylabel('ESE count')
    plt.xlabel('exon expression counts')
    
    plt.plot(x,y, color = 'darkorange', label='%s [8,10]'%region_key,linestyle='dashed')

    for plot_pair  in ESE_plot_dict.items():
        plt.plot(plot_pair[1][0],plot_pair[1][1], color = color_dict[plot_pair[0]], label=plot_pair[0])
        
    ax.legend()
    plt.xticks(range(0,6),[insert_commas(10**x)for x in range(0,6)])
    plt.xlim(2,5)
    plt.ylim(17.5,37.5)
    plt.tight_layout()
    pdf_ESE.savefig(fig)




fig,ax = plt.subplots()
plt.title('Intergenic - median ESE count vs exon expression\n%s'%(PESE_name))
plt.ylabel('ESE count')
plt.xlabel('exon expression counts')
    
for key in main_figure_line_intergenic['Intergenic']:
    x,y=main_figure_line_intergenic['Intergenic'][key]
    plt.plot(x,y, label=key)

ax.legend()
plt.xticks(range(0,6),[insert_commas(10**x)for x in range(0,6)])
plt.xlim(2,5)
plt.ylim(17.5,35)
plt.tight_layout()
pdf_ESE.savefig(fig)








pdf_ESE.close()









