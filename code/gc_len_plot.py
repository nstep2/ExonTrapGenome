






outdir = exp_output_path.gc_len_subplots 
pdf_plots = PdfPages(outdir+'gc_length_subplots.pdf')







exon_id_list=pc_middle_exon_id_list
unique_exon_id_list, duplicate_seq_pc_middle_exon_id_list = el.scan_seq_to_exon_id_list_dict(exon_id_list,aggregate_exon_dict, genome_fasta)

pc_middle_length_array = np.zeros(501)
for exon_id in unique_exon_id_list:
    ex = el.exon_id_values(exon_id)
    if ex.length <= 500:
        pc_middle_length_array[ex.length] += 1
    

recovered_exon_plus_fuzzy_list = [x[0] for x in pc_highly_overlapping.values()] + pc_exact_list
recovered_exon_plus_fuzzy_list = set(recovered_exon_plus_fuzzy_list).difference(duplicate_seq_pc_middle_exon_id_list)

ET_pc_middle_length_array = np.zeros(501)
for exon_id in el.exon_id_intersection(unique_exon_id_list, recovered_exon_plus_fuzzy_list):
    ex = el.exon_id_values(exon_id)
    if ex.length <= 500:
        ET_pc_middle_length_array[ex.length] += 1









def plot_len_recovery(pc_middle_length_list,ET_pc_middle_length_list, **kwargs):
    
    
    temp_color_dict = dict()
    temp_color_dict[0] = '#00B9F1'
    temp_color_dict[1] = '#DA6FAB'
    temp_color_dict[2] = '#0072BC'
    temp_color_dict[3] = '#F15A22'
    temp_color_dict[4] = '#00A875'
    temp_color_dict[5] = '#ECDE38'
    
    
    
    fig, ax1 = plt.subplots()
    
    plt.title('recovered annotated exons vs all annotated exons by size\nmoving average window = 9')
    
    if 'gc_split' not in kwargs:
        plt.plot(range(4,len(ET_pc_middle_length_list)-4),moving_average(ET_pc_middle_length_list/pc_middle_length_list, 9),color='blue')
    
    percent_list=list()
    if 'gc_split' in kwargs:
        for ii, gc in enumerate(kwargs['gc_split']):
            percent, et_lengths, an_lengths = gc
            plt.plot(range(4,len(et_lengths)-4),moving_average(et_lengths/an_lengths, 9), color=temp_color_dict[ii]) #,color='blue')
            percent_list.append(percent)
    
    
    if 'gc_split' in kwargs:
        plt.legend(percent_list)    
    else:
        plt.legend(['Exon trapping recovered','all known interal mRNA exons'])
    
    
    
    plt.plot([0,0],[0,0],color='red')
    plt.xlabel('length (bp)')
    plt.ylabel('proportion recovered by exon trapping')
    plt.ylim(0,1)
    ax1_twin=ax1.twinx()
    
    if 'gc_split' in kwargs:
        plt.xlim(0,650)
    else:    
        plt.xlim(0,500)
        
    
    frac_an=moving_average(pc_middle_length_list,9)/max(moving_average(pc_middle_length_list,9))
    max_frac_an = max(frac_an/sum(frac_an))
    
    max_ex_count = max(moving_average(pc_middle_length_list,9)) 
    
    top = int(max_ex_count)
    mid = top/2
    
    if 'gc_split' not in kwargs:
    
        ax1_twin.plot(range(4,len(ET_pc_middle_length_list)-4),frac_an/sum(frac_an),color='red')
    
    
    if 'gc_split' in kwargs:
        for gc in kwargs['gc_split']:
            percent, et_lengths, an_lengths = gc
            sum_gc_split_an_lengths = np.zeros(len(an_lengths))
            
        for gc in kwargs['gc_split']:
            percent, et_lengths, an_lengths = gc
            sum_gc_split_an_lengths += an_lengths
            
        gc_frac_an_sum = moving_average(sum_gc_split_an_lengths,9)/(max(moving_average(sum_gc_split_an_lengths,9) ) )
        
        for ii, gc in enumerate(kwargs['gc_split']):
            percent, et_lengths, an_lengths = gc
            gc_frac_an=moving_average(an_lengths,9)/(max(moving_average(an_lengths,9) ) )

            proportion_an_gc = sum(an_lengths)/sum(sum_gc_split_an_lengths)
            
            ax1_twin.plot(range(4,len(an_lengths)-4),gc_frac_an/sum(gc_frac_an) *proportion_an_gc, color=temp_color_dict[ii] )   #,color='red')
    
    
    
    
    plt.ylim(0,max_frac_an)

    plt.yticks([0,.5*max_frac_an, max_frac_an], ['0%%','%d' % ( mid ), '%d' % (top)])
    plt.ylabel('Exon count')
    

    plt.tight_layout()
    
    pdf_plots.savefig(fig)








def plot_len_recovery_subgroups(pc_middle_length_list,ET_pc_middle_length_list, **kwargs):
    
    
    temp_color_dict = dict()
    temp_color_dict[0] = '#00B9F1'
    temp_color_dict[1] = '#DA6FAB'
    temp_color_dict[2] = '#0072BC'
    temp_color_dict[3] = '#F15A22'
    temp_color_dict[4] = '#00A875'
    temp_color_dict[5] = '#ECDE38'
    
    
    
    fig, ax1 = plt.subplots(2,1,figsize=(6,8))
    
    plt.title('recovered annotated exons vs all annotated exons by size\nmoving average window = 9')
    
    if 'gc_split' not in kwargs:
        ax1[0].plot(range(4,len(ET_pc_middle_length_list)-4),moving_average(ET_pc_middle_length_list/pc_middle_length_list, 9),color='blue')
    
    percent_list=list()
    if 'gc_split' in kwargs:
        for ii, gc in enumerate(kwargs['gc_split']):
            percent, et_lengths, an_lengths = gc
            ax1[0].plot(range(4,len(et_lengths)-4),moving_average(et_lengths/an_lengths, 9), color=temp_color_dict[ii]) #,color='blue')
            percent_list.append(percent)
    
    
    
    
    #plt.plot([0,0],[0,0],color='red')
    plt.xlabel('length (bp)')
    plt.ylabel('proportion recovered by exon trapping')
    plt.ylim(0,1)
    #ax1_twin=ax1.twinx()
    
    if 'gc_split' in kwargs:
        plt.xlim(0,650)
    else:    
        plt.xlim(0,500)
        
    
    frac_an=moving_average(pc_middle_length_list,9)/max(moving_average(pc_middle_length_list,9))
    max_frac_an = max(frac_an/sum(frac_an))*0.3984962406015038  #*.5
    
    max_ex_count = max(moving_average(pc_middle_length_list,9)) 
    
    top = int(max_ex_count*0.3984962406015038)   #*0.37593984962406013
    mid = top/2
    
    if 'gc_split' not in kwargs:
        ax1_twin.plot(range(4,len(ET_pc_middle_length_list)-4),frac_an/sum(frac_an),color='red')
    
    
    if 'gc_split' in kwargs:
        for gc in kwargs['gc_split']:
            percent, et_lengths, an_lengths = gc
            sum_gc_split_an_lengths = np.zeros(len(an_lengths))
            
        for gc in kwargs['gc_split']:
            percent, et_lengths, an_lengths = gc
            sum_gc_split_an_lengths += an_lengths
            
        gc_frac_an_sum = moving_average(sum_gc_split_an_lengths,9)/(max(moving_average(sum_gc_split_an_lengths,9) ) )
        
        for ii, gc in enumerate(kwargs['gc_split']):
            percent, et_lengths, an_lengths = gc
            gc_frac_an=moving_average(an_lengths,9)/(max(moving_average(an_lengths,9) ) )

            proportion_an_gc = sum(an_lengths)/sum(sum_gc_split_an_lengths)
            #ax1[1].plot(range(4,len(an_lengths)-4),gc_frac_an/sum(gc_frac_an) *proportion_an_gc, linestyle='--'  , color=temp_color_dict[ii] )   #,color='red')
            ax1[1].plot(range(4,len(an_lengths)-4),gc_frac_an/sum(gc_frac_an) *proportion_an_gc, color=temp_color_dict[ii] )   #,color='red')
    
    
    
    if 'gc_split' in kwargs:
        ax1[1].legend(percent_list)    
    else:
        ax1[1].legend(['Exon trapping recovered','all known interal mRNA exons'])
    
    
    plt.ylim(0,max_frac_an)
    
    plt.yticks([0,.5*max_frac_an, max_frac_an], ['0','%d' % ( mid ), '%d' % (top)])
    plt.ylabel('Exon count')
    ax1[0].set_ylim(0,1)
    

    plt.tight_layout()
    
    pdf_plots.savefig(fig)












all_middle_length_array = np.zeros(501)
for exon_id in  primary_3ss_exon_id_set:
    ex = el.exon_id_values(exon_id)
    if ex.length <= 500:
        all_middle_length_array[ex.length] += 1
     
     

plot_len_recovery(pc_middle_length_array,ET_pc_middle_length_array)





gc_bin_list = [30, 40, 50, 60, 70]
gc_an_results = el.get_exon_id_list_binned_GC(el.size_exon_id_list(unique_exon_id_list,50,500), gc_bin_list, genome_fasta)

gc_et_results = el.get_exon_id_list_binned_GC(el.size_exon_id_list(el.exon_id_intersection(unique_exon_id_list, recovered_exon_plus_fuzzy_list),50,500), gc_bin_list, genome_fasta)


gc_extra_list = list()
for ii, key in enumerate(gc_an_results):
    
    
    pc_length_array = np.zeros(501)
    for exon_id in gc_an_results[key]:
        ex = el.exon_id_values(exon_id)
        if ex.length <= 500:
            pc_length_array[ex.length] += 1
    
    et_length_array = np.zeros(501)
    for exon_id in gc_et_results[key]:
        ex = el.exon_id_values(exon_id)
        if ex.length <= 500:
            et_length_array[ex.length] += 1
            
    gc_extra_list.append([key, et_length_array, pc_length_array])



















plot_len_recovery_subgroups(pc_middle_length_array,ET_pc_middle_length_array, gc_split=gc_extra_list)













pdf_plots.close()












