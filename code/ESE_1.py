
PESE_name = 'chasin_ESEseq'


#6mers
if PESE_name == 'chasin_ESEseq':
    ESE_seq = PESE_class(ESE_seq_score,0 , 'ESEseq')    
    PESE=ESE_seq



pdf_stats = PdfPages(exp_output_path.ESE_scan+'ESE_scan_ESE_database_%s_exon_count_thresh_%d_%d_reads.pdf' % (PESE_name,exon_count_threshold,exon_count_build))



pickle_path = exp_output_path.Dfam_pickle + "Dfam_pickle.pickle"
with open(pickle_path, "rb") as input_file:
    exon_ids_overlapping_repeat_list=pickle.load(input_file)


















region_exon_id_sets_list = list()
region_exon_id_sets_list.append([pc_middle_exon_id_list, 'mRNA'])
region_exon_id_sets_list.append([lncRNA_middle_exon_ids, 'lncRNA'])
region_exon_id_sets_list.append([no_overlap_set, 'Intergenic'])
region_exon_id_sets_list.append([antisense_transcript_set, 'Antisense'])
region_exon_id_sets_list.append([intron_interior_set, 'Intronic'])



for pair in region_exon_id_sets_list:
     pair[0] = list(set(pair[0]).difference(exon_ids_overlapping_repeat_list))






ESE_exon_len_vs_count_list=dict()
ESE_exon_ids=dict()
ESE_100_bp=dict()
ESE_inclusion_vs_count_list=dict()
exon_fraction_covered_list = dict()
for pair in region_exon_id_sets_list:
    exon_id_list = pair[0]
    key = pair[1]
    ESE_exon_ids[key], ESE_100_bp[key], ESE_exon_len_vs_count_list[key], ESE_inclusion_vs_count_list[key], exon_fraction_covered_list[key] =  get_exon_id_ESE_scores(exon_id_list, PESE, genome_fasta,aggregate_exon_dict, exon_count_threshold)








def generate_250_bp_offset_exon_ids(exon_id_list, offset):
    offset_exon_id_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        new_ex = "%s:%d-%d:%s" % (ex.chrom,ex.start+offset,ex.end+offset,ex.strand)
        offset_exon_id_list.append(new_ex)
    return offset_exon_id_list







offset_250_exon_ids_dict = dict()
for key in ESE_exon_ids:
    offset_250_exon_ids_dict[key] = generate_250_bp_offset_exon_ids(ESE_exon_ids[key],250)
    


median_ESE_offset_dict = dict()
offset_exon_len_vs_count_dict = dict()
offset_fraction_covered_dict=dict()
for key in offset_250_exon_ids_dict:
    exon_id_list = offset_250_exon_ids_dict[key]
    median_ESE_offset_dict[key], offset_exon_len_vs_count_dict[key], offset_fraction_covered_dict[key]  = get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta)




    
ESE_plot_dict=dict()
for key in ESE_inclusion_vs_count_list:
    inclusion_vs_count = ESE_inclusion_vs_count_list[key]
    
    
    x,y = zip(*inclusion_vs_count)
    xlabel=''
    ylabel=''
    title = key
    x_lower_lim=0
    y_lower_lim=0
    num_bins=31
    ESE_plot_dict[key] =  plot_2d_hist_ESE(x,y,xlabel,ylabel,title,x_lower_lim,y_lower_lim,num_bins,pdf_stats,PESE_name)




with open(exp_output_path.ESE_scan + '4C.txt', 'w') as f:
    key_list = list()
    list_of_list = list()
    for key in ESE_plot_dict:
        key_list.append(key)
        plot_pair = ESE_plot_dict[key]
        
        x_list = 10**np.array(plot_pair[0])
        
        if len(list_of_list) == 0:
            list_of_list.append(x_list)
        
        list_of_list.append(plot_pair[1])
        
    
    
    zip_list = zip(list_of_list)
    
    f.write('\t'.join( ['exon_read_counts'] + key_list) + '\n')
    matrix = np.array(list(zip_list))
    for row in matrix.T:
        row_str = "\t".join(str(element) for element in row[0])
        f.write(row_str + "\n")





    

    
fig,ax = plt.subplots()
plt.title('median ESE count vs exon expression\n%s'%(PESE_name))
plt.ylabel('ESE count')
plt.xlabel('exon expression counts')
for plot_pair  in ESE_plot_dict.items():
    plt.plot(plot_pair[1][0],plot_pair[1][1], color = color_dict[plot_pair[0]], label=plot_pair[0])

ax.legend()

plt.xticks(range(0,6),[insert_commas(10**x)for x in range(0,6)])
plt.xlim(2,5)
plt.ylim(12.5,37.5)
plt.tight_layout()
pdf_stats.savefig(fig)

    












bar_25_75_list=list()
for ii, key in enumerate(ESE_exon_len_vs_count_list):
    len_vs_count = ESE_exon_len_vs_count_list[key]
    
    dummy,m = zip(*len_vs_count)    
    bar_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75),key])



bar_offset_25_75_list=list()
for ii, key in enumerate(ESE_exon_len_vs_count_list):
    len_vs_count = offset_exon_len_vs_count_dict[key]
    
    dummy,m = zip(*len_vs_count)
    bar_offset_25_75_list.append([np.percentile(m,25),np.percentile(m,50),np.percentile(m,75),key])




with open(exp_output_path.ESE_scan + '4E.txt', 'w') as f:
    f.write('\t'.join(['dataset','exon_25_percentile','exon_50_percentile','exon_75_percentile','offset_25_percentile','offset_50_percentile','offset_75_percentile','\n']))
    for ii, key in enumerate(ESE_exon_len_vs_count_list):
        out_list = [key] + bar_25_75_list[ii][:-1] + bar_offset_25_75_list[ii][:-1]
        f.write('\t'.join(str(element) for element in out_list)+'\n')





m1,m2,m3, k_exon = zip(*bar_25_75_list)
n1,n2,n3, key_offset = zip(*bar_offset_25_75_list)

m1,m2,m3 = np.array(m1),np.array(m2),np.array(m3)
n1,n2,n3 = np.array(n1),np.array(n2),np.array(n3)


x = np.arange(len(bar_offset_25_75_list))  # the label locations
width = 0.35  # the width of the bars



fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, m2, width, label='exon',yerr = [m2-m1,m3-m2])
rects2 = ax.bar(x + width/2, n2, width, label='offset',yerr = [n2-n1,n3-n2])
plt.xticks(range(len(k_exon)),k_exon)
plt.ylabel('PESE counts per 100 bp\n%s octomers with threshold %.2f' % (insert_commas(len(PESE.PESE)),PESE.threshold))
plt.title('offset right 250 bases & threshold %s PESE counts\n%s, exon_thresh: %s' % (insert_commas(exon_count_threshold),PESE_name,insert_commas(exon_count_threshold)))
plt.ylim(0,max(m3)+5)
ax.legend()
plt.tight_layout()
pdf_stats.savefig(fig)






m1,m2,m3, k_exon = zip(*bar_25_75_list)
n1,n2,n3, key_offset = zip(*bar_offset_25_75_list)

m1,m2,m3 = np.array(m1),np.array(m2),np.array(m3)
n1,n2,n3 = np.array(n1),np.array(n2),np.array(n3)


x = np.arange(len(bar_offset_25_75_list))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, m2, width, label='exon')
rects2 = ax.bar(x + width/2, n2, width, label='offset')
plt.xticks(range(len(k_exon)),k_exon)
plt.ylabel('PESE counts per 100 bp\n%s octomers with threshold %.2f' % (insert_commas(len(PESE.PESE)),PESE.threshold))
plt.title('offset right 250 bases & threshold %s PESE counts\n%s, exon_thresh: %s' % (insert_commas(exon_count_threshold),PESE_name,insert_commas(exon_count_threshold)))
plt.ylim(0,max(m2)+5)
ax.legend()
plt.tight_layout()
pdf_stats.savefig(fig)










pdf_stats.close()

























