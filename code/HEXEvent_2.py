
"""

Generates HEVEVENT figure panel. Creates supplemental data file. Saves data to pickle.

"""




exon_id_cat_dict = region_exon_id_sets_dict



outdir = exp_output_path.HEXEvent 
pdf_plots = PdfPages(outdir+'HEXEvent_%d_reads.pdf' % (exon_count_build))



def above_threshold_exon_id_set(threshold,exon_id_list, aggregate_exon_dict):
    exon_id_list = set(exon_id_list).intersection( aggregate_exon_dict.keys())
    new_exon_id_list = list()
    for exon_id in exon_id_list:
        exon = aggregate_exon_dict[exon_id]
        if exon['count'] >= threshold:
            new_exon_id_list.append(exon_id)
    return new_exon_id_list
    






threshold=100

#gc_bins=[25,40,50,60,75]
gc_bins=[100]



def make_exon_id_list_data(exon_id_list, primary_3ss_exon_id_set, aggregate_exon_dict, gc_bins, genome_fasta):
    
        
    exon_id_list = set(exon_id_list).intersection(primary_3ss_exon_id_set)
    
    exon_id_list = above_threshold_exon_id_set(threshold,exon_id_list, aggregate_exon_dict)
    
    #this isn't really doing its intended job since the gc_bin includes all the data.
    gc_result = el.get_exon_id_list_binned_GC(exon_id_list, gc_bins, genome_fasta)
    
    new_exon_list = list()
    for key in gc_result:
        counts = el.get_exon_dict_counts(gc_result[key],aggregate_exon_dict)
        new_exon_list.append(counts)
    
    return new_exon_list
    
    
    
    
    



exon_id_list = HEXEvent_exon_id_list_cassette
HEXEvent_cassette_list = make_exon_id_list_data(exon_id_list, primary_3ss_exon_id_set, aggregate_exon_dict, gc_bins, genome_fasta)


exon_id_list = HEXEvent_exon_id_list_constitutive
HEXEvent_constitutive_list = make_exon_id_list_data(exon_id_list, primary_3ss_exon_id_set, aggregate_exon_dict, gc_bins, genome_fasta)


exon_id_list = HEXEvent_exon_id_list
HEXEvent_list = make_exon_id_list_data(exon_id_list, primary_3ss_exon_id_set, aggregate_exon_dict, gc_bins, genome_fasta)



exon_id_list = exon_id_cat_dict['mRNA']
pc_middle_list = make_exon_id_list_data(exon_id_list, primary_3ss_exon_id_set, aggregate_exon_dict, gc_bins, genome_fasta)

exon_id_list = exon_id_cat_dict['lncRNA']
lncRNA_middle_list = make_exon_id_list_data(exon_id_list, primary_3ss_exon_id_set, aggregate_exon_dict, gc_bins, genome_fasta)

exon_id_list = exon_id_cat_dict['Intergenic']
no_overlap_list = make_exon_id_list_data(exon_id_list, primary_3ss_exon_id_set, aggregate_exon_dict, gc_bins, genome_fasta)

exon_id_list = exon_id_cat_dict['Intronic']
intronic_list = make_exon_id_list_data(exon_id_list, primary_3ss_exon_id_set, aggregate_exon_dict, gc_bins, genome_fasta)




## output supplemental data
with open(exp_output_path.out_supplemental+'3C.txt','w') as f:
    outstring = 'chrom\tstart\tend\tstrand\texon_id\tgenomic_category\tread_counts\n'
    f.write(outstring)
    
    for exon_id_list_name in exon_id_cat_dict:
        exon_id_list = exon_id_cat_dict[exon_id_list_name]
        exon_id_list = el.exon_id_intersection(exon_id_list, primary_3ss_exon_id_set)
        
        exon_id_list = above_threshold_exon_id_set(threshold,exon_id_list, aggregate_exon_dict)
        
        for exon_id in exon_id_list:
            ex=el.exon_id_values(exon_id)
            
            outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id,exon_id_list_name)
            
            val = aggregate_exon_dict[exon_id]['count']
            outstring = '{:}\t{:}\n'.format(outstring, val)
            f.write(outstring)

    



fig,ax = plt.subplots()
flierprops = dict(marker='o', markersize=1)

box4=plt.boxplot(np.log10(pc_middle_list[0]), positions=[0], widths=.6, patch_artist=True, flierprops=flierprops, showfliers=False,whis=[10,90])
for patch in box4['boxes']:
    patch.set_facecolor(color_dict['mRNA'])

box3=plt.boxplot(np.log10(HEXEvent_cassette_list[0]), positions=[1], widths=0.6, patch_artist=True, flierprops=flierprops, showfliers=False,whis=[10,90])
for patch in box3['boxes']:
    patch.set_facecolor(color_dict['Alternative'])
    
box5=plt.boxplot(np.log10(lncRNA_middle_list[0]), positions=[2], widths=0.6, patch_artist=True, flierprops=flierprops, showfliers=False,whis=[10,90])
for patch in box5['boxes']:
    patch.set_facecolor(color_dict['lncRNA'])
    
box6=plt.boxplot(np.log10(no_overlap_list[0]), positions=[3], widths=0.6, patch_artist=True, flierprops=flierprops, showfliers=False,whis=[10,90])
for patch in box6['boxes']:
    patch.set_facecolor(color_dict['Intergenic'])

box7=plt.boxplot(np.log10(intronic_list[0]), positions=[4], widths=0.6, patch_artist=True, flierprops=flierprops, showfliers=False,whis=[10,90])
for patch in box7['boxes']:
    patch.set_facecolor(color_dict['Intronic'])
    
plt.xticks(range((5)),['mRNA','Alternative','lncRNA','Intergenic', 'Intronic'],rotation=35)
plt.xlabel('gc percent range')
plt.ylabel('splicing inclusion')
plt.title('Threshold: %d reads\nconstitutive_threshold fraction: %.3f' % (threshold,constitutive_threshold))
plt.legend([box4["boxes"][0],box3["boxes"][0],box5["boxes"][0],box6["boxes"][0],box7["boxes"][0]], ['mRNA','Alternative','lncRNA','Intergenic','Intronic'], loc='upper right')
plt.xlim(-.5,4.5)
y_t = [100,1000,10000,100000]
plt.yticks(range(2,6),["{:,}".format(a) for a in y_t])
plt.show()
plt.tight_layout()
pdf_plots.savefig(fig)




pdf_plots.close()

