

exon_id_cat_dict = region_exon_id_sets_dict


outdir = exp_output_path.exon_id_conservations 
pdf_plots = PdfPages(outdir+'exon_id_conservations.pdf')






runfile('bigwig_lib.py', wdir='/home/pdf/repositories/exon_def/Process_ET', current_namespace=True)




def insert_commas(value):
    return f'{value:,}'



'''
bigwig_path = exp_output_path.bigwig_path+'phastcons_17way/hg38.phastCons17way.bw'
bw_data_phastcons_17 = bigwig_cons(bigwig_path)
'''

bigwig_path = exp_output_path.bigwig_path+'phastcons_7way/hg38.phastCons7way.bw'
bw_data_phastcons_7 = bigwig_cons(bigwig_path)


bigwig_path = exp_output_path.bigwig_path+'phylop_7way/hg38.phyloP7way.bw'
bw_data_phylop_7 = bigwig_cons(bigwig_path)


'''
bigwig_path = exp_output_path.bigwig_path+'phylop_17way/hg38.phyloP17way.bw'
bw_data_phylop_17 = bigwig_cons(bigwig_path)
'''



selected_cons_metric = bw_data_phylop_7
selected_cons_metric_name = 'bw_data_phylop_7'




f_supp = open(exp_output_path.out_supplemental+'7A.txt','w')

outstring = 'chrom\tstart\tend\tstrand\texon_id\tgenomic_category\tcons_score_array_3\tcons_score_array_5\n'
f_supp.write(outstring)



def get_mean_conservation_array(check_id_list, check_distance, aggregate_exon_dict, bw_data):
    
    
    average_result_array_3ss = np.zeros(2*check_distance+1)
    average_result_count_3ss = 0
    
    average_result_array_5ss = np.zeros(2*check_distance+1)
    average_result_count_5ss = 0
    
    tenths = int(len(check_id_list)/10)
    for ii, exon_id in enumerate(check_id_list):
        
        if ii % tenths == 0:
            print("%d%%"%(ii/tenths*10), ii)
        
        ex = el.exon_id_values(exon_id)
        
        pos_3ss = ex.pos_3ss
        pos_5ss = ex.pos_5ss
        
        
        try:
            
            score_array3 = bw_data.get_base_scores(ex.chrom, pos_3ss-check_distance, pos_3ss+check_distance)
            
            if bw_data.check_if_none(score_array3) == False:
                if ex.strand =='-':
                    score_array3 = np.flip(score_array3)
                average_result_array_3ss = average_result_array_3ss + np.array(score_array3)
                average_result_count_3ss += 1
            else:
                score_array3 = '.'
            
                
            
            score_array5 = bw_data.get_base_scores(ex.chrom, pos_5ss-check_distance, pos_5ss+check_distance)
            
            if bw_data.check_if_none(score_array5) == False:
                if ex.strand =='-':
                    score_array5 = np.flip(score_array5)
                average_result_array_5ss = average_result_array_5ss + np.array(score_array5)
                average_result_count_5ss += 1
            else:
                score_array5 = '.'
               
            
            
            outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id,exon_id_list_name)
                
            val_3 = ','.join(map(str, score_array3[1:]))
            val_5 = ','.join(map(str, score_array5[1:]))
            outstring = '{:}\t{:}\n'.format(outstring, val_3, val_5)
            f_supp.write(outstring)

        
        
        
        except RuntimeError:
            
            pass
        
        
        
    return (average_result_array_3ss/average_result_count_3ss), (average_result_array_5ss/average_result_count_5ss)






def conservation_figure(data, label, text, file_path):
    fig = plt.figure()
    plt.plot(data)
    plt.ylim(.0,.9)
    plt.xticks([1,25,51,76,101],[-50,-25,0,25,50])
    plt.yticks([0,.2,.4,.6,.8],[0,.2,.4,.6,.8])
    plt.xlabel('bp from %s splice site' % (label))
    plt.title('conservation around %s using \n%s' % (label,text))
    plt.ylabel('conservation score')
    plt.tight_layout()
    pdf_plots.savefig(fig)







average_result_array_3ss_dict = dict()
average_result_array_5ss_dict = dict()
for ii, key in enumerate(['mRNA', 'lncRNA', 'Intergenic', 'Intronic']):

    exon_id_list_name = key
    check_id_list = exon_id_cat_dict[key]

    average_result_array_3ss, average_result_array_5ss = get_mean_conservation_array(check_id_list, 50, aggregate_exon_dict, selected_cons_metric)
    average_result_array_3ss_dict[key] = average_result_array_3ss
    average_result_array_5ss_dict[key] = average_result_array_5ss




outdir = exp_output_path.exon_id_conservations 
pdf_plots = PdfPages(outdir+'exon_id_conservations.pdf')






fig = plt.figure()

label='3ss'
text=selected_cons_metric_name
for ii, key in enumerate(average_result_array_3ss_dict):
    plt.plot(average_result_array_3ss_dict[key], color=color_dict[key])


plt.ylim(.0,.9)
plt.yticks([0,.2,.4,.6,.8],[0,.2,.4,.6,.8])
plt.xticks([1,25,51,76,101],[-50,-25,0,25,50])
plt.xlabel('bp from %s splice site' % (label))
plt.title('conservation around %s using \n%s' % (label,text))
plt.ylabel('conservation score')
plt.tight_layout()
pdf_plots.savefig(fig)


fig = plt.figure()


label='5ss'
text=selected_cons_metric_name
for ii, key in enumerate(average_result_array_5ss_dict):
    plt.plot(average_result_array_5ss_dict[key], color=color_dict[key],label=key)



plt.ylim(.0,.9)
plt.yticks([0,.2,.4,.6,.8],[0,.2,.4,.6,.8])
plt.xticks([1,25,51,76,101],[-50,-25,0,25,50])
plt.xlabel('bp from %s splice site' % (label))
plt.title('conservation around %s using \n%s' % (label,text))
plt.ylabel('conservation score')
plt.legend(['mRNA','lncRNA','intronic','intergenic'])
plt.tight_layout()
pdf_plots.savefig(fig)




exon_id_list_name='mRNA'
check_id_list = exon_id_cat_dict['mRNA']
average_result_array_3ss, average_result_array_5ss = get_mean_conservation_array(check_id_list, 50, aggregate_exon_dict, selected_cons_metric)


average_result_array_3ss_pc_phylop_17, average_result_array_5ss_pc_phylop_17 = (average_result_array_3ss, average_result_array_5ss)





'''
data = average_result_array_5ss
label='5ss'
text='pc_middle_exon bw_data_phylop_17'
#conservation_figure(data, label, text, file_path)

data = average_result_array_3ss
label='3ss'
text='pc_middle_exon bw_data_phylop_17'
#conservation_figure(data, label, text, file_path)
'''











exon_id_list_name='lncRNA'
check_id_list = lncRNA_middle_exon_ids
average_result_array_3ss, average_result_array_5ss = get_mean_conservation_array(check_id_list, 50, aggregate_exon_dict, selected_cons_metric)

average_result_array_3ss_lncRNA_phylop_17, average_result_array_5ss_lncRNA_phylop_17 = (average_result_array_3ss, average_result_array_5ss)

'''
data = average_result_array_5ss
label='5ss'
text='lncRNA_middle_exon bw_data_phylop_17'
#conservation_figure(data, label, text, file_path)

data = average_result_array_3ss
label='3ss'
text='lncRNA_middle_exon bw_data_phylop_17'
#conservation_figure(data, label, text, file_path)

'''




















































'''

try:
    exon_id_cat_dict['Intergenic']
    intron_interior_set_flag = True
except:
    intron_interior_set_flag = False
    
    

if intron_interior_set_flag == True:
    '''
    


exon_id_list_name='Intergenic'
check_id_list = exon_id_cat_dict['Intergenic']
average_result_array_3ss, average_result_array_5ss = get_mean_conservation_array(check_id_list, 50, aggregate_exon_dict, selected_cons_metric)

average_result_array_3ss_intergenic_phylop_17, average_result_array_5ss_intergenic_phylop_17 = (average_result_array_3ss, average_result_array_5ss)




exon_id_list_name=['Intronic']
check_id_list = exon_id_cat_dict['Intronic']
#check_id_list = mRNA_introns #difference didn't matter.
average_result_array_3ss, average_result_array_5ss = get_mean_conservation_array(check_id_list, 50, aggregate_exon_dict, selected_cons_metric)

average_result_array_3ss_intronic_phylop_17, average_result_array_5ss_intronic_phylop_17 = (average_result_array_3ss, average_result_array_5ss)
   



#def conservation_figure(data, label, text, file_path):
fig = plt.figure()
data0 = average_result_array_3ss_intronic_phylop_17
data1 = average_result_array_3ss_intergenic_phylop_17
data2 = average_result_array_3ss_pc_phylop_17
data3 = average_result_array_3ss_lncRNA_phylop_17

label='3ss'
text='phylop_17'
plt.plot(data2, color=color_dict['mRNA'])
plt.plot(data3, color=color_dict['lncRNA'])
plt.plot(data0, color=color_dict['Intronic'])
plt.plot(data1, color=color_dict['Intergenic'])


plt.ylim(.0,.9)
plt.yticks([0,.2,.4,.6,.8],[0,.2,.4,.6,.8])
plt.xticks([1,25,51,76,101],[-50,-25,0,25,50])
plt.xlabel('bp from %s splice site' % (label))
plt.title('conservation around %s using \n%s' % (label,text))
plt.ylabel('conservation score')
plt.legend(['mRNA','lncRNA', 'intronic','intergenic'])
plt.tight_layout()
pdf_plots.savefig(fig)



#def conservation_




#figure(data, label, text, file_path)
fig = plt.figure()
data0 = average_result_array_5ss_intronic_phylop_17
data1 = average_result_array_5ss_intergenic_phylop_17
data2 = average_result_array_5ss_pc_phylop_17
data3 = average_result_array_5ss_lncRNA_phylop_17
label='5ss'
text='phylop_17'
plt.plot(data2, color=color_dict['mRNA'])
plt.plot(data3, color=color_dict['lncRNA'])
plt.plot(data0, color=color_dict['Intronic'])
plt.plot(data1, color=color_dict['Intergenic'])


plt.ylim(.0,.9)
plt.yticks([0,.2,.4,.6,.8],[0,.2,.4,.6,.8])
plt.xticks([1,25,51,76,101],[-50,-25,0,25,50])
plt.xlabel('bp from %s splice site' % (label))
plt.title('conservation around %s using \n%s' % (label,text))
plt.ylabel('conservation score')
plt.legend(['mRNA','lncRNA','intronic','intergenic'])
plt.tight_layout()
pdf_plots.savefig(fig)






f_supp.close()
pdf_plots.close()







