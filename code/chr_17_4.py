





from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
matrix5 = load_matrix5()
matrix3 = load_matrix3()







D = {'annotated':set(chr17_pc_middle_exons_10_200+chr17_lncRNA_middle_exons_10_200), 'ET':set(ET_chr_17_keys), 'SpliceAI':set(splice_AI_exon_id_list), 'MaxEntScan':set(pair_exon_id_list)}


ET_only = D['ET'].difference(D['annotated']).difference(D['SpliceAI']).difference(D['MaxEntScan'])

SAI_only = D['SpliceAI'].difference(D['annotated']).difference(D['ET']).difference(D['MaxEntScan'])

MES_only = D['MaxEntScan'].difference(D['annotated']).difference(D['ET']).difference(D['SpliceAI'])


ALL_exon_intersect = D['MaxEntScan'].intersection(D['annotated']).intersection(D['ET']).intersection(D['SpliceAI'])

ALL_exon_finder_intersect = D['MaxEntScan'].intersection(D['ET']).intersection(D['SpliceAI']).difference(D['annotated'])




pickle_path = exp_output_path.pickle_merged + "ET_SAI_MES_exon_sets_%d.pickle" % (exon_count_build)
    
with open(pickle_path, "wb") as output_file:
    pickle.dump(D, output_file)
    pickle.dump(ET_only, output_file)
    pickle.dump(SAI_only, output_file)
    pickle.dump(MES_only, output_file)
    pickle.dump(ALL_exon_intersect, output_file)
    pickle.dump(ALL_exon_finder_intersect, output_file)










#el.get_splice
def get_genome_5ss_list(exon_id_list, genome_fasta):
    score_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        seq = genome_fasta[ex.chrom][ex.pos_5ss-4:ex.pos_5ss+5]
        seq=str(seq)
        if seq.find('N')>=0:
            continue
        score = maxent.score5(seq, matrix = matrix5)
        score_list.append(score)
    #print(np.median(score_list),np.mean(score_list))
    return score_list

def get_genome_3ss_list(exon_id_list, genome_fasta):
    score_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        seq = genome_fasta[ex.chrom][ex.pos_3ss-21:ex.pos_3ss+2]
        seq=str(seq)
        if seq.find('N')>=0:
            continue
        score = maxent.score3(seq, matrix = matrix3)
        score_list.append(score)
    #print(np.median(score_list),np.mean(score_list))
    return score_list

def get_genome_gc_list(exon_id_list, genome_fasta):
    score_list = list()
    for exon_id in exon_id_list:
        score = el.get_gc_for_region(exon_id,genome_fasta)
        score_list.append(score)
    #print(np.median(score_list),np.mean(score_list))
    return score_list

def get_genome_5ss_gc_ratio_list(exon_id_list, genome_fasta):
    score_list = list()
    for exon_id in exon_id_list:
        score = el.get_gc_ratio_5ss(exon_id,genome_fasta)
        #print(score[0])
        #break
        score_list.append(score[0])
    #print(np.median(score_list),np.mean(score_list))
    return score_list

def get_genome_3ss_gc_ratio_list(exon_id_list, genome_fasta):
    score_list = list()
    for exon_id in exon_id_list:
        score = el.get_gc_ratio_3ss(exon_id,genome_fasta)
        #print(score[0])
        #break
        score_list.append(score[0])
    #print(np.median(score_list),np.mean(score_list))
    return score_list

def get_ET_5ss_gc_ratio_vs_count_list(exon_id_list, genome_fasta,aggregate_exon_dict):
    score_list = list()
    for exon_id in exon_id_list:
        score = el.get_gc_ratio_5ss(exon_id,genome_fasta)
        count = aggregate_exon_dict[exon_id]['count']
        #print(score[0])
        #break
        score_list.append([score[0], count])
    #print(np.median(score_list),np.mean(score_list))
    return score_list


def get_ET_3ss_gc_ratio_vs_count_list(exon_id_list, genome_fasta,aggregate_exon_dict):
    score_list = list()
    for exon_id in exon_id_list:
        score = el.get_gc_ratio_3ss(exon_id,genome_fasta)
        count = aggregate_exon_dict[exon_id]['count']
        #print(score[0])
        #break
        score_list.append([score[0], count])
    #print(np.median(score_list),np.mean(score_list))
    return score_list

ET_5ss_gc_vs_count_list = get_ET_5ss_gc_ratio_vs_count_list(el.exon_id_intersection(D['ET'],D['ET']), genome_fasta,aggregate_exon_dict)

ET_3ss_gc_vs_count_list = get_ET_3ss_gc_ratio_vs_count_list(el.exon_id_intersection(D['ET'],D['ET']), genome_fasta,aggregate_exon_dict)






def get_seq_5ss_list(exon_id_list,genome_fasta):
    list_5ss_seq = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        if ex.strand == '+':
            seq_5ss = str(genome_fasta[ex.chrom][ex.pos_5ss-3:ex.pos_5ss+6])
            #print(len(seq_5ss))
            if len(seq_5ss) == 9:
                
            #    continue
                list_5ss_seq.append(str(seq_5ss).upper())
        else:
            int('function does not handle - strand!')
    return list_5ss_seq





outdir = exp_output_path.ET_spliceAI_venn_2
pdf_stats = PdfPages(outdir+'SAI_ET_MES_stats.pdf')

ET_only_5ss = get_genome_5ss_list(ET_only,genome_fasta)
SAI_only_5ss = get_genome_5ss_list(SAI_only,genome_fasta)
MES_only_5ss = get_genome_5ss_list(MES_only,genome_fasta)
ALL_5ss = get_genome_5ss_list(ALL_exon_intersect,genome_fasta)
ALL_finder_5ss = get_genome_5ss_list(ALL_exon_finder_intersect,genome_fasta)
annotated_5ss = get_genome_5ss_list(D['annotated'],genome_fasta)

fig = plt.figure()

plt.boxplot([ET_only_5ss],positions = [1],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#EBE837', color='#EBE837'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([SAI_only_5ss],positions = [2],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#31BBED', color='#31BBED'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([MES_only_5ss],positions = [3],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#FBAA28', color='#FBAA28'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([annotated_5ss],positions = [4],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#D6D6D6', color='#D6D6D6'),widths=[0.5],medianprops=dict(color='red'))
plt.xticks(range(1,5),['ET only','SAI only','MES only','annotated'],rotation = 30)
plt.ylabel(' 5ss score (maxentscan)')
plt.ylim(0,15)
pdf_stats.savefig(fig)

ET_only_3ss = get_genome_3ss_list(ET_only,genome_fasta)
SAI_only_3ss = get_genome_3ss_list(SAI_only,genome_fasta)
MES_only_3ss = get_genome_3ss_list(MES_only,genome_fasta)
ALL_3ss = get_genome_3ss_list(ALL_exon_intersect,genome_fasta)
ALL_finder_3ss = get_genome_3ss_list(ALL_exon_finder_intersect,genome_fasta)
annotated_3ss = get_genome_3ss_list(D['annotated'],genome_fasta)


fig = plt.figure()

plt.boxplot([ET_only_3ss],positions = [1],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#EBE837', color='#EBE837'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([SAI_only_3ss],positions = [2],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#31BBED', color='#31BBED'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([MES_only_3ss],positions = [3],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#FBAA28', color='#FBAA28'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([annotated_3ss],positions = [4],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#D6D6D6', color='#D6D6D6'),widths=[0.5],medianprops=dict(color='red'))
plt.xticks(range(1,5),['ET only','SAI only','MES only','annotated'],rotation = 30)
plt.ylabel(' 5ss score (maxentscan)')
plt.ylim(0,15)


pdf_stats.savefig(fig)



ET_only_gc = get_genome_gc_list(ET_only,genome_fasta)
SAI_only_gc = get_genome_gc_list(SAI_only,genome_fasta)
MES_only_gc = get_genome_gc_list(MES_only,genome_fasta)
ALL_gc = get_genome_gc_list(ALL_exon_intersect,genome_fasta)
ALL_finder_gc = get_genome_gc_list(ALL_exon_finder_intersect,genome_fasta)
annotated_gc = get_genome_gc_list(D['annotated'],genome_fasta)


fig = plt.figure()

plt.boxplot([ET_only_gc],positions = [1],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#EBE837', color='#EBE837'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([SAI_only_gc],positions = [2],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#31BBED', color='#31BBED'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([MES_only_gc],positions = [3],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#FBAA28', color='#FBAA28'),widths=[0.5],medianprops=dict(color='red'))

plt.boxplot([annotated_gc],positions = [4],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#D6D6D6', color='#D6D6D6'),widths=[0.5],medianprops=dict(color='red'))
plt.xticks(range(1,5),['ET only','SAI only','MES only','annotated'],rotation = 30)
plt.ylabel(' 5ss score (maxentscan)')
plt.ylim(0,15)

plt.xticks(range(1,5),['ET only','SAI only','MES only','annotated'],rotation=30)
plt.ylabel('exon gc ')
plt.ylim(.25,.75)
pdf_stats.savefig(fig)


def get_exon_id_lengths_list(exon_id_list):
    exon_length_list = list()
    for exon_id in exon_id_list:
        exon_length_list.append(el.exon_id_values(exon_id).length)
    return exon_length_list
        

ET_only_length_list = get_exon_id_lengths_list(ET_only)
SAI_only_length_list = get_exon_id_lengths_list(SAI_only)
MES_only_length_list = get_exon_id_lengths_list(MES_only)
ALL_length_list = get_exon_id_lengths_list(ALL_exon_intersect)
ALL_finder_length_list = get_exon_id_lengths_list(ALL_exon_finder_intersect)
annotated_length_list = get_exon_id_lengths_list(D['annotated'])


fig = plt.figure()

plt.boxplot([ET_only_length_list],positions = [1],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#EBE837', color='#EBE837'),widths=[0.5],medianprops=dict(color='red'))
plt.boxplot([SAI_only_length_list],positions = [2],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#31BBED', color='#31BBED'),widths=[0.5],medianprops=dict(color='red'))
plt.boxplot([MES_only_length_list],positions = [3],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#FBAA28', color='#FBAA28'),widths=[0.5],medianprops=dict(color='red'))
plt.boxplot([annotated_length_list],positions = [4],showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='#D6D6D6', color='#D6D6D6'),widths=[0.5],medianprops=dict(color='red'))
plt.xticks(range(1,5),['ET only','SAI only','MES only','annotated'],rotation = 30)
plt.ylabel(' 5ss score (maxentscan)')


plt.xticks(range(1,5),['ET only','SAI only','MES only','annotated'],rotation=30)
plt.text(1, 140, '%d' % np.median(ET_only_length_list), horizontalalignment='center', size=10,rotation=90)
plt.text(2, 140, '%d' % np.median(SAI_only_length_list), horizontalalignment='center', size=10,rotation=90)
plt.text(3, 140, '%d' % np.median(MES_only_length_list), horizontalalignment='center', size=10,rotation=90)

plt.text(4, 140, '%d' % np.median(annotated_length_list), horizontalalignment='center', size=10,rotation=90)
plt.ylabel('exon lengths')

pdf_stats.savefig(fig)





pdf_stats.close()



import intervaltree

chr17_IT = dict()
chr17_IT['annotated']   = intervaltree.IntervalTree()
chr17_IT['ET']          = intervaltree.IntervalTree()
chr17_IT['SpliceAI']    = intervaltree.IntervalTree()
chr17_IT['MaxEntScan']  = intervaltree.IntervalTree()

count_0_len = 0
for key in D:
    for exon_id in D[key]:
        ex=el.exon_id_values(exon_id)
        if ex.end-ex.start <= 0:
            count_0_len += 1
            continue
        chr17_IT[key][ex.start:ex.end]=exon_id
    

a = set(D['ET'] ).difference(D['SpliceAI']).difference(D['MaxEntScan']).difference(D['annotated'])
len(a)

a_without_ET = set().difference(D['SpliceAI']).difference(D['MaxEntScan']).difference(D['annotated'])
len(a)
len(a_without_ET)

a_without_main = dict()
for main_key in D:
    a_without_main[main_key] = set()
    for key in set(D.keys()).difference([main_key]):
        a_without_main[main_key] = a_without_main[main_key].union(D[key])
len(a_without_main['MaxEntScan'])
len(a_without_main['ET'])



pair_mat = np.zeros([4,8])

k2i = dict([(k,i) for i,k in enumerate(D.keys())])
i2k = dict([(i,k) for i,k in enumerate(D.keys())])


for main_key in D.keys():
    print('main_key: %s' % (main_key))
    all_exons_found = list()    
    for key in set(D.keys()).difference([main_key]):
        count_found = 0
        exons_found = list()
        for exon_id in D[main_key]:
            ex=el.exon_id_values(exon_id)
            found=chr17_IT[key][ex.start:ex.end]
            for interval in found:
                exons = interval[2]
                exons_found.append(exon_id)

                
        b = set(exons_found).difference(a_without_main[main_key])
        print('overlap', key, insert_commas(len(b)))
        print('intersection', insert_commas(len(set(D[main_key]).intersection(D[key]))))
        all_exons_found += b
        all_exons_found = list(set(all_exons_found))
        
        pair_mat[k2i[main_key]][2*k2i[key]]= len(b)
        pair_mat[k2i[main_key]][2*k2i[key]+1]= len(set(D[main_key]).intersection(D[key]))
    
    main_exact_others = set(D[main_key] ).intersection(a_without_main[main_key])
    print('')
    print('All %s: %s' % (main_key, insert_commas(len(D[main_key]))))
    
    print('overlap the others, but not exact', insert_commas(len(set(all_exons_found))))
    
    print('overlap none, and not exact', insert_commas(len( set(D[main_key] ).difference(main_exact_others ).difference( all_exons_found ) ) ) )
    
    print('exact with others', insert_commas(len( main_exact_others) ) )
    print('\n')

print('cols',[(key,k2i[key]*2) for key in k2i])
print(pair_mat)
print('rows', [(key,k2i[key]) for key in k2i])

#len(set(exons_found).difference(a))
len(set(D['ET']).intersection(D['MaxEntScan']))
len(set(D['ET']).intersection(D['SpliceAI']))
len(set(D['ET']).intersection(D['annotated']))

for key in D:
    print('%s has an exon every %s bases' % (key,insert_commas(int(len_chr_17/len(D[key])))))
    
len_chr_17/len(D['ET'])
len_chr_17/len(D['MaxEntScan'])
len_chr_17/len(D['SpliceAI'])
len_chr_17/len(D['annotated'])
    
    
    

    
    
    
    


##### shares one splice site


close_ss_thresh = 6
D_exon_overlaps_dict = dict()
#for main_key in D:
for main_key in D.keys():
    print('main_key: %s' % (main_key))
    
    D_exon_overlaps_dict[main_key]=dict()
    
    all_exons_found = list()    
    for key in set(D.keys()).difference([main_key]):
        D_exon_overlaps_dict[main_key][key] = dict()
        count_found = 0
        exons_found = list()
        exons_close_overlap_5ss = list()
        exons_alternate_5ss = list()
        exons_close_overlap_3ss = list()
        exons_alternate_3ss = list()
        exons_both_alternate_ss = list()
        exons_no_overlap = list()
        
        key_exons_found = list()
        key_exons_close_overlap_5ss = list()
        key_exons_alternate_5ss = list()
        key_exons_close_overlap_3ss = list()
        key_exons_alternate_3ss = list()
        key_exons_both_alternate_ss = list()
        key_exons_no_overlap = list()
        
        for exon_id in D[main_key]:
            ex=el.exon_id_values(exon_id)
            found=chr17_IT[key][ex.start:ex.end]
            
            if len(found)==0:
                    exons_no_overlap.append(exon_id)
                    #key_exons_no_overlap.append([])
            
            for interval in found:
                exons = interval[2]
                ex_found = el.exon_id_values(exons)
                
                
                if exon_id == exons:
                    exons_found.append(exon_id)
                    key_exons_found.append(exon_id)
                
                elif ex.pos_5ss == ex_found.pos_5ss and abs(ex.length-ex_found.length) != 0 and abs(ex.length-ex_found.length) <= close_ss_thresh:
                    exons_close_overlap_3ss.append(exon_id)
                    exons_close_overlap_3ss.append(exons)
                    
                elif ex.pos_5ss == ex_found.pos_5ss and abs(ex.length-ex_found.length) != 0 and abs(ex.length-ex_found.length) > close_ss_thresh:
                    exons_alternate_3ss.append(exon_id)
                    key_exons_alternate_3ss.append(exons)
                    
                elif ex.pos_3ss == ex_found.pos_3ss and abs(ex.length-ex_found.length) != 0 and abs(ex.length-ex_found.length) <= close_ss_thresh:
                    exons_close_overlap_5ss.append(exon_id)
                    key_exons_close_overlap_5ss.append(exons)
                elif ex.pos_3ss == ex_found.pos_3ss and abs(ex.length-ex_found.length) != 0 and abs(ex.length-ex_found.length) > close_ss_thresh:
                    exons_alternate_5ss.append(exon_id) 
                    key_exons_alternate_5ss.append(exons) 
                    
                elif ex.pos_3ss != ex_found.pos_3ss and ex.pos_5ss != ex_found.pos_5ss :
                    exons_both_alternate_ss.append(exon_id)
                    key_exons_both_alternate_ss.append(exons)
                
                
                    
                    
                if main_key == 'ET' and key == 'annotated':
                    1
                    #if exon_id != exons:
                    #    print(exon_id,exons,'\n')
        
        D_exon_overlaps_dict[main_key][key]['exons_found_exactly'] = [exons_found,key_exons_found]
        D_exon_overlaps_dict[main_key][key]['exons_close_overlap_5ss'] =         [exons_close_overlap_5ss,key_exons_close_overlap_5ss]
        D_exon_overlaps_dict[main_key][key]['exons_alternate_5ss'] =         [exons_alternate_5ss,key_exons_alternate_5ss]
        D_exon_overlaps_dict[main_key][key]['exons_close_overlap_3ss'] =         [exons_close_overlap_3ss,key_exons_close_overlap_3ss]
        D_exon_overlaps_dict[main_key][key]['exons_alternate_3ss'] =         [exons_alternate_3ss,key_exons_alternate_3ss]
        D_exon_overlaps_dict[main_key][key]['exons_both_alternate_ss'] =         [exons_both_alternate_ss,key_exons_both_alternate_ss]
        D_exon_overlaps_dict[main_key][key]['exons_no_overlap'] =         [exons_no_overlap,key_exons_no_overlap]
         
        
        b = set(exons_found).difference(a_without_main[main_key])
        print('overlap', key, insert_commas(len(b)))
        print('intersection', insert_commas(len(set(D[main_key]).intersection(D[key]))))
        all_exons_found += b
        all_exons_found = list(set(all_exons_found))
    
    
    main_exact_others = set(D[main_key] ).intersection(a_without_main[main_key])
    print('')
    print('All %s: %s' % (main_key, insert_commas(len(D[main_key]))))
    
    print('overlap the others, but not exact', insert_commas(len(set(all_exons_found))))
    
    print('overlap none, and not exact', insert_commas(len( set(D[main_key] ).difference(main_exact_others ).difference( all_exons_found ) ) ) )
    
    print('exact with others', insert_commas(len( main_exact_others) ) )
    print('\n')







from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


outdir = exp_output_path.ET_spliceAI_venn_2

pdf_5ss = PdfPages(outdir+'pdf_5ss.pdf')
pdf_3ss = PdfPages(outdir+'pdf_3ss.pdf')
pdf_gc = PdfPages(outdir+'pdf_gc.pdf')






text_file_handle = open(outdir + 'overlap_types.txt', 'wt')
for ii, main_key in enumerate(D_exon_overlaps_dict):
    

    def section_page(text, pdf_handle):
        
        page = plt.figure()
        page.clf()
        page.text(0.5,0.5,text, transform=page.transFigure, size=24, ha="center")
        pdf_handle.savefig(page)
        plt.close()
    
    section_page('%s' % main_key, pdf_5ss)
    section_page('%s' % main_key, pdf_3ss)
    section_page('%s' % main_key, pdf_gc)
        
    
    mat = np.zeros((3,7))
    print_string = ''
    
    print_string = print_string + ('main %s:\n' % main_key)
    for jj, key in enumerate(D_exon_overlaps_dict[main_key]):
        print_string = print_string + ('\t%s vs %s:\n' % (main_key,key))
        all_exons_list = list()
        
        
        labels_list = list()
        boxplot_data_5ss = list()
        boxplot_data_3ss = list()
        boxplot_data_gc = list()
        for kk, cat in enumerate(D_exon_overlaps_dict[main_key][key]):
            val_main = list(set(D_exon_overlaps_dict[main_key][key][cat][0]))
            val_key = list(set(D_exon_overlaps_dict[main_key][key][cat][1]))
            
            if cat != 'exons_no_overlap' and cat != 'exons_no_overlap':
                all_exons_list += val_key
            mat[jj][kk] = len(val_main)
            
            exon_id_list_main = el.exon_id_intersection(val_main, aggregate_exon_dict.keys())
            main_5ss_scores = el.get_exon_dict_5ss_scores(exon_id_list_main, aggregate_exon_dict)
            main_5ss_gc = list()
                        
            main_3ss_scores = el.get_exon_dict_3ss_scores(exon_id_list_main, aggregate_exon_dict)
            main_gc_content = list()
            for eid in exon_id_list_main:
                main_gc_content.append(el.get_gc_for_region(eid, genome_fasta))
            
            labels_list.append(cat)
            boxplot_data_5ss.append((main_5ss_scores))
            boxplot_data_3ss.append((main_3ss_scores))
            boxplot_data_gc.append((main_gc_content))
            
            print_string = print_string + ('\t%s:\t%s\t%.1f\t%.1f\t%.1f\n' % (cat,insert_commas(len(val_main)), np.mean(main_3ss_scores), np.mean(main_5ss_scores), 100*np.mean(main_gc_content)))
        all_exons_list = list(set(all_exons_list))
        #print(print_string)
        
        print_string = print_string + ('\t%s:\t%s\n\n' % ('all in %s with at least 1 shared ss' % key,insert_commas( len(all_exons_list) )))
        print(print_string)
        print('')
        #print(mat)
        
        fig = plt.figure()
        plt.boxplot(boxplot_data_5ss,showfliers=False,whis=[10,90])
        plt.title('5ss for %s vs %s' % (main_key,key))
        plt.xticks(range(len(labels_list)),labels_list, rotation = 35)
        plt.ylim(0,15)
        plt.tight_layout()
        pdf_5ss.savefig(fig)
        plt.close()
        
        
        fig = plt.figure()
        plt.boxplot(boxplot_data_3ss,showfliers=False,whis=[10,90])
        plt.title('3ss for %s vs %s' % (main_key,key))
        plt.xticks(range(len(labels_list)),labels_list, rotation = 35)
        plt.ylim(0,15)
        plt.tight_layout()
        pdf_3ss.savefig(fig)
        plt.close()
        
        fig = plt.figure()
        plt.boxplot(boxplot_data_gc,showfliers=False,whis=[10,90])
        plt.title('gc for %s vs %s' % (main_key,key))
        plt.xticks(range(len(labels_list)),labels_list, rotation = 35)
        plt.ylim(0.2,0.8)
        plt.tight_layout()
        pdf_gc.savefig(fig)
        plt.close()
        

    text_file_handle.write(print_string)



pdf_5ss.close()
pdf_3ss.close()
pdf_gc.close()

text_file_handle.close()












exon_id_list = splice_AI_exon_id_list


def get_exon_ids_overlapping_region(exon_id_list, exon_IT):
    overlapping_exon_ids = list()
    antisense_overlapping_exon_ids = list()
    anti_strand = {'+':'-','-':'+'}
    total_exon_IT_bases = list()
    exon_id=list(exon_id_list)[0]
    ex = el.exon_id_values(exon_id)
    for interval in exon_IT[ex.chrom][ex.strand]:
        total_exon_IT_bases += list(range(interval[0],interval[1]))
    
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        
        intervals = exon_IT[ex.chrom][ex.strand][ex.start:ex.end]
        for interval in intervals:
            overlapping_exon_ids.append(exon_id)
        
        intervals = exon_IT[ex.chrom][anti_strand[ex.strand]][ex.start:ex.end]
        for interval in intervals:
            antisense_overlapping_exon_ids.append(exon_id)
            
    overlapping_exon_ids=list(set(overlapping_exon_ids))
    antisense_overlapping_exon_ids=list(set(antisense_overlapping_exon_ids))
    
    print('total # bases exon_IT = %s' % (insert_commas(len(set(total_exon_IT_bases)))))
    
    return overlapping_exon_ids, antisense_overlapping_exon_ids



'''
    #supplemental 5C
    with open(exp_output_path.out_supplemental+'5C.txt','w') as f:
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
'''


def get_exon_overlap_categories(D, exon_id_list, gencode_pc_exon_IT, gencode_lncRNA_exon_IT, gencode_transcript_IT):
    
    #supplemental 5C
    f_5C = open(exp_output_path.out_supplemental+'5C.txt','w') 
    outstring = 'genomeic_category\tmRNA\tIntronic\tIntergenic\tAntisense\n'
    f_5C.write(outstring)
    
    
    
    outpath = exp_output_path.ET_spliceAI_venn_2
    pdf_ratio = PdfPages(outpath+'pdf_ratios.pdf')
    
    
    D_regions = dict()
    for key in D:
        D_regions[key]={'intronic':list(),'antisense':list(),'intergenic':list(),'exon_overlapping':list()}
   
    
    region_ratio_list = list()
        
    for key in D:
        exon_id_list = D[key]
        
        print('gencode_pc_exon_IT')
        mRNA_exon_overlapping_exon_ids, dummy  = get_exon_ids_overlapping_region(exon_id_list, gencode_pc_exon_IT)
        
        print('gencode_lncRNA_exon_IT')
        lncRNA_exon_overlapping_exon_ids, dummy  = get_exon_ids_overlapping_region(exon_id_list, gencode_lncRNA_exon_IT)
        
        print('gencode_transcript_IT')
        transcript_overlapping_exon_ids, anti_sense_transcript_overlapping_exon_ids  = get_exon_ids_overlapping_region(exon_id_list, gencode_transcript_IT)
        
        exon_overlapping = list(set(mRNA_exon_overlapping_exon_ids+lncRNA_exon_overlapping_exon_ids))
        D_regions[key]['exon_overlapping']=exon_overlapping
        
        intronic = set(transcript_overlapping_exon_ids).difference(exon_overlapping)
        D_regions[key]['intronic']=intronic
        
        #intron_overlapping
        antisense_only = el.exon_id_difference(anti_sense_transcript_overlapping_exon_ids, transcript_overlapping_exon_ids)
        D_regions[key]['antisense']=antisense_only
        
        no_overlap = set(exon_id_list).difference(exon_overlapping).difference(intronic).difference(anti_sense_transcript_overlapping_exon_ids).difference(transcript_overlapping_exon_ids)
        D_regions[key]['intergenic']=no_overlap
        
        denom = len(exon_id_list)
        exon_num = len(exon_overlapping)
        intron_num = len(intronic)
        antisense_only_num = len(antisense_only)
        no_overlap_num = len(no_overlap)
        transcript_num = len(transcript_overlapping_exon_ids)
        
        
        a=exon_num/denom
        b=intron_num/denom
        c=no_overlap_num/denom
        d=transcript_num/denom
        e=antisense_only_num/denom
        
        region_ratio_list.append([key,a,b,c,d,e])
        
        #write figure 5C line
        outstring = '{:}\t{:}\t{:}\t{:}\t{:}\n'.format(key,a,b,c,e)
        f_5C.write(outstring)
    
    #close figure 5C
    f_5C.close()
        
    key_vals=list()
    a_vals=np.zeros(len(region_ratio_list))
    b_vals=np.zeros(len(region_ratio_list))
    c_vals=np.zeros(len(region_ratio_list))
    d_vals=np.zeros(len(region_ratio_list))
    e_vals=np.zeros(len(region_ratio_list))
    
    


    
    for ii, region in enumerate(region_ratio_list):
        key = region[0]
        vals = np.array(region[1:6])
        key_vals.append(key)
        a_vals[ii]=vals[0] #exon
        b_vals[ii]=vals[1] #intron
        c_vals[ii]=vals[2] #intergenic
        d_vals[ii]=vals[3] #sense transcript (lncRNA + mRNA)
        e_vals[ii]=vals[4] #antisense
    
    
    
    
    fig = plt.figure()
    plt.title('ratio unannotated intron to antisense exons')
    plt.bar(range(len(a_vals)),b_vals/e_vals,color='red',width=.3)
    plt.xticks(range(len(a_vals)), key_vals)
    plt.plot([4,0],[1,1])
    outpath = exp_output_path.ET_spliceAI_venn_2
    name='fraction_ratio_intron_antisense_ET_SAI_MES.eps'
    plt.savefig(outpath+name)
    pdf_ratio.savefig(fig)
    
        
    #supplemental S5F
    with open(exp_output_path.out_supplemental+'S5F.txt','w') as f:
        outstring = 'exon_finder\tratio\n'
        f.write(outstring)
        for ii, val in enumerate(key_vals):
            if ii == 0:
                continue
            outstring = '{:}\t{:}\n'.format(val, b_vals/e_vals[ii])
            f.write(outstring)
    
    
    
    
    
    pdf_ratio.close()
    return D_regions
    
    
    
D_regions = get_exon_overlap_categories(D, exon_id_list, gencode_pc_exon_IT, gencode_lncRNA_exon_IT, gencode_transcript_IT)




print("D_regions['ET'] keys:", D_regions['ET'].keys())

D_regions['ET']['intronic']
















