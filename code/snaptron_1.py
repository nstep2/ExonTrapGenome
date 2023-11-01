'''

run 
old /home/pdf/repositories/exon_def/Process_ET/intronic_exon_conservation.py

old /home/pdf/repositories/exon_def/Process_ET/Parse_ENSEMBLE_GENCODE_exons_load_figures_2.py

data from other scripts:
lncRNA_database_exons
aaa

'''




'''
#run this once to create an index to speed searches up (if needed)
CREATE INDEX exon_id ON introns (chromosome, start, end, strand)
'''




#exon_sql_path


outdir = exp_output_path.snaptron_exon 
pdf_plots = PdfPages(outdir+'snaptron_exon.pdf')




import sqlite3
import time

#exp_output_path.exon_sql_path = '/mnt/hgfs/Linux/snaptron/exons.sqlite'


#sexp_output_path.exon_sql_path = '/mnt/hgfs/main_ssd/snaptron/srav3h/redownload_230702/exons.sqlite'  #not valid database




conn = sqlite3.connect(exp_output_path.exon_sql_path)
conn.row_factory = sqlite3.Row   #a sqlite list/dict (both col# and col key access column values)
c = conn.cursor()


'''
res = conn.execute("SELECT name FROM sqlite_master WHERE type='table';")
for name in res.fetchall():
    print(name[0])
c.execute("PRAGMA table_info(intron);")
column_names = [column[1] for column in c.fetchall()]
'''

#run this once to create an index to speed searches up (if needed)
#c.execute("CREATE INDEX exon_id ON intron (chrom, start, end, strand)")


selected = c.execute('select * from intron')



col_list=[k[0] for ii, k in enumerate(selected.description) if k[0] != 'samples']
selected_string = ', '.join(col_list)


c.execute("SELECT count(*) FROM sqlite_master WHERE type='index';")




'''
result = c.execute("SELECT name FROM sqlite_master WHERE type='table';")
print(list(result.fetchall()))

result = c.execute("SELECT name FROM sqlite_master;")
print(list(result.fetchall()))

c.execute("SELECT name FROM sqlite_master WHERE type='table';")
print(c.fetchall())
print(c.fetchone())
'''
#a=c.fetchone()









def db_search_exon_id(c, exon_id, window):
    ex = el.exon_id_values(exon_id)
    
    start = ex.start
    end   = ex.end-1
    cmd = "SELECT {:} FROM intron WHERE chrom='{:}' AND start BETWEEN {:} AND {:} AND end BETWEEN {:} AND {:} AND strand='{:}';".format(selected_string, ex.chrom,start-window,start+window,end-window,end+window,ex.strand)
    
    
    r = c.execute(cmd)
    
    return r





def get_snaptron_exon_id_hits(exon_id_list, c, window):

    count_found = 0
    count_checked = 0
    found_r = dict()
    check_len = len(exon_id_list)
    for ii, exon_id in enumerate(exon_id_list):
        #if ii % ( int(check_len/40)) == 0:
        #    print('percent checked = {:.1%}'.format(ii/check_len))
        count_checked += 1
        result = db_search_exon_id(c, exon_id, window)
        found_hit = False
        r = result.fetchall()
        for entry in r:
            
            if len(r) > 0:
                found_hit = True
                found_r[exon_id] = r
            
            if found_hit == True:
                count_found += 1
            
        
    return found_r, count_found, count_checked







def db_search_inside_exon_id(c, exon_id, window):
    ex = el.exon_id_values(exon_id)
    
    start = ex.start
    end   = ex.end
    cmd = "SELECT {:} FROM intron WHERE chrom='{:}' AND start > {:} AND end < {:} AND strand='{:}';".format(selected_string, ex.chrom,start-window,end+window,ex.strand)
    
    #print(cmd)
    
    r = c.execute(cmd)
    
    return r



def get_snaptron_inside_exon_id_hits(exon_id_list, c, window):

    count_found = 0
    count_checked = 0
    found_r = dict()
    check_len = len(exon_id_list)
    for ii, exon_id in enumerate(exon_id_list):
        #if ii % ( int(check_len/40)) == 0:
        #    print('percent checked = {:.1%}'.format(ii/check_len))
        count_checked += 1
        result = db_search_inside_exon_id(c, exon_id, window)
        found_hit = False
        r = result.fetchall()
        for entry in r:
            
            if len(r) > 0:
                found_hit = True
                found_r[exon_id] = r
            
            if found_hit == True:
                count_found += 1
            
        
    return found_r, count_found, count_checked





def exon_id_from_pop(q1):
    return "%s:%d-%d:%s" % (q1['chrom'], q1['start'], q1['end'], q1['strand'])











'''



def exon_id_lengths(exon_id_list, genome_fasta):
    length_list = list()
    score_5ss_list = list()
    score_3ss_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        score_3ss_list.append(ex.score_3ss(genome_fasta)[0])
        score_5ss_list.append(ex.score_5ss(genome_fasta)[0])
        length_list.append(ex.length)
    return length_list, score_3ss_list, score_5ss_list

'''






import matplotlib.pyplot as plt





'''

def plot_found_missing_stats(found_exon_id_list, missing_exon_id_list, exon_id_label):
    
    exon_id_list = missing_exon_id_list
    missing_lengths, missing_3ss_list, missing_5ss_list = exon_id_lengths(exon_id_list, genome_fasta)
    
    exon_id_list = found_exon_id_list
    found_lengths, found_3ss_list, found_5ss_list  = exon_id_lengths(exon_id_list, genome_fasta)
    
    
    
    
    fig, ax = plt.subplots()
    plt.boxplot([found_lengths], positions=[0], labels=['found'], showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='y', color='y', alpha = 0.5))
    
    plt.boxplot([missing_lengths], positions=[1], labels=['missing'], showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='cyan', color='cyan', alpha = 0.5))
    
    plt.ylabel('lengths')
    plt.ylim(50,350)
    plt.title(exon_id_label)
    plt.tight_layout()
    
    
    
    
    fig, ax = plt.subplots()
    plt.boxplot([found_3ss_list], positions=[0], labels=['found'], showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='y', color='y', alpha = 0.5))
    
    plt.boxplot([missing_3ss_list], positions=[1], labels=['missing'], showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='cyan', color='cyan', alpha = 0.5))
    
    plt.ylabel('maxent 3ss')
    plt.ylim(4,12)
    plt.title(exon_id_label)
    plt.tight_layout()
    
    
    
    
    
    fig, ax = plt.subplots()
    plt.boxplot([found_5ss_list], positions=[0], labels=['found'], showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='y', color='y', alpha = 0.5))
    
    plt.boxplot([missing_5ss_list], positions=[1], labels=['missing'], showfliers=False,whis=[10,90], patch_artist=True, boxprops=dict(facecolor='cyan', color='cyan', alpha = 0.5))
    
    plt.ylabel('maxent 5ss')
    plt.ylim(4,12)
    plt.title(exon_id_label)
    plt.tight_layout()


'''




def get_snaptron_gene_stats(transcript_IT, GENCODE_gene_dict,exon_id_list,c):
    pre_mRNA_snaptron_hits_count = 0
    window = 3
    
    found_exon_id_list = list()
    
    found_tid_dict = dict()
    found_gene_dict = dict()
    
    count_transcript_not_found = 0
    
    aaa={gene_id:GENCODE_gene_dict[gene_id] for gene_id in GENCODE_gene_dict if GENCODE_gene_dict[gene_id]['tags']['gene_type']=='protein_coding'}
#len(aaa)

    def make_IT(aaa, transcript_IT):
        IT = dict()
        for chrom in transcript_IT:
            IT[chrom] = {'+':intervaltree.IntervalTree(),'-':intervaltree.IntervalTree()}
        
        for gene_id in aaa:
            start = aaa[gene_id]['start']
            end = aaa[gene_id]['end']
            chrom = aaa[gene_id]['chrom']
            strand = aaa[gene_id]['strand']
            
            new_inteval = intervaltree.Interval(start, end, gene_id)
            IT[chrom][strand].add(new_inteval)
    
        return IT
            
    IT = make_IT(aaa, transcript_IT)
    
    all_gene_id_dict = dict() 
    for ii, exon_id in enumerate(exon_id_list):
        
        ex=el.exon_id_values(exon_id)
        found_intervals = IT[ex.chrom][ex.strand][ex.start:ex.end]
        count_found=0
        for interval in found_intervals:
            
            found_gene = interval[2]
            all_gene_id_dict[ found_gene ] = 0
            
            found_r, count_found, count_checked =get_snaptron_exon_id_hits([exon_id], c, window)

            if count_found > 0:
                found_exon_id_list.append(exon_id)
                pre_mRNA_snaptron_hits_count += 1
                found_gene_dict[ found_gene  ] = 1

            
    search_count = len(exon_id_list)
    found_rate = len(found_gene_dict)/len(all_gene_id_dict)
    all_gene_rate = len(found_gene_dict)/len(aaa)
    
    found_exon_id_list = list(set(found_exon_id_list))
    
    return pre_mRNA_snaptron_hits_count, found_rate, all_gene_rate, found_exon_id_list
        


pre_mRNA_snaptron_hits_count, found_rate, all_gene_rate, found_exon_id_list =  get_snaptron_gene_stats(gencode_pc_transcript_IT, GENCODE_gene_dict,el.threshold_exon_ids( intron_interior_set, 10000, aggregate_exon_dict ),c) 
#print(pre_mRNA_snaptron_hits_count, found_rate, all_gene_rate )




threshold=10000
pre_mRNA_snaptron_hits_count, found_rate, all_gene_rate, found_exon_id_list =  get_snaptron_gene_stats(gencode_pc_transcript_IT, GENCODE_gene_dict,el.threshold_exon_ids( intron_interior_set, threshold, aggregate_exon_dict ),c) 




intronic_snaptron_found_exon_id_list = list()
genes_with_intron = list()
for ii, threshold in enumerate([100,300,700,1000,3000,7000,10000,15000]):
    pre_mRNA_snaptron_hits_count, found_rate, all_gene_rate, found_exon_id_list =  get_snaptron_gene_stats(gencode_pc_transcript_IT, GENCODE_gene_dict,el.threshold_exon_ids( intron_interior_set, threshold, aggregate_exon_dict ),c) 
    genes_with_intron.append([threshold, pre_mRNA_snaptron_hits_count, found_rate, all_gene_rate ])
    intronic_snaptron_found_exon_id_list += found_exon_id_list
    
    print('%d/%d found' % (len(found_exon_id_list),len(el.threshold_exon_ids( intron_interior_set, threshold, aggregate_exon_dict ))))

intronic_snaptron_found_exon_id_list  = list(set(intronic_snaptron_found_exon_id_list ))

len(intronic_snaptron_found_exon_id_list)












int('unneeded??')


































'''
x,y,z,z1 = zip(*genes_with_intron)

plt.figure()
plt.plot(x,z)
plt.title("percent of genes with intronic introns that are also \n snaptron verified intron")
plt.ylim(0,.25)
plt.xscale('log')
plt.xticks([100,1000,10000],["100","1000","10000"])
plt.yticks([0,.125,.25],["0%","12.5%","25%"])
plt.xlabel('Sequencing reads threshold')
plt.ylabel('percent found')


plt.figure()
plt.plot(x,z1)
plt.title("percent of genes with snaptron verified intronic intron")
plt.ylim(0,.15)
plt.xscale('log')
plt.xticks([100,1000,10000],["100","1000","10000"])
plt.yticks([0,.05,.1,.15],["0%","5%","10%","15%"])
plt.xlabel('Sequencing reads threshold')
plt.ylabel('percent found')

'''









import random
random.seed(2222)

def random_exon_discovery_rate(exon_id_list, c, genome_fasta):
    new_exon_id_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        new_start = random.randrange( len(genome_fasta[ex.chrom]) - 20000 ) + 10000
        new_end = new_start + ex.length
        new_exon_id = "%s:%d-%d:%s" % (ex.chrom, new_start, new_end, ex.strand)
        new_exon_id_list.append(new_exon_id)
        
        
    window = 3
    found_r, count_found, count_checked =get_snaptron_exon_id_hits(new_exon_id_list, c, window)

    return count_found


random.seed(2222)
random_count_found = random_exon_discovery_rate(intron_interior_set, c, genome_fasta)      
window = 3  
        
found_r, intronic_count_found, count_checked =get_snaptron_exon_id_hits(intron_interior_set, c, window)

print("{:,} intron interio exons found and {:,} random exons found with an equal number of exons of the same length at a random place on the same chromosome and strand.".format(intronic_count_found, random_count_found))











    




'''
threshold_found_dict={'Intergenic':list(),'Intronic':list()}
for key in region_exon_id_sets_dict:
    threshold_found_dict[key] = list()
    
threshold_list = [10**x for x in np.arange(2,4.5,.25)]
for key in region_exon_id_sets_dict:
    
    for ii, threshold in enumerate(threshold_list):
        
        #if key not in ['Intergenic','Intronic']:
        #    continue
        
        exon_id_list = region_exon_id_sets_dict[key]   
        exon_id_list = el.exon_id_intersection(exon_id_list, primary_3ss_exon_id_set)
        exon_id_list = el.threshold_exon_ids(exon_id_list, threshold, aggregate_exon_dict)
        if key == 'Intergenic':
            exon_id_list = set(exon_id_list).difference(lncRNA_database_exons) # requires gencode #2 script
        
        
        window = 3
        start_time = time.time()
        found_r, count_found, count_checked = get_snaptron_exon_id_hits(exon_id_list, c, window)
        end_time = time.time()-start_time
        print("{:} search completed in {:.2f} seconds".format(key, end_time))
        
        
        found_keys = list(found_r.keys())
        found_fraction = count_found/count_checked
        
        print('threshold {:.0f} \twith {:.3f} found'.format(threshold, found_fraction))
        
        threshold_found_dict[key].append([threshold, found_fraction, count_found])
'''
        

'''




termini_found_fraction_dict = dict()
for key in region_exon_id_sets_dict:
    #exon_id_list = el.threshold_exon_ids(region_exon_id_sets_dict[key], 1000, aggregate_exon_dict)
    
    exon_id_list = region_exon_id_sets_dict[key]    
    
    
    try:
        exon_id_list = el.threshold_exon_ids(region_exon_id_sets_dict[key], 100, aggregate_exon_dict)
    except:
        exon_id_list = region_exon_id_sets_dict[key]    
    
    
    
    
    exon_id_list = el.exon_id_intersection(exon_id_list, primary_3ss_exon_id_set)
    
    window = 3
    start_time = time.time()
    found_r, count_found, count_checked = get_snaptron_exon_id_hits(exon_id_list, c, window)
    end_time = time.time()-start_time
    print("{:} search completed in {:.2} seconds".format(key, end_time))
    
    
    found_keys = list(found_r.keys())
    missing_exon_ids = list( set(exon_id_list).difference(found_keys))
    
    
    
    len(found_r)
    found_fraction = count_found/count_checked
    
    print("{:}/{:} = {:.4f} exon_ids".format(len(found_r),len(exon_id_list),found_fraction))
    print("{:}/{:} = {:.4f} found_exons\n\n".format(count_found,len(found_r),count_found/len(found_r)))
    
    termini_found_fraction_dict[key]= {'found_r':found_r, 'count_found':count_found, 'count_checked':count_checked}
    termini_found_fraction_dict[key]['missing_exon_ids']=missing_exon_ids    
    termini_found_fraction_dict[key]['found_fraction']=found_fraction 
    termini_found_fraction_dict[key]['found_keys']=found_keys 
    
    exon_id_label=key
    #plot_found_missing_stats(found_keys, missing_exon_ids, exon_id_label)






found_rate_list = list()
found_rate_keys = list()
for ii, key in enumerate(termini_found_fraction_dict):
    found_rate = 100*termini_found_fraction_dict[key]['found_fraction']
    found_rate_keys.append(key)
    
    found_rate_list.append(found_rate)
    #found_rate_keys.append(key)



termini_found_fraction_dict['Intronic']['found_keys'][:20]





found_label_list = list()
found_c_list = list()
missing_c_list = list()
m_color_list = list()
f_color_list = list()
f_positions = list()
m_positions = list()

found_3ss_list = list()
missing_3ss_list = list()

found_5ss_list = list()
missing_5ss_list = list()

found_len_list = list()
missing_len_list = list()


for ii, key in enumerate(termini_found_fraction_dict.keys()):
    f = termini_found_fraction_dict[key]['found_keys']
    m = termini_found_fraction_dict[key]['missing_exon_ids']
    
    f=el.exon_id_intersection(f, aggregate_exon_dict.keys())
    m=el.exon_id_intersection(m, aggregate_exon_dict.keys())
    
    found_c = el.get_exon_dict_counts(f,aggregate_exon_dict)
    missing_c = el.get_exon_dict_counts(m,aggregate_exon_dict)
    f_color_list.append('yellow')
    m_color_list.append('cyan')
    found_label_list.append(key)
    
    found_c_list.append(found_c)
    missing_c_list.append(missing_c)
    f_positions.append(ii+1-.15)
    m_positions.append(ii+1+.15)
    
    found_3ss = el.get_exon_dict_3ss_scores(f,aggregate_exon_dict)
    missing_3ss = el.get_exon_dict_3ss_scores(m,aggregate_exon_dict)
    
    found_5ss = el.get_exon_dict_5ss_scores(f,aggregate_exon_dict)
    missing_5ss = el.get_exon_dict_5ss_scores(m,aggregate_exon_dict)
    
    found_len = el.get_exon_dict_lengths(f,aggregate_exon_dict)
    missing_len = el.get_exon_dict_lengths(m,aggregate_exon_dict)
    
    found_3ss_list.append(found_3ss)
    missing_3ss_list.append(missing_3ss)

    found_5ss_list.append(found_5ss)
    missing_5ss_list.append(missing_5ss)

    found_len_list.append(found_len)
    missing_len_list.append(missing_len)
    








'''
































pdf_plots.close()






c.close()
conn.close()



