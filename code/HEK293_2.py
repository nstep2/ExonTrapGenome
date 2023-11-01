




"""

pip install pyensembl

#pip install pyEntrezId




"""

import requests
import pyensembl 


from pyensembl import EnsemblRelease
# specify the desired Ensembl release
#Run pyensembl install --release 103 --species homo_sapiens
ensembl_data = EnsemblRelease(103)  #GENCODE version 37 corresponds to Ensembl 103.







outdir = exp_output_path.HEK293_recovery +'HEK293.pdf'
pdf_plots = PdfPages(outdir)









#alt_splice_input_file = '/home/dude/work/HEK293_alt_splicing/GSE221838_d0_v_d1000_SE.MATS.JCEC.csv'

alt_splice_input_file = exp_output_path.HEK293_PSI_input #+ 'GSE221838_d0_v_d1000_SE.MATS.JCEC.csv'




exon_percents_dict = dict()
NA_skip_count = 0
with open(alt_splice_input_file) as f:
    header = next(f)
    for jj, line in enumerate(f):
        line_split = line.strip().split('\t')
        percents_string = line_split[20]
        percents_string_split = percents_string.split(',')
        if 'NA' in percents_string_split:
            NA_skip_count += 1
            continue
        percents = [float(x) for ii,x in enumerate(percents_string_split)]
        
        exon_coords = "%s:%s-%s:%s" % (line_split[3],line_split[5],line_split[6],line_split[4])
        exon_percents_dict[exon_coords] = {'percents':percents,'gene_id':line_split[1], 'gene_name':line_split[2]}







multi_percent_ids_count = 0
multi_exon_skip_count = 0
alt_chrom_count = 0
check_keys = gencode_pc_exon_IT.keys()
exon_id_to_exon_percent = dict()
new_percents_list = list()
exons_not_in_pc_dict = dict()
exons_in_pc_dict     = dict()
multi_HEK293_exon_list = list()
for key in exon_percents_dict:
    ex = el.exon_id_values(key)
    start  = ex.start
    end    = ex.end
    chrom  = ex.chrom
    strand = ex.strand
    
    new_exon_id = "%s:%d-%d:%s" % (chrom, start+1, end+1, strand) 
    
    if chrom not in check_keys:
        alt_chrom_count += 1
        continue
    intervals = gencode_pc_exon_IT[chrom][strand][start:end]
    if len(intervals) > 1:
        multi_exon_skip_count += 1
        #continue
    if len(intervals) == 0:
        exons_not_in_pc_dict[key] =  exon_percents_dict[key]
        continue
    
    new_percents_list.append(np.mean(exon_percents_dict[key]['percents']))
    for interval in intervals:   #only one interval at this point, use counts of the max one, input file has several exons off by a few base pairs
        if interval[2] not in exons_in_pc_dict:
            exons_in_pc_dict[interval[2]] = [key, np.mean(exon_percents_dict[key]['percents']), exon_percents_dict[key], new_exon_id == interval[2]]
        else:
            exons_in_pc_dict[interval[2]][1] = max(np.mean(exon_percents_dict[key]['percents']), exons_in_pc_dict[interval[2]][1])
            multi_HEK293_exon_list.append(interval[2])
            multi_percent_ids_count += 1




means = list()
for key in exon_percents_dict:
    means.append(np.mean(exon_percents_dict[key]['percents']))    





PSI_vs_exon_reads = list()

ET_found_count  = 0
ET_missed_count = 0

found_exon_percents  = list()
missed_exon_percents = list()
for key in exons_in_pc_dict:
    if key in aggregate_exon_dict:
        if exons_in_pc_dict[key][3] == False:
            continue
        
        ET_found_count += 1
        found_exon_percents.append(exons_in_pc_dict[key][1])
        
        PSI_vs_exon_reads.append([exons_in_pc_dict[key][1], aggregate_exon_dict[key]['count']])
    else:
        ET_missed_count += 1
        missed_exon_percents.append(exons_in_pc_dict[key][1])
    




'''
fig = plt.figure()
plt.title('Density - exons PSI for exons either found or no found in HEK293 ')
plt.hist([found_exon_percents, missed_exon_percents], density=True)
plt.legend(['in exon trapping data', 'not in exon trapping data'])
plt.ylabel('count')
plt.xlabel('Pecent Spliced Index')
pdf_plots.savefig(fig)
'''


fig = plt.figure(figsize=(6,4))
plt.title('Absolute - exons PSI for exons either found or no found in HEK293 ')
plt.hist([found_exon_percents, missed_exon_percents], density=False, orientation='horizontal')
plt.legend(['in exon trapping data', 'not in exon trapping data'])
plt.ylabel('count')
plt.xlabel('Pecent Spliced Index')
pdf_plots.savefig(fig)


with open(exp_output_path.HEK293_recovery + 'S4C.txt', 'w') as f:
    f.write('{:}\t{:}\n'.format('ET_found','ET_missing'))
    
    for ii, val in enumerate(found_exon_percents):
        f.write('{:}\t{:}\n'.format(found_exon_percents[ii], missed_exon_percents[ii]))
        


from matplotlib.colors import LogNorm




import scipy.stats as stats

y,x = zip(*PSI_vs_exon_reads)

# specify bin edges
bin_edges = np.arange(0, 1.1, 0.1)

x=np.array(x)
y=np.array(y)

# bin value1 into groups with specified bin edges
bin_means, _, binnumber = stats.binned_statistic(y, x, statistic='mean', bins=bin_edges)

fig = plt.figure()
# create boxplot of value2 for each group
plt.boxplot([x[binnumber == i] for i in range(1, len(bin_edges))], showfliers=False)
plt.yscale('log')
plt.ylabel('exon trapping sequencing read count')
plt.xlabel('Pecent Spliced In')
pdf_plots.savefig(fig)





pdf_plots.close()



with open(exp_output_path.HEK293_recovery + 'S4B.txt', 'w') as f:
    f.write('{:}\t{:}\n'.format('ET_read_count','PSI'))
    
    for ii, val in enumerate(x):
        f.write('{:}\t{:}\n'.format(x[ii], y[ii]))
        





















