

pairwise_lib_table = np.zeros([24,24])

for exon_id in aggregate_exon_dict:
    exon = aggregate_exon_dict[exon_id]
    for ii in range(1,24):
        for jj in range(1,24):
            if exon['lib_array'][ii] >0 and exon['lib_array'][jj] > 0:
                pairwise_lib_table[ii][jj] += 1
    
pairwise_lib_table[:,1:]


path =  exp_output_path.library_scatter + 'pairwise_lib_counts.txt'

with open(path,'w') as f:
    for ii in range(1,24):
        f.write('\t%d' % (ii) )    
    f.write('\n')
    
    for ii in range(1,24):
        f.write('%d' % (ii) )    
        for jj in range(1,24):
            f.write('\t%d' % pairwise_lib_table[ii][jj])
        f.write('\n')



for exon_id in aggregate_exon_dict:
    exon = aggregate_exon_dict[exon_id]
    for ii in range(1,24):
        for jj in range(1,24):
            if ii in exon['lib_array'] and jj in exon['lib_array']:
                pairwise_lib_table[ii][jj] += 1











library_codes = dict()
code_to_libraries = dict()

library_codes[1]='1'
library_codes[2]='1'
library_codes[3]='1'
library_codes[4]='1'
library_codes[5]='1'
code_to_libraries['1'] = [1,2,3,4,5]

library_codes[6]='2'
library_codes[7]='2'
library_codes[21]='2'
library_codes[22]='2'
code_to_libraries['2'] = [6,7,21,22]

library_codes[8]='3'
library_codes[9]='3'
library_codes[10]='3'
library_codes[11]='3'
library_codes[23]='3'
code_to_libraries['3'] = [8,9,10,11,23]

library_codes[12]='4'
library_codes[13]='4'
library_codes[14]='4'
library_codes[15]='4'
code_to_libraries['4'] = [12,13,14,15]

library_codes[16]='5'
library_codes[17]='5'
library_codes[18]='5'
library_codes[19]='5'
library_codes[20]='5'
code_to_libraries['5'] = [16,17,18,19,20]





pairwise_lib_group_table = np.zeros([6,6])

for exon_id in aggregate_exon_dict:
    lib_group_counts_dict = {str(ii):0 for ii in range(1,6)}
    lib_group_counts_sum_dict = {str(ii):0 for ii in range(1,6)}
    exon = aggregate_exon_dict[exon_id]
    
    for ii in exon['lib_count']:
        lib_group_counts_dict[library_codes[ii]]=1

    
    for ii in lib_group_counts_dict:
        for jj in lib_group_counts_dict:
            if lib_group_counts_dict[ii] >= 1 and lib_group_counts_dict[jj] >= 1:
                pairwise_lib_group_table[int(ii)][int(jj)] += 1
        

#pairwise_lib_group_count_table


path =  exp_output_path.library_scatter + 'pairwise_lib_group_counts.txt'

with open(path,'w') as f:
    for ii in range(1,6):
        f.write('\t%d' % (ii) )    
    f.write('\n')
    
    for ii in range(1,6):
        f.write('%d' % (ii) )    
        for jj in range(1,6):
            f.write('\t%d' % pairwise_lib_group_table[ii][jj])
        f.write('\n')








outdir = exp_output_path.library_scatter
pdf_plots = PdfPages(outdir+'libraries_per_exon.pdf')








exon_library_count_list = list()
for exon_id in aggregate_exon_dict:
    exon = aggregate_exon_dict[exon_id]
    exon_library_count_list.append( (np.array([exon['lib_count'][k] for k in exon['lib_count']]) > exon['count']/10000).sum() )



pc_exon_library_count_list = list()
for exon_id in set(pc_middle_exon_id_list).intersection(aggregate_exon_dict.keys()):
    exon = aggregate_exon_dict[exon_id]
    pc_exon_library_count_list.append( ( np.array( [exon['lib_count'][k] for k in exon['lib_count']] ) > (exon['count']/10000) ).sum() )
    




fig = plt.figure()
plt.hist(exon_library_count_list, bins=23)
plt.title('number of libraries per exon')
plt.ylabel('count')
plt.xlabel('number of libraries')
pdf_plots.savefig(fig)



fig = plt.figure()
plt.hist(pc_exon_library_count_list, bins=23)
plt.title('number of libraries per mRNA exon')
plt.ylabel('count')
plt.xlabel('number of libraries')
pdf_plots.savefig(fig)



pdf_plots.close()




















