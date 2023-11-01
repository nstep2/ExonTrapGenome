import scipy 
import statistics

pickle_path = exp_output_path.pickle_merged + "ET_SAI_MES_exon_sets_%d.pickle" % (exon_count_build)


exon_id_cat_dict = region_exon_id_sets_dict



finder_unique=dict()
with open(pickle_path, "rb") as input_file:
    D = pickle.load(input_file)
    
    finder_unique['ET_only'] = pickle.load(input_file)
    finder_unique['SAI_only'] = pickle.load(input_file)
    finder_unique['MES_only'] = pickle.load(input_file)

    ALL_exon_intersect = pickle.load(input_file)
    ALL_exon_finder_intersect = pickle.load(input_file)

finder_unique['annotated']=exon_id_cat_dict['mRNA']


with open(exp_output_path.pickle_main+'chasin_ESEseq.pickle', 'rb') as infile:
    ESE_scorer = pickle.load(infile)



def kmer_hash(kmer):
    kmer = kmer.upper()
    if kmer.find('N') >= 0:
        return -1
    
    base_to_num = {'A':0,'C':1,'G':2,'T':3}
    kmer_hash = 0 
    for ii in range(len(kmer)):
        kmer_hash = 4*kmer_hash
        kmer_hash += base_to_num[kmer[ii]]
    return kmer_hash


def num_hash_kmer(hash_num, kmer_len):
    kmer_str = ''
    num_to_base = {0:'A',1:'C',2:'G',3:'T'}
    ii=hash_num
    for jj in range(kmer_len):
        base = num_to_base[int(hash_num/4**(kmer_len-jj-1))%4]
        kmer_str = '%s%s' % (kmer_str, base)
    return kmer_str
    


def fraction_covered(found_PESE,exon_len):
    exon_range = np.zeros(exon_len)
    
    for hit in found_PESE:
        exon_range[hit[0]:hit[0]+len(hit[2])] += 1
    
    fraction_covered = sum(exon_range>0)/exon_len
    return fraction_covered





def get_PESE_per_100_bp_exon_id_list(exon_id_list, PESE, genome_fasta):
    count_PESE = 0
    total_seq_len = 0
    exon_len_vs_count_list = list()
    exon_fraction_covered_list = list()
    #defined_found_PESE = False
    
    kmer_array = np.zeros(4**PESE.motif_len)
    
    found_ESE_list = list()
    
    for exon_id in exon_id_list:
        if el.exon_id_values(exon_id).chrom == 'chrM':
            continue
        
            
        seq = el.get_seq_from_exon_id(exon_id, genome_fasta)
        found_PESE = PESE.score_seq(seq)
        
        if found_PESE != -1:
            
            for entry in found_PESE:
                if entry[1] > PESE.threshold:
                    index_hash = kmer_hash(entry[2])
                    kmer_array[index_hash] += 1
            
            total_seq_len+=len(seq)
            count_PESE += len(found_PESE)
            found_ESE_list.append(len(found_PESE)/len(seq))
            exon_len_vs_count_list.append([len(seq),len(found_PESE)])
            frac_covered = fraction_covered(found_PESE,len(seq))
            exon_fraction_covered_list.append(frac_covered)
        
        if total_seq_len != 0:
            found_rate_100_bp = 100*count_PESE/total_seq_len
        else:
            found_rate_100_bp =0
        
    mean  = np.mean(found_ESE_list)
    stdev = statistics.stdev(found_ESE_list)
        
    return kmer_array, mean, stdev, exon_len_vs_count_list, exon_fraction_covered_list





def generate_250_bp_offset_exon_ids(exon_id_list, offset):
    offset_exon_id_list = list()
    for exon_id in exon_id_list:
        ex = el.exon_id_values(exon_id)
        new_ex = "%s:%d-%d:%s" % (ex.chrom,ex.start+offset,ex.end+offset,ex.strand)
        offset_exon_id_list.append(new_ex)
    return offset_exon_id_list







ESE_finder_unique_vals_dict = dict()

for key in finder_unique:
    ESE_finder_unique_vals_dict[key]=dict()
    
    summed_found_PESE, found_rate_100_bp, stdev, exon_len_vs_count_list, exon_fraction_covered_list = get_PESE_per_100_bp_exon_id_list(finder_unique[key], ESE_scorer, genome_fasta)
    
    total_ESE_bp = [x[0]-5 for x in exon_len_vs_count_list]
    
    ESE_finder_unique_vals_dict[key]['ESE_rate'] = found_rate_100_bp
    ESE_finder_unique_vals_dict[key]['ESE_stdev'] = stdev
    ESE_finder_unique_vals_dict[key]['summed_found_PESE'] = summed_found_PESE/sum(total_ESE_bp)
    
    
    offset_exon_ids=generate_250_bp_offset_exon_ids(finder_unique[key], 150)
    summed_found_PESE, found_rate_100_bp, stdev, exon_len_vs_count_list, exon_fraction_covered_list = get_PESE_per_100_bp_exon_id_list(offset_exon_ids, ESE_scorer, genome_fasta)
    ESE_finder_unique_vals_dict[key]['offset_ESE_rate'] = found_rate_100_bp
    ESE_finder_unique_vals_dict[key]['offset_ESE_stdev'] = stdev












ESE_finder_all_vals_dict = dict()

for key in D:
    ESE_finder_all_vals_dict[key]=dict()
    
    summed_found_PESE, found_rate_100_bp, stdev, exon_len_vs_count_list, exon_fraction_covered_list = get_PESE_per_100_bp_exon_id_list(D[key], ESE_scorer, genome_fasta)
    
    total_ESE_bp = [x[0]-5 for x in exon_len_vs_count_list]
    
    ESE_finder_all_vals_dict[key]['ESE_rate'] = found_rate_100_bp
    ESE_finder_all_vals_dict[key]['summed_found_PESE'] = summed_found_PESE/sum(total_ESE_bp)
    
    
    
    









outdir = exp_output_path.ESE_scan 
pdf_plots = PdfPages(outdir+'ET_SAI_MES_exon_ESE.pdf')





ESE_finder_all_vals_dict


colors = {'annotated':'#D6D6D6','ET_only':'#31BBED','SAI_only':'#FBAA28','MES_only':'#EBE837'}


fig, ax = plt.subplots()
xlabels = list()
offset_labels=list()

for ii, key in enumerate(ESE_finder_unique_vals_dict):
    vals = 100*ESE_finder_unique_vals_dict[key]['ESE_rate']
    stdev = 100*ESE_finder_unique_vals_dict[key]['ESE_stdev']
    vals_offset = 100*ESE_finder_unique_vals_dict[key]['offset_ESE_rate']
    stdev_offset = 100*ESE_finder_unique_vals_dict[key]['offset_ESE_stdev']
    
    plt.bar([ii-.15], vals, yerr=stdev, label=[key], color=colors[key], width=[.3])
    plt.bar([ii+.15], vals_offset, yerr=stdev_offset, label=[key], color=colors[key], alpha=.5, width=[.3])
    xlabels.append( key )
    offset_labels.append( key+'_offset' )

plt.xticks(list(np.arange(-.15,len(xlabels)-.15,1)) + list(np.arange(.15,len(xlabels)+.15,1)), xlabels+offset_labels, rotation=30)

plt.ylabel('ESE count')
plt.title('exon finder unique exons')
plt.ylim(0,40)
pdf_plots.savefig(fig)


with open(exp_output_path.out_supplemental+'5B.txt','w') as f:
    f.write('category\tmedian\tstddev\toffset_median\toffset_stddev\n')
    
    for ii, key in enumerate(ESE_finder_unique_vals_dict):
        vals = 100*ESE_finder_unique_vals_dict[key]['ESE_rate']
        stdev = 100*ESE_finder_unique_vals_dict[key]['ESE_stdev']
        vals_offset = 100*ESE_finder_unique_vals_dict[key]['offset_ESE_rate']
        stdev_offset = 100*ESE_finder_unique_vals_dict[key]['offset_ESE_stdev']
        
        out_string = '{:}\t{:}\t{:}\t{:}\t{:}\n'.format(key,vals,stdev,vals_offset,stdev_offset)
        f.write(out_string)
        





pdf_plots.close()






