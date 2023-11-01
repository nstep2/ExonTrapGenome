



import intervaltree

import exon_id_library.exon_id_lib as el
import numpy as np



outdir = exp_output_path.calculate_chromosome_coverage
pdf_plots = PdfPages(outdir+'calculate_chromosome_coverage_%d_reads.pdf'  % (exon_count_build))







threshold = 100

average_coverage_list = list()
for chrom in aggregate_exon_IT:
    chrom_bases_list = list()
    for strand in aggregate_exon_IT[chrom]:
        strand_chrom_bases_list = list()
        for interval in aggregate_exon_IT[chrom][strand]:
            for exon_id in interval[2]:
                if threshold <= aggregate_exon_dict[exon_id]['count']:
                    ex = el.exon_id_values(exon_id)
                    bases = range(ex.start,ex.end)
                    chrom_bases_list += bases
                    strand_chrom_bases_list += bases
                    
                
        genome_len = len(genome_fasta[chrom])
        genome_coverage = len(set(strand_chrom_bases_list))/genome_len
        print('\t',chrom, strand, genome_coverage)
        average_coverage_list.append(genome_coverage)
    
    genome_len = len(genome_fasta[chrom])
    genome_coverage = len(set(chrom_bases_list))/genome_len
    print(chrom, genome_coverage)
    
average_coverage = np.mean(average_coverage_list)
print("Average chromosome stranded coverage: %.1f%%" % (100*average_coverage))
      


       
lncRNA_transcript_exon_overlap = list()
pc_transcript_exon_overlap = list()

lncRNA_antisense_transcript_exon_overlap = list()
pc_antisense_transcript_exon_overlap = list()

no_overlap_set = list()
intron_interior_set = list()
sense_transcript_set = list()
antisense_transcript_set = list()
exon_overlap_set = list()
lncRNA_overlap_set = list()
anti_strand = {'+':'-','-':'+'}

exon_id_list = list(primary_3ss_exon_id_set)
for ii, exon_id in enumerate(exon_id_list):
    ex = el.exon_id_values(exon_id)
    
    
    
    
    intervals = gencode_exon_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) > 0:
        exon_overlap_set.append(exon_id)
    
    intervals = gencode_lncRNA_transcript_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) > 0:
        lncRNA_overlap_set.append(exon_id)
    
    intervals = gencode_transcript_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) > 0:
        sense_transcript_set.append(exon_id)
    
    intervals = gencode_pc_transcript_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) > 0:
        pc_transcript_exon_overlap.append(exon_id)
    
    intervals = gencode_lncRNA_transcript_IT[ex.chrom][ex.strand][ex.start:ex.end]
    if len(intervals) > 0:
        lncRNA_transcript_exon_overlap.append(exon_id)
    
    intervals = gencode_pc_transcript_IT[ex.chrom][anti_strand[ex.strand]][ex.start:ex.end]
    if len(intervals) > 0:
        pc_antisense_transcript_exon_overlap.append(exon_id)
    
    intervals = gencode_lncRNA_transcript_IT[ex.chrom][anti_strand[ex.strand]][ex.start:ex.end]
    if len(intervals) > 0:
        lncRNA_antisense_transcript_exon_overlap.append(exon_id)
    
    intervals = gencode_transcript_IT[ex.chrom][anti_strand[ex.strand]][ex.start:ex.end]
    if len(intervals) > 0:
        antisense_transcript_set.append(exon_id)








lncRNA_transcript_exon_overlap = set(lncRNA_transcript_exon_overlap)
pc_transcript_exon_overlap = set(pc_transcript_exon_overlap)

lncRNA_antisense_transcript_exon_overlap = set(lncRNA_antisense_transcript_exon_overlap)
pc_antisense_transcript_exon_overlap = set(pc_antisense_transcript_exon_overlap)

no_overlap_set          = set(no_overlap_set)
intron_interior_set     = set(intron_interior_set)
sense_transcript_set    = set(sense_transcript_set)
antisense_transcript_set = set(antisense_transcript_set)
exon_overlap_set        = set(exon_overlap_set)
lncRNA_overlap_set      = set(lncRNA_overlap_set)


lncRNA_overlap_set = lncRNA_overlap_set.difference(exon_overlap_set)


pc_antisense_transcript_exon_overlap = pc_antisense_transcript_exon_overlap.difference(sense_transcript_set)
len(pc_antisense_transcript_exon_overlap)

lncRNA_antisense_transcript_exon_overlap = lncRNA_antisense_transcript_exon_overlap.difference(sense_transcript_set)
len(lncRNA_antisense_transcript_exon_overlap)





other_antisense_transcript_set = antisense_transcript_set.difference(sense_transcript_set).difference(lncRNA_antisense_transcript_exon_overlap).difference(pc_antisense_transcript_exon_overlap)


antisense_transcript_set = antisense_transcript_set.difference(sense_transcript_set).difference(other_antisense_transcript_set)




intron_interior_set     = set(pc_transcript_exon_overlap).union(lncRNA_transcript_exon_overlap).difference(exon_overlap_set)
intron_interior_pc_set  = intron_interior_set.intersection( pc_transcript_exon_overlap).difference(lncRNA_transcript_exon_overlap)
intron_interior_lncRNA_set  = intron_interior_set.intersection( lncRNA_transcript_exon_overlap).difference(intron_interior_pc_set)
other_intron_interior_set   = intron_interior_set.difference(intron_interior_pc_set).difference(intron_interior_lncRNA_set)


no_overlap_set = set(primary_3ss_exon_id_set)
no_overlap_set = no_overlap_set.difference(sense_transcript_set)
no_overlap_set = no_overlap_set.difference(antisense_transcript_set)


mRNA_exon_overlap   = pc_transcript_exon_overlap.intersection(exon_overlap_set)
lncRNA_exon_overlap = lncRNA_transcript_exon_overlap.intersection(exon_overlap_set)


labels = ['mRNA exon overlapping','mRNA intron','mRNA antisense','lncRNA exon overlapping','lncRNA intron','lncRNA antisense','Intergenic']

sizes  = [mRNA_exon_overlap,intron_interior_pc_set,pc_antisense_transcript_exon_overlap, lncRNA_exon_overlap,intron_interior_lncRNA_set,lncRNA_antisense_transcript_exon_overlap,no_overlap_set]

with open(exp_output_path.out_supplemental+'3A.txt','w') as f:
    outstring = 'chrom\tstart\tend\tstrand\texon_id\tgenomic_category\n'
    f.write(outstring)
    
    for ii, val in enumerate(labels):
        exon_id_list = sizes[ii]
        for exon_id in exon_id_list:
            ex=el.exon_id_values(exon_id)
            
            outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id)
            
            outstring = '{:}\t{:}\n'.format(outstring, val)
            f.write(outstring)
        


colors = list()
colors.append([color_dict['mRNA'],1])
colors.append([color_dict['Intronic'],1])
colors.append([color_dict['Antisense'],1])
colors.append([color_dict['lncRNA'],1])
colors.append([color_dict['Intronic'],.7])
colors.append([color_dict['Antisense'],.7])
colors.append([color_dict['Intergenic'],1])

c,a = zip(*colors)

sizes_len=list()
for ii, exon_set in enumerate(sizes):
    exon_set=el.threshold_exon_ids(exon_set,100, aggregate_exon_dict)
    sizes_len.append(len(exon_set))

fig, ax1 = plt.subplots()
wedg, txt, atxt = ax1.pie(sizes_len, labels=labels, autopct='%1.1f%%',
        shadow=False, startangle=90, counterclock=False, colors =c )
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
for ii, w in enumerate(wedg):
    w.set_alpha(a[ii])
plt.show()
plt.tight_layout()
pdf_plots.savefig(fig)







labels = ['mRNA', 'lncRNA','mRNA antisense','lncRNA antisense','mRNA intron','lncRNA intron','Intergenic']
colors = list()
colors.append([color_dict['mRNA'],1])
colors.append([color_dict['lncRNA'],1])
colors.append([color_dict['Antisense'],1])
colors.append([color_dict['Antisense'],.7])
colors.append([color_dict['Intronic'],1])
colors.append([color_dict['Intronic'],.7])
colors.append([color_dict['Intergenic'],1])

sizes  = [(exon_overlap_set), lncRNA_overlap_set,(pc_antisense_transcript_exon_overlap),(lncRNA_antisense_transcript_exon_overlap),(intron_interior_pc_set),(intron_interior_lncRNA_set),(no_overlap_set)]
sizes_len=list()
for ii, exon_set in enumerate(sizes):
    exon_set=el.threshold_exon_ids(exon_set,100, aggregate_exon_dict)
    sizes_len.append(len(exon_set))

c,a = zip(*colors)
fig, ax1 = plt.subplots(figsize=(6,6))
wedg, txt, atxt = ax1.pie(sizes_len, labels=labels, colors =c , autopct='%1.1f%%',
        shadow=False, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
for ii, w in enumerate(wedg):
    w.set_alpha(a[ii])
    
for ii, at in enumerate(atxt):
    val = at.get_text()
    at.set_text("{0}\n{1:,}".format(val,sizes_len[ii]))

plt.title('proportion of sequencing reads in different regions')

plt.show()
plt.tight_layout()
pdf_plots.savefig(fig)








above_100_primary_3ss_no_overlap = no_overlap_set.intersection(primary_3ss_exon_id_set)

exon_length_sum = 0
for exon_id in above_100_primary_3ss_no_overlap:
    ex = el.exon_id_values(exon_id)
    exon_length_sum += ex.length
    

print('primary_3ss_and_above_100_no_overlap_set: %s' % insert_commas(len(no_overlap_set.intersection(primary_3ss_exon_id_set))))

print('bases in primary_3ss_and_above_100_no_overlap_set : %s' % insert_commas(exon_length_sum))







###
### Calculate the intron portion of each transcript
### Also calc. exon, transcript bases

import time

def calculate_intronic_bases(gencode_pc_exon_IT,gencode_pc_transcript_IT):
    
    total_number_GENCODE_pc_intron_bases = 0
    total_number_GENCODE_pc_exon_bases = 0
    total_number_GENCODE_pc_transcript_bases = 0
    
    
    for chrom in gencode_pc_transcript_IT:
        for strand in gencode_pc_transcript_IT[chrom]:
            start_time = time.time()
            print(chrom,strand)
            
            transcript_bases = set()
            exon_bases = set()
            transcript_bases_list = list()
            exon_bases_list = list()
            
            for ii, interval in enumerate(gencode_pc_transcript_IT[chrom][strand]):
                if ii % 16 == 0: #break the lists into smaller sets to decrease processing time
                    transcript_bases_list.append(transcript_bases)
                    exon_bases_list.append(exon_bases)
                    transcript_bases = set()
                    exon_bases = set()
                if ii % 1000==0:
                    print('%d/%d processed in %.1f seconds' % (ii,len(gencode_pc_transcript_IT[chrom][strand]),time.time() - start_time))
                
                start = interval[0]
                end   = interval[1]
                gencode_exons = gencode_pc_exon_IT[chrom][strand][start:end] 
                transcript_bases=transcript_bases.union(range(start,end))
                #transcript_bases = set(range(start,end))
                for exon_interval in gencode_exons:
                    exon_id = exon_interval[2]

                    ex = el.exon_id_values(exon_id)
                    exon_bases = exon_bases.union(range(ex.start,ex.end))
    
            transcript_bases_list.append(transcript_bases)
            exon_bases_list.append(exon_bases)
            
            from itertools import chain
            union_transcript_bases_list =  set(chain.from_iterable(transcript_bases_list))
            union_exon_bases_list = set(chain.from_iterable(exon_bases_list))
    
            total_number_GENCODE_pc_intron_bases += len(transcript_bases.difference(exon_bases))    
            total_number_GENCODE_pc_exon_bases += len(exon_bases)
            total_number_GENCODE_pc_transcript_bases += len(transcript_bases)
    
    transcript_bases = set() #to release memory
    exon_bases = set() #to release memory
    
    union_transcript_bases_list=set()
    union_exon_bases_list=set()
    
    print( "total intron bases (ignoring overlapping lncRNA): %s" %  insert_commas(total_number_GENCODE_pc_intron_bases))
    print( "total exon bases (ignoring overlapping lncRNA): %s" %  insert_commas(total_number_GENCODE_pc_exon_bases))
    print( "total transcript bases (ignoring overlapping lncRNA): %s" %  insert_commas(total_number_GENCODE_pc_transcript_bases))

calculate_intronic_bases(gencode_lncRNA_exon_IT,gencode_lncRNA_transcript_IT)


































def above_threshold_exon_id(threshold, exon_id_set, aggregate_exon_dict):
    count_above_threshold = 0
    for exon_id in exon_id_set:
        if aggregate_exon_dict[exon_id]['count'] >= threshold:
            count_above_threshold += 1
            #print(exon_id)
    return count_above_threshold

def above_threshold_exon_id_set(threshold, exon_id_set, aggregate_exon_dict):
    count_above_threshold_set = list()
    for exon_id in exon_id_set:
        if aggregate_exon_dict[exon_id]['count'] >= threshold:
            count_above_threshold_set.append(exon_id)
            #print(exon_id)
    return set(count_above_threshold_set)



















region_exon_id_sets_list = list()
region_exon_id_sets_list.append([set(pc_middle_exon_id_list), 'mRNA'])
region_exon_id_sets_list.append([set(lncRNA_middle_exon_ids), 'lncRNA'])
region_exon_id_sets_list.append([no_overlap_set, 'Intergenic'])
region_exon_id_sets_list.append([antisense_transcript_set, 'Antisense'])
region_exon_id_sets_list.append([intron_interior_set, 'Intronic'])


region_exon_id_sets_dict = dict()
region_exon_id_sets_dict['mRNA'] = set(pc_middle_exon_id_list)
region_exon_id_sets_dict['lncRNA'] = set(lncRNA_middle_exon_ids)
region_exon_id_sets_dict['Intergenic'] = no_overlap_set
region_exon_id_sets_dict['Antisense'] = antisense_transcript_set
region_exon_id_sets_dict['Intronic'] = intron_interior_set














def gc_list_len(gc_result):
    total_len = 0
    for ii, key in enumerate(gc_result):
        total_len += len(gc_result[key])
    fraction_list = list()
    key_list=list()
    for ii, key in enumerate(gc_result):
        print(key, len(gc_result[key])/total_len)
        fraction_list.append(len(gc_result[key])/total_len)
        key_list.append(key)
    return fraction_list, key_list






threshold = 1
gc_bins=[25,40,50,60,75]

region_list_gc_dict = dict()
for key in region_exon_id_sets_dict:
    
    exon_id_list = region_exon_id_sets_dict[key]
    #exon_id_set_result=above_threshold_exon_id_set(threshold,intron_interior_set, aggregate_exon_dict)
    
    gc_result = el.get_exon_id_list_binned_GC(exon_id_list, gc_bins, genome_fasta)
    
    region_list_gc_dict[key], key_list = gc_list_len(gc_result)
    








fig=plt.figure()
plt.bar(np.arange(0,len(region_list_gc_dict['mRNA']),1),region_list_gc_dict['mRNA'], width=0.2)
plt.bar(np.arange(0.2,len(region_list_gc_dict['mRNA'])+.2,1),region_list_gc_dict['Intronic'], width=0.2)
plt.bar(np.arange(0.4,len(region_list_gc_dict['mRNA'])+.4,1),region_list_gc_dict['Intergenic'], width=0.2)
plt.bar(np.arange(0.6,len(region_list_gc_dict['mRNA'])+.6,1), region_list_gc_dict['Antisense'], width=0.2)
plt.xticks(range(len(region_list_gc_dict['mRNA'])),key_list)
plt.legend(['pc_middle_exons','intron_interior','no_overlap_set','antisense_transcript_set_gc'])
plt.title('Threshold: %d reads' % (threshold))
plt.xlabel('GC content range')
plt.tight_layout()
pdf_plots.savefig(fig)









exon_counts_list_dict = dict()

region_list_gc_dict = dict()
for lib_key in region_exon_id_sets_dict:
    
    exon_id_list = region_exon_id_sets_dict[lib_key]    
    exon_id_list = el.exon_id_intersection(exon_id_list,primary_3ss_exon_id_set)
    exon_id_list = above_threshold_exon_id_set(threshold,exon_id_list, aggregate_exon_dict)
    
    gc_result = el.get_exon_id_list_binned_GC(exon_id_list, gc_bins, genome_fasta)
    exon_counts_list_dict[lib_key]  = list()
    for key in gc_result:
        counts = el.get_exon_dict_counts(gc_result[key],aggregate_exon_dict)
        exon_counts_list_dict[lib_key].append(counts)






fig= plt.figure()
box1=plt.boxplot(exon_counts_list_dict['mRNA'], positions=np.arange(0,len(exon_counts_list_dict['mRNA'])), widths=0.2, patch_artist=True,whis=[10,90],showfliers=False)

for patch in box1['boxes']:
    patch.set_facecolor(color_dict['mRNA'])
    
box2=plt.boxplot(exon_counts_list_dict['Intronic'], positions=np.arange(0.25,len(exon_counts_list_dict['mRNA'])+.25), widths=0.2, patch_artist=True,whis=[10,90],showfliers=False)
for patch in box2['boxes']:
    patch.set_facecolor(color_dict['Intronic'])
    
box3=plt.boxplot(exon_counts_list_dict['Intergenic'], positions=np.arange(0.5,len(exon_counts_list_dict['mRNA'])+.5), widths=0.2, patch_artist=True,whis=[10,90],showfliers=False)
for patch in box3['boxes']:
    patch.set_facecolor(color_dict['Intergenic'])
    
box4=plt.boxplot(exon_counts_list_dict['Antisense'], positions=np.arange(0.75,len(exon_counts_list_dict['mRNA'])+.75), widths=0.2, patch_artist=True,whis=[10,90],showfliers=False)
for patch in box4['boxes']:
    patch.set_facecolor(color_dict['Antisense'])
    
plt.yscale('log')
plt.xticks(range(len(exon_counts_list_dict['mRNA'])),key_list,rotation=35)
plt.xlabel('gc percent range')
plt.ylabel('splicing inclusion')
plt.title('Threshold: %d reads' % (threshold))
#plt.legend(['pc middle', 'intron interior'])
plt.legend([box1["boxes"][0],box2["boxes"][0],box3["boxes"][0],box4["boxes"][0]], ['pc middle', 'intron interior','no_overlap_set','antisense_transcript'], loc='upper left')
plt.tight_layout()
pdf_plots.savefig(fig)



width_val= .15
fig= plt.figure()
box1=plt.boxplot(exon_counts_list_dict['mRNA'], positions=np.arange(0,len(exon_counts_list_dict['mRNA'])), widths=width_val, patch_artist=True,whis=[10,90],showfliers=False)
for patch in box1['boxes']:
    patch.set_facecolor(color_dict['mRNA'])
    
box2=plt.boxplot(exon_counts_list_dict['lncRNA'], positions=np.arange(width_val*1.25,len(exon_counts_list_dict['mRNA'])+width_val*1.25), widths=width_val, patch_artist=True,whis=[10,90],showfliers=False)
for patch in box2['boxes']:
    patch.set_facecolor(color_dict['lncRNA'])
    
box3=plt.boxplot(exon_counts_list_dict['Intronic'], positions=np.arange(width_val*2.5,len(exon_counts_list_dict['mRNA'])+width_val*2.5), widths=width_val, patch_artist=True,whis=[10,90],showfliers=False)
for patch in box3['boxes']:
    patch.set_facecolor(color_dict['Intronic'])
    
box4=plt.boxplot(exon_counts_list_dict['Intergenic'], positions=np.arange(width_val*3.75,len(exon_counts_list_dict['mRNA'])+width_val*3.75), widths=width_val, patch_artist=True,whis=[10,90],showfliers=False)
for patch in box4['boxes']:
    patch.set_facecolor(color_dict['Intergenic'])
    
box5=plt.boxplot(exon_counts_list_dict['Antisense'], positions=np.arange(width_val*5,len(exon_counts_list_dict['mRNA'])+width_val*5), widths=width_val, patch_artist=True,whis=[10,90],showfliers=False)
for patch in box5['boxes']:
    patch.set_facecolor(color_dict['Antisense'])
    
plt.yscale('log')
plt.xticks(range(len(exon_counts_list_dict['mRNA'])),key_list,rotation=35)
plt.xlabel('gc percent range')
plt.ylabel('splicing inclusion')
plt.title('Threshold: %d reads' % (threshold))

plt.legend([box1["boxes"][0],box2["boxes"][0],box3["boxes"][0],box4["boxes"][0],box5["boxes"][0]], ['pc middle', 'lncRNA middle', 'intron interior','no_overlap_set','antisense_transcript'], loc='upper left')
plt.tight_layout()
pdf_plots.savefig(fig)























#make supplemental 2D
with open(exp_output_path.out_supplemental+"2D.txt", 'w') as f:
    
    outstring = 'chrom\tstart\tend\tstrand\texon_id\tgenomic_category\tgc_category\tcount\n'
    f.write(outstring)
    for key in region_exon_id_sets_dict:
        exon_id_list = region_exon_id_sets_dict[key]
        exon_id_list = el.exon_id_intersection(exon_id_list, primary_3ss_exon_id_set)
        
        gc_result = el.get_exon_id_list_binned_GC(exon_id_list, gc_bins, genome_fasta)
        
        for gc_key in gc_result:
            for exon_id in gc_result[gc_key]:
                ex=el.exon_id_values(exon_id)
                
                outstring = '{:}\t{:}\t{:}\t{:}\t{:}\t'.format(ex.chrom,ex.start,ex.end,ex.strand,exon_id)
                outstring = '{:}\t{:}\t{:}\t{:}\n'.format(outstring, key, gc_key, aggregate_exon_dict[exon_id]['count'])
                f.write(outstring)
                
    
    
    
    



















pdf_plots.close()









#/home/pdf/repositories/exon_def/Process_ET/calculate_chromosome_coverage.py

pickle_path = exp_output_path.pickle_merged + "calculate_chromosome_coverage__%d.pickle" % (exon_count_build)

with open(pickle_path, "wb") as output_file:
    pickle.dump(no_overlap_set, output_file)
    pickle.dump(intron_interior_set, output_file)
    pickle.dump(intron_interior_pc_set, output_file)
    pickle.dump(intron_interior_lncRNA_set, output_file)
    pickle.dump(antisense_transcript_set, output_file)
    
    pickle.dump(exon_overlap_set, output_file)
    pickle.dump(sense_transcript_set, output_file)
    pickle.dump(pc_middle_exon_id_list, output_file)
    pickle.dump(lncRNA_middle_exon_ids, output_file)
    pickle.dump(region_exon_id_sets_list, output_file)
    pickle.dump(region_exon_id_sets_dict, output_file)
    
    pickle.dump(lncRNA_overlap_set, output_file)
    pickle.dump(pc_antisense_transcript_exon_overlap, output_file)
    pickle.dump(lncRNA_antisense_transcript_exon_overlap, output_file)
    
    
    

    
    print('calculate_chromosome_coverage pickle saved to %s' % (pickle_path))
    
    








'''
exon_id_list = no_overlap_set
exon_id_list_name = 'intergenic_exon_ids'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)



exon_id_list = antisense_transcript_set
exon_id_list_name = 'antisense_exon_ids'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)



exon_id_list = intron_interior_set
exon_id_list_name = 'intronic_exon_ids'
out_file_path = ('{:}/{:}_count_{:}.bed').format(exp_output_path.out_bed_main, exon_id_list_name, exon_count_build)
el.export_exon_id_list_to_bed(exon_id_list, out_file_path)

'''















pickle_path = exp_output_path.pickle_merged + "calculate_chromosome_coverage__%d.pickle" % (exon_count_build)

'''
with open(pickle_path, "rb") as input_file:
    no_overlap_set = pickle.load(input_file)
    intron_interior_set = pickle.load(input_file)
    intron_interior_pc_set = pickle.load(input_file)
    intron_interior_lncRNA_set = pickle.load(input_file)
    antisense_transcript_set = pickle.load(input_file)
    exon_overlap_set = pickle.load(input_file)
    sense_transcript_set = pickle.load(input_file)
    pc_middle_exon_id_list = pickle.load(input_file)
    lncRNA_middle_exon_ids = pickle.load(input_file)
    #pickle.dump(aggregate_exon_dict, output_file)
    region_exon_id_sets_list = pickle.load(input_file)
    region_exon_id_sets_dict = pickle.load(input_file)
    
    
    print('calculate_chromosome_coverage pickle saved to %s' % (pickle_path))
    
    '''
    















