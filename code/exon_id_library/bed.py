#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 15:28:09 2021

@author: pdf
"""

import exon_id_library.exon_id_library as el


class simple_bed():
    def __init__(self, bed_file, aggregate_exon_IT):
        
        
        #chromsomes_list = ['chrX',  'chrY',  'chrMT', 'chrM']
        chromsomes_list = ['chrX',  'chrY',  'chrMT']
        for i in range(1,23):
            chromsomes_list.append( "chr%d" % (i) )
        chromsomes_list = set(chromsomes_list)
        
        strand_set = set()
        strand_set.add('+')
        strand_set.add('-')
        
        self.middle_exon_id_list = list()
        self.first_last_exon_id_list = list()
        self.name_to_exon_id = dict()
        self.exon_id_to_name = dict()
        
        self.fuzzy_overlap_middle_exon_ids = list()
        
        
        with open(bed_file,'r') as f:
            for line in f:
                line_split = line.strip().split('\t')
                
                number_exons = int(line_split[9])
                
                if number_exons < 2:
                    #print(line)
                    continue
                
                exon_len_list = line_split[10].split(',')
                exon_pos_list = line_split[11].split(',')
                
                exon_len_list = [int(x) for x in exon_len_list[:-1]]
                exon_pos_list = [int(x) for x in exon_pos_list[:-1]]
                
                chrom = line_split[0]
                if chrom not in chromsomes_list:
                    continue
                
                strand = line_split[5]
                if strand not in strand_set:
                    continue
                
                genome_start = int(line_split[1])
                
                name = line_split[3]                
                self.name_to_exon_id[name] = list()
                #self.exon_id_to_name[exon_id] = list()
                for ii, exon_len in enumerate(exon_len_list):
                    exon_pos = exon_pos_list[ii]
                    
                    start = genome_start + exon_pos + 1
                    end = start + exon_len 
                    
                    if end - start < 0:
                        continue
                    
                    
                    exon_id = '%s:%d-%d:%s' % (chrom, start,end, strand)
                    
                    self.name_to_exon_id[name].append(exon_id)
                    self.exon_id_to_name[exon_id] = name
                    
                    if ii == 0 or ii == len(exon_len_list) - 1:
                        self.first_last_exon_id_list.append(exon_id)
                    else:
                        self.middle_exon_id_list.append(exon_id)
                    
                    
                        
                        intervals = aggregate_exon_IT[chrom][strand][start:end]
                        if len(intervals) > 0:
                            #exons_with_interval += 1 #.append(len(intervals))
        
                            for interval in intervals:
                                for overlap_id in interval[2]:
                                   overlap, len_A, len_B = el.get_overlap_exon_A_with_B(exon_id, overlap_id)
                                   if abs(overlap-len_A) <= 3 and abs(overlap-len_B) <= 3:
                                       self.fuzzy_overlap_middle_exon_ids.append(overlap_id)
                                       #int('nnnnonnnooooooooo')
       
        self.middle_exon_id_list = list(set(self.middle_exon_id_list))
        self.first_last_exon_id_list = list(set(self.first_last_exon_id_list))
        self.fuzzy_overlap_middle_exon_ids = list(set(self.fuzzy_overlap_middle_exon_ids))
        
           
    
 
