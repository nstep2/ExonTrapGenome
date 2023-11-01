#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 04:15:45 2022

@author: pdf
"""

#Dfam_volcano_calcs.py





avg_exon_exp_bases = rep_genome_stats_df['avg_exon_exp_bases']
avg_exon_exp = rep_genome_stats_df['avg_exon_exp']

e_genome_bases = rep_genome_stats_df['exon_bases'] 
r_genome_bases = rep_genome_stats_df['genome_bases']
genome_len=3088286401
for key in repeat_stat_rows_dict:
    entry = repeat_stat_rows_dict[key]
    
    r_repeat_bases = entry['rep_genome_bases']
    e_repeat_bases = entry['exon_rep_intersect_bases']
    e_exon_bases   = entry['exon_genome_bases']
    
    
    
    
    #over representation
    test_hyper_parial_sf = scipy.stats.hypergeom(M=genome_len, n=(e_genome_bases/2), N=r_repeat_bases).logsf(e_repeat_bases-1) 
    
    #under representation
    test_hyper_partial_cdf = scipy.stats.hypergeom(M=genome_len, n=(e_genome_bases/2), N=r_repeat_bases).logcdf(e_repeat_bases) #should not have -1?
    
    
    test_hyper_3_ratio= (e_genome_bases/r_repeat_bases)/((e_genome_bases/2)/(genome_len))    



    #count all exon bases, not just overlap
    test_hyper_all_sf = scipy.stats.hypergeom(M=genome_len, n=(e_genome_bases/2), N=r_repeat_bases).logsf(e_exon_bases-1) 
    
    #count all exon bases, not just overlap
    #under representation
    test_hyper_partial_cdf = scipy.stats.hypergeom(M=genome_len, n=(e_genome_bases/2), N=r_repeat_bases).logcdf(e_exon_bases)


















