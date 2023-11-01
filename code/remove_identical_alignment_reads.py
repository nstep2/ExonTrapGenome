#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 14:00:01 2021

@author: pdf
"""
import sys

if len(sys.argv) > 1:
    #fasta_file_location = experiment_location()
    input_file = sys.argv[1]
    #mode = sys.argv[1]
else:
    print('missing input arg. Rerun with file path')
    
import gzip, os, subprocess

with gzip.open(input_file,'rt') as f:
    
    #tmp_path = os.path.basepath(input_file)
    tmp_name = input_file[:-7] + '_collapsed' + '.txt.gz'

    selected_alignments_name_handle = gzip.open( tmp_name, 'wt',compresslevel=3) 
    
    
    
    line_1 = ''
    line_2 = ''
    
    line_1_count = 0
    line_2_count = 0
    
    line_1_stuff = ['','','']
    line_2_stuff = ['','','']
    
    for ii, line in enumerate(f):
        if line[0] == '@':
            selected_alignments_name_handle.write(line)
            continue
        
        else:
            #
            line_2 = line
            line_2_split = line_2.split('\t')
            line_2_stuff = [line_2_split[2],line_2_split[3],line_2_split[7]]
            
            #check if line 1 has not been initialized
            if line_1 == '':
                line_1 = line
                line_1_split = line_1.split('\t')
                line_1_count = 0
                line_1_stuff = [line_1_split[2],line_1_split[3],line_1_split[7]]
        
        
        
        # check if read alignments are identical
        if line_1_stuff[0] == line_2_stuff[0] and line_1_stuff[1] == line_2_stuff[1] and line_1_stuff[2] == line_2_stuff[2]:
            line_1_count += 1
            
        else:
            #write out the old line 1
            out_string = line_1.strip()
            out_string = out_string + '\t%d\n' % (line_1_count)
            selected_alignments_name_handle.write(out_string)
            #now reset line_1 to line_2
            line_1 = line_2
            line_1_count = 1
            line_1_stuff = line_2_stuff
    
    
    #write the final read
    out_string = line_1.strip()
    out_string = out_string + '\t%d\n' % (line_1_count)
    selected_alignments_name_handle.write(out_string)
    
    
    #close the read handle
    selected_alignments_name_handle.close()
        
    
    
print('finished collapsing %s' % (input_file))   
    
    
    
    
    
    
    
    











