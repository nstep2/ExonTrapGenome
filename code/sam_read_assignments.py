



#Test data
#input_file = '/media/8TB_1/ETF_data/20200202/reads/trimmed/4_ETF_cleaned_S11_L001_R1_001.fastq.gz_all_trimmed_simple_paired.sam.gz'
#backbone_type = 'T3'



from termcolor import colored
#print(colored('Hello, World!', 'red'))
#print(colored('Hello, World!', 'cyan', on_color='on_red'))


from experiment_paths.experiment_paths import *

exp_output_path.ssd_tmp 


import pyfaidx
import datetime
import sys, time, subprocess, os

if len(sys.argv) > 1:
    
    input_file = sys.argv[2]
    print(colored('\n\n' + sys.argv[2] + '\n\n', 'red'))
    backbone_type = sys.argv[1]
else:
    print('missing input arg. Rerun with file path')
    
import time



total_start_time = time.time()

start_time=total_start_time

datetime.datetime.now().strftime('%H:%M:%S')
print('selecting')
input_file_path = input_file
extract_sam_command ='python ' + exp_output_path.git_repo +'make_sam_hits_data_file.py %s %s %s %s' % (input_file_path,'dummy',0,backbone_type)
output_text = subprocess.Popen(extract_sam_command, shell=True).wait()
print(output_text)
extract_output_file = input_file_path[:-7]+'_selected.txt.gz'

current_time = time.time()
start_time=current_time
print('section time elapsed: %d minutes' % ((current_time-start_time)/60))
print('total elapsed time 1: %d minutes (%s)' % ((current_time-total_start_time)/60,datetime.datetime.now().strftime('%H:%M:%S')))
start_time=current_time

'''
input_path = input_file[:-7]+'_selected_alignments.txt.gz'
sort (linux, < 1 hour)
'_sorted.txt.gz'
'''

print('sorting')
input_file_path=extract_output_file
output_path = input_file_path[:-7]+'_sorted.txt.gz'
#sort_command_1 = "GZIP_OPT=-3 LC_ALL=C gunzip -c %s | sort -T /media/2TB_Samsung_1/ssd_scratch --compress-program=gzip  -k3,3 -k4,4n -k8,8n  - | gzip -f -3 > %s" %(input_file_path, output_path)
sort_command_1 = "GZIP_OPT=-3 LC_ALL=C gunzip -c %s | sort -T %s --compress-program=gzip  -k3,3 -k4,4n -k8,8n  - | gzip -f -3 > %s" %(input_file_path, exp_output_path.ssd_tmp, output_path)

subprocess.Popen(sort_command_1, shell=True).wait()
sorted_output_path_1=output_path

current_time = time.time()
print('section time elapsed: %d minutes' % ((current_time-start_time)/60))
print('total elapsed time 2: %d minutes (%s)' % ((current_time-total_start_time)/60,datetime.datetime.now().strftime('%H:%M:%S')))
start_time=current_time



print('collapsing')
input_file_path = sorted_output_path_1

collapse_command = 'python ' + exp_output_path.git_repo +'remove_identical_alignment_reads.py %s' % (input_file_path)
output_text = subprocess.Popen(collapse_command, shell=True).wait()
print(output_text)
collape_output_path_1 = input_file_path[:-7]+'_collapsed.txt.gz'


current_time = time.time()
print('section time elapsed: %d minutes' % ((current_time-start_time)/60))
print('total elapsed time 3: %d minutes (%s)' % ((current_time-total_start_time)/60,datetime.datetime.now().strftime('%H:%M:%S')))
start_time=current_time







import pyfaidx
input_fasta = exp_output_path.hg38_genome_fasta 
genome_fasta = pyfaidx.Fasta(input_fasta)


from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
matrix5 = load_matrix5()
matrix3 = load_matrix3()


import exon_id_library.exon_id_lib as el

import sys
import re
import gzip
import subprocess










if backbone_type == 'T3':
    darkcycle_offset = 4
if backbone_type == 'T4':
    darkcycle_offset = 3
if backbone_type == 'T5':
    darkcycle_offset = 3
    
    
print('reshaping')
input_file_path = collape_output_path_1
input_file_handle = gzip.open(input_file_path,'rt')

output_file = collape_output_path_1[:-7] + '_reshape.txt.gz'

output_file_handle = gzip.open(output_file,'wt',compresslevel=3)
    


check_window_len = 3


        
count_lagging_soft_clipped_reads = 0
count_leading_soft_clipped_reads = 0
count_complicated_cigar_string   = 0

count_altered_5ss = 0
count_altered_3ss = 0

last_good_line = ''
for ii, line in enumerate(input_file_handle):
    #if ii == 50:
    #    break
    
    if line[0] == '@':
        continue
    
    
        
        
    
    
    
    split_line = line.split('\t')
    
    
    #print('line', line)
    #int('a')
    
    
    cigar_string = split_line[5]
    cigar_split = re.findall(r'[A-Za-z]+|\d+', cigar_string)
    
    
    if len(cigar_string) > 1 and (cigar_split[1] == 'S' or cigar_split[-1] == 'S') :
        #if int(cigar_split[0]) > 0 or True:
        if cigar_split[1] == 'S':
            count_leading_soft_clipped_reads += 1
        if cigar_split[-1] == 'S':
            count_lagging_soft_clipped_reads += 1
        continue
    
    if len(cigar_split) >= 6:
        
        found_N = False
        for cigar_ii in range(3,len(cigar_split), 2):
            if cigar_split[cigar_ii] == 'N':
                found_N = True
        
        if found_N == True:
            count_complicated_cigar_string += 1
            
            continue
    
    
    cigar_string = split_line[-3]
    cigar_split = re.findall(r'[A-Za-z]+|\d+', cigar_string)
    
    
    if len(cigar_string) > 1 and (cigar_split[1] == 'S' or cigar_split[-1] == 'S') :
        if int(cigar_split[0]) > 0 or True:
            if cigar_split[1] == 'S':
                count_leading_soft_clipped_reads += 1
            if cigar_split[-1] == 'S':
                count_lagging_soft_clipped_reads += 1
                
                
            continue
    
    
    if len(cigar_split) >= 6:
        
        found_N = False
        for cigar_ii in range(3,len(cigar_split), 2):
            if cigar_split[cigar_ii] == 'N':
                found_N = True
        
        #if cigar_split[3] == 'N' or cigar_split[5] == 'N':
        if found_N == True:
            
            count_complicated_cigar_string += 1
            continue
    

    

    ### determine start/end and strand
    
    if int(split_line[1]) & 64 == 64:
        mate_read_1_flag = True
    else:
        mate_read_1_flag = False
    
    if mate_read_1_flag == True and abs(int(split_line[8]) ) < 1000 and int(split_line[8]) != 0 and split_line[2] != '*':
        
        last_good_line = line
        
        #adjust for the backbone read dark cycled bases
        if (int(split_line[8])) > 0: # read orientation positive or neg strand
    
            start_position_estimate = int(split_line[3])
            end_position_estimate   = int(split_line[3]) + int(split_line[8])
        else:   
            
            start_position_estimate = int(split_line[7])
            end_position_estimate   = int(split_line[7]) + abs(int(split_line[8]))
        
        
        
        if int(split_line[1]) & 16 != 16:
            strand_flag = '-'
            strand = '-'
            start_position_estimate = start_position_estimate - darkcycle_offset
            
        else:
            strand_flag = '+'
            strand = '+'
            end_position_estimate = end_position_estimate + darkcycle_offset
            
    
    
        
        
        exon_id = "%s:%d-%d:%s" % (split_line[2], start_position_estimate, end_position_estimate, strand)
        
        
        ex = el.exon_id_values(exon_id)
        best_5ss_pos = ex.pos_5ss
        best_3ss_pos = ex.pos_3ss
        
        improved_score_5ss = False
        if ex.strand == '+':
            up_5ss   = ex.pos_5ss-4
            down_5ss = ex.pos_5ss+5
            
            seq_5ss = str(genome_fasta[ex.chrom][up_5ss:down_5ss].seq)
            if seq_5ss.find('N') >= 0 or seq_5ss.find('n') >= 0:
                score_5ss = -5000
            else:
                try:
                    score_5ss = maxent.score5(seq_5ss, matrix=matrix5)
                except:
                    int('error 1')
        
        
        
            up_5ss   = ex.pos_5ss-4
            down_5ss = ex.pos_5ss+5
            
            
            seq_5ss = str(genome_fasta[ex.chrom][up_5ss:down_5ss])
            if seq_5ss.find('N') >= 0 or seq_5ss.find('n') >= 0:
                score_5ss = -5000
            else:
                try:
                    score_5ss = maxent.score5(seq_5ss, matrix=matrix5)
                except:
                    int('error 2')
    
            
            best_5ss_pos = ex.pos_5ss
            best_5ss_score = score_5ss

            if score_5ss > -5000:
                
                for check_5ss_pos in range(ex.pos_5ss-check_window_len, ex.pos_5ss + check_window_len+1):
                    
                    check_seq_5ss = str(genome_fasta[ex.chrom][check_5ss_pos-4:check_5ss_pos+5])
                    
                    if check_seq_5ss.find('N') >= 0 or check_seq_5ss.find('n') >= 0 or len(check_seq_5ss) != 9:
                        continue
                    
                    try:
                        check_score_5ss =  maxent.score5(check_seq_5ss, matrix=matrix5)
                    except:
                        int('error 5')
            
                    if check_score_5ss > best_5ss_score:
                        best_5ss_pos = check_5ss_pos
                        best_5ss_score = check_score_5ss
                if ex.pos_5ss != best_5ss_pos:
                    count_altered_5ss += 1
                    improved_score_5ss = True
            
        
        
        
        if ex.strand == '-':
            up_5ss   = ex.pos_5ss-7
            down_5ss = ex.pos_5ss+2
            
            
            seq_5ss = str(genome_fasta[ex.chrom][up_5ss:down_5ss].complement.reverse)
            if seq_5ss.find('N') >= 0 or seq_5ss.find('n') >= 0:
                score_5ss = -5000
            else:
                try:
                    score_5ss = maxent.score5(seq_5ss, matrix=matrix5)
                except:
                    int('error 3')
        
            
            best_5ss_pos = ex.pos_5ss
            best_5ss_score = score_5ss
            
            
            if score_5ss > -5000:
                    
                for check_5ss_pos in range(ex.pos_5ss-check_window_len, ex.pos_5ss + check_window_len+1):
                    
                    
                    try:  #to skip a chromosome end issue
                        check_seq_5ss = str(genome_fasta[ex.chrom][check_5ss_pos-7:check_5ss_pos+2].complement.reverse)
                    except:
                        print('ex.chrom:', ex.chrom)
                        print('check_5ss_pos-7:', check_5ss_pos-7)
                        print('check_5ss_pos+2:', check_5ss_pos+2)
                        #int('')
                        if ex.chrom != 'chrM':
                            check_seq_5ss = str(genome_fasta[ex.chrom][check_5ss_pos-7:check_5ss_pos+2].complement.reverse)
                        else:
                            continue  #I'm ok skipping a circular overlap M chromosome
                            
                    
                    check_seq_5ss = str(genome_fasta[ex.chrom][check_5ss_pos-7:check_5ss_pos+2].complement.reverse)
                    
                    if check_seq_5ss.find('N') >= 0 or check_seq_5ss.find('n') >= 0 or len(check_seq_5ss) != 9:
                        continue
                    
                    try:
                        check_score_5ss =  maxent.score5(check_seq_5ss, matrix=matrix5)
                    except:
                        int('error 6')
            
                    if check_score_5ss > best_5ss_score:
                        best_5ss_pos = check_5ss_pos
                        best_5ss_score = check_score_5ss
                if ex.pos_5ss != best_5ss_pos:
                    count_altered_5ss += 1
                    improved_score_5ss = True
            
            

        
        # 3ss
        improved_score_3ss = False
        if ex.strand == '+':
            up_3ss   = ex.pos_3ss-21
            down_3ss = ex.pos_3ss+2
            
            
            try:
                seq_3ss = str(genome_fasta[ex.chrom][up_3ss:down_3ss].seq)
            except:
                print('ex.chrom:', ex.chrom)
                print('up_3ss:', up_3ss)
                print('down_3ss:', down_3ss)
                #int('')
                if ex.chrom != 'chrM':
                    seq_3ss = str(genome_fasta[ex.chrom][up_3ss:down_3ss].seq)
                else:
                    continue  #I'm ok skipping a circular overlap M chromosome
                
            seq_3ss = str(genome_fasta[ex.chrom][up_3ss:down_3ss].seq)
            if seq_3ss.find('N') >= 0 or seq_3ss.find('n') >= 0 or len(seq_3ss) != 23:
                score_3ss = -5000
            else:
                try:
                    score_3ss = maxent.score3(seq_3ss, matrix=matrix3)
                except:
                    print('seq_3ss:', seq_3ss)
                    print(exon_id, up_3ss, down_3ss)
                    int('error 7')
        
            
            best_3ss_pos = ex.pos_3ss
            best_3ss_score = score_3ss
            
            if score_3ss > -5000:
                
                for check_3ss_pos in range(ex.pos_3ss-check_window_len, ex.pos_3ss + check_window_len+1):
                    if best_3ss_score == -5000:
                        break
                        
                        
                    check_seq_3ss = str(genome_fasta[ex.chrom][check_3ss_pos-21:check_3ss_pos+2])
                    
                    
                    if check_seq_3ss.find('N') >= 0 or check_seq_3ss.find('n') >= 0 or len(check_seq_3ss) != 23:
                        continue
                    
                    try:
                        check_score_3ss =  maxent.score3(check_seq_3ss, matrix=matrix3)
                    except:
                        int('error 8')
            
            
                    #print(check_seq_3ss, "%.2f" %check_score_3ss)
                    if check_score_3ss >= best_3ss_score:
                        best_3ss_pos = check_3ss_pos
                        best_3ss_score = check_score_3ss
                        #print(ii,best_3ss_pos)
                if ex.pos_3ss != best_3ss_pos:
                    count_altered_3ss += 1
                    improved_score_3ss = True
            
            
            
            
        
        
        #3ss
        
        
        if ex.strand == '-':
            up_3ss   = ex.pos_3ss-4
            down_3ss = ex.pos_3ss+19
            
            
            seq_3ss = str(genome_fasta[ex.chrom][up_3ss:down_3ss].complement.reverse)
            if seq_3ss.find('N') >= 0 or seq_3ss.find('n') >= 0:
                score_3ss = -5000
                
            else:
                try:
                    score_3ss = maxent.score3(seq_3ss, matrix=matrix3)
                except:
                    int('error 4')
        
            
            #break
            best_3ss_pos = ex.pos_3ss
            best_3ss_score = score_3ss
            
            if score_3ss > -5000:
                
                for check_3ss_pos in range(ex.pos_3ss-check_window_len, ex.pos_3ss + check_window_len+1):
                    
                    check_seq_3ss = str(genome_fasta[ex.chrom][check_3ss_pos-4:check_3ss_pos+19].complement.reverse)
                    
                    if check_seq_3ss.find('N') >= 0 or check_seq_3ss.find('n') >= 0 or len(check_seq_3ss) != 23:
                        continue
                    try:
                        check_score_3ss =  maxent.score3(check_seq_3ss, matrix=matrix3)
                    except:
                        int('error 9')
            
                    if check_score_3ss > best_3ss_score:
                        best_3ss_pos = check_3ss_pos
                        best_3ss_score = check_score_3ss
                if ex.pos_3ss != best_3ss_pos:
                    count_altered_3ss += 1
                    improved_score_3ss = True
                
                
                
                
        
        new_start_position_estimate = start_position_estimate
        new_end_position_estimate = end_position_estimate
        if improved_score_5ss == True:
            
            if strand =='+':
                new_end_position_estimate = best_5ss_pos
            else:
                new_start_position_estimate = best_5ss_pos
                    
                
        if improved_score_3ss == True:
            
            if strand =='+':
                new_start_position_estimate = best_3ss_pos
            else:
                new_end_position_estimate = best_3ss_pos




        
        new_exon_id = "%s:%d-%d:%s" % (split_line[2], new_start_position_estimate, new_end_position_estimate, strand)
        
        #outline = ''
        outline = '%s' %(split_line[0])
        outline = '%s\t%s' %(outline,exon_id)
        outline = '%s\t%s' %(outline,new_exon_id)
        outline = '%s\t%s' %(outline,ex.chrom)
        outline = '%s\t%d' %(outline,new_start_position_estimate)
        outline = '%s\t%d' %(outline,new_end_position_estimate)
        outline = '%s\t%s' %(outline,strand)
        outline = '%s\t%s' %(outline,split_line[-1].strip())
        outline = '%s\t%s' %(outline,split_line[5])
        outline = '%s\t%s\n' %(outline,split_line[-3])
        #outline = '%s\t%' %(outline,)
        
        output_file_handle.write(outline)
     
        
        
        
        
        
        
        
        
        
input_file_handle.close()
output_file_handle.close()


current_time = time.time()
print('section time elapsed: %d minutes' % ((current_time-start_time)/60))
print('total elapsed time 4: %d minutes (%s)' % ((current_time-total_start_time)/60,datetime.datetime.now().strftime('%H:%M:%S')))
start_time=current_time





print('window %d has %0.2f%% altered 5ss pos' %(check_window_len, count_altered_5ss/ii*100))


print('window %d has %0.2f%% altered 3ss pos' %(check_window_len, count_altered_3ss/ii*100))



print('sorting')
in_name = output_file
sorted_name = in_name[:-7]+'_sorted.txt.gz'
#sort_command_string = "GZIP_OPT=-3 LC_ALL=C gunzip -c %s | sort -T /media/2TB_Samsung_1/ssd_scratch --compress-program=gzip  -k4,4 -k5,5n -k6,6n -k8,8n - | gzip -f -3 > %s" % (in_name,sorted_name)
sort_command_string = "GZIP_OPT=-3 LC_ALL=C gunzip -c %s | sort -T %s --compress-program=gzip  -k4,4 -k5,5n -k6,6n -k8,8n - | gzip -f -3 > %s" % (in_name,exp_output_path.ssd_tmp,sorted_name)

print("BEGIN command: %s\n" % (sort_command_string))
subprocess.Popen(sort_command_string, shell=True).wait()
#print("\nFINISH command: %s\n" % (sort_command_string))



current_time = time.time()
print('section time elapsed: %d minutes' % ((current_time-start_time)/60))
print('total elapsed time 5: %d minutes (%s)' % ((current_time-total_start_time)/60,datetime.datetime.now().strftime('%H:%M:%S')))
start_time=current_time









print('collapsing')
input_file = sorted_name


with gzip.open(input_file,'rt') as f:
    
    tmp_name = input_file[:-7] + '_collapse' + '.txt.gz'

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
            line_2_stuff = [line_2_split[3],line_2_split[4],line_2_split[5],line_2_split[6]]
            
            #check if line 1 has not been initialized
            if line_1 == '':
                line_1 = line
                line_1_split = line_1.split('\t')
                line_1_count = 0
                line_1_exon_id_list = ''
                line_1_stuff = [line_1_split[3],line_1_split[4],line_1_split[5],line_1_split[6]]
        
        
        
        # check if read alignments are identical
        if line_1_stuff[0] == line_2_stuff[0] and line_1_stuff[1] == line_2_stuff[1] and line_1_stuff[2] == line_2_stuff[2] and line_1_stuff[3] == line_2_stuff[3]:
            line_1_count += int(line_2_split[7])
            line_1_exon_id_list = line_1_exon_id_list+'%s_%s;'%(line_2_split[2],line_2_split[7])
            
        else:
            #write out the old line 1
            #out_string = line_1.strip()
            
            #out_string = out_string + '\t%d\n' % (line_1_count)
            #make some int values into str
            for ii, val in enumerate(line_1_stuff):
                line_1_stuff[ii] = str(line_1_stuff[ii])
            out_string = line_2_split[3]  #this can be the wrong chromosme. Also, the chromosome is in line_1_stuff. I think I meant this to be an exon_id of some sort.. not sure
            out_string = "%s\t%s" % (out_string, '\t'.join(line_1_stuff))
            out_string = "%s\t%d" % (out_string,line_1_count)
            out_string = "%s\t%s\n" % (out_string, line_1_exon_id_list)
            selected_alignments_name_handle.write(out_string)
            #now reset line_1 to line_2
            line_1 = line_2
            line_1_count = int(line_2_split[7]) 
            line_1_stuff = line_2_stuff
            line_1_exon_id_list = '%s_%s;'%(line_2_split[2],line_2_split[7])
    
    
    #write the final read
    for ii, val in enumerate(line_1_stuff):
        line_1_stuff[ii] = str(line_1_stuff[ii])
    out_string = line_2_split[3]
    out_string = "%s\t%s" % (out_string, '\t'.join(line_1_stuff))
    out_string = "%s\t%d" % (out_string,line_1_count)
    out_string = "%s\t%s\n" % (out_string, line_1_exon_id_list)
    selected_alignments_name_handle.write(out_string)
    
    #close the read handle
    selected_alignments_name_handle.close()
        
    
    
print('finished collapsing %s' % (input_file))   
    
    
    



current_time = time.time()
print('section time elapsed: %d minutes' % ((current_time-start_time)/60))
print('total elapsed time 6: %d minutes (%s)' % ((current_time-total_start_time)/60,datetime.datetime.now().strftime('%H:%M:%S')))
start_time=current_time




print(colored('\n\n ***Finished***: \n%s  \n\n' % (os.path.basename(input_file)), 'green')) 






