#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 21:41:36 2021

@author: pdf
"""


chromsomes_list = ['chrX',  'chrY',  'chrMT', 'chrM']
for i in range(1,23):
    chromsomes_list.append( "chr%d" % (i) )




def get_annotated_exons_from_GENCODE_GFF3(input_file, version, genome_fasta ):
    
    exon_data_dict = dict()
    exon_data_dict_25_to_400 = dict()
    exon_data_dict_50_to_200 = dict()
    
    #gene_to_transcript_dict = dict()
    transcript_to_gene_dict = dict()
    
    transcript_list = list()
    transcript_id_count = dict()
    
    seq_to_exon_region_dict = dict()
    duplicate_exon_seq_to_exon_region_dict = dict()
    
    exon_id_sets_dict = dict()
    exon_id_transcript_sets_dict = dict()
    UTR_annotation = {'5_utr':dict(), '3_utr':dict()}
    
    GENCODE_exon_list=list()
    GENCODE_exon_dict=dict()
    
    GENCODE_gene_dict=dict()
    
    with open(input_file, "r") as handle:
        for line in handle:
            if line[0] == '#':
                continue
            line_split = line.split('\t')
            
            
            if line_split[2] == 'five_prime_UTR' or line_split[2] == 'three_prime_UTR':
                chrom = line_split[0]
                start = int(line_split[3])
                end   = int(line_split[4])
                strand = line_split[6]
                
                tmp = line_split[8].split('Parent=')[1]
                transcript_id = tmp.split(';')[0]
                
                custom_chrom_region = "%s:%d-%d:%s" % (chrom, start, end+1, strand)
                if line_split[2] == 'five_prime_UTR':
                    UTR_annotation['5_utr'][custom_chrom_region] = transcript_id
                if line_split[2] == 'three_prime_UTR':
                    UTR_annotation['3_utr'][custom_chrom_region] = transcript_id
                    
            if version == 'basic':
                tags = line.split('tag=')
                if len(tags) > 1:
                    tags = tags[1].split(';')[0]
                    if 'basic' not in tags.split(','):
                    #if tag == 'basic':
                        continue
                else:
                    continue  #since there is no tag to indicate that it is basic

            if line_split[2] == 'gene':
                gene_id=line_split[8].split('gene_id=')[1].split(';')[0]
                
                if gene_id in GENCODE_gene_dict:
                    gene_count = 1
                    while("%s_%d"%(gene_id, gene_count) in GENCODE_gene_dict):
                        print(gene_id)
                        gene_count += 1
                    gene_id = "%s_%d"%(gene_id, gene_count)
                  
                    
                GENCODE_gene_dict[gene_id] = dict()
                GENCODE_gene_dict[gene_id]['start'] = int(line_split[3])
                GENCODE_gene_dict[gene_id]['end'] = int(line_split[4])
                GENCODE_gene_dict[gene_id]['strand'] = line_split[6]
                GENCODE_gene_dict[gene_id]['chrom'] = line_split[0]
                
                
                GENCODE_gene_dict[gene_id]['tags'] = dict()
                GENCODE_gene_dict[gene_id]['tags']['keyword'] = list()
                tags = line_split[8].split('tag=')[0].split(';')
                
                for tag in tags:
                    tag_split=tag.split('=')
                    if len(tag_split) == 2:
                        GENCODE_gene_dict[gene_id]['tags'][tag_split[0]] = tag_split[1]
                    else:
                        GENCODE_gene_dict[gene_id]['tags']['keyword'].append(tag_split[0])
                        
                
                
                
                
                #GENCODE_gene_dict[gene_id]['tags']
                

                    

            if line_split[2] != 'exon':
                continue
            
            if line_split[8].find('transcript_support_level') >= 0:
                
                TSL = line_split[8].split('transcript_support_level=')[1].split(';')[0]
                if TSL == 'NA': #NA means it was not analyzed, b/c pseudogene, HLA, immunoglobin, T-cell, single exon
                    TSL = 6
                TSL = int(TSL)
                if TSL > 6:     #???
                    continue
                
            #Level cutoff
            if line_split[8].find(';level=') >= 0:
                
                level = line_split[8].split(';level=')[1].split(';')[0]
                if level == 'NA': #NA means it was not analyzed, b/c pseudogene, HLA, immunoglobin, T-cell, single exon
                    level = 6
                level = int(level)
                if level > 3:
                    continue
            
            
            
            
            
            GENCODE_exon = GENCODE_exon_class()
            GENCODE_exon.build_from_gff3_entry(line)
            GENCODE_exon_dict[GENCODE_exon.gencode_exon_id]=GENCODE_exon
            
            
            
            chrom = line_split[0]
            start = int(line_split[3])
            end   = int(line_split[4])
            strand = line_split[6]
            
            codes_tab = line.split('\t')[8]
            gene_type_list = codes_tab.split('gene_type=')[1].split(';')[0]
            #print(gene_type_list)
            if gene_type_list.find(',') > 0:
                print('line:')
                print(line)
                print('codes_tabl:')
                print(codes_tab)
                print('gene_type_list')
                print(gene_type_list)
                print("codes_tab.split('gene_type=')[1].split(';')[0]")
                print(codes_tab.split('gene_type=')[1].split(';')[0])
                int('error. exception.')
            
            transcript_type = codes_tab.split('transcript_type=')[1].split(';')[0]
            if transcript_type not in exon_id_transcript_sets_dict:
                exon_id_transcript_sets_dict[transcript_type] = set()
            
            if transcript_type.find(',') > 0:
                print(line)
                int('unholy cow')
            
            #if line.find('nonsense_mediated_decay') > 0:
            #    print(line)
            #    int('unholy cow')
                
            #for gene_type in gene_type_list:
            gene_type = gene_type_list
            if gene_type not in exon_id_sets_dict:
                exon_id_sets_dict[gene_type] = set()
                
            if chrom not in chromsomes_list:
                continue
            
            
            
            tmp = line_split[8].split('ID=exon:')[1]
            tmp_exon = tmp.split(';')[0]
            transcript_id = tmp_exon.split(':')[0]
            exon_id = transcript_id + '_' + str(int(tmp_exon.split(':')[1] ) - 1)
            
            
            gene_id=line_split[8].split('gene_id=')[1].split(';')[0]
            
            region_tid = "%s:%d-%d:%s" % (chrom,start,end,strand)
            transcript_to_gene_dict[region_tid]={'start':start,'end':end,'transcript_id':transcript_id,'strand':strand,'chrom':chrom,'gene_id':gene_id}
            
            
            try:
                seq = genome_fasta[chrom][start-1:end]
            except:
                print(chrom,start,end)
                print('line', line)
                int('MOOON')
            if strand == '-':
                seq = seq.reverse.complement
            seq = str(seq)
            
            if len(seq) < 20:
                continue
            
            
            custom_chrom_region = "%s:%d-%d:%s" % (chrom, start, end+1, strand)
            
            #for gene_type in gene_type_list:
            exon_id_sets_dict[gene_type].add(custom_chrom_region)
            exon_id_transcript_sets_dict[transcript_type].add(custom_chrom_region)
                
            #make a dictionary with every exon sequence
            if seq not in seq_to_exon_region_dict:
                seq_to_exon_region_dict[seq] = list()
            if custom_chrom_region not in seq_to_exon_region_dict[seq]:
                seq_to_exon_region_dict[seq].append(custom_chrom_region)
                
                
            if len(seq_to_exon_region_dict[seq]) > 1:
                if seq not in duplicate_exon_seq_to_exon_region_dict:
                    duplicate_exon_seq_to_exon_region_dict[seq] = seq_to_exon_region_dict[seq]
                if custom_chrom_region not in duplicate_exon_seq_to_exon_region_dict[seq]:
                    duplicate_exon_seq_to_exon_region_dict[seq].append(custom_chrom_region)
                    
            transcript_list.append(transcript_id)
            
            if transcript_id not in transcript_id_count:
                transcript_id_count[transcript_id] = [ [exon_id, custom_chrom_region] ]
            else:
                transcript_id_count[transcript_id].append([exon_id, custom_chrom_region])
                
                
                
            if custom_chrom_region not in exon_data_dict:
                exon_data_dict[custom_chrom_region] = [chrom, start, end, abs(start-end), transcript_id, seq]
                
            if custom_chrom_region not in exon_data_dict_25_to_400 and abs(start-end) >= 25 and abs(start-end) <= 400:
                exon_data_dict_25_to_400[custom_chrom_region] = [chrom, start, end, abs(start-end), transcript_id, seq]
                
            if custom_chrom_region not in exon_data_dict_50_to_200 and abs(start-end) >= 50 and abs(start-end) <= 200:
                exon_data_dict_50_to_200[custom_chrom_region] = [chrom, start, end, abs(start-end), transcript_id, seq]
                
    exon_annotation_datasets = dict()
    exon_annotation_datasets['all'] = exon_data_dict
    exon_annotation_datasets['50_200'] = exon_data_dict_50_to_200
    exon_annotation_datasets['25_400'] = exon_data_dict_25_to_400
    
    dup_exon_dict = dict()
    dup_exon_dict['seq_to_exon_region_dict'] = seq_to_exon_region_dict
    dup_exon_dict['duplicate_exon_seq_to_exon_region_dict'] = duplicate_exon_seq_to_exon_region_dict
    
    if version == 'comprehensive':
        return exon_annotation_datasets, transcript_id_count, dup_exon_dict, UTR_annotation, exon_id_sets_dict,GENCODE_exon_dict,exon_id_transcript_sets_dict, transcript_to_gene_dict, GENCODE_gene_dict
    else:
        return exon_annotation_datasets, transcript_id_count, dup_exon_dict, UTR_annotation, exon_id_sets_dict,GENCODE_exon_dict,exon_id_transcript_sets_dict, transcript_to_gene_dict, GENCODE_gene_dict








TSL_string = "transcript_support_level"










class GENCODE_exon_class:
    def build_from_gff3_entry(self, line):
        line_split = line.strip('\n').split('\t')
        #name = line_split
        self.start = int(line_split[3])
        self.end = int(line_split[4]) + 1
        self.chrom = line_split[0]
        self.source_type = line_split[1]
        self.feature_type = line_split[2]
        self.score = line_split[5]
        self.strand = line_split[6]
        self.phase = line_split[7]
        self.key_dict = dict()
        meta = line_split[8]
        meta_split = meta.split(';')
        for entry in meta_split:
            entry_split = entry.split('=')
            key = entry_split[0]
            val = entry_split[1]
            if val.isnumeric() == True:
                self.key_dict[key]=int(val)
            else:
                self.key_dict[key]=val
        self.exon_id = "%s:%d-%d:%s" % (self.chrom,self.start,self.end,self.strand)
        self.transcript_id = self.key_dict['transcript_id']
        self.id=self.key_dict['transcript_id']+':'+str(self.key_dict['exon_number'])
        self.gencode_exon_id = self.key_dict['exon_id']+'_'+self.key_dict['ID']
        self.length=abs(self.start-self.end)
        self.exon_number=int(self.key_dict['exon_number'])


