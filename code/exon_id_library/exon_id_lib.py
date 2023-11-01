"""

This is a set of functions meant to operate of exon_id(s)

exon_id(s) are genomic coordinates of the format:
    chrom:start-end:strand



"""



'''
try:
    int('skip')
    genome_fasta = 1
except:
    1
'''




from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
matrix5 = load_matrix5()
matrix3 = load_matrix3()









def write_exon_ids_to_bed(exon_id_list, f_path):
    with open(f_path,'w') as f:
        header = 'chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n'
        f.write(header)
        for exon_id in exon_id_list:
            ex=el.exon_id_values(exon_id)
            outline = '%s\t%d\t%d\t%s\t0\t%s\n' % (ex.chrom, ex.start, ex.end, exon_id,ex.strand)
            f.write(outline)
            


def write_exon_ids_and_counts_to_bed(exon_id_list, f_path, exon_dict):
    with open(f_path,'w') as f:
        header = 'chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n'
        f.write(header)
        for exon_id in exon_id_list:
            ex=exon_id_values(exon_id)
            outline = '%s\t%d\t%d\t%s\t%d\t%s\n' % (ex.chrom, ex.start, ex.end, exon_id,exon_dict[exon_id]['count'],ex.strand)
            f.write(outline)
            
            











def get_exon_dict_pair_3ss_5ss_gc_len_scores(exon_id_list,exon_dict,genome_fasta):
    score_list = list()
    for exon_id in exon_id_list:
        entry_list = list()
        score = exon_dict[exon_id]['3ss_score']
        if score < 50 and score > -90:
            entry_list.append(score)
        else:
            continue
        
        score = exon_dict[exon_id]['5ss_score']
        if score < 50 and score > -90:
            entry_list.append(score)
        else:
            continue
        
        gc = get_gc_for_region(exon_id,genome_fasta)
        entry_list.append(gc)
        
        
        length = exon_id_values(exon_id).length
        entry_list.append(length)
        
        
        score_list.append(entry_list)

    return score_list



def get_exon_dict_3ss_scores(exon_id_list,exon_dict):
    score_list = list()
    for exon_id in exon_id_list:
        score = exon_dict[exon_id]['3ss_score']
        if score < 50 and score > -90:
            score_list.append(exon_dict[exon_id]['3ss_score'])

    return score_list



def get_exon_dict_5ss_scores(exon_id_list,exon_dict):
    score_list = list()
    for exon_id in exon_id_list:
        score = exon_dict[exon_id]['5ss_score']
        if score < 50 and score > -90:
            score_list.append(exon_dict[exon_id]['5ss_score'])

    return score_list



def get_exon_dict_counts(exon_id_list,exon_dict):
    counts_list = list()
    for exon_id in exon_id_list:
        counts_list.append(exon_dict[exon_id]['count'])

    return counts_list

def get_exon_dict_lengths(exon_id_list,exon_dict):
    counts_list = list()
    for exon_id in exon_id_list:
        counts_list.append(exon_dict[exon_id]['length'])

    return counts_list


def size_exon_id_list(exon_id_list,lower, upper):
    sized_exon_id_list = [exon_id for ii, exon_id in enumerate(exon_id_list) if exon_id_values(exon_id).length > lower and exon_id_values(exon_id).length < upper]
    return sized_exon_id_list


def exon_id_intersection(exon_id_list_1,exon_id_list_2):
    return list(set(exon_id_list_1).intersection(set(exon_id_list_2)))

def exon_id_difference(exon_id_list_1,exon_id_list_2):
    return list(set(exon_id_list_1).difference(set(exon_id_list_2)))






def get_lncRNA_transcript_ids(GENCODE_exon_dict_basic):
    lncRNA_exon_tid_list = list()
    for key in GENCODE_exon_dict_basic:
        exon = GENCODE_exon_dict_basic[key]
        
        if exon.key_dict['transcript_type']=='lncRNA':
            lncRNA_exon_tid_list.append(exon.transcript_id)
    return list(set(lncRNA_exon_tid_list))


def get_protein_coding_transcript_ids(GENCODE_exon_dict_basic):
    pc_exon_tid_list = list()
    for key in GENCODE_exon_dict_basic:
        exon = GENCODE_exon_dict_basic[key]
        
        if exon.key_dict['transcript_type']=='protein_coding':
            pc_exon_tid_list.append(exon.transcript_id)
    return list(set(pc_exon_tid_list))



def get_lncRNA_exon_ids(GENCODE_exon_dict_basic):
    lncRNA_exon_id_list = list()
    for key in GENCODE_exon_dict_basic:
        exon = GENCODE_exon_dict_basic[key]
        if exon.key_dict['transcript_type']=='lncRNA':
            lncRNA_exon_id_list.append(exon.exon_id)
    return lncRNA_exon_id_list




class exon_id_values:
    def __init__(self, exon_id):
        split = exon_id.split(':')
        coord=split[1].split('-')
        self.chrom = split[0]
        self.start=int(coord[0])
        self.end =int(coord[1])
        self.strand=split[2]
        self.exon_id=exon_id
        self.length=abs(self.start-self.end)

        if self.strand == '+':
            self.id_3ss = "%s:%d:%s" % (self.chrom, self.start, self.strand)
            self.id_5ss = "%s:%d:%s" % (self.chrom,self.end,self.strand)
            self.pos_5ss = self.end
            self.pos_3ss = self.start
        else:  # strand == '-'
            self.id_5ss = "%s:%d:%s" % (self.chrom,self.start,self.strand)
            self.id_3ss = "%s:%d:%s" % (self.chrom,self.end,self.strand)
            self.pos_3ss = self.end
            self.pos_5ss = self.start
    def score_5ss(self, genome_fasta):
        
        if self.strand == '+':
            seq = genome_fasta[self.chrom][self.end-4:self.end+5]
            seq = str(seq).upper()
            if seq.find('N') >= 0:
                return -5000, seq,  'found N'
            score = maxent.score5(str(seq), matrix = matrix5)
            
        if self.strand == '-':
            seq = genome_fasta[self.chrom][self.start-7:self.start+2].reverse.complement
            seq = str(seq).upper()
            if seq.find('N') >= 0:
                return -5000, seq,  'found N'
            seq = str(seq)
            score = maxent.score5(str(seq), matrix = matrix5)
                
        return score, seq, 'no N found'
    
        
    def score_3ss(self, genome_fasta):
        
        if self.strand == '+':
            seq = genome_fasta[self.chrom][self.start-21:self.start+2]
            seq = str(seq).upper()
            if seq.find('N') >= 0:
                return -5000, seq,  'found N'
            score = maxent.score3(str(seq), matrix = matrix3)
            
        if self.strand == '-':
            seq = genome_fasta[self.chrom][self.end-4:self.end+19].reverse.complement
            seq = str(seq).upper()
            if seq.find('N') >= 0:
                return -5000, seq,  'found N'
            seq = str(seq)
            score = maxent.score3(str(seq), matrix = matrix3)
                      
        return score, seq, 'no N found'
     





def exon_id_share_5ss(exon_id_1,exon_id_2):
    ex_1 = exon_id_values(exon_id_1)
    ex_2 = exon_id_values(exon_id_2)

    result=dict()
    result['share_5ss'] = (ex_1.id_5ss == ex_2.id_5ss)
    result['share_3ss'] = (ex_1.id_3ss == ex_2.id_3ss)
    result['ex_1']=ex_1
    result['ex_2']=ex_2

    return result

def get_overlap_exon_A_with_B(exon_id_A,exon_id_B):
    A = exon_id_values(exon_id_A)
    B = exon_id_values(exon_id_B)

    if A.strand == B.strand:
        range_A=set(range(A.start,A.end))
        range_B=set(range(B.start,B.end))
        overlap = len(range_A.intersection(range_B))
        len_A = len(range_A)
        len_B = len(range_B)
        return overlap, len_A, len_B
    else:

        return 0, len(range(A.start,A.end)),len(range(B.start,B.end))


#returns list of exon ids that 1. overlap with 2. provided exon interval tree
#also return list of exon_ids tested and whether found True/False

def query_overlap_with_recovered_exons(exon_id, IT):
    ex = exon_id_values(exon_id)

    found_intervals = IT[ex.chrom][ex.strand].overlap(ex.start,ex.end)
    if len(found_intervals) > 0:
        found_overlap_flag = True
    else:
        found_overlap_flag = False
    found_list = list()
    for interval in found_intervals:
        for val in interval[2]:
            found_list.append(val)

    result=dict()
    result['found']=found_list
    result['found_overlap_flag']=found_overlap_flag
    return result

#return query_result_dict['found'] which is dict of exon_ids with list of exon_ids in intervals found
#query_result_dict['not_found'] which is list of exon_ids that do not find overalpping intervals
def query_IT_with_exon_id_list(exon_id_list, IT):
    query_result_dict = dict()
    query_result_dict['found']=dict()
    query_result_dict['not_found']=list()
    #found_as_dict = dict()
    #query_result_dict['found_as_dict']=dict()
    for exon_id in exon_id_list:
        result=query_overlap_with_recovered_exons(exon_id, IT)
        if result['found_overlap_flag'] == True:
            query_result_dict['found'][exon_id]=result['found']
        else:
            query_result_dict['not_found'].append(exon_id)
    return query_result_dict





def build_single_ss_dict_from_exon_id_list(exon_id_list):
    exon_5ss_id_dict = dict()
    exon_3ss_id_dict = dict()
    for exon_id in exon_id_list:
        ex = exon_id_values(exon_id)
        if ex.id_5ss not in exon_5ss_id_dict:
            exon_5ss_id_dict[ex.id_5ss]=[exon_id]
        else:
            exon_5ss_id_dict[ex.id_5ss].append(exon_id)

        if ex.id_3ss not in exon_3ss_id_dict:
            exon_3ss_id_dict[ex.id_3ss]=[exon_id]
        else:
            exon_3ss_id_dict[ex.id_3ss].append(exon_id)

        #exon_3ss_id_dict[ex.id_3ss]
    return exon_5ss_id_dict, exon_3ss_id_dict




















def get_highly_overlapping_non_exact_exon_dict(exon_id_list, exon_IT, exon_dict):


    query_result_dict = query_IT_with_exon_id_list(exon_id_list, exon_IT)
    #query_result_dict['found'][exon_id]=result['found']   #hashed is_found search
    #query_result_dict['not_found'].append(exon_id) #nothing in the data overlaps this set

    has_overlapping = list(query_result_dict['found'].keys())
    #has_overlapping = list(annotated_overlap_ET_exon_id['found].keys())

    exact_in_overlapping = list( set(has_overlapping).intersection(exon_dict.keys()) )
    only_overlapping = list(set(has_overlapping).difference(set(exact_in_overlapping)))

    best_overlapping_annotated_not_exact=dict()
    for annotated_entry in only_overlapping:

        for ii, overlapping_ET_entry in enumerate(query_result_dict['found'][annotated_entry]):
            exon_id_A = annotated_entry
            exon_id_B = overlapping_ET_entry
            len_overlap, len_A, len_B = get_overlap_exon_A_with_B(exon_id_A,exon_id_B)
            if abs(len_B-len_A) < 5 and abs(len_overlap-len_A) < 5:
                best_overlapping_annotated_not_exact[exon_id_A]=[exon_id_A,abs(len_overlap-len_A)]

    exact_list = exact_in_overlapping
    highly_overlapping=best_overlapping_annotated_not_exact

    return highly_overlapping, exact_list, has_overlapping



class exon_assignment_class:
    def __init__(self, check_list, exon_IT, exon_dict):
        ET_data_keys = set(exon_dict.keys())
        check_list = set(check_list)

        exact        = check_list.intersection(ET_data_keys)
        missed_exact = check_list.difference(ET_data_keys)

        highly_overlapping, exact_list, has_overlapping = get_highly_overlapping_non_exact_exon_dict(missed_exact, exon_IT, exon_dict)

        fuzzy = set(highly_overlapping)
        exon_id_list = fuzzy
        query_IT_result_dict = query_IT_with_exon_id_list(exon_id_list, exon_IT)

        share_5ss_exon_id_list,share_3ss_exon_id_list, recovering_5ss_exon_id_list, recovering_3ss_exon_id_list, exact_match_to_alternat_dict, alternate_to_exact_mach_dict = get_exon_ids_that_share_5ss_3ss_in_ET_data(query_IT_result_dict)

        list_exon_id_captured_by_dual_exons, list_dual_exon_id, dual_exon_id_to_pair_individual_dict = get_individual_exons_from_dual(exact_match_to_alternat_dict, alternate_to_exact_mach_dict)

        list_exon_id_captured_by_dual_exons=set(list_exon_id_captured_by_dual_exons)

        self.missed_exact_list         = list(check_list.difference(exact))
        self.missed_fuzzy_list         = list(check_list.difference(fuzzy))
        self.missed_dual_list          = list(check_list.difference(list_exon_id_captured_by_dual_exons))
        self.missed_all_list           = list(check_list.difference(exact).difference(fuzzy).difference(list_exon_id_captured_by_dual_exons))

        self.exon_id_list                        = list(check_list)
        self.recovered_list                      = list(exact.union(fuzzy).union(list_exon_id_captured_by_dual_exons))
        self.recovered_exact_list                = list(exact)
        self.recovered_fuzzy_list                = list(fuzzy)
        self.recovered_part_of_dual_exon_list    = list(list_exon_id_captured_by_dual_exons)
        self.dual_exon_list                      = list(list_dual_exon_id)

        self.len_exon_id_list                   = len(check_list)
        self.len_recovered_list                 = len(self.recovered_list)
        self.len_recovered_exact_list           = len(self.recovered_exact_list)
        self.len_recovered_fuzzy_list           = len(self.recovered_fuzzy_list)
        self.len_recovered_part_of_dual_exon_list = len(self.recovered_part_of_dual_exon_list)
        self.len_dual_exon_list                 = len(self.dual_exon_list)

        if self.len_exon_id_list == 0:
            self.recovery_ratio = -1
        else:
            self.recovery_ratio = self.len_recovered_list/self.len_exon_id_list

        self.dual_exon_id_to_pair_individual_dict = dual_exon_id_to_pair_individual_dict
        self.share_5ss_exon_id_list               = share_5ss_exon_id_list
        self.share_3ss_exon_id_list               = share_3ss_exon_id_list
        self.recovering_5ss_exon_id_list          = recovering_5ss_exon_id_list
        self.recovering_3ss_exon_id_list          = recovering_3ss_exon_id_list
        self.exact_match_to_alternat_dict         = exact_match_to_alternat_dict
        self.alternate_to_exact_mach_dict         = alternate_to_exact_mach_dict
        self.has_overlapping                      = has_overlapping
    def print_ratio(self):
        print("recovery_ratio: %.2f" % self.recovery_ratio)
        print("len_exon_id_list: ", self.len_exon_id_list)
        print("len_recovered_list: ", self.len_recovered_list)
        print("len_recovered_exact_list: ", self.len_recovered_exact_list)
        print("len_recovered_fuzzy_list: ", self.len_recovered_fuzzy_list)
        print("len_recovered_part_of_dual_exon_list: ", self.len_recovered_part_of_dual_exon_list)
        print("len_dual_exon_list: ", self.len_dual_exon_list)






# share_5ss_exon_id_list - input list (often annotated) exon_ids that share a 5ss with an overlapping ET exon_id
# share_3ss_exon_id_list - input list (often annotated) exon_ids that share a 3ss with an overlapping ET exon_id
# recovering_5ss_exon_id_list - ET exon_ids that recover an annotated 5'ss
# recovering_3ss_exon_id_list - ET exon_ids that recover an annotated 3'ss
# note that this data does not allow for sharing both the 5ss and 3ss between the input list and the ET data dict
def get_exon_ids_that_share_5ss_3ss_in_ET_data(query_IT_result_dict):
    share_5ss_exon_id_list=list()
    share_3ss_exon_id_list=list()
    recovering_5ss_exon_id_list=list()
    recovering_3ss_exon_id_list=list()
    exact_match_to_alternat_dict = dict()
    alternate_to_exact_mach_dict = dict()
    for exon_id_1 in query_IT_result_dict['found']:   #these represent the checked exon_ids, often annotated
        #result_dict = dict()
        share_5ss = False
        share_3ss = False
        for exon_id_2 in query_IT_result_dict['found'][exon_id_1]:   #these represent the  exon_ids found in a cluster overlapping the (often annotated) checked exon_id
            #for exon_id_2 in interval[2]:
            result = exon_id_share_5ss(exon_id_1,exon_id_2)
            if result['share_5ss'] == True and result['share_3ss']  == False:
                share_5ss_exon_id_list.append(exon_id_1)
                recovering_5ss_exon_id_list.append(exon_id_2)
            if result['share_3ss'] == True and result['share_5ss']  == False:
                share_3ss_exon_id_list.append(exon_id_1)
                recovering_3ss_exon_id_list.append(exon_id_2)


            #exon_id_1 is the (often annotated) exon
            #exon_id_2 is a found overlapping ET data exon
            if result['share_3ss'] == True or result['share_5ss']  == True:
                if exon_id_2 == exon_id_1: #make sure this is an alternate isoform
                    continue               #if they match skip
                if exon_id_1 not in exact_match_to_alternat_dict:
                    exact_match_to_alternat_dict[exon_id_1]=list()

                #for exon_id_3 in query_IT_result_dict['found'][exon_id_1]:  # we want to make sure that overlapping exon_id shares a splice site since clusters can consist of a group of exon_ids overlapping a large region

                #result_2 = exon_id_share_5ss(exon_id_1,exon_id_3)
                #if result_2['share_3ss'] == True or result_2['share_5ss']  == True:  #make sure this shares a splice site
                exact_match_to_alternat_dict[exon_id_1].append(exon_id_2)
                #exact_match_to_alternat_dict[exon_id_1]=set(exact_match_to_alternat_dict[exon_id_1])
                #exact_match_to_alternat_dict[exon_id_1]=exact_match_to_alternat_dict[exon_id_1].difference(set([exon_id_1]))
                exact_match_to_alternat_dict[exon_id_1]=list(set(exact_match_to_alternat_dict[exon_id_1]))

                if exon_id_2 not in alternate_to_exact_mach_dict:
                    alternate_to_exact_mach_dict[exon_id_2]=list()
                alternate_to_exact_mach_dict[exon_id_2].append(exon_id_1)
                alternate_to_exact_mach_dict[exon_id_2] = list(set(alternate_to_exact_mach_dict[exon_id_2]))

    share_5ss_exon_id_list=list(set(share_5ss_exon_id_list))
    share_3ss_exon_id_list=list(set(share_3ss_exon_id_list))
    recovering_3ss_exon_id_list=list(set(recovering_3ss_exon_id_list))
    recovering_5ss_exon_id_list=list(set(recovering_5ss_exon_id_list))

    return share_5ss_exon_id_list,share_3ss_exon_id_list, recovering_5ss_exon_id_list, recovering_3ss_exon_id_list,exact_match_to_alternat_dict, alternate_to_exact_mach_dict


#result_1 = exon_id_share_5ss(exon_id_1,exon_id_2)
#if result['share_3ss'] == True or result['share_5ss']  == True:
#uses exact_match_to_alternat_dict returned by another function
def get_individual_exons_from_dual(exact_match_to_alternat_dict, alternate_to_exact_mach_dict):
    list_exon_id_captured_by_dual_exons = list()
    list_dual_exon_id=list()
    dual_exon_id_to_pair_individual_dict = dict()
    for alternate_id_1 in alternate_to_exact_mach_dict:

        check_list = alternate_to_exact_mach_dict[alternate_id_1] #list of (often annotated) exon_ids associated with this alternate exon

        for exon_id_1 in check_list:
            for exon_id_2 in check_list:
                if exon_id_1 == exon_id_2:
                    continue
                if alternate_id_1 in exact_match_to_alternat_dict[exon_id_1] and alternate_id_1 in exact_match_to_alternat_dict[exon_id_2]:
                    1
                else:
                    continue

                result_1 = exon_id_share_5ss(alternate_id_1,exon_id_1)
                result_2 = exon_id_share_5ss(alternate_id_1,exon_id_2)
                either_5ss_3ss = result_1['share_5ss'] and result_2['share_3ss']
                either_3ss_5ss = result_1['share_3ss'] and result_2['share_5ss']
                
                if either_5ss_3ss == True or either_3ss_5ss == True:
                    list_exon_id_captured_by_dual_exons.append(exon_id_1)
                    list_exon_id_captured_by_dual_exons.append(exon_id_2)
                    list_dual_exon_id.append(alternate_id_1)
                    dual_exon_id_to_pair_individual_dict[alternate_id_1]=[exon_id_1,exon_id_2]

    list_exon_id_captured_by_dual_exons = list(set(list_exon_id_captured_by_dual_exons))
    list_dual_exon_id = list(set(list_dual_exon_id))

    return list_exon_id_captured_by_dual_exons, list_dual_exon_id, dual_exon_id_to_pair_individual_dict






def get_seq_from_exon_id(exon_id, genome_fasta):
    ex = exon_id_values(exon_id)
    
    #if ex.strand == '+':
    seq = genome_fasta[ex.chrom][ex.start-1:ex.end-1]
    if ex.strand == '-':
        #seq = genome_fasta[ex.chrom][ex.start-1:ex.end-1]
        seq = seq.reverse.complement
    seq = str(seq)

    return seq





def get_gc_for_region(exon_id,genome_fasta):
    seq = get_seq_from_exon_id(exon_id,genome_fasta)
    gc=0
    for base in seq:
        if base == 'c' or base == 'C' or base == 'g' or base == 'G':
                gc += 1
    return gc/len(seq)



def get_gc_ratio_5ss(exon_id,genome_fasta):
    ex = exon_id_values(exon_id)
    
    if ex.strand == '+':
        exon_last_gc_region = "%s:%d-%d:%s" % (ex.chrom,ex.end-50,ex.end,ex.strand)
        exon_downstream_gc_region = "%s:%d-%d:%s" % (ex.chrom,ex.end,ex.end+50,ex.strand)
        exon_gc = get_gc_for_region(exon_last_gc_region,genome_fasta)
        intron_gc = get_gc_for_region(exon_downstream_gc_region,genome_fasta)
        ratio = exon_gc/intron_gc
    
    if ex.strand == '-':
        exon_last_gc_region = "%s:%d-%d:%s" % (ex.chrom,ex.start,ex.start+50,ex.strand)
        exon_downstream_gc_region = "%s:%d-%d:%s" % (ex.chrom,ex.start-50,ex.start,ex.strand)
        exon_gc = get_gc_for_region(exon_last_gc_region,genome_fasta)
        intron_gc = get_gc_for_region(exon_downstream_gc_region,genome_fasta)
        ratio = exon_gc/intron_gc
    
    return ratio, exon_gc, intron_gc


text_exon_id = 'chr21:41070537-41070817:+'
text_exon_id = 'chr4:1457317-1457463:-'
#get_gc_ratio_5ss(text_exon_id,genome_fasta)



def get_gc_ratio_3ss(exon_id,genome_fasta):
    ex = exon_id_values(exon_id)
    
    if ex.strand == '+':
        exon_gc_region = "%s:%d-%d:%s" % (ex.chrom,ex.start,ex.start+50,ex.strand)
        intron_gc_region = "%s:%d-%d:%s" % (ex.chrom,ex.start-50,ex.start,ex.strand)
        exon_gc = get_gc_for_region(exon_gc_region,genome_fasta)
        intron_gc = get_gc_for_region(intron_gc_region,genome_fasta)
        ratio = exon_gc/intron_gc
    
    if ex.strand == '-':
        exon_gc_region = "%s:%d-%d:%s" % (ex.chrom,ex.end-50,ex.end,ex.strand)
        intron_gc_region = "%s:%d-%d:%s" % (ex.chrom,ex.end,ex.end+50,ex.strand)
        exon_gc = get_gc_for_region(exon_gc_region,genome_fasta)
        intron_gc = get_gc_for_region(intron_gc_region,genome_fasta)
        ratio = exon_gc/intron_gc
    
    return ratio, exon_gc, intron_gc


text_exon_id = 'chr21:41070537-41070817:+'
#get_gc_ratio_3ss(text_exon_id,genome_fasta)

text_exon_id = 'chr4:1457317-1457463:-'
#get_gc_ratio_3ss(text_exon_id,genome_fasta)







def get_exon_id_list_binned_GC(exon_id_list, gc_bins, genome_fasta):
    gc_bins = sorted(gc_bins)
    
    if gc_bins[0]==0:
        1
    else:
        gc_bins = sorted([0] + gc_bins)
        
    if gc_bins[-1]==1:
        1
    else:
        gc_bins.append(100)
        
    gc_results = dict()
    gc_bin_pairs = dict()
    print(gc_bins)
    for ii in range(len(gc_bins)-1):
        bin_range = [ gc_bins[ii], gc_bins[ii+1] ]
        bin_id = "%d-%d" % (bin_range[0], bin_range[1])
        gc_results[bin_id]=list()
        gc_bin_pairs[bin_id]=bin_range
        
    for ii, exon_id in enumerate(exon_id_list):
        gc = 0
        
        seq = get_seq_from_exon_id(exon_id, genome_fasta)
        
        for base in seq:
            if base == 'c' or base == 'C' or base == 'g' or base == 'G':
                gc += 1
        ratio = gc/len(seq)*100
        for bin_id in gc_bin_pairs:
            bin_range = gc_bin_pairs[bin_id]
            if  bin_range[0] < ratio and ratio <= bin_range[1]:
                gc_results[bin_id].append(exon_id)
    return gc_results










def scan_seq_to_exon_id_list_dict(exon_id_list,exon_dict,genome_fasta):
    seq_to_exon_id_list_dict=dict()

    for exon_id in exon_id_list:
        #exon_seq = exon_dict[exon_id]['seq']
        ex = exon_id_values(exon_id)
        exon_seq = genome_fasta[ex.chrom][ex.start:ex.end]
        if ex.strand == '-':
            exon_seq = exon_seq.reverse.complement
        exon_seq=str(exon_seq)

        if exon_seq not in seq_to_exon_id_list_dict:
            seq_to_exon_id_list_dict[exon_seq]=list()
        seq_to_exon_id_list_dict[exon_seq].append(exon_id)

    print(len(seq_to_exon_id_list_dict))
    unique_seq_exons_ids = list()
    for x in seq_to_exon_id_list_dict:
        if len(seq_to_exon_id_list_dict[x]) == 1:
            unique_seq_exons_ids.append(seq_to_exon_id_list_dict[x][0])
    #[seq_to_exon_id_list_dict[x][0] for x in seq_to_exon_id_list_dict if len(seq_to_exon_id_list_dict[x]) == 1]


    duplicate_seq_exon_ids = list()
    duplicate_seq = [x for x in seq_to_exon_id_list_dict if len(seq_to_exon_id_list_dict[x]) > 1]
    for x in duplicate_seq:
        duplicate_seq_exon_ids += seq_to_exon_id_list_dict[x]

    return unique_seq_exons_ids, duplicate_seq_exon_ids




def threshold_exon_ids(exon_id_list, count_threshold, aggregate_exon_dict):
    threshold_exon_id_list = set()
    for exon_id in exon_id_list:
        ex = aggregate_exon_dict[exon_id]
        if ex['count'] >= count_threshold:
            threshold_exon_id_list.add(exon_id)
        
    return list(threshold_exon_id_list)


def threshold_5ss_exon_ids(exon_id_list, score_threshold, aggregate_exon_dict):
    threshold_exon_id_list = set()
    for exon_id in exon_id_list:
        ex = aggregate_exon_dict[exon_id]
        if ex['5ss_score'] >= score_threshold:
            threshold_exon_id_list.add(exon_id)
        
    return list(threshold_exon_id_list)


def threshold_3ss_exon_ids(exon_id_list, score_threshold, aggregate_exon_dict):
    threshold_exon_id_list = set()
    for exon_id in exon_id_list:
        ex = aggregate_exon_dict[exon_id]
        if ex['3ss_score'] >= score_threshold:
            threshold_exon_id_list.add(exon_id)
        
    return list(threshold_exon_id_list)






def get_max_count_exon_id_in_list(exon_id_list, aggregate_exon_dict):
    max_count = 0
    max_count_exon_id = ''
    for exon_id in exon_id_list:
        if aggregate_exon_dict[exon_id]['count'] > max_count:
            max_count = aggregate_exon_dict[exon_id]['count']
            max_count_exon_id = exon_id
    return max_count_exon_id





def recovery_percent(exon_id_list, aggregate_exon_IT, aggregate_exon_dict, genome_fasta):
    
    exon_id_list=size_exon_id_list(exon_id_list,50,500)
    
    unique_exon_id_list, duplicate_seq_pc_middle_exon_id_list = scan_seq_to_exon_id_list_dict(exon_id_list,aggregate_exon_dict, genome_fasta)
    
    z=get_highly_overlapping_non_exact_exon_dict(unique_exon_id_list, aggregate_exon_IT, aggregate_exon_dict)
    highly_overlapping, exact_list, has_overlapping=z


    recovery_percentage = (len(highly_overlapping) + len(exact_list))/len(size_exon_id_list(unique_exon_id_list,50,500))
    
    return recovery_percentage, highly_overlapping, exact_list, has_overlapping





def exon_id_to_bed_line(exon_id):
    ex = exon_id_values(exon_id)
    line = '{:}\t{:}\t{:}'.format(ex.chrom, ex.start, ex.end)
    line = '{:}\t{:}\t0\t{:}'.format(line, exon_id, ex.strand)
    return line


def exon_id_with_counts_to_bed_line(exon_id, aggregate_exon_dict):
    ex = exon_id_values(exon_id)
    line = '{:}\t{:}\t{:}'.format(ex.chrom, ex.start, ex.end)
    line = '{:}\t{:}\t0\t{:}\t'.format(line, exon_id, aggregate_exon_dict[exon_id]['count'], ex.strand)
    
    library_counts_string = ','.join([str(x) for ii, x in enumerate(aggregate_exon_dict[exon_id]['lib_array']) ])
        
    
    line="{:}{:}".format(line, library_counts_string)
    
    return line


def bed_line_to_exon_id(line):
    line_split = line.split('\t')
    exon_id = "{:}:{:}-{:}:{:}".format(line_split[0],line_split[1],line_split[2],line_split[5])
    return exon_id




def export_exon_id_list_to_bed(exon_id_list, file_name_path):
    with open(file_name_path, 'w') as f:
        for exon_id in exon_id_list:
            f.write(exon_id_to_bed_line(exon_id) + '\n')


def load_bed_to_exon_id_list(exon_id_list, file_name_path):
    exon_id_list = list()
    with open(file_name_path,'r') as f:
        for line in f:
            exon_id = bed_line_to_exon_id(line)
            exon_id_list.append(exon_id)



def export_exon_id_list_with_counts_to_bed(exon_id_list, aggregate_exon_dict, file_name_path):
    with open(file_name_path, 'w') as f:
        for exon_id in exon_id_list:
            f.write(exon_id_with_counts_to_bed_line(exon_id, aggregate_exon_dict) + '\n')










import subprocess






