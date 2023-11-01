






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
        else:
            self.id_5ss = "%s:%d:%s" % (self.chrom,self.start,self.strand)
            self.id_3ss = "%s:%d:%s" % (self.chrom,self.end,self.strand)


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













import numpy as np
def calculate_mer_in_bit_seq(bit_seq,mer_len):
    print("convert_bitarray_to_seq(bit_seq)", convert_bitarray_to_seq(bit_seq))
    mer_array = np.zeros(4**mer_len)
    #print(range(0,int(len(bit_seq)/2)-mer_len))
    for ii in range(0,len(bit_seq)-2*mer_len+1,2):
        print('iteration', ii)
        print(bit_seq[ii:ii+2*mer_len])
        mer_id = bit_seq[ii:ii+2*mer_len]
        mer_id = util.ba2int(mer_id)
        print("mer_id",mer_id)
        mer_array[mer_id] += 1
        mer_id
        print("int2ba",util.int2ba(mer_id,endian="big",length=2*mer_len))
        print(convert_bitarray_to_seq( util.int2ba(mer_id,endian="big",length=2*mer_len) ))

    return mer_array

#util.int2ba(mer_id,endian="big",length=2*mer_len)

def convert_mer_id_to_ba(mer_id, mer_len):
    return util.int2ba(mer_id,endian="big",length=2*mer_len)

def convert_mer_id_to_seq(mer_id, mer_len):
    return convert_bitarray_to_seq(util.int2ba(mer_id,endian="big",length=2*mer_len))































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


            
            if result['share_3ss'] == True or result['share_5ss']  == True:
                if exon_id_2 == exon_id_1: #make sure this is an alternate isoform
                    continue               #if they match skip
                if exon_id_1 not in exact_match_to_alternat_dict:
                    exact_match_to_alternat_dict[exon_id_1]=list()

                
                exact_match_to_alternat_dict[exon_id_1].append(exon_id_2)
                
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

    seq = genome_fasta[ex.chrom][ex.start:ex.end]
    if ex.strand == '-':
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








