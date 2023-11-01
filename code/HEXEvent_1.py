


HEXEvent_file = exp_output_path.HEXEVENT_input_folder + 'HEXEvent_constitutive_cassette.txt'

#9 constitLevel
#10 inclLevel

HEXEvent_exon_id_list = list()
HEXEvent_exon_id_list_constitutive = list()
HEXEvent_exon_id_list_cassette = list()

higher_constitutive_threshold = 0.99
lower_constitutive_threshold  = 0.9
constitutive_threshold = higher_constitutive_threshold
with open(HEXEvent_file) as f:
    line = next(f)
    print(line)
    for ii, col in enumerate(line.split('\t')):
        print(ii, col)
    for line in f:
        if line == '\n':
            break
        line_split = line.split('\t')
        exon_id = "%s:%d-%d:%s" % (line_split[0],int(line_split[2])+1,int(line_split[3])+1,line_split[1])
        HEXEvent_exon_id_list.append(exon_id)
        
        if float(line_split[9]) >= higher_constitutive_threshold:
            HEXEvent_exon_id_list_constitutive.append(exon_id)
        elif float(line_split[9]) <= lower_constitutive_threshold:
            HEXEvent_exon_id_list_cassette.append(exon_id)






import pickle 


pickle_path = exp_output_path.pickle_main
with open(exp_output_path.pickle_main+'HEXEvent__%d.pickle'%exon_count_build, 'wb') as outfile:
    pickle.dump(HEXEvent_exon_id_list_cassette,outfile)
    pickle.dump(HEXEvent_exon_id_list_constitutive,outfile)
    pickle.dump(HEXEvent_exon_id_list,outfile)




alt = el.threshold_exon_ids( el.exon_id_intersection(HEXEvent_exon_id_list_cassette, aggregate_exon_dict.keys() ), 1, aggregate_exon_dict)

len( el.exon_id_intersection(HEXEvent_exon_id_list_cassette, aggregate_exon_dict.keys()) )/len(HEXEvent_exon_id_list_cassette)
len(el.threshold_exon_ids( aggregate_exon_dict.keys() , 100, aggregate_exon_dict))


