




from sys import argv
if 'exon_count_build' in argv:
    exon_count_build = int(argv[2])
else:
    exon_count_build = 100
    #exon_count_build = 1










pickle_path = exp_output_path.pickle_merged + "calculate_chromosome_coverage__%d.pickle" % (exon_count_build)

print("start loading calculate_chromosome_coverage.pickle")

opened_pickle = open(pickle_path,'rb')

no_overlap_set = pickle.load(opened_pickle)
intron_interior_set = pickle.load(opened_pickle)
intron_interior_pc_set = pickle.load(opened_pickle)
intron_interior_lncRNA_set = pickle.load(opened_pickle)
antisense_transcript_set = pickle.load(opened_pickle)
exon_overlap_set = pickle.load(opened_pickle)
sense_transcript_set = pickle.load(opened_pickle)
pc_middle_exon_id_list = pickle.load(opened_pickle)
lncRNA_middle_exon_ids = pickle.load(opened_pickle)
region_exon_id_sets_list = pickle.load(opened_pickle)
region_exon_id_sets_dict = pickle.load(opened_pickle)


opened_pickle.close()


print("finnished loading calculate_chromosome_coverage.pickle")


