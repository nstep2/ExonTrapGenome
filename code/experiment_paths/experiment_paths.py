"""
File paths to input folders/files and output locations

This also loads a genome file handle and a python library to calculate maxentscan scores

This library also create a dictionary of ggplot suggested colors 

"""





import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import datetime
import pyfaidx    #genome_fasta created below



#MaxEntScan python wrappers
from maxentpy import maxent
from maxentpy.maxent import load_matrix5, load_matrix3
matrix5 = load_matrix5()
matrix3 = load_matrix3()


#Dictionary of plot colors
color_dict = dict()
color_dict['lncRNA'] = '#00B9F1'
color_dict['Antisense'] = '#DA6FAB'
color_dict['mRNA'] = '#0072BC'
color_dict['Intergenic'] = '#F15A22'
color_dict['Intronic'] = '#00A875'
color_dict['Alternative'] = '#ECDE38'


color_key_dict = dict()
color_key_dict['lncRNA_middle_exon_ids'] = color_dict['lncRNA']
color_key_dict['antisense_transcript_set'] = color_dict['Antisense']
color_key_dict['pc_middle_exon_id_list'] = color_dict['mRNA']
color_key_dict['no_overlap_set'] = color_dict['Intergenic']
color_key_dict['intron_interior_set'] = color_dict['Intronic']








'''
This checks if a path exists and if it does not then makes it. If the path didn't exist then also print the path
This function is called at the end of the script to generate the subfolder pathse.

'''
def make_dir_path(dir_path):
    import os

    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print("Created directory path:" , dir_path )
    else:    
        if True == False:
            print("Path already exists:" , dir_path )    




''' DEFINE PATHS TO CREATE FOR OUTPUTS + DEFINE SOME PARAMETERS'''

class exp_output_path_class():
    def __init__(self):
        self



#main project class declaration
exp_output_path = exp_output_path_class()



#MAIN script location. It is assumed all scripts are in this folder. Some classes and shell scripts exist in folders.

#exp_output_path.git_repo = '/home/dude/repositories/exon_def/et_main/'


# script_repository
exp_output_path.git_repo = '/mnt/hgfs/main_ssd/slim_test/code/exon_def/et_main/submit_aug_2023/'

# topmost level output paths
exp_output_path.main_dir = '/mnt/hgfs/main_ssd/et_test/'








### parameters
exp_output_path.hisat2_threads_count = 3  # alignment threads
exp_output_path.bbduk_threads_count = 3   #

exp_output_path.limt_reads_processed = 500000  # parameter > 0 limits run to parameter reads, otherwise all data processed



#models can be obtained from an installed SplicAI package, e.g. a conda install or a pip install
exp_output_path.spliceAI_model_path = exp_output_path.main_dir + '/input/SpliceAI_release_1.3.1/SpliceAI-1.3.1/spliceai/models/'





# input file paths

exp_output_path.adapter_fastq_file_path             = exp_output_path.main_dir + 'support_files/bbduk_files/fasta/trim_fasta_15mer_both_plasmid_forward.fa'

#exp_output_path.adapter_fastq_file_path_absolute    = exp_output_path.main_dir + 'support_files/bbduk_files/fasta/'
exp_output_path.adapter_fastq_file_path_absolute = exp_output_path.git_repo


#exp_output_path.initial_fastq_input_files           = '/mnt/hgfs/main_ssd/Sequencing_data/' + 'L8YVI8X/STE15968.20210202/210127_A00987_0132_AHTLTVDSXY/'

exp_output_path.initial_fastq_input_files           = '/mnt/hgfs/aux_ssd/test_download_GEO/bioproject_PRJNA878769/'



#exp_output_path.hg38_genome_fasta   =  exp_output_path.main_dir + 'input/genome/Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'

exp_output_path.hg38_genome_fasta   =  exp_output_path.main_dir + 'input/UCSC_hg_38/hg38.fa'


exp_output_path.GENCODE_GFF3        =  exp_output_path.main_dir + 'input/gencode_v37/gencode.v37.annotation.gff3'

exp_output_path.exon_sql_path       =  exp_output_path.main_dir + 'input/snaptron/exons.sqlite'   #snaptron database

exp_output_path.exon_database_bed   =  exp_output_path.main_dir + 'input/exon_database_bed/'

exp_output_path.lncNRA_database_bed =  '/mnt/hgfs/main_ssd/Linux_2TB/lncRNA_databases/'
#exp_output_path.lncNRA_database_bed =  exp_output_path.main_dir + 'input/lncRNA_databases/'

exp_output_path.HEK293_PSI_input    =  exp_output_path.main_dir + 'input/HEK293_alt_splicing/GSE221838_d0_v_d1000_SE.MATS.JCEC.csv'

exp_output_path.housekeeping_list   =  exp_output_path.main_dir + 'input/HEK293_exp/housekeeping/gkaa609_supplemental_files/Supplementary_Table1.txt'

exp_output_path.HEK293_exp_file     =  exp_output_path.main_dir + 'input/HEK293/GSE235387_HEK293.xlsx'

exp_output_path.embl_Dfam_consensus =  exp_output_path.main_dir + 'input/Dfam/Dfam.embl'

exp_output_path.dfam_hg38_hits      =  exp_output_path.main_dir + 'input/Dfam/hg38_dfam.nrph.hits'


#exp_output_path.bigwig_path = '/mnt/hgfs/Linux_pc/bigwig/'

exp_output_path.bigwig_path = '/mnt/hgfs/main_ssd/et_test/input/phyloP/'



###############    CREATE HANDLE TO FASTA USING pyfaidx

genome_fasta = pyfaidx.Fasta(exp_output_path.hg38_genome_fasta)

###############    





# top level output paths

#exp_output_path.main_dir                    =              '/mnt/hgfs/main_ssd/et_main/'

exp_output_path.initial_fastq_main_files    = exp_output_path.initial_fastq_input_files  #not actually fastq location?

exp_output_path.trimmed_fastq_input_files   = exp_output_path.main_dir + 'trimmed/'

exp_output_path.exons_data                  = exp_output_path.main_dir + 'main_data/'    # many subfolders are built on this location

exp_output_path.pickle_main                 = exp_output_path.main_dir + 'pickles/'   # intermediate data is stored as pickles here

exp_output_path.ssd_tmp                     = exp_output_path.main_dir + '/scratch_ssd/'

exp_output_path.out_bed_main                = exp_output_path.pickle_main + 'bed_output_files/'





# custom input files - modifications to existing files

#non-standard chromosomes removed


exp_output_path.hisat2_index_std_chroms = exp_output_path.main_dir + 'input/UCSC_hg_38/hg38_standard'





###other_files

'''
exp_output_path.Tra2B_pre_align_index ='/mnt/hgfs/main_ssd/work/human_genome/indexes/2bit/Tra2B_reporter_intron_2bit/Tra2B_reporter_intron_index_2bit'

exp_output_path.log_odds_framework_path = '/home/pineapple/repositories/exon_def/tools/log_odds_framework/'

exp_output_path.hg38_genome_fasta ='/home/pineapple/work/indexes/downloaded/Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'
'''





###project_files



exp_output_path.bam_dir                  = exp_output_path.trimmed_fastq_input_files + 'bam_merged/'
exp_output_path.bam_dir_sort             = exp_output_path.bam_dir + 'sort_dir/'

exp_output_path.trimmed_fastq_input_files_figures = exp_output_path.trimmed_fastq_input_files + 'figures/'

exp_output_path.trimmed_fastq_SAM_files     = exp_output_path.main_dir + 'SAM/'
exp_output_path.trimmed_fastq_SAM_files_md5 = exp_output_path.main_dir + 'SAM/md5/'
exp_output_path.trimmed_fastq_SAM_files_synthetic = exp_output_path.trimmed_fastq_SAM_files + 'synthetic/'


exp_output_path.main_database_dir           = exp_output_path.exons_data+'sqlite3/'   #???





exp_output_path.exons_data_merge_plots                      = exp_output_path.exons_data + 'plots/'
exp_output_path.out_supplemental                            = exp_output_path.exons_data + 'supplemental/'
exp_output_path.exons_data_merge_plots_weblog_3ss           = exp_output_path.exons_data_merge_plots + 'weblogo/3ss/aggregate/'
exp_output_path.exons_data_merge_plots_weblog_5ss           = exp_output_path.exons_data_merge_plots + 'weblogo/5ss/aggregate/'
exp_output_path.exons_data_merge_dinucleotide_box_plots     = exp_output_path.exons_data_merge_plots + 'dinucleotide_box_plots/'
exp_output_path.exons_data_merge_minor_location_histogram   = exp_output_path.exons_data_merge_plots + 'minor_location_histogram/'
exp_output_path.exons_data_merge_score_histograms           = exp_output_path.exons_data_merge_plots + 'score_histograms/'
exp_output_path.dfam_plots                                  = exp_output_path.exons_data + 'dfam/plots/'


exp_output_path.GENCODE_dir = exp_output_path.exons_data + 'GENCODE/'
exp_output_path.GENCODE_text = exp_output_path.GENCODE_dir +'text/'
exp_output_path.GENCODE_plot = exp_output_path.GENCODE_dir +'plot/'







exp_output_path.exons_data_merge_exp_an_un_plots = exp_output_path.exons_data_merge_plots + 'expression/an_vs_un/'



exp_output_path.pickle_individual = exp_output_path.pickle_main + 'libraries/'
exp_output_path.pickle_pre_merge = exp_output_path.pickle_main + 'libraries/pre_merge/'

exp_output_path.pickle_merged = exp_output_path.pickle_main + 'merged/'
exp_output_path.SpliceAI_ET_exons = exp_output_path.pickle_main + 'SpliceAI_ET_exons/'
exp_output_path.Dfam_pickle = exp_output_path.pickle_main + 'Dfam_pickle/'







exp_output_path.ESE_folder              = exp_output_path.main_dir + 'input/tmp/ESE/'
exp_output_path.ESE_input_folder        = exp_output_path.main_dir + 'input/tmp/ESE/'
exp_output_path.HEXEVENT_input_folder   = exp_output_path.main_dir + 'input/tmp/HEXEvent/'





#exp_output_path.chr17_folder = '/mnt/hgfs/main_ssd/Linux_2TB/Splice_AI_computation/chr17/'
#exp_output_path.chr17_shuffle_folder = '/mnt/hgfs/aux_ssd/aux_ssd/Splice_AI_computation/chr17_shuffle/'
#exp_output_path.chr17_shuffle_input_folder = '/mnt/hgfs/main_ssd/work/genome/SpliceAI/'



exp_output_path.chr17_folder = exp_output_path.main_dir + 'Splice_AI_computation/chr17/'
exp_output_path.chr17_shuffle_folder = exp_output_path.main_dir + 'Splice_AI_computation/chr17_shuffle/'
#exp_output_path.chr17_shuffle_input_folder = '/mnt/hgfs/main_ssd/work/genome/SpliceAI/'
exp_output_path.chr17_shuffle_input_folder = '/mnt/hgfs/main_ssd/work/genome/SpliceAI/'







exp_output_path.reshape_sam_read_assignments = exp_output_path.exons_data_merge_plots + 'reshape_sam_read_assignments/'

exp_output_path.ET_spliceAI_venn_2 = exp_output_path.exons_data_merge_plots + "ET_spliceAI_venn_2/"

exp_output_path.spliceAI_chr_17_analysis = exp_output_path.exons_data_merge_plots + "spliceAI_chr_17_analysis/"

exp_output_path.load_merged_aggregate_exon_dict = exp_output_path.exons_data_merge_plots + "load_merged_aggregate_exon_dict/"


exp_output_path.Parse_ENSEMBLE_GENCODE_exons_load = exp_output_path.exons_data_merge_plots + "Parse_ENSEMBLE_GENCODE_exons_load/"

exp_output_path.exon_number_expectation = exp_output_path.exons_data_merge_plots + "exon_number_expectation/"

exp_output_path.library_specifics = exp_output_path.exons_data_merge_plots + "library_specifics/"


exp_output_path.Parse_ENSEMBLE_GENCODE_exons_load = exp_output_path.exons_data_merge_plots + "Parse_ENSEMBLE_GENCODE_exons_load/"



exp_output_path.calculate_chromosome_coverage = exp_output_path.exons_data_merge_plots + "calculate_chromosome_coverage/"



exp_output_path.Dfam_repeats = exp_output_path.exons_data_merge_plots + "Dfam_repeats/"



exp_output_path.HEXEvent = exp_output_path.exons_data_merge_plots + "HEXEvent/"

exp_output_path.exon_id_conservations = exp_output_path.exons_data_merge_plots + "exon_id_conservations/"

exp_output_path.weblogo = exp_output_path.exons_data_merge_plots + "weblogo/"

exp_output_path.Parse_ENSEMBLE_GENCODE_exons_load_figures_2 = exp_output_path.exons_data_merge_plots + "Parse_ENSEMBLE_GENCODE_exons_load_figures_2/"

exp_output_path.ESE_scan = exp_output_path.exons_data_merge_plots + "ESE_scan/"


exp_output_path.old_gene = exp_output_path.exons_data_merge_plots + "old_gene/"

exp_output_path.rewrite_exon_number_expectation = exp_output_path.exons_data_merge_plots + "rewrite_exon_number_expectation/"


exp_output_path.Compute_spliceAI_exon_id_scratch = exp_output_path.exons_data_merge_plots + "Compute_spliceAI_exon_id_scratch/"

exp_output_path.proportion_bases_different_annotations = exp_output_path.exons_data_merge_plots + "proportion_bases_different_annotations/"



exp_output_path.load_merged_aggregate_exon_dict_figures_2 = exp_output_path.exons_data_merge_plots + "load_merged_aggregate_exon_dict_figures_2/"





exp_output_path.region_scores = exp_output_path.exons_data_merge_plots + "region_scores/"

exp_output_path.Compute_spliceAI_exon_id_scratch = exp_output_path.exons_data_merge_plots + "Compute_spliceAI_exon_id_scratch/"

exp_output_path.libraries_per_exon = exp_output_path.exons_data_merge_plots + "libraries_per_exon/"

exp_output_path.recovery_ratio_lncRNA_pc = exp_output_path.exons_data_merge_plots + "recovery_ratio_lncRNA_pc/"

exp_output_path.snaptron_exon = exp_output_path.exons_data_merge_plots + "snaptron_exon/"






#review folders GR
exp_output_path.GR_review_response      = exp_output_path.exons_data_merge_plots + "GR_review_response/"
exp_output_path.library_scatter         = exp_output_path.GR_review_response + "library_scatter/"
exp_output_path.review_response_ESE     = exp_output_path.GR_review_response + "review_response_ESE/"
exp_output_path.gc_len_subplots         = exp_output_path.GR_review_response + "gc_len_subplots/"
exp_output_path.dual_exons              = exp_output_path.GR_review_response + "dual_exons/"
exp_output_path.HEK293_recovery         = exp_output_path.GR_review_response + 'HEK293_recovery/'









a = vars(exp_output_path)
for path in vars(exp_output_path):
    if a[path] == exp_output_path.git_repo:
        continue
    if isinstance(a[path], str): 
        make_dir_path(a[path])


print('Finished loading experiment paths')









#setup custom run parameters ( threads and parallel instance numbers)
# alignment threads
# bbduk
# spliceAI




# parallel count
#   alignment
#   build exons
#   build vector bam
#   build chr17 bam?





