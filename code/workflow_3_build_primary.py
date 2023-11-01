
from experiment_paths.experiment_paths import *

repo_path = ''


#runfile(repo_path + 'build_exon_interval_from_SAM_1.py', wdir = exp_output_path.git_repo, current_namespace=True)

runfile(repo_path + 'build_exon_interval_from_SAM_2.py', wdir = exp_output_path.git_repo, current_namespace=True)


#synthetic make SAM
runfile(repo_path + 'synthetic_make_sam.py', wdir = exp_output_path.git_repo, current_namespace=True)
runfile(repo_path + 'workflow_8_synthetic_collapse.py', wdir = exp_output_path.git_repo, current_namespace=True)


#synthetic evaluate
runfile(repo_path + 'build_exon_interval_from_SAM_3.py', wdir = exp_output_path.git_repo, current_namespace=True)



runfile(repo_path + 'build_exon_interval_from_SAM_4.py', wdir = exp_output_path.git_repo, current_namespace=True)















