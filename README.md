# Introduction

This repository contains scripts associated with the publication "The human genome contains over a million autonomous exons." The paper focuses on querying the human genome for splicing potential using deep sequencing of randomly fragmented genomic DNA assayed in an exon trapping reporter assay. The scripts reproduce the analyses presented in the publication.

## Getting Started

### System Requirements

- The analysis was performed on an Ubuntu PC workstation with 64 GB of RAM. The runtime of the scripts is several days.

### Acquiring Sequencing Data

- The associated sequencing data can be acquired from NCBI GEO Series GSE213006. This series requires approximately 700 GB of disk space for storage. Ensure you have sufficient disk space to accommodate this dataset.

### Dependency Setup

- The scripts have dependencies on external resources. Most of these resources can be downloaded using a provided script, `download_inputs.py`. However, some resources may require manual downloading and parsing. Be sure to follow the instructions provided in the repository.

### Running the Scripts

- The scripts are primarily designed to run in a Spyder IDE session using a specific Conda environment as the Python interpreter. Make sure you have this environment set up.
- Enable the 'run in console's namespace' setting in Spyder.
- Perform the analyses by running the scripts with file names starting with "workflow_." These scripts should be run sequentially, with 'workflow_1.py' being run first.

### Configuration

- Before running the scripts, you will need to manually edit a Python script to specify the paths of the repository, the location of the sequencing reads, and the analysis directory.

## Setup

Run the following commands to set up the `ET_env` Conda environment and activate it:

```bash
conda env create -f supporting_files/ET_env.yml
conda activate ET_env


### Conda Environment

- These analyses depend on the Conda package management system. Refer to the README file in the repository for detailed instructions on creating the necessary Conda environment.
- Additionally, some packages are installed using pip or apt-get, so please follow the installation instructions provided in the repository.


