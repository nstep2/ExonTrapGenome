# Introduction

This repository contains scripts associated with the publication "The human genome contains over a million autonomous exons." The work focuses on querying the human genome for splicing potential using deep sequencing of randomly fragmented genomic DNA assayed in an exon trapping reporter assay. The scripts reproduce the analyses presented in the publication.

## Getting Started
**Getting Started**

To set up and run the molecular biology research project scripts, follow these steps:

1. **Install Conda:**
   - If not already installed, follow the instructions to install Conda on your system.

2. **Download the Script:**
   - Download the script `install_ET.py` to a directory where you want to install the scripts.

3. **Run the Installation Script:**
   - Execute `install_ET.py` to create the `ET_directory`. This will also download the necessary analysis scripts.

4. **Navigate to the Installation Scripts Directory:**
   - Change your working directory to the `install_scripts` folder.

5. **Create Conda Environment:**
   - Run `ET_create_env.sh` to create a Conda environment named `ET_env`.

6. **Download Sequence Data:**
   - Execute `ET_download_fastq.sh` to download the required data files into the `ET_seq_reads` folder. Note that this may require approximately 700 GB of available disk space and may take several hours, as data will be fetched from ENA (European Nucleotide Archive).

7. **Download Additional Resources:**
   - Run `ET_download_resources.sh` to download many of the cited resources needed by the scripts.

8. **Activate Conda Environment:**
   - Activate the Conda environment `ET_env`.

9. **Initialize the Project:**
   - Run the Python script `ET_initialize.py` to set up the project for analysis.

Following these steps will ensure that you have all the necessary files and environments in place to run the analyses for your molecular biology research project.



### System Requirements

- The analysis was performed on an Ubuntu PC workstation with 64 GB of RAM. The runtime of the scripts is several days.

### Acquiring Sequencing Data

- The associated sequencing data can be acquired from NCBI GEO Series GSE213006. This series requires approximately 700 GB of disk space for storage. 

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
```


### Conda Environment

## Setup

Run the following commands to set up the `ET_env` Conda environment and activate it:

```bash
conda env create -f supporting_files/ET_env.yml
conda activate ET_env

- Additionally, some packages are installed using pip or apt-get, so please follow the installation instructions provided in the repository.


