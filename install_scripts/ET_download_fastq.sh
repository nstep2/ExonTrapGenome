#!/bin/bash


conda create -n ET_fastq-dl -c conda-forge -c bioconda fastq-dl
conda activate ET_fastq-dl 
conda install pysradb


# Read GSM accession numbers from a file
mapfile -t GSMs < GSM_Acc_List.txt

# Function to download a single GSM accession
download() {
    gsm=$1
    counter=$2
    # Convert GSM ID to SRP ID
    srp=$(pysradb gsm-to-srp $gsm | tail -n 1 | cut -f 2)
    echo -e "\033[34mProcessing GSM id $gsm (number $counter of ${#GSMs[@]}), corresponding to SRP id $srp\033[0m"
    fastq-dl --cpus 1 --accession $srp --provider ENA --outdir ../bioproject_PRJNA878769
}

# Loop over GSM accession numbers
for i in "${!GSMs[@]}"; do
    # Start download in the background
    download "${GSMs[$i]}" "$((i+1))" &
    
    # If we have 4 background jobs running, wait until one finishes
    if (( $(jobs -p | wc -l) >= 4 )); then
        wait -n
    fi
done

# Wait for any remaining background jobs to finish
wait

