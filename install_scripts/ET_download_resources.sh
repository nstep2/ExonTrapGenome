#!/bin/bash -i

#scipt uses pigz to avoid error message for Homo_sapiens.GRCh38.bed.gz


conda env remove --name ET_downloader


# Create the conda environment
conda create --name ET_downloader --yes

# Activate the environment
conda activate ET_downloader

# Check if pigz is installed
which pigz &> /dev/null



if [ $? -eq 0 ]; then
    echo "pigz is already installed."
else
    echo "pigz is not installed. Installing pigz."

    # Install pigz from conda-forge channel
    conda install -c conda-forge pigz --yes
fi

# Check if curl is installed
which curl &> /dev/null

if [ $? -eq 0 ]; then
    echo "curl is already installed."
else
    echo "curl is not installed. Installing curl."

    # Install curl from conda-forge channel
    conda install -c conda-forge curl --yes
fi



# Define the data in an array
data=(
"SpliceAI_release_1.3.1	https://github.com/Illumina/SpliceAI/archive/refs/tags/v1.3.1.zip"
"Dfam		https://dfam.org/releases/Dfam_3.4/annotations/hg38/hg38_dfam.nrph.hits.gz"
"Dfam		https://dfam.org/releases/Dfam_3.4/families/Dfam.embl.gz"
"snaptron		https://snaptron.cs.jhu.edu/data/srav2/exons.sqlite"
"housekeeping_gene_list		https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkaa609#supplementary-data"
"UCSC_hg_38		rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
"hg_38_2bit		rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit"
"lncRNA_databases		https://lncipedia.org/downloads/lncipedia_5_2/high-confidence-set/lncipedia_5_2_hc_hg38.bed"
"lncRNA_databases		https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/9.0/genome_coordinates/Homo_sapiens.GRCh38.bed.gz"
"lncRNA_databases		http://www.noncode.org/datadownload/NONCODEv6_human.fa.gz"
"lncRNA_databases		https://public-docs.crg.es/rguigo/CLS/data/gencodePlus/hs.GENCODE+.lncRNA.bed"
"phyloP		https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP7way/hg38.phyloP7way.bw"
"phyloP	https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons30way/)"
"CHESS		https://github.com/chess-genome/chess/releases/download/v2.2/chess2.2.gff.gz"
"HEK293_inclusion		https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221838/suppl/GSE221838_d0_v_d1000_SE.MATS.JCEC.txt.gz"
"HEK293_expression		https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235387/suppl/GSE235387_HEK293.xlsx"
"gencode_v37		https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gff3.gz"
"IGV		https://data.broadinstitute.org/igv/projects/downloads/2.16/IGV_Linux_2.16.1_WithJava.zip"
"UCSC_hg_38_hisat2_index	https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz"
)


# Function to print processing message
function print_processing_message() {
    echo "File processed: $1"
    echo ""
}

# Function to convert seconds to hh:mm:ss format
function format_time() {
    local seconds="$1"
    local formatted_time=$(date -u -d @"$seconds" +'%H:%M:%S')
    echo "$formatted_time"
}

# Initialize total time counter
total_time=0

# Loop through the data array
for line in "${data[@]}"; do
    # Split the line by tabs
    IFS=$'\t' read -ra parts <<< "$line"
    folder="${parts[0]}"
    link="${parts[1]}"

    # Create folder for column 1
    mkdir -p "$folder"

    # Download the link from column 3 into that folder
    cd "$folder" || exit
    filename=$(basename "$link")

    # Print the filename being downloaded
    echo "Downloading: $filename"

    # Remove the file if it already exists
    if [ -f "$filename" ]; then
        rm "$filename"
    fi

    # Measure the start time for this iteration
    start_time=$(date +%s)

    if [[ "$link" == rsync* ]]; then
        # Download the file using rsync
        rsync -q "$link" "$filename"
    elif [[ "$link" == http* ]]; then
        # Download the file using curl
        curl -L -s -o "$filename" "$link"
    else
        # Download the file using wget
        wget -q -O "$filename" "$link"
    fi

    
    
    # If the link ends in .gz, use pigz to decompress it
    if [[ "$filename" == *.gz ]]; then
        pigz -d -f "$filename"
        echo "unzipping.."
    elif [[ "$filename" == *.zip ]]; then
    	unzip "$filename"
        echo "unzipping.."
    else
    	echo "not unzipping"
    fi
    
    
    
    
    

    # Print the file has been processed
    print_processing_message "$filename"

    # Measure the end time for this iteration
    end_time=$(date +%s)

    # Calculate the time taken for this iteration
    iteration_time=$((end_time - start_time))

    # Add the iteration time to the total time
    total_time=$((total_time + iteration_time))

    # Print current time taken for the entire script
    echo "Current time taken: $(format_time $total_time)"

    cd .. || exit
done

# Print total time taken for the entire script
echo "Total time taken: $(format_time $total_time)"



#cd UCSC_hg_38
#samtools faidx hg38.fa $(seq -f 'chr%.0f' 1 22) chrX chrY chrM > hg38_standard.fa
#hisat2-build hg38_standard.fa hg38_standard









conda deactivate 

conda env remove --name ET_downloader













