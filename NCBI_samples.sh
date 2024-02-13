#!/bin/sh
# Slurm directives:
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G
#SBATCH --time 1-5:10
#SBATCH --output=logs/NCBIdownload.out
#SBATCH --error=logs/NCBIdownload.err
#SBATCH --job-name=NCBIdownloading
# Define the input file with the names of RNAseq samples
SAMPLES=~/RNA_samples_names.txt  # txt file with the names of the RNA sample accessions from NCBI

# Download and convert SRA data to FASTQ format
while read -r RNA_samples
do
    ${HOME}/sratoolkit.3.0.7-centos_linux64/bin/fastq-dump --split-files "$RNA_samples"
done < "$SAMPLES"

# Create directory for FASTQ files if it doesn't exist
mkdir -p Fastq
# Move downloaded files to Fastq directory
mv SRR* Fastq/
