#!/bin/sh
# Slurm directives:
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G
#SBATCH --time 1-5:10
#SBATCH --output=logs/salmon.out
#SBATCH --error=logs/salmon.err
#SBATCH --job-name=salmon
# Quantify transcript abundance using Salmon
mkdir -p salmon_quant
cd salmon_quant || exit  # Change directory or exit if it fails
while IFS= read -r RNA_samples
do
    salmon quant -t ~/Genome/transcripts.apple.fa \
                 -l A \
                 -a ~/sorted_bam/"${RNA_samples}noncanonical"Aligned.toTranscriptome.out.bam \
                 -o "${RNA_samples}noncanonical"
done < "$SAMPLES"

# Merge Salmon quantification results
salmon quantmerge --quants SRR*noncanonical -o salmon_merge.txt
# salmon_merge.txt is the file that will be used as input for WGCNA in R
