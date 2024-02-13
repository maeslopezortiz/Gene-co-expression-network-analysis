#!/bin/sh
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
