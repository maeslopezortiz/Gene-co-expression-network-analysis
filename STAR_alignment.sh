#!/bin/sh
# Slurm directives:
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G
#SBATCH --time 1-5:10
#SBATCH --output=logs/STAR.out
#SBATCH --error=logs/STAR.err
#SBATCH --job-name=STAR
# Generate the index genome. Make sure you have downloaded all the necessary files from the repository. e.g. GFF files
mkdir -p Genome/index
cd Genome || exit  # Change directory or exit if it fails
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/Genome/index --genomeFastaFiles Genome/GalaChrs.fasta --sjdbGTFfile Genome/Gala_haploid_v2.gff --sjdbOverhang 100 --genomeSAindexNbases 12
cd ..

# Align FASTQ files using STAR
mkdir -p sorted_bam
cd sorted_bam || exit  # Change directory or exit if it fails
while IFS= read -r RNA_samples
do
    STAR --genomeDir ~/Genome/index \
         --runThreadN 8 \
         --sjdbGTFfile ~/Genome/Gala_haploid_v2.gtf \
         --readFilesIn ~/Fastq/"$RNA_samples"_1.fastq ~/Fastq/"$RNA_samples"_2.fastq \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMstrandField intronMotif \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI AS NM MD \
         --outFilterIntronMotifs RemoveNoncanonical \
         --outFileNamePrefix "${RNA_samples}noncanonical"
done < "$SAMPLES"

# Erase unnecessary files
rm *.sam
cd ..
