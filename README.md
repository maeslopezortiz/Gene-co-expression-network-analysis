# Gene co-expression network analysis: Example in Apple
## Description
In this project, our aim is the creation of co-expression networks for gene correlation analysis using public apple RNA datasets. We will examine modular relationships among genes across various functional domains and tissue types.

## Downloading a genome reference
First, we need to select and download the apple reference genome in which we will do the mapping of the raw reads.
> I downloaded the apple genome published by Sun, _et al_.,2020.
> 
>  We need the genome FASTA file, annotation GFF file, and  blast2go file for functional analysis.
```sh
mkdir Genome
wget http://bioinfo.bti.cornell.edu/ftp/Apple_genome/genome/haploid/Gala/Gala_haploid_v2.chr.fa.gz
wget http://bioinfo.bti.cornell.edu/ftp/Apple_genome/genome/haploid/Gala/Gala_haploid_v2.blast2go.gz
wget http://bioinfo.bti.cornell.edu/ftp/Apple_genome/genome/haploid/Gala/Gala_haploid_v2.gff.gz
```
Use gunzip to unzip the files.

It is important to index the genome according the mapper that it will be used. Here, we use STAR.
```sh
#!/bin/sh
# Slurm directives:
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G
#SBATCH --time 1-5:10
#SBATCH --output=logs/STAR.out
#SBATCH --error=logs/STAR.err
#SBATCH --job-name=STAR
# Generate the index genome. Make sure you have downloaded all the necessary files from the repository.
mkdir -p Genome/index
cd Genome || exit  # Change directory or exit if it fails
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/Genome/index --genomeFastaFiles Genome/GalaChrs.fasta --sjdbGTFfile Genome/Gala_haploid_v2.gff --sjdbOverhang 100 --genomeSAindexNbases 12
cd ..
```
## Downloading RNAseq samples from NCBI
You can select which type of date is of your interest. For this example, I select samples from different apple tissues. I create an txt file with the number accession [RNA_samples.txt](https://github.com/maeslopezortiz/Gene-co-expression-network-analysis/blob/main/RNA_samples.txt). You can find the description of each sample in [Samples_NCBI_info.xlsx](https://github.com/maeslopezortiz/Gene-co-expression-network-analysis/blob/main/Samples_NCBI_info.xlsx).

Dowbloading the samples using a loop script:
```sh
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
```

Create directory for FASTQ files if it doesn't exist
```sh
mkdir -p Fastq
```

Move downloaded files to Fastq directory
```sh
mv SRR* Fastq/
```

## STAR alignments
Now we will align FASTQ files using STAR.

```sh
SAMPLES=~/RNA_samples_names.txt 
mkdir -p sorted_bam #create the output directory for the bam files
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
```
It is important to erase unnecessary files to save space.

```sh
rm *.sam
```
The output files that we will use for the transcript quantification are the sorted bam files. 

e. g. **SRR26729830noncanonicalAligned.toTranscriptome.out.bam**

## Salmon quantification
Now, we are going to quantify transcript abundance using Salmon.

```sh
#!/bin/sh
# Slurm directives:
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G
#SBATCH --time 1-5:10
#SBATCH --output=logs/salmon.out
#SBATCH --error=logs/salmon.err
#SBATCH --job-name=salmon
SAMPLES=~/RNA_samples_names.txt 
mkdir -p salmon_quant
cd salmon_quant || exit  # Change directory or exit if it fails
while IFS= read -r RNA_samples
do
    salmon quant -t ~/Genome/transcripts.apple.fa \
                 -l A \
                 -a ~/sorted_bam/"${RNA_samples}noncanonical"Aligned.toTranscriptome.out.bam \
                 -o "${RNA_samples}noncanonical"
done < "$SAMPLES"
```

We need to merge all the quantification files per accession into one file [salmon_merge.txt](https://github.com/maeslopezortiz/Gene-co-expression-network-analysis/blob/main/salmon_merge_all.tsv.zip) to be used as input for WGCNA in R.

```sh
salmon quantmerge --quants SRR*noncanonical -o salmon_merge.txt
```
# WGCNA in R
To generate the modules and cluster of genes, we will work in R with the package [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html). 

This example was created under the example of jashapiro Fix WGCNA url (#477).  Here the link of the original script: https://github.com/AlexsLemonade/refinebio-examples/blob/staging/04-advanced-topics/network-analysis_rnaseq_01_wgcna.Rmd

First download all the libraries:
```sh
require(devtools)
library(dplyr)
library(tidyverse)     
library(magrittr)      
library(WGCNA)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(patchwork)
library(DESeq2)
library(igraph)
```
Then we are using the based on an example of WGCNA. Run the script [apple_WGCNA.R](https://github.com/maeslopezortiz/Gene-co-expression-network-analysis/blob/main/apple_WGCNA.R)

Here is an example of the outputs of the WGCNA analysis.
>From the main WGCNA results, we have a file **apple_wgcna_gene_allgenes_to_module.tsv** from where we can find all the modules that were created. We can choose a specific Module for a downstream analysis and plots.
>

Let’s make a plot of module 27.
![](https://github.com/maeslopezortiz/Gene-co-expression-network-analysis/blob/main/boxplot_me27.png)

Then, we can create a summary heatmap of a given module, here module 27. It is using the full gene expression matrix and by default is `normalized_counts'. 

![](https://github.com/maeslopezortiz/Gene-co-expression-network-analysis/blob/main/heatmap_module27.png)

We can extract this information and the correlation of the important genes which where clasified in the "module 27"
```sh

LISTOFGENES <- gene_module_key[gene_module_key$module == "ME27", 1] #Select genes from gene_module_key where the module is "ME27".
SUBNET <- normalized_counts[, colnames(normalized_counts) %in% LISTOFGENES$gene] #Create a subset of normalized_counts based on the selected genes
CR <- cor((SUBNET)) #Calculate the correlation matrix
DF <- as.data.frame(as.table(CR)) #Convert the correlation matrix into a data frame
M27_Correlations <- DF[DF$Freq > 0.7 | DF$Freq < -0.7,] #Select rows from the data frame where the absolute correlation coefficient is greater than 0.7
write.table(M27_Correlations, file = "FM27corrleations", row.names = FALSE)

```
Here we use a strong filter of a correlation  -0.7 ; 0.7.  

```sh
Var1	Var2	Freq
Mdg_02g014990-mRNA1	Mdg_15g019970-mRNA1	-0.865905603
Mdg_03g005210-mRNA1	Mdg_02g014990-mRNA1	-0.83285841
Mdg_06g019020-mRNA1	Mdg_02g014990-mRNA1	-0.830533544
Mdg_02g014990-mRNA1	Mdg_02g014440-mRNA1	-0.819903242
Mdg_10g028000-mRNA1	Mdg_02g014990-mRNA1	-0.817650655
Mdg_02g014990-mRNA1	Mdg_10g011740-mRNA1	-0.811750223
Mdg_02g014990-mRNA1	Mdg_09g009450-mRNA1	-0.81066779
Mdg_02g014990-mRNA1	Mdg_13g007380-mRNA1	-0.804955294
Mdg_02g014990-mRNA1	Mdg_12g007560-mRNA1	-0.802932795
```
 This new file, FM27corrleations.txt, we can be used as input for Cytoescape.

## Cytoescape visualization

Here we have the visualization network of module 27 created in Cytoescape.

![](https://github.com/maeslopezortiz/Gene-co-expression-network-analysis/blob/main/net_me27.png)
 
## **References**

  Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). https://doi.org/10.1186/1471-2105-9-559
  
  Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
  
  Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T. Cytoscape: a software environment for integrated models of biomolecular interaction networks
  
  Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2024 https://www.ncbi.nlm.nih.gov/pubmed/23104886. Alex Dobin, dobin@cshl.edu. https://github.com/alexdobin/STAR/issues 

  Sun, X., Jiao, C., Schwaninger, H. et al. Phased diploid genome assemblies and pan-genomes provide insights into the genetic history of apple domestication. Nat Genet 52, 1423–1432 (2020). https://doi.org/10.1038/s41588-020-00723-9

 Alex's Lemonade Stand Foundation.Childhood Cancer Data Lab of ALSF.2020. Link: https://github.com/AlexsLemonade
