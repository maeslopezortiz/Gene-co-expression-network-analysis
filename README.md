# Gene co-expression network analysis: Example in Apple
## Description
In this project, our aim is the creation of co-expression networks for gene correlation analysis using public apple RNA datasets. We will examinate modular relationships among genes across various functional domains and tissue types.

## Downloading a genome reference
First, we need to select and download the apple reference geneome in which we will do the mapping of the raw reads.
> I downloaded the apple genome published by Sun, _et al_.,2020.
> 
>  We need the genome FASTA file, annotation gff file, and  blast2go file for functional analysis.
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

## **References**

  Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). https://doi.org/10.1186/1471-2105-9-559
  
  Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
  
  Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T. Cytoscape: a software environment for integrated models of biomolecular interaction networks
  
  Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2024 https://www.ncbi.nlm.nih.gov/pubmed/23104886. Alex Dobin, dobin@cshl.edu. https://github.com/alexdobin/STAR/issues 

  Sun, X., Jiao, C., Schwaninger, H. et al. Phased diploid genome assemblies and pan-genomes provide insights into the genetic history of apple domestication. Nat Genet 52, 1423–1432 (2020). https://doi.org/10.1038/s41588-020-00723-9
