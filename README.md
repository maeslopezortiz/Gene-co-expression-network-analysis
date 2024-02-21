# Gene co-expression network analysis: Example in Apple
## Description
In this project, our aim is the creation of co-expression networks for gene correlation analysis using public apple RNA datasets. We will examinate modular relationships among genes across various functional domains and tissue types.

## Downloading a genome reference
First, we need to select and download the apple reference geneome in which we will do the mapping of the raw reads.
> I downloaded the apple genome published by Sun, _et al_.,2020. 
We need the genome FASTA file, annotation gff file, and  blast2go file for functional analysis.
```sh
mkdir Genome
wget http://bioinfo.bti.cornell.edu/ftp/Apple_genome/genome/haploid/Gala/Gala_haploid_v2.chr.fa.gz
wget http://bioinfo.bti.cornell.edu/ftp/Apple_genome/genome/haploid/Gala/Gala_haploid_v2.blast2go.gz
wget http://bioinfo.bti.cornell.edu/ftp/Apple_genome/genome/haploid/Gala/Gala_haploid_v2.gff.gz
```

## **References**

  Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). https://doi.org/10.1186/1471-2105-9-559
  
  Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods.
  
  Shannon P, Markiel A, Ozier O, Baliga NS, Wang JT, Ramage D, Amin N, Schwikowski B, Ideker T. Cytoscape: a software environment for integrated models of biomolecular interaction networks
  
  Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2024 https://www.ncbi.nlm.nih.gov/pubmed/23104886. Alex Dobin, dobin@cshl.edu. https://github.com/alexdobin/STAR/issues 

  Sun, X., Jiao, C., Schwaninger, H. et al. Phased diploid genome assemblies and pan-genomes provide insights into the genetic history of apple domestication. Nat Genet 52, 1423–1432 (2020). https://doi.org/10.1038/s41588-020-00723-9
