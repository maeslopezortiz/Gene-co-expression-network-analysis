
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

# using the cowplot theme for ggplot
theme_set(theme_cowplot())
# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"
# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
# Define the file path to the data directory
# Replace with the path of the folder the files will be in
data_dir <- file.path("data", "apple")
# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file
data_file <- file.path(data_dir, "salmon_merge_all.tsv")
data_file
# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
# Replace with the path to your metadata file
metadata_file <- file.path(data_dir, "RNAappleinfo_2.tsv")
metadata_file
# Check if the gene expression matrix file is at the path stored in `data_file`
file.exists(data_file)
# Check if the metadata file is at the file path stored in `metadata_file`
file.exists(metadata_file)


# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)
# Read in data TSV file
df <- readr::read_tsv(data_file) %>%
  # Here we are going to store the gene IDs as row names so that we can have a numeric matrix to perform calculations on later
  tibble::column_to_rownames("Name")
# Make the data in the order of the metadata
df <- df %>%
  dplyr::select(metadata$Run)
# Check if this is in the same order
all.equal(colnames(df), metadata$Run)

# The next DESeq2 functions need the values to be converted to integers
df <- round(df) %>%
  # The next steps require a data frame and round() returns a matrix
  as.data.frame() %>%
  # Only keep rows that have total counts above the cutoff
  dplyr::filter(rowSums(.) >= 100)

#metadata <- metadata %>%
#  dplyr::mutate(
#    time_point = dplyr::case_when(
#      # Create our new variable based on refinebio_title containing AV/CV
#      stringr::str_detect(refinebio_title, "_AV_") ~ "acute illness",
#      stringr::str_detect(refinebio_title, "_CV_") ~ "recovering"
#    ),
#    # It's easier for future items if this is already set up as a factor
#    time_point = as.factor(time_point)
#  )
#levels(metadata$time_point)
#levels(metadata$time_point)
#levels(metadata$Tissue)

# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = df, # Our prepped data frame with counts
  colData = metadata, # Data frame with annotation for our samples
  design = ~1 # Here we are not specifying a model
)

# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)
#check for outliers here because it can affect WGCNA results.
#this example skip this step (there are no obvious outliers).

# Retrieve the normalized data from the `DESeqDataSet`
normalized_counts <- assay(dds_norm) %>%
  t() # Transpose this data

sft <- pickSoftThreshold(normalized_counts,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)



sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

##So using the plot above, going with a power soft-threshold of 7!
###A 4GB standard desktop or a laptop may handle up to 8000-10000 probes, depending on operating system and other running programs.
#Run WGCNA!
cor=WGCNA::cor
bwnet <- blockwiseModules(normalized_counts,
                          maxBlockSize = 20000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = 12, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 1234, # there's some randomness associated with this calculation
                          # so we should set a seed
)

#Write main WGCNA results object to file
readr::write_rds(bwnet,
                 file = file.path("results", "apple1_wgcna_results_all.RDS")
)
#Explore our WGCNA results
module_eigengenes <- bwnet$MEs

# Print out a preview
head(module_eigengenes)
#Which modules have biggest differences across treatment groups?

all.equal(metadata$Run, rownames(module_eigengenes))

# Create the design matrix from the `tissue` variable
des_mat <- model.matrix(~ metadata$Tissue)


# lmFit() needs a transposed version of the matrix
fit <- limma::lmFit(t(module_eigengenes), design = des_mat)

# Apply empirical Bayes to smooth standard errors
fit <- limma::eBayes(fit)

# Apply multiple testing correction and obtain stats
stats_df <- limma::topTable(fit, number = ncol(module_eigengenes)) %>%
  tibble::rownames_to_column("module")
head(stats_df)


#Let’s make plot of module 26
module_df <- module_eigengenes %>%
  tibble::rownames_to_column("accession_code") %>%
  # Here we are performing an inner join with a subset of metadata
  dplyr::inner_join(metadata %>%
                      dplyr::select(Run, Tissue),
                    by = c("accession_code" = "Run")
  )

#Now we are ready for plotting.
ggplot(
  module_df,
  aes(
    x = Tissue,
    y = ME27,
    color = Tissue
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()+ theme(axis.text.x=element_text(color = "black", size=8, angle=90))



#What genes are a part of module 0?
gene_module_key <- tibble::enframe(bwnet$colors, name = "gene", value = "module") %>%
  # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
  dplyr::mutate(module = paste0("ME", module))

gene_module_key %>%
  dplyr::filter(module == "ME27")

readr::write_tsv(gene_module_key, file = file.path("results", "apple_wgcna_gene_allgenes_to_module.tsv")
)
readr::write_tsv(gene_module_key %>%
                   dplyr::filter(module == "ME27"), file = file.path("results", "apple1_M27_wgcna_gene_all_to_module.tsv")
)


#Let’s save this gene to module key to a TSV file for future use.
readr::write_tsv(gene_module_key,
                 file = file.path("results", "apple1_wgcna_gene_to_module26.tsv")
)

#Make a custom heatmap function
make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = metadata,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with refinebio_accession_code and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its refinebio_accession_code
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("Run")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    # Only select the treatment and sample ID columns
    dplyr::select(Run, Tissue) %>%
    # Add on the eigengene expression by joining with sample IDs
    dplyr::inner_join(module_eigengene, by = "Run") %>%
    # Arrange by patient and time point
    dplyr::arrange(Tissue) %>%
    # Store sample
    tibble::column_to_rownames("Run")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    # Supply treatment labels
    Tissue = col_annot_df$Tissue,
    # Add annotation barplot
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name))
    # Pick colors for each experimental group in time_point
    #col = list(Tissue = c("anthers" = "#f1a340", "apple_fruit" = "#998ec3","filaments" = "#408ef1","ovaries" = "#f14b40","petals" = "#40f13a","pollen" = "#e7f140","receptacles" = "#564ec3","sepals" = "#f15b40","apple_roots_infect"="#f15","apple_roots_mock"="#f99","bud_stage"="#408","flower" ="#500","fruit_skin"="#998", "leaf"="#100ec2", "leaf_N"="#990ec2", "leaf_phospho"="#99ec22","root"="#881","root_N"="#212","root_phospho"="#111ec2","stigmas"="#900" ))
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    # Only keep genes from this module
    dplyr::filter(rownames(.) %in% module_genes) %>%
    # Order the samples to match col_annot_df
    dplyr::select(rownames(col_annot_df)) %>%
    # Data needs to be a matrix
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    # Scale can work on matrices, but it does it by column so we will need to
    # transpose first
    t() %>%
    scale() %>%
    # And now we need to transpose back
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     # Supply color function
                                     col = color_func,
                                     # Supply column annotation
                                     bottom_annotation = col_annot,
                                     # We don't want to cluster samples
                                     cluster_columns = FALSE,
                                     # We don't need to show sample or gene labels
                                     show_row_names = FALSE,
                                     show_column_names = FALSE
  )
  
  # Return heatmap
  return(heatmap)
}

#Make module heatmaps
mod_0_heatmap <- make_module_heatmap(module_name = "ME0")
mod_25_heatmap <- make_module_heatmap(module_name = "ME25")
# Print out the plot
mod_0_heatmap
print(mod_27_heatmap)
#We can save this plot to PNG.
pdf(file.path("results", "apple_module_27_heatmap.pdf"))
mod_27_heatmap
dev.off()

#For comparison, let’s try out the custom heatmap function with a different, not differentially expressed module

mod_46_heatmap <- make_module_heatmap(module_name = "ME46")

# Print out the plot
mod_25_heatmap
png(file.path("results", "SRP140558_module_25_heatmap."))
mod_25_heatmap
dev.off()
#Session info
# Print session info
sessioninfo::session_info()

bwnet

# Convert labels to colors for plotting
#mergedColors = labels2colors(bwnet$colors)
# Plot the dendrogram and the module colors underneath
#plotDendroAndColors(
#  bwnet$dendrograms[[1]],
#  mergedColors[bwnet$blockGenes[[1]]],
#  "Module colors",
#  dendroLabels = FALSE,
#  hang = 0.03,
#  addGuide = TRUE,
#  guideHang = 0.05 )


#module_df <- data.frame(
#  gene_id = names(bwnet$colors),
#  colors = labels2colors(bwnet$colors)
#)


module_df[1:5,]
write_delim(module_df,
            file = "gene_modules_apple_F_all.txt",
            delim = "\t")

# pick out a few modules of interest here
modules_of_interest = c("green", "turquoise", "tan")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
#normalized_counts[1:5,1:10]

#genes_of_interest = module_df %>%
#  subset(colors %in% modules_of_interest)
#expr_of_interest = normalized_counts[genes_of_interest$gene_id,]
#expr_of_interest[1:5,1:5]

#table(gene_module_key$module)
#gene_module_key[gene_module_key$module == "ME0",]
#gene_module_key[gene_module_key$module == "ME63",]
#dim(gene_module_key[gene_module_key$module == "ME63",])
#(gene_module_key[gene_module_key$module == "ME0",1])

LISTOFGENES <- (gene_module_key[gene_module_key$module == "ME27",1])

str(normalized_counts)
head(normalized_counts)[1:5]
normalized_counts[1:5,1:5]
normalized_counts[,colnames(normalized_counts) %in% LISTOFGENES]
head(LISTOFGENES)
head(colnames(normalized_counts))
#colnames(normalized_counts) %in% LISTOFGENES

#any(colnames(normalized_counts) %in% LISTOFGENES$gene)
normalized_counts[,colnames(normalized_counts) %in% LISTOFGENES$gene]
SUBNET <- normalized_counts[,colnames(normalized_counts) %in% LISTOFGENES$gene]

CR <- cor((SUBNET))
head(CR)
hist(CR)
DF <- as.data.frame(as.table(CR))
M27_Correlations <- DF[DF$Freq > 0.7 | DF$Freq < -0.7,]
#DFOUT <- DF[DF$Freq > 0.5 | DF$Freq < -0.5,]
write.table(M27_Correlations,file="FM27corrleations",row.names=F)









