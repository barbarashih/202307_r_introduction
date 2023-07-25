#### Section 0: Settings
# !!! Requirements:
# --- Gene count data requirement:
# - gene_id is ensembl gene ID (if not supplying gene annotation)
# - The ensembl gene IDs do not have version number (for example, they do not end with .1 or .2 or .3 and so on) (if not supplying gene annotation)
# - Other than the first row (sample names) and the first column (gene IDs), all other values are count data/numbers 
# --- Sample annotation file requirements: 
# - The first row is your annotation names (i.e. sample_id, tissue, timepoint.. so on), and that these names don't have spaces or special characters (underscores are fine)
# - The first column is sample_id (that matches count file), and subsequent columns are sample annotations.
# - Make sure there is a column named "group" that you would use for quality control/filtering lowly expressed genes. 
# --- Count file requirements: 
# - The first row is sample_id (that matches the sample annotation)
# - The first column is gene_id 
# - Make sure that other than the first column, 
# - You will need to change design and contrast definitions to match your study
# --- Sections where you need to change are labelled with !!! CHANGE !!!
# --- Gene annotation (optional):
# - The first column is the gene id that matches those in count file

# Note that if a file is in a subfolder in your working directory, you can use / to indicate the subfolder
# For instance, if my count_mx.csv file is in a folder "data" in my working directory, 
# then I would refer to it as "data/count_mx.csv"
working_directory <- "C:/Users/shihb/AppData/Local/r_working"     # !!! CHANGE !!! Change to the folder you want to work from (where your files are)
species <- "mouse"												                        # !!! CHANGE !!! "mouse" or "human"
count_filepath <- "count_mx.csv" 					                  	    # !!! CHANGE !!! The gene count matrix (first column gene_id, first row sample_id, everything else numbers)
sample_annotation_filepath <-  "sample_annotation.csv"         		# !!! CHANGE !!! 
gene_annotation_filepath <- ""									                  # !!! CHANGE !!! Leave this as "" if you don't have genome annotation
gene_id_type <- "entrezgene"                                      # !!! CHANGE !!! choose between entrezgene or ensembl

setwd(working_directory)										# set your working directory



#### Section 1: Download library
# install packagws if you dont already have them
# Install BiocManager
if (!require("BiocManager", quietly = TRUE)){ install.packages("BiocManager") }
# Install limma
if (!require("limma", quietly = TRUE)){ BiocManager::install("limma", update = TRUE) }
# Install Glimma
if (!require("Glimma", quietly = TRUE)){ BiocManager::install("Glimma") }
# Install edgeR
if (!require("edgeR", quietly = TRUE)){ BiocManager::install("edgeR") }
# Install gplots for plotting heatmap
if (!require("gplots", quietly = TRUE)){ install.packages("gplots") }	
# Install biomaRt for gene annotation
if (!require("biomaRt", quietly = TRUE)){ BiocManager::install("biomaRt") }	
# Install R.utils
if (!require("R.utils", quietly = TRUE)){ install.packages("R.utils") }	


#### Section 2: Load library
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(R.utils)
library(RColorBrewer)
library(gplots)
library(biomaRt)


#### Section 3: Read in csv files for your study and make a DEGList object
# Read in files
count_mx <- read.csv(count_filepath, row.names=1)
sample_annotation <- read.csv(sample_annotation_filepath)

# Download gene annotation if you don't already have one
if(nchar(gene_annotation_filepath) == 0){
	biomart_dataset <- ifelse(tolower(species) == "mouse", "mmusculus_gene_ensembl", "hsapiens_gene_ensembl")
	ensembl <- useEnsembl(biomart = "genes", dataset = biomart_dataset)
	if(gene_id_type == "entrezgene"){
	  gene_annotation <- getBM(attributes = c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol', 'description', 'chromosome_name'), mart = ensembl)
	} else {
	  gene_annotation <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description', 'chromosome_name'), mart = ensembl)
	}
	colnames(gene_annotation)[1] <- "gene_id"
# if there is a input gene_annotation_filepath, read in the file and set the first column name as gene_id
} else {
	gene_annotation <- read.csv(gene_annotation_filepath)
	colnames(gene_annotation)[1] <- "gene_id"
}
gene_annotation <- merge(data.frame(gene_id = row.names(count_mx)), gene_annotation, by.x="gene_id", by.y=colnames(gene_annotation)[1], all.x=TRUE)
gene_annotation <- gene_annotation[!duplicated(gene_annotation$gene_id),]
  
# Make DGElist object
x <- DGEList(counts = count_mx, genes = gene_annotation, group = sample_annotation$group, samples= sample_annotation)


#### Section 4: Remove lowly expressed genes
## Keep a record of the log(count per million) data before filtering the genes
lcpm <- cpm(x, log=TRUE)


## Filter out lowly expressed genes
keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]

## Look at data distribution before and after filtering
# organise plot data
nsamples <- ncol(x)							# Total number of samples
L <- mean(x$samples$lib.size) * 1e-6  		# Mean library size (sum of all read counts in each sample)
M <- median(x$samples$lib.size) * 1e-6  	# Median library size (sum of all read counts in each sample)
lcpm_filtered <- cpm(x, log=TRUE)			# Store the log(count per million) for filtered data
lcpm_cutoff <- log2(10/M + 2/L) 			# Roughly around the cutoff for lowly expressed genes; for drawing a vertical line in the graph

# Sort out plotting colours/parameters
col_vec <- brewer.pal(12, "Paired") 										# generate 12 colours
col_vec_plot <- rep(col_vec, ceiling(length(nsamples)/length(col_vec)))		# Repeat col_vec until col_list_plot has a length longer than the number of samples
col_vec_plot <- col_vec_plot[1:nsamples]									# Trim this vector back down so it's the same length as number of samples

par(mfrow=c(1,2)) 																		# Split the plotting window into 1 row and 2 columns

## Density plot for the raw data
plot(density(lcpm[,1]), col=col_vec_plot[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")		# Plot the first sample to create a graph
title(main="A. Raw data", xlab="Log-cpm")												# Label the plot
abline(v=lcpm_cutoff, lty=3)															# Add a vertical line
# Use a loop to plot the remaining samples (i.e. 2nd to last sample)
for (i in 2:nsamples){
den <- density(lcpm[,i])																# Generate data for density plot
lines(den$x, den$y, col=col_vec_plot[i], lwd=2)											# Plot the second to last sample
}
legend("topright", colnames(x), text.col=col_vec_plot, bty="n")							# Add a legend for the samples

## Density plot for after removing lowly expressed genes (similar to above, but change lcpm to lcpm_filtered)
plot(density(lcpm_filtered[,1]), col=col_vec_plot[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")				# Plot the first sample to create a graph
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm_cutoff, lty=3)
for (i in 2:nsamples){
den <- density(lcpm_filtered[,i])
lines(den$x, den$y, col=col_vec_plot[i], lwd=2)
}
legend("topright", colnames(x), text.col=col_vec_plot, bty="n")



#### Section 5: Normalising gene distribution
# In the tutorial, it created x2 to generate a modified dataset where there are two unusual samples that have wildly different distribution
# You don't need to do that for actually running your analysis
# Just visualise your normalised data and make sure the samples have roughly similar distributions after normalisation
# (note this will probably not be the case when you're analysing metabolomic or proteomic data due to the large number of missing values)
x <- calcNormFactors(x, method = "TMM")				# Normalise data
par(mfrow=c(1,1)) 									# Split the plotting window into 1 row and 2 columns
boxplot(lcpm_filtered, las=2, col=col_vec_plot)
title(main="B. Example: Normalised data",ylab="Log-cpm")


#### Section 6: Unsupervised clustering
# Use Glimma to interactively explore through the MD plots
glMDSPlot(lcpm_filtered, labels=colnames(x),  
	groups=x$samples[, colnames(x$samples)[!colnames(x$samples) %in% c("files", "norm.factors")]], 
	launch=TRUE) 


#### Section 7: Differential gene expression 
# !!! You need to edit this section according to your study
# Carefully read through the comments to make sure you the relevant variables in each line
design <- model.matrix(~0+group+lane)						# !!!! CHANGE !!!! Change group and lane to your variables of interest (i.e. tissue, timepoint). Colon indicates interaction term
#design <- model.matrix(~0+group+lane+group:lane)			# Colon indicates interaction term and might be relevent if for example you want to compare brain Day7 vs brain Day0
colnames(design) <- gsub("group", "", colnames(design))		# !!!! CHANGE !!!! Do this for the variables you have included above
colnames(design) <- gsub("lane", "", colnames(design))		# !!!! CHANGE !!!! Do this for the variables you have included above
head(design)												# !!!! CHANGE !!!! Use the column names in design for making contrast fit (below)

contrast.matrix <- makeContrasts(							# !!!! CHANGE !!!! change this to what your comparisons are
   BasalvsLP = Basal-LP, 									# !!!! CHANGE !!!! BasalvsLP can be changed to a name that you like. Change Basal - LP to the conditions you are interested in. Make sure these correspond to the column names in design
   BasalvsML = Basal - ML, 									# !!!! CHANGE !!!! Same as above. You can add on as many additional comparisons as you want
   LPvsML = LP - ML, 										# !!!! CHANGE !!!! Same as above. You can add on as many additional comparisons as you want
   levels = colnames(design))
contrast.matrix

# Voom converts raw counts to log-CPM while taking into account of the library size and variances
# Variances are then rescaled to quarter-root variances (or square-root of standard deviations) and plotted against the average log2 count for each gene
# This is so that limma linear modelling (statistics) can then be performed 
# This conversion is necessary because the limma statistical method assumes normally distributed log-CPM regardless the level of gene expression, 
# but typically in RNAseq data, lower expressing genes show a higher variance
v <- voom(x, design, plot=TRUE)

# These graphs would help you determine if the level of lowly-expressed gene filtering in the earlier step is sufficient 
# If the lowly expressed gene filtering is insufficient, you would typtically see a drop in variance level at the low end of the expression scale
# You would typically want to see a gradually decreasing trend that flattens at high log2(count size + 0.5) for the "voom: Mean-variance trend"
# If you see a "n" shaped graph in "voom: Mean-variance trend", you probably didn't filter out enough lowly expressed genes
# You would typically want to see a flat-ish trend for the "Final model: Mean-variance trend" graph
par(mfrow=c(1,2))
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrast.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")


#### Section 8: Explore the numbers of differentially expressed genes
# this sums up the total number of significantly differentially expressed genes
summary(decideTests(efit))						
# this sums up the total number of significantly differentially expressed genes, but is a bit more stringent
# treat is used to calculate p values with a minimum log-fold-change requirement. The log in limma refer to log2
tfit <- treat(vfit, lfc=1)									# !!! CHANGE !!! change lfc (log-fold-change) to other values if you want a slightly different threshold for log fold change
summary(decideTests(efit))						
dt <- decideTests(tfit)

# Find the proportion of commonly differentially expressed genes across your comparisons
# Here is an example comparing the first and second comparisons you stated in contrast.matrix
de.common <- which(dt[,1]!=0 & dt[,2]!=0)					# !!!! CHANGE !!!!  dt[,1] and dt[,2] refers to the first and second comparisons you stated in contrast.matrix. Change accordingly
length(de.common)											# Number of overlapping genes
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))	# !!!! CHANGE !!!!  dt[,1:2] refers to the first and second comparisons you stated in contrast.matrix. Change accordingly


#### Section 9: Explore the differentially expressed genes
statsResutls1 <- topTreat(tfit, coef=1, n=Inf)				# Get differentially expressed genes for the first comparison (coef refers to the comparison number in contrast matrix
statsResutls2 <- topTreat(tfit, coef=2, n=Inf)				# Get differentially expressed genes for the second comparison (coef refers to the comparison number in contrast matrix
head(statsResutls1)											# Look at the top of the stats table
head(statsResutls2)											# Look at the top of the stats table

# Save stats table
write.csv(statsResutls1, "stats_table1.csv")          # repeat the same for each coefficient

#### Explore the differentially expressed genes in MDplot
gene_id_name <- "gene_id"									# you can change this to a different column name in your gene_annotation
contrast_mx_comp_coef <- 1									# !!! CHANGE !!! change contrast_mx_comp_coef to the comparison ID from contrast.matrix (i.e. 1 for first comparison, 2 for second... so on)
glMDPlot(tfit, coef=contrast_mx_comp_coef, status=dt, main=colnames(tfit)[contrast_mx_comp_coef],
         side.main=gene_id_name, counts=lcpm_filtered, groups=x$samples$group, launch=TRUE)
		 

#### Section 9: Heatmap for top significant genes
gene_plot_name <- "SYMBOL"										# !!! CHANGE !!! change this to your the column name of gene_annotation that you want to use for the plot
comparison1_top_genes <- row.names(statsResutls1)[1:100]		# Change 100 to the number of top genes you want to plot
i <- which(v$genes$gene_id %in% comparison1_top_genes)			# find out the indexes for the genes of interest 
mycol <- colorpanel(1000,"blue","white","red")					# generate plot colours
heatmap.2(lcpm[i,], scale="row",								# Plot heatmap 
   labRow=v$genes[[gene_plot_name]][i], labCol=x$samples$group, 
   col=mycol, trace="none", density.info="none", 
   margin=c(8,6), lhei=c(2,10), dendrogram="column")
   
   