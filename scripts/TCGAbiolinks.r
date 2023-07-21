#### Script for calculating the significance of gene expression on survival in RNAseq data from TCGA
# Other notes (things to think about)
# - Include Option to colour samples using markers published by other papers 
#   (see clinical_annotation for the columns starting with paper_) ? (this is for PCA)
# - Duplicated samples (patients having multiple samples)
# - Patient treatment (splitting the survival analysis into two to take this into account?)

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

setwd("\\\\lancs/homes/43/shihb/My Documents/r_working")

#### Define query
project <- "TCGA-SKCM"

query_RNASeq <- GDCquery(project = project,          
	data.category = "Transcriptome Profiling" ,          
	data.type = "Gene Expression Quantification" ,          
	workflow.type = "STAR - Counts",
	experimental.strategy = "RNA-Seq" )
RNAseq_samples <- getResults(query_RNASeq)

#### Prepare data from query
GDCdownload(query_RNASeq)
data_RNAseq <- GDCprepare(query_RNASeq)
data_clinic <- GDCquery_clinic(project, "clinical")

# use SummarizedExperiment library to extract relevant information
# https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html#GDCdownload
# tpm_mx <- assay(data_RNAseq, "tpm_unstrand") # tpm matrix; don't need this
count_mx <- assay(data_RNAseq, "unstranded") # gene count matrix
sample_annotation <- colData(data_RNAseq) # clinical annotation data
gene_annotation <- rowData(data_RNAseq) # gene annotation data

# Normalise count matrix 
dgelist <- DGEList(counts=count_mx, genes=gene_annotation)
dgelist <- calcNormFactors(dgelist)
normcounts<-cpm(dgelist, normalized.lib.sizes = T)


### Note (not done yet)
# You might want to consider:
# - Removing control samples (or note down control samples)
# - Removing samples of poor quality
# - Removing duplicated samples from the same patient
tmp_reorganised_log2 <- log2(normcounts + 1)

data_clinic_filt <- data_clinic[, c("bcr_patient_barcode", "vital_status", "days_to_death", "days_to_last_follow_up")]
data_clinic_filt <- data_clinic_filt[data_clinic_filt$bcr_patient_barcode %in% substr(colnames(tmp_reorganised_log2),1,12),]
tmp_reorganised_log2 <- data.matrix(tmp_reorganised_log2[,substr(colnames(tmp_reorganised_log2),1,12) %in% data_clinic_filt$bcr_patient_barcode])

# Finding control samples
group1 <- TCGAquery_SampleTypes(colnames(tmp_reorganised_log2), typesample = c("NT", "NB", "NBC", "NEBV", "NBM"))
group2 <- TCGAquery_SampleTypes(colnames(tmp_reorganised_log2), typesample = c("TP", "TR", "TB", "TAP"))
group3 <- TCGAquery_SampleTypes(colnames(tmp_reorganised_log2), typesample = c("TRBM", "TM", "TB", "TAP")) 



#### Calculate the survival for all genes
tabSurvKMcomplete <- NULL
tokenStop <- 1
for( i in 1: round(nrow(tmp_reorganised_log2)/100)){
    message( paste( i, "of ", round(nrow(tmp_reorganised_log2)/100)))
    # run the analysis on batches of 100 genes
    tokenStart <- tokenStop
	tokenStop <- 100 * i
    tabSurvKM <- TCGAanalyze_SurvivalKM(
        data_clinic_filt,
        tmp_reorganised_log2,
        Genelist = rownames(tmp_reorganised_log2)[tokenStart:tokenStop],
        Survresult = F,
		p.cut = 0.2,
        ThreshTop = 0.67,
        ThreshDown = 0.33
    )
    
    tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)
}

tabSurvKMcomplete <- tabSurvKMcomplete[tabSurvKMcomplete$pvalue < 0.01,]
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]


#### Snippets of code that might be useful
# This currently doesn't work.
# input_gene_list would be the list of genes of interest, 
# but depending on the user input (ensembl ID or gene name), this would need to be processed
tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
    rownames(tabSurvKMcomplete) %in% input_gene_list,
]


# Example of how to analyse survival using clinical annotation/annotation 
# i.e. you can add the annotation for low, medium, high for gene expression onto data_clinic, 
# then perform survial analysis
test <- TCGAanalyze_survival(
data = data_clinic,
clusterCol = "ajcc_pathologic_stage",
main = "TCGA Set\n GBM",
height = 10,
width=10
)


library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)

for(gene in c("GLIS1", "GLIS2", "GLIS3")){
#for(gene in "RNF169"){
	gene_ensembl <- gene_annotation$gene_id[gene_annotation$gene_name == gene]

	plot_df <- data.frame(gene= as.numeric(as.data.frame(normcounts)[gene_ensembl,]), 
		sample_type = sample_annotation$sample_type, 
		ajcc_pathologic_stage = sample_annotation$ajcc_pathologic_stage,
		ajcc_pathologic_t = sample_annotation$ajcc_pathologic_t,
		ajcc_pathologic_m = sample_annotation$ajcc_pathologic_m,
		ajcc_pathologic_n = sample_annotation$ajcc_pathologic_n
		)
	if(TRUE){
		p1 <- ggplot(data = plot_df, aes(y=gene, x=sample_type)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw() + scale_y_log10()
		p2 <- ggplot(data = plot_df, aes(y=gene, x=ajcc_pathologic_stage)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw()  + scale_y_log10() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		p3 <- ggplot(data = plot_df, aes(y=gene, x=ajcc_pathologic_t)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw()  + scale_y_log10()
		p4 <- ggplot(data = plot_df, aes(y=gene, x=ajcc_pathologic_m)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw()  + scale_y_log10()
		p5 <- ggplot(data = plot_df, aes(y=gene, x=ajcc_pathologic_n)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw()  + scale_y_log10()
	} else {
		p1 <- ggplot(data = plot_df, aes(y=gene, x=sample_type)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw()
		p2 <- ggplot(data = plot_df, aes(y=gene, x=ajcc_pathologic_stage)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		p3 <- ggplot(data = plot_df, aes(y=gene, x=ajcc_pathologic_t)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw()  
		p4 <- ggplot(data = plot_df, aes(y=gene, x=ajcc_pathologic_m)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw() 
		p5 <- ggplot(data = plot_df, aes(y=gene, x=ajcc_pathologic_n)) + geom_boxplot(outlier.shape=NA) + geom_point(position=position_jitter(w=0.2)) + theme_bw()  
	}
	g <- arrangeGrob(p1, p2, p3, p4,p5,
				 ncol = 2, nrow = 3, top=text_grob(gene))
	ggsave(paste0(gene, ".png"), g)
	
}