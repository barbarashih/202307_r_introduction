#### Script for calculating the significance of gene expression on survival in RNAseq data from TCGA
# Other notes
# Include Option to colour samples using markers published by other papers 
# (see clinical_annotation for the columns starting with paper_) ? (this is for PCA)

library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:/Users/bshih/Documents/tmp/nwcr")

project <- "TCGA-SKCM"

query_RNASeq <- GDCquery(project = project,          
	data.category = "Transcriptome Profiling" ,          
	data.type = "Gene Expression Quantification" ,          
	workflow.type = "STAR - Counts",
	experimental.strategy = "RNA-Seq" )
RNAseq_samples <- getResults(query_RNASeq)

data_RNAseq <- GDCprepare(query_RNASeq)
data_clinic <- GDCquery_clinic(project, "clinical")


# use SummarizedExperiment to look at the data
# https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html#GDCdownload
count_mx <- assay(data_RNAseq, "unstranded") # count matrix
tpm_mx <- assay(data_RNAseq, "tpm_unstrand") # tpm matrix
clinical_annotation <- colData(data_RNAseq)
gene_annotation <- rowData(data_RNAseq)

# Normalise count matrix 
dgelist <- DGEList(counts=count_mx, genes=gene_annotation)
dgelist <- calcNormFactors(dgelist)
normcounts<-cpm(dgelist, normalized.lib.sizes = T)


### Note
# Remove control samples (or note down control samples)
tmp_reorganised_log2 <- log2(normcounts + 1)

data_clinic_filt <- data_clinic[, c("bcr_patient_barcode", "vital_status", "days_to_death", "days_to_last_follow_up")]
data_clinic_filt <- data_clinic_filt[data_clinic_filt$bcr_patient_barcode %in% substr(colnames(tmp_reorganised),1,12),]
tmp_reorganised_log2 <- data.matrix(tmp_reorganised_log2[,substr(colnames(tmp_reorganised_log2),1,12) %in% data_clinic_filt$bcr_patient_barcode])

# Finding control samples
group1 <- TCGAquery_SampleTypes(colnames(tmp_reorganised_log2), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(tmp_reorganised_log2), typesample = c("TP"))



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

#
tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
    rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
]




# Example of how to analyse survival using clinical annotation/annotation 
# i.e. you can add the annotation for low, medium, high for gene expression onto data_clinic, 
# then perform survial analysis
test <- TCGAanalyze_survival(
data = data_clinic,
clusterCol = "gender",
main = "TCGA Set\n GBM",
height = 10,
width=10
)


