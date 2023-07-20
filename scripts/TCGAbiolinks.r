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
# Option to split by other markers? (more for PCA)


tmp_reorganised_log2 <- log2(normcounts + 1)

data_clinic_filt <- data_clinic[, c("bcr_patient_barcode", "vital_status", "days_to_death", "days_to_last_follow_up")]
data_clinic_filt <- data_clinic_filt[data_clinic_filt$bcr_patient_barcode %in% substr(colnames(tmp_reorganised),1,12),]
tmp_reorganised_log2 <- data.matrix(tmp_reorganised_log2[,substr(colnames(tmp_reorganised_log2),1,12) %in% data_clinic_filt$bcr_patient_barcode])

group1 <- TCGAquery_SampleTypes(colnames(tmp_reorganised_log2), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(tmp_reorganised_log2), typesample = c("TP"))

tokenStop <- 1

tabSurvKMcomplete <- NULL

for( i in 1: round(nrow(tmp_reorganised_log2)/100)){
    message( paste( i, "of ", round(nrow(tmp_reorganised_log2)/100)))
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

tabSurvKMcompleteDEGs <- tabSurvKMcomplete[
    rownames(tabSurvKMcomplete) %in% dataDEGsFiltLevel$mRNA,
]









test <- TCGAanalyze_survival(
data = data_clinic,
clusterCol = "gender",
main = "TCGA Set\n GBM",
height = 10,
width=10
)




GDCdownload(query = query_RNASeq, method = "api")


####
## Load the required library
library('TCGAbiolinks')
project_name <- "TCGA-ACC"

## Defines the query to the GDC
query <- GDCquery(project = project_name,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts")

## Get metadata matrix
metadata <- query[[1]][[1]]

## Download data using api
GDCdownload(query, method = "api")

## Get main directory where data is stored
main_dir <- file.path("GDCdata", project_name)
## Get file list of downloaded files
file_list <- file.path("GDCdata", project_name,list.files(main_dir,recursive = TRUE)) 

## Read first downloaded to get gene names
test_tab <- read.table(file = file_list[1], sep = '\t', header = TRUE)
## Delete header lines that don't contain usefull information
test_tab <- test_tab[-c(1:4),]
## STAR counts and tpm datasets
tpm_data_frame <- data.frame(test_tab[,1])
count_data_frame <- data.frame(test_tab[,1])

## Append cycle to get the complete matrix
for (i in c(1:length(file_list))) {
  ## Read table
  test_tab <- read.table(file = file_list[i], sep = '\t', header = TRUE)
  ## Delete unwanted lines
  test_tab <- test_tab[-c(1:4),]
  ## Column bind of tpm and counts data
  tpm_data_frame <- cbind(tpm_data_frame, test_tab[,7])
  count_data_frame <- cbind(count_data_frame, test_tab[,4])
  ## Print progres from 0 to 1
  print(i/length(file_list))
}