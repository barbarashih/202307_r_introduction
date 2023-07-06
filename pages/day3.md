# Introduction to R (Day 3)

## 1. List + function. Apply, sapply, lapply
### Task 1.1 Download data from GTEx
Go to https://gtexportal.org/home/datasets. Under GTEx Analysis V8 > Annotations, download:
1) GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt 

2) GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

Copy the files to your working directory.

### Task 1.2. Read in these files
```r
# set working directory
working_dir <- "C:/Users/shihb/OneDrive - Lancaster University/work/teaching/workshop/202307_r_introduction"
setwd(working_dir)

# Read in metadata
sample_annotation_fp <- "data/day2/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
subject_annotation_fp <- "data/day2/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
sample_annotation <- read.delim(sample_annotation_fp)
subject_annotation <- read.delim(subject_annotation_fp)

# Look at the top few rows for each dataframe
head(sample_annotation)
head(subject_annotation)

# In order to merge the two columns, We need to have a column in sample_annotation that matches to SUBJID in subject_annotation
# Splitting a string to generate list of character vectors
sampleID_split <- strsplit(sample_annotation$SAMPID, "-")
head(sampleID_split)

# For loop
# declare an empty vector
combined <- vector()
for(idx in 1:length(sampleID_split)){
    current_strings <- sampleID_split[[idx]]
    combined[idx] <- paste0(current_strings[1], "-", current_strings[2])
}
```

### Task 1.3
```r
# sapply (takes vector or lists as inputs and returns vectors)
combined <- sapply(sampleID_split, FUN=function(x)paste0(x[1], "-", x[2]))
combined <- sapply(strsplit(sample_annotation$SAMPID, "-"), FUN=function(x)paste0(x[1], "-", x[2]))

# lapply (takes vector or lists as inputs returns lists)
combined_lapply <- lapply(sampleID_split, FUN=function(x)paste0(x[1], "-", x[2]))

# apply (takes data frame or matrix as inputs and return vectors, lists or array)
head(sample_annotation[,1:5])
combined_apply <- apply(sample_annotation, MARGIN=1, FUN=function(x)paste0(strsplit(x[1], "-")[[1]][1:2], collapse="-"))
head(combined_apply)

```

### Task 1.4. Merge sample_annotation with subject_annotation
```r
# Now you can add a SUBJID column to sample_annotation 
sample_annotation$SUBJID <- combined_apply
str(sample_annotation)
sample_annotation <- merge(sample_annotation, subject_annotation, by.x="SUBJID", by.y="SUBJID")
str(sample_annotation)

```


<br> </br>


## 2. Libraries
### Task 2.1. Install a library
```r
# Other people have made libraries of useful functions
# Installing libraries
install.packages("ggplot2")
install.packages("gridExtra")
```

### Task 2.2. Use ggplot2 library
```r
# Importing libraries
library(ggplot2)
library(grid)
library(gridExtra)

plot_df <- data.frame(Gene1=1:10, Gene2=10:1, Group=c(rep("Group1", 5), rep("Group2", 5)))
plot_df

p1 <- ggplot(data = plot_df, aes(x=Group, y=Gene1)) + geom_point()
p1

p2 <- ggplot(data = plot_df, aes(x=Group, y=Gene2)) + geom_point() + ggtitle("title") + theme_bw()
p2



g <- arrangeGrob(p1, p2, nrow=1)
grid.draw(g)
ggsave("output/plot_d2_t5-2.pdf", g)

```


## 3. Putting it all together
## Aim to do
1. Import the GTEx gene expression data
1. Annotate the gene expression data with sample and subject/gender annotation
1. Create a function that can be used to plot a gene
1. Plot 10 genes in a loop (TERT,ELOVL3,FADS1,KRT79,ACO1,MGST1,PLAUR,CSF3R,MAPK10,DKK3)
### Task 3.1. Download trimmed GTEx gene expression data
https://raw.githubusercontent.com/barbarashih/202307_r_introduction/main/data/day3/gtex_gene_exprs.csv

### Task 3.2. Import and reorganise the gene expression data 
```r
  # Read in the file
  in_fp <- "data/day2/gtex_gene_exprs.csv"
  gene_exprs <- read.table(in_fp, header=TRUE, sep=',')
  # Rotate the table
  gene_exprs_num <- gene_exprs[,grep("GTEX", colnames(gene_exprs))]
  gene_exprs_num <- t(gene_exprs_num)
  gene_exprs_num <- as.data.frame(gene_exprs_num)
  colnames(gene_exprs_num) <- gene_exprs$Description
  gene_exprs_num$SAMPID <- rownames(gene_exprs_num)
  # R automatically replace - with . for column names
  # You need to turn it back so it can match the sample annotation
  # you can use gsub to replace characters, but 
  # this code here doesn't really work
  test <- gsub(".", "-", gene_exprs_num$SAMPID)
  head(test)
  # This is because . has special meaning for gsub and grep (it means any character), \\ is used to escape that special meaning
  gene_exprs_num$SAMPID <- gsub("\\.", "-", gene_exprs_num$SAMPID)
  head(gene_exprs_num)

  # Filter sample annotation 
  sample_annotation_filt <- sample_annotation[,c("SUBJID", "SAMPID", "AGE", "SEX", "SMTSD")]

  # Annotate the gene expression data
  gene_exprs_annotated <- merge(sample_annotation_filt, gene_exprs_num, by.x="SAMPID", by.y="SAMPID")

```
<details><summary>Challenge 3.1 </summary>
    
```r
# Try out different uses of . in gsub 
gsub("a.p", "12345", c("apple", "banana", "guava", "avocado"))
gsub("a..", "12345", c("apple", "banana", "guava", "avocado"))
gsub("a.", "12345", c("apple", "banana", "guava", "avocado"))
gsub("a.*a", "12345", c("apple", "banana", "guava", "avocado"))
```

</details>


### Task 3.3. Make a plotting function
```r
	myFun_plot <- function(in_df, x, y){
		plot_df <- data.frame(x = in_df[[x]], y = in_df[[y]])
		p1 <- ggplot(data=plot_df, aes(x=x, y=y)) + 
          geom_point() + 
          ggtitle(y) + 
          xlab(x) +
          ylab(y) +
          theme_bw()
        print(p1)
	}
```

### Task 3.4. Loop through the genes
```r
	# Organise your genes
	genes <- "TERT,ELOVL3,KRT79"
	genes <- strsplit(genes, ",")[[1]]
	# Loop through the genes
	for(current_gene in genes){
		myFun_plot(gene_exprs_annotated, x="AGE", y=current_gene)
	}
```

### Task 3.5. Nested loop through genes and annotations
```r
	# Loop through the genes
	for(current_gene in genes){
        for(current_annotation in c("AGE", "SMTSD", "SEX")){
		    myFun_plot(gene_exprs_annotated, x=current_annotation, y=current_gene)
        }
    }
    
```

## 4. If statements
### Task 4.1
```r
# Annotate age groups
unique(sample_annotation_filt$AGE)
ageGroup <- vector()
for(idx in 1:length(sample_annotation_filt$AGE )){
    current_age <- sample_annotation_filt$AGE[idx]
    if(current_age %in% c("60-69", "70-79")){
        ageGroup[idx] <- "Old"
    } else if(current_age %in% c("20-29", "30-39")){
        ageGroup[idx] <- "Young"
    } else {
        ageGroup[idx] <- "Middle"
    }
}
sample_annotation_filt$ageGroup <- ageGroup
head(sample_annotation_filt)
```
<details><summary>Challenge 4.1</summary>

```r
# Short way to do the same thing
# Check out what ifelse does
?ifelse
sample_annotation_filt$ageGroup <- "Middle"
sample_annotation_filt$ageGroup <- ifelse(sample_annotation_filt$AGE %in% c("60-69", "70-79"), "Old", sample_annotation_filt$ageGroup)
sample_annotation_filt$ageGroup <- ifelse(sample_annotation_filt$AGE %in% c("20-29", "30-39"), "Young", sample_annotation_filt$ageGroup)
```
</details>

<br> </br>



## 5. More on ggplot
### Task 5.1. Make a basic graph
```r
# Import library
library(ggplot2)

# Annotate the gene expression data with the sample_annotation with age group
plot_df <- merge(sample_annotation_filt, gene_exprs_num, by.x="SAMPID", by.y="SAMPID")

# plot
p <- ggplot(data = plot_df, aes(x=SMTSD, y=KRT79)) + geom_point()
p

```

### Task 5.2. Build on existing graph
```r
# Split by group
p + facet_wrap(~ageGroup)

# Change the theme
p + facet_wrap(~ageGroup) + theme_bw()
```

### Task 5.3. Group by categorical data
```r
# Spread the points so they don't overlap each other
p <- ggplot(data = plot_df, aes(x=SMTSD, y=KRT79)) + 
  geom_point(position = position_jitter(w = 0.2, h = 0)) + 
  theme_bw()
# Colour points by AGE
p <- ggplot(data = plot_df, aes(x=SMTSD, y=KRT79, colour=AGE)) + 
  geom_point(size = 5, position = position_jitter(w = 0.2, h = 0)) + 
  theme_bw()
p
# Change point shapes
p <- ggplot(data = plot_df, aes(x=SMTSD, y=KRT79, colour=AGE, shape=SEX)) + 
  geom_point(size = 5, position = position_jitter(w = 0.2, h = 0)) + 
  theme_bw()
p
head(plot_df)
# if TRUE, do something, otherwise something else
plot_df$SEX <- ifelse(plot_df$SEX == 1, "M", "F")
```

### Task 5.4. Dot shape, colour and fill
```r
# Change all shapes into a specified shpae
p <- ggplot(data = plot_df, aes(x=SMTSD, y=KRT79)) + 
  geom_point(size = 5, shape=21, position = position_jitter(w = 0.2, h = 0)) + 
  theme_bw() + scale_shape_manual(values=c(21, 25))
p

# Change shape by the values of a column
p <- ggplot(data = plot_df, aes(x=SMTSD, y=KRT79, shape=SEX)) + 
  geom_point(size = 5, position = position_jitter(w = 0.2, h = 0)) + 
  theme_bw() + scale_shape_manual(values=c(21, 25))
p

# Change colour by the values of a column
p <- ggplot(data = plot_df, aes(x=SMTSD, y=KRT79, colour=AGE, shape=SEX)) + 
  geom_point(size = 5, position = position_jitter(w = 0.2, h = 0)) + 
  theme_bw() + scale_shape_manual(values=c(21, 25))
p

# Change fill by the values of a column
p <- ggplot(data = plot_df, aes(x=SMTSD, y=KRT79, fill=AGE, shape=SEX)) + 
  geom_point(size = 5, position = position_jitter(w = 0.2, h = 0)) + 
  theme_bw() + scale_shape_manual(values=c(21, 25))
p

```

### Task 5.5. Factors
```r
# The plots are automatically sorted by alphabetical order 
# You can specify the order by using factors
plot_df$AGE_f <- factor(plot_df$AGE, levels = c("Young", "Middle", "Old"))

plot_df$AGE_f2 <- plot_df$AGE_f 
levels(plot_df$AGE_f2) <- c("Y", "M", "O")
head(plot_df[, c("AGE", "AGE_f", "AGE_f2")])

```


### Task 5.6. Change default colours
```r
# Change colour schemes
# Install library
install.packages("viridis")
library(viridis)

# Overwrite colour scheme
# https://sjmgarnier.github.io/viridis/reference/scale_viridis.html
p <- p + scale_color_viridis()
p <- p + scale_color_viridis(option="A")

# You can also specify 
```


