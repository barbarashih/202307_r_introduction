# Introduction to R (Day 2)
## 1. List
Lists are R objects that can contain elements of different types.
### Task 1.1
```r
# Create a vector with a mixture of characters and numbers
myMeals <- c(breakfast="toast", lunch="sandwich", dinner=1.2)
# The numbers are automatically converted into characters (note the quotation marks around 1.2)
myMeals

# List allows a mixutre of characters and numbers
myMeals <- list(breakfast="toast", lunch="sandwich", dinner=1.2)
myMeals

# You can have multiple vectors in a list
myMeals <- list(breakfast=c("toast", "egg", "coffee"), 
                lunch=c("sandwich", "crisps"), 
                dinner=c("pie", "apple"))

# You can have a mixture of vectors and dataframes in a list
myMeals <- list(breakfast=data.frame(food=c("toast", "egg", "coffee"), cost=c(1,3,3)), 
                lunch=c("sandwich", "crisps"), 
                dinner=c("pie", "apple"))

# Look at the structure of the list
str(myMeals)
```

### Task 1.2
```r
# You can use double square brackets or dollar sign to refer to an element within the list
# double square bracket is a more accurate way of referencing
# Dollar sign looks neater and is more commonly used, but is less accurate and occassionally doesn't work for some R packages
myMeals[["breakfast"]]
myMeals$breakfast
# You can also chain the referencing
# (i.e. refer to a dataframe element in the list, then a column in the dataframe)
myMeals$breakfast$food
myMeals[["breakfast"]][["food"]]
myMeals[["breakfast"]]$food


```

<details>
<summary>Challenge 1.1</summary>

### Challenge 1.1
```r
# The difference between double/single square brackets 
doubleBracket <- myMeals[["breakfast"]]
singleBracket <- myMeals["breakfast"]
class(doubleBracket)
class(singleBracket)
myMeals[["breakfast"]][["food"]]
myMeals[["breakfast"]]["food"]

```
Q. What is the difference between double and single square brackets?
<details><summary>Answer</summary>
Single square bracket is used to access a subset of the dataframe/list, and the data type remain the same. For example, if you use single square brackets to access a column in a dataframe, you would get a dataframe. If you use double square brackets, you would extract the element and get a vector.
</details>

```r
# The difference between $ and [[]] (double square brackets)
# Double square brackets require exact match, whereas $ sign allows partial matching
myMeals$b # $ sign called the element name most similar to "b" (in this case, breakfast)
myMeals$b$fo # $ sign called the element name most similar to "fo" (in this case, food)
myMeals[["b"]] # This would give an error because there isn't an element named "b"
```
Q. Make a vector object containing food items from breakfast and lunch for the list below.
```r
myMeals <- list(breakfast=data.frame(food=c("toast", "egg", "coffee"), cost=c(1,3,3)), 
                lunch=c("sandwich", "crisps"), 
                dinner=c("pie", "apple"))
```

<details>
<summary>Answer</summary>

```r
	# Take a look at the data structure
	myMeals
	breakfast_and_lunch <- c(myMeals[["breakfast"]]$food, myMeals[["lunch"]])
	breakfast_and_lunch
```
</details>
</details>

### Task 1.3. Split a dataframe into a list
```r
# Set working directory
working_dir <- "C:/Users/shihb/OneDrive - Lancaster University/work/teaching/workshop/202307_r_introduction"
setwd(working_dir)

# Read in file
sample_annotation <- read.delim("data/day2/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
# Keep rows where the SMTSD column that has "Skin - Not Sun Exposed (Suprapubic)" or "Skin - Sun Exposed (Lower leg)"
keep <- sample_annotation$SMTSD %in% c("Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)")
keep <- grep("Skin", sample_annotation$SMTSD)
sample_annotation_filt <- sample_annotation[ keep , ]


# Split the dataframe by the column SMTSD
sample_annotation_list <- split(sample_annotation_filt, sample_annotation_filt$SMTSD)
# Explore what the generated list looks like
class(sample_annotation_list)
head(sample_annotation_list)
names(sample_annotation_list)
```

### Task 1.4. Split a character or a character vector into a list
```r
# You can split characters 
a_sentence <- "A white cat is fighting a black cat."
strsplit(a_sentence, " ")
# You can do strsplit with words
strsplit(a_sentence, "cat")

# You can split a character vector
another_sentence <- "It is sunny outside."
multiple_sentences <- c(a_sentence, another_sentence)
multiple_sentences
sentence_strsplit <- strsplit(multiple_sentences, " ")
sentence_strsplit

```
Q. Why do you think strsplit returns a list instead of a vector or a dataframe?
<details>
 <summary>Answer</summary>

```r
	# Look at the struture of sentence_strsplit
	str(sentence_strsplit)
	
	# Note that the first and the second sentence have different number of spaces. If strsplit output returns a dataframe, the second row will need to be padded with empty values.

```
</details>


<br> </br>

## 2. Loops
### Task 2.1. Basic loop
```r
# Loops are very useful for repeating a series of processes 
# Go from 1 to 5
for(idx in 1:5){
	print(idx)
}
# You can refer to each element of a vector or list in loops
items <- c("banana", "orange", "apple")
for(current_item in items){
	print_text <- paste0("This is: ", current_item)
	print(print_text)
}

```

### Task 2.2. Loop through a vector/list through its index
```r
# You can use indexes in loops to refer to different parts of a vector or list
for(idx in 1:length(items)){
	current_item <- items[idx]
	print_text <- paste0("This is: ", current_item)
	print(print_text)
}

```

### Task 2.3
```r
# Loop through the list we generate earlier
# Print the number of rows for each element
for(idx in 1:length(sample_annotation_list)){
	current_df <- sample_annotation_list[[idx]]
	current_sampleType <- names(sample_annotation_list)[idx]
	num_row <- nrow(current_df)
	print_text <- paste0(current_sampleType, "has a total of ", num_row, "rows.")
	print(print_text)	
}

# Save a copy of each element
# Create an output folder
dir.create("output", showWarnings = FALSE)
for(idx in 1:length(sample_annotation_list)){
	current_df <- sample_annotation_list[[idx]]
	current_sampleType <- names(sample_annotation_list)[idx]
	num_row <- nrow(current_df)

	# save the data frame
	out_fp <- paste0("output/", current_sampleType, ".csv")
	write.csv(current_df, out_fp, row.names=FALSE)
}
```

Q. Can you find the output folder?
<details>
 <summary>Answer</summary>
 
 ```r
 # It would be in your working directory
 # If you're not sure where it is currently, you can find it by 
 getwd()
 ```
</details>


### Task 2.4. Nested loop

```r
# Loop 1 (go through each element in the list)
for(idx1 in 1:length(sample_annotation_list)){
	current_df <- sample_annotation_list[[idx1]]
	current_sampleType <- names(sample_annotation_list)[idx1]
	# Loop 2 (print the first 5 elements in the SAMPID column in the dataframe)
	for(idx2 in 1:5){
		current_sampleName <- current_df$SAMPID[idx2]
		print_text <- paste0(current_sampleType, "sample ", idx2, ": ", current_sampleName)
		print(print_text)
	}
}
```


<br> </br>

## 3. Functions
### Task 3.1
```r
# We have already used a few in-built functions
class(read.delim)
class(class)
```
### Task 3.2
```r
# We can make our own functions to save time 
myFun_colour <- function(favCol){
	out_line <- paste0("My favourite colour is ", favCol)
	return(out_line)
}
# Try out the function
myFun_colour("red")
myFun_colour("green")

# More than one input
myFun_colour <- function(favCol, leastFavCol){
	out_line <- paste0("My favourite colour is ", favCol, ", and my least favourite colour is ", leastFavCol)
	return(out_line)
}
myFun_colour("red", "black")
myFun_colour("red")
```
<details><summary>Challenge 3.1</summary>

### Challenge 3.1
Write a function that takes in degree Celsius and returns Fahrenheit.

</details>

### Task 3.3
```r
# Setting a default
myFun_colour <- function(favCol, leastFavCol = "black"){
	out_line <- paste0("My favourite colour is ", favCol, ", and my least favourite colour is ", leastFavCol)
	return(out_line)
}
myFun_colour("red")
myFun_colour("red", "green")
```

Q. Look at some of the in-built functions. Do they have default values?
```r
?paste
```
<details><summary>Answer</summary>
Yes. By default, paste uses space for its separator (sep = " "). 
</details>

### Task 3.4
```r
# If you don't specify the return value
# the last evaluated value in your function will be returned
myFun_colour <- function(favCol, leastFavCol = "black"){
	paste0("My favourite colour is ", favCol, ", and my least favourite colour is ", leastFavCol)
	123
}
myFun_colour("red")
myFun_colour("black")

```
<br> </br>

## 4. List + function. Apply, sapply, lapply
### Task 4.1 Download data from GTEx
Go to https://gtexportal.org/home/datasets. Under GTEx Analysis V8 > Annotations, download:
1) GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt 

2) GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt

Copy the files to your working directory.

### Task 4.2. Read in these files
```r
# set working directory
working_dir <- "C:/Users/shihb/OneDrive - Lancaster University/work/teaching/workshop/202307_r_introduction"
setwd(working_dir)

# Read in metadata
sample_annotation_fp <- "data/day2/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
subject_annotation_fp <- "data/day2/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
sample_annotation <- read.delim(sample_annotation_fp)
subject_annotation <- read.delim(subject_annotation_fp)

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

### Task 4.3
```r
# sapply (takes vector or lists as inputs and returns vectors)
combined <- sapply(sampleID_split, FUN=function(x)paste0(x[1], "-", x[2]))
combined <- sapply(strsplit(sample_annotation$SAMPID, "-"), FUN=function(x)paste0(x[1], "-", x[2]))

# lapply (takes vector or lists as inputs returns lists)
combined_lapply <- lapply(sampleID_split, FUN=function(x)paste0(x[1], "-", x[2]))

# apply (takes data frame or matrix as inputs and return vectors, lists or array)
head(sample_annotation[,1:5])
combined_apply <- apply(sample_annotation, MARGIN=1, FUN=function(x)paste0(strsplit(x[1], "-")[[1]][1:2], collapse="-"))

```

### Task 4.4. Merge sample_annotation with subject_annotation
```r
# Now you can add a SUBJID column to sample_annotation 
sample_annotation$SUBJID <- combined_apply
sample_annotation <- merge(sample_annotation, subject_annotation, by.x="SUBJID", by.y="SUBJID")

```


<br> </br>


## 5. Libraries
### Task 5.1. Install a library
```r
# Other people have made libraries of useful functions
# Installing libraries
install.packages("ggplot2")
install.packages("gridExtra")
```

### Task 5.2. Use ggplot2 library
```r
# Importing libraries
library(ggplot2)
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


## 6. Putting it all together
## Aim to do
1. Import the GTEx gene expression data
1. Annotate the gene expression data with sample and subject/gender annotation
1. Create a function that can be used to plot a gene
1. Plot 10 genes in a loop (TERT,ELOVL3,FADS1,KRT79,ACO1,MGST1,PLAUR,CSF3R,MAPK10,DKK3)
### Task 6.1. Download trimmed GTEx gene expression data

### Task 6.2. Import and reorganise the gene expression data 
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
  gene_exprs_num$SAMPID <- gsub("\\.", "-", gene_exprs_num$SAMPID)

  # Filter sample annotation 
  sample_annotation_filt <- sample_annotation[,c("SUBJID", "SAMPID", "AGE", "SEX", "SMTSD")]

  # Annotate the gene expression data
  gene_exprs_annotated <- merge(sample_annotation_filt, gene_exprs_num, by.x="SAMPID", by.y="SAMPID")

```

### Task 6.3. Make a plotting function
```r
	myFun_plot <- function(in_df, x1, x2, y){
		plot_df <- data.frame(x1 = in_df[[x1]], x2 = in_df[[x2]], y = in_df[[y]])
		p1 <- ggplot(data=plot_df, aes(x=x1, y=y)) + geom_point() + ggtitle(y) + theme_bw()
		p2 <- ggplot(data=plot_df, aes(x=x2, y=y)) + geom_point() + ggtitle(y) + theme_bw()
		g <- arrangeGrob(p1, p2, nrow=1)
		
	}
```

### Task 6.4. Loop through the genes
```r
	# Organise your genes
	genes <- "TERT,ELOVL3,FADS1,KRT79,ACO1,MGST1,PLAUR,CSF3R,MAPK10,DKK3"
	genes <- strsplit(genes, ",")[[1]]
	# Loop through the genes
	for(gene in genes){
		myFun_plot(gene_exprs_annotated, x1="AGE", x2="SMTSD", y=gene)
	}
```


	myFun_plot2 <- function(in_df, x1, x2, y){
		plot_df <- data.frame(x1 = in_df[[x1]], x2 = in_df[[x2]], x3 = in_df[[x3]], y = in_df[[y]])
		p1 <- ggplot(data=plot_df, aes(x=x1, y=y)) + geom_point() + ggtitle(y) + theme_bw()
		p2 <- ggplot(data=plot_df, aes(x=x2, y=y)) + geom_point() + ggtitle(y) + theme_bw()
		p3 <- ggplot(data=plot_df, aes(x=x2, y=y)) + geom_point() + ggtitle(y) + theme_bw()
		grid.arrange(p1, p2, nrow=1)
	}