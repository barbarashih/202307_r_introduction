#### Principal Component Analysis walkthrough
# Using data from https://fdc.nal.usda.gov/download-datasets.html#bkmk-1 to explain PCA
# https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/Dimension%20Reduction.pdf
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
# https://blog.bioturing.com/2018/06/14/principal-component-analysis-explained-simply/
# varying details on PCA
# https://www.youtube.com/watch?v=HMOI_lkzW08		# short video
# https://www.youtube.com/watch?v=_UVHneBUBW0&t=0s  # a bit more into the theory
# https://www.youtube.com/watch?v=0Jp4gsfOLMs		# more on the R practical sidev
# https://setosa.io/ev/principal-component-analysis/	# text explaination with interactive animation

#### Section 0: Settings
food_filepath <- "data/pca/food.csv"
nutrient_filepath <- "data/pca/nutrient.csv"
food_nutrient_filepath <- "data/pca/food_nutrient.csv"
working_directory <- "C:/Users/shihb/OneDrive - Lancaster University/work/teaching/workshop/202307_r_introduction/"

setwd(working_directory) # set working directory

#### Load library
if (!require("tidyverse", quietly = TRUE)) { install.packages("tidyverse") } 	 # this library contains stringr
if (!require("ggcorrplot", quietly = TRUE)) { install.packages("ggcorrplot") } 
if (!require("factoextra", quietly = TRUE)) { install.packages("factoextra") } 
if (!require("ggfortify", quietly = TRUE)) { install.packages("ggfortify") } 
if (!require("umap", quietly = TRUE)) { install.packages("umap") } 
if (!require("plotly", quietly = TRUE)) { install.packages("umap") } 
library(stringr)
library(ggcorrplot)
library(factoextra)
library(ggfortify)
library(umap)
library(plotly)


#### Section 1: Import data
food <- read.csv(food_filepath)
nutrient <- read.csv(nutrient_filepath)
food_nutrient <- read.csv(food_nutrient_filepath)
# manually added a few foods by their id


#### Section 2: Select num_per_category number from each food_category_id (so there is less data to look at)
num_per_category <- 5
food$description_modified <- gsub("raw, ", "", food$description)
food$description_modified <- gsub("cooked, ", "", food$description_modified)
food$food_firstword <- sapply(strsplit(food$description_modified, " "), function(x2)gsub(",|\\.", "", tolower(x2[1])))  # make a column with the first word of description (turned lower case) (removing any full stop or commas)
food$food_firsttwoword <- sapply(strsplit(food$description_modified, " "), function(x2)gsub(",|\\.", "", tolower(paste(x2[1], x2[2]))))  # make a column with the first word of description (turned lower case) (removing any full stop or commas)
food <- food[!duplicated(food$food_firsttwoword), ]											# remove any records that are duplicated for the "food_firstword" column
food <- food[!is.na(food$food_firsttwoword), ]													# remove any records that have NA (not available) for food description
food <- food[!food$food_firsttwoword %in% c("liquid", "solid"), ]								# remove any records that have "liquid" or "solid" as the first word
food_list <- split(food, food$food_category_id) 											# split food dataframe by food_category_id
food_list <- lapply(food_list, function(x)if(nrow(x)>num_per_category){x[1:num_per_category,]}else{x[1:nrow(x),]})		# keep the first 5 enteries of each food_category_id 
food_reduced <- do.call(rbind, food_list)													# turn the list back to a dataframe



#### Section 3: Merge the reduced food dataframe with nutrient information and turn it from a long dataframe into a wide dataframe
# Merge food_nutrient information 
food_reduced <- merge(food_reduced, food_nutrient, by = "fdc_id")
# Annotate nutrient names
colnames(nutrient) <- gsub("^name$", "nutrient_name", colnames(nutrient))					# rename the column name "name" as "nutrient_name" so it's more understandable ^ indicate the start of the word, $ indicate the end of the word
food_reduced <- merge(food_reduced, nutrient, by.x="nutrient_id", by.y="id", all.x=TRUE)
# Turn the dataframe from a long format to a wide format 
# There will be warning messages because of missing values
food_reduced_wide <- reshape(data = food_reduced[,c("nutrient_name", "amount", "fdc_id")],
                    idvar= "nutrient_name",
                    v.names= c("amount"),
                    timevar= "fdc_id",
                    direction = "wide")
					
# The word "amount." in column names is automatically added in front of each of the column name. They mean the amount of the nutrient in each food id - replacing "amount" to avoid confusion)
colnames(food_reduced_wide) <- gsub("amount\\.", "food.", colnames(food_reduced_wide))



#### Section 4: Reorganise the wide dataframe and food annotation
# PCA analysis takes in a numeric data matrix, so setting the nutrients as row.names
# Remove all special characters and spaces in food_reduced_wide$nutrient so it can be used as row.names
food_reduced_wide$nutrient_name <- str_replace_all(food_reduced_wide$nutrient_name, "[[:punct:]]", " ")
food_reduced_wide$nutrient_name <- gsub(",", ".", food_reduced_wide$nutrient_name)
food_reduced_wide$nutrient_name <- gsub(" ", ".", food_reduced_wide$nutrient_name)
# Remove nutrient name that is NA
food_reduced_wide <- food_reduced_wide[!is.na(food_reduced_wide$nutrient_name),]
row.names(food_reduced_wide) <- food_reduced_wide$nutrient_name
food_reduced_wide$nutrient_name <- NULL

	
#### Section 5: Data pre-processing 
# Take the top 50 nutrients
food_reduced_wide_50nutrients <- food_reduced_wide[1:50,]
# Replace NA with zeros
food_reduced_wide_50nutrients[is.na(food_reduced_wide_50nutrients)] <- 0
# Remove nutrients and food that are zeros for all foods
food_reduced_wide_50nutrients <- food_reduced_wide_50nutrients[rowSums(food_reduced_wide_50nutrients) > 0,]
food_reduced_wide_50nutrients <- food_reduced_wide_50nutrients[, colSums(food_reduced_wide_50nutrients) > 0]
# Remove nutrients that have zero variance
# PCA input uses column as the values so we need to transpose the data
# Data normalisation; this ensure that a more abundent nutrient that have larger values wouldn't have more contribution to the PCA analysis
# Using log scale to attempt to minimise the strong influence of extreme values. Add a small amount (the smallest non-zero value divided by 10) before logging the data because you cannot log zero
food_reduced_wide_50nutrients <- data.matrix(food_reduced_wide_50nutrients)



#### Section 6: Create non-duplicated food annotation 
# Create a food annotation in the same order as the food columns in wide - this is so that we can colour and label easier later
food_reduced$wide_df_colnames_match <- paste0("food.", food_reduced$fdc_id)
food_reduced_annotation <- food_reduced[food_reduced$wide_df_colnames_match %in% colnames(food_reduced_wide_50nutrients),]		# keep only foods that are still in food_reduced_wide after the data processing
food_reduced_annotation <- food_reduced_annotation[, c("wide_df_colnames_match", "description", "food_category_id", "food_firsttwoword")]
food_reduced_annotation <- food_reduced_annotation[!duplicated(food_reduced_annotation$wide_df_colnames_match),]
food_reduced_annotation <- food_reduced_annotation[match(colnames(food_reduced_wide_50nutrients), food_reduced_annotation$wide_df_colnames_match),]	# match the value order between colnames(food_reduced_wide) and food_reduced$wide_df_colnames_match

# add nutrient information for each food to the dataframe
food_reduced_wide_t <- as.data.frame(t(food_reduced_wide_50nutrients_log2)) 				# transpose nutrient information
food_reduced_annotation$food_category_id <- factor(food_reduced_annotation$food_category_id) 	# food group is catogerial data rather than numeric, therefore changing this to factor
food_reduced_annotation_with_nutrients <- merge(food_reduced_annotation	, food_reduced_wide_t, by.x="wide_df_colnames_match", by.y=0)
food_reduced_annotation_with_nutrients$description_short <- sapply(strsplit(food_reduced_annotation_with_nutrients$description, " "), function(x){x_len <- ifelse(length(x)<3, length(x), 2); paste(x[1:x_len], collapse=" ")})

#### Output files that will be used as input for PCA
write.csv(food_reduced_annotation_with_nutrients, "data/pca/food_annotated.csv", row.names=FALSE)  # manually added a column called food_name and renamed the file as pca_input
write.csv(food_reduced_wide_50nutrients, "data/pca/pca_input_food_nutrient_matrix.csv")
