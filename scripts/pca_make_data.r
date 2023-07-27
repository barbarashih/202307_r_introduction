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


#### Section 7: PCA
food_annotation <- read.csv("data/pca/pca_input_food_annotated.csv")				# This would be your sample annotation
pca_input <- read.csv("data/pca/pca_input_food_nutrient_matrix.csv", row.names=1)	# This would be your count matrix
pca_input <- data.matrix(pca_input)
pca_input_log <- log2(pca_input + min(pca_input[pca_input>0])/10)	# a small amount is added to all values in the data to remove zeros (zeros cannot be logged)

res_pca <- prcomp(t(pca_input), scale=TRUE)		
fviz_eig(res_pca, addlabels = TRUE)
nutrient_contribution <- res_pca$rotation  		# with strength and direction. Number of PC equal to number of variables in the input
coordinates <- res_pca$x
pca_summary <- summary(res_pca)$importance[2,] * 100


# this block uses a slightly different method to do PCA
#factoextra_res <- PCA(pca_input_log, scale=TRUE)	# !!! NOTE !!! scale normalises the data - you might want to set to false if your data is already normalised
# % variation explained per PC 
#loadings<-sweep(factoextra_res$var$coord,2,sqrt(factoextra_res$eig[1:5,1]),FUN="/")			# Contributions with positive/negative directions
#eig.val <- get_eigenvalue(factoextra_res)
# Contributions for foods
#res_var <- get_pca_var(factoextra_res)
#res_var$coord          # Coordinates.  loadings * the component standard deviations
#res_var$contrib        # % Contributions for each nutrient to the PCs (just the strength of contribution, no negative/positive direction). The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
#res_var$cos2           # Quality of representation. var.coord^2
# Contributions for nutrients
#res_ind <- get_pca_ind(factoextra_res)
#res_ind$coord          # Coordinates.  loadings * the component standard deviations
#res_ind$contrib        # % Contributions for each sample to the PCs (just the strength of contribution, no negative/positive  direction). The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
#res_ind$cos2           # Quality of representation. var.coord^2




# food relationships
pca_plot_df <- merge(food_annotation, res.ind$coord, by.x="wide_df_colnames_match", by.y=0)

#### Plotting every possible combination of dimentions and save the plots
dir.create("pca_output", showWarnings = FALSE)
for(PC_x in 1:ncol(coordinates)){
	for(PC_y in 1:ncol(coordinates)){
		# the if statement is to remove duplicate plots (i.e. it will just plot PC3 vs PC4, and not PC4 vs PC3)
		if(PC_x < PC_y){
			plotting_dimentions <- c(PC_x, PC_y)
			# Find top contributing nutriends
			top_contrib_nutrients <- data.frame(contribution = as.data.frame(nutrient_contribution$contrib)[[PC_x]], 
							nutrient=row.names(res_var$contrib))
			print(paste0("Top contributing nutrients for ", PC_x, ":"))
			tail(top_contrib_nutrients$nutrient[order(top_contrib_nutrients$contribution) ])


			current_pca_plot_df <- pca_plot_df[,c(paste0("Dim", PC_x), paste0("Dim", PC_x))]
			colnames(current_pca_plot_df)[1:2] <- c("PC_x", "PC_y")
			p <- ggplot(data=pca_plot_df, aes(x=PC_x, y=PC_y)) + 
				geom_point(aes(colour=Lactose)) +
				theme_bw() + geom_text(aes(label=food_name)) +
				xlab(paste0(PC_x, " (", pca_summary[PC_x], ")")) + ylab(paste0(PC_y, " (", pca_summary[PC_y], ")")) + 
				scale_colour_viridis()
			# ggplotly(p) # remove the comment if you want a interactive plot
			ggsave(paste0("pca_output/PC", PC_x, "_PC", PC_y, ".png"))
		}
	}
}




#### Function that takes in 
# - food_nutrient.matrix
# - food_annotation
# - nutrient of interest
## outputs
# plot PCA plot using the 2 PC this nutrient has highest contributing in
pca_plot_function <- function(in_mx, in_annotation, query_nutrient){
	# log transform the data
	in_mx <- data.matrix(in_mx)
	in_mxlog <- log2(in_mx + min(in_mx[in_mx>0])/10)
	res_pca <- prcomp(t(in_mxlog), rank = 20, scale=TRUE)	# !!! NOTE !!! scale normalises the data - you might want to set to false if your data is already normalised
	nutrient_contribution <- as.data.frame(res_pca$rotation )					# use absolute to make all values positive

	# Identify the top ranking PCs
	#res_nutrient <- get_pca_var(res_pca) # contribution for gene each PC
	# convert res_nutrient into rank
	nutrient_contribution_rank <- apply(abs(nutrient_contribution), 2, FUN=rank)
	nutrient_contribution_rank <- as.data.frame(nutrient_contribution_rank)
	nutrient_contribution_rank_query <- t(nutrient_contribution_rank)[, query_nutrient]
	nutrient_contribution_rank_query <- nutrient_contribution_rank_query[order(nutrient_contribution_rank_query)]
	PC_x <- names(nutrient_contribution_rank_query)[1]			# PC dispalyed in the x axis. This is the PC where the query shows the strongest contribution
	PC_y <- names(nutrient_contribution_rank_query)[2]			# PC dispalyed in the y axis. This is the PC where the query shows the second strongest contribution
	# Find the % contribution for PC_x and PC_y
	pca_summary <- summary(res_pca)$importance * 100
	# plot PCA
	plot_df <- merge( in_annotation, res_pca$x, by.x="wide_df_colnames_match", by.y=0)
	colnames(plot_df)[c(2,3)] <- c("PC_x", "PC_y")
	p <- ggplot(data=plot_df, aes(x=PC_x, y=PC_y)) + 
		geom_point(aes(colour=query_nutrient)), size=2)) +
		geom_text(aes(label=food_name)) +
		xlab(paste0(PC_x, " (", format(pca_summary[PC_x], digit=1), ")")) +
		xlab(paste0(PC_y, " (", format(pca_summary[PC_y], digit=1), ")")) +
		scale_colour_viridis()
	# Also, print a dataframe with other high contributing nutrient in these PCs
	other_nutrients_PC_x <- data.frame(nutrient.x=row.names(nutrient_contribution),
										contribution.x = as.data.frame(nutrient_contribution)[[PC_x]],
									contribution_rank.x = nutrient_contribution_rank[[PC_x]]) # percentile for the rank
	other_nutrients_PC_y <- data.frame(nutrient.y=row.names(nutrient_contribution),
									contribution.y = nutrient_contribution[[PC_y]], 
									contribution_rank.y = (nutrient_contribution_rank[[PC_y]])) # percentile for the rank
	other_nutrients_PC_x <- other_nutrients_PC_x[order(other_nutrients_PC_x$contribution_rank), ][1:5,]				
	other_nutrients_PC_y <- other_nutrients_PC_y[order(other_nutrients_PC_y$contribution_rank), ][1:5,]
	output_df <- cbind(other_nutrients_PC_x, other_nutrients_PC_y)
	print(output_df)

	return(p)
}

# Plot the two PCs where a given nutrient has the strongest contribution to
pca_plot_function(pca_input, food_annotation, "Energy")
pca_plot_function(pca_input, food_annotation, "Protein")
pca_plot_function(pca_input, food_annotation, "Glucose")
pca_plot_function(pca_input, food_annotation, "Calcium..Ca")


