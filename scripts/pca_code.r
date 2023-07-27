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
if (!require("plotly", quietly = TRUE)) { install.packages("plotly") } 
if (!require("viridis", quietly = TRUE)) { install.packages("viridis") } 
library(stringr)
library(ggcorrplot)
library(factoextra)
library(ggfortify)
library(umap)
library(plotly)
library(viridis)


#### Section 1: PCA
food_annotation <- read.csv("data/pca/pca_input_food_annotated.csv")				# This would be your sample annotation
pca_input <- read.csv("data/pca/pca_input_food_nutrient_matrix.csv", row.names=1)	# This would be your count matrix
pca_input <- data.matrix(pca_input)
pca_input_log <- log2(pca_input + min(pca_input[pca_input>0])/10)	# a small amount is added to all values in the data to remove zeros (zeros cannot be logged)

## PCA
res_pca <- prcomp(t(pca_input), scale=TRUE)		
fviz_eig(res_pca, addlabels = TRUE)


#### Section 2: reorganise PCA output for plotting
# reorganise to work out the contribution of each nutrient ( strength of contribution and the percentile-rank)
nutrient_contribution <- as.data.frame(res_pca$rotation)  		# with strength and direction. Number of PC equal to number of variables in the input
nutrient_contribution_abs <- abs(nutrient_contribution)			# use abs to make all values positive
nutrient_contribution_rank <- apply(nutrient_contribution_abs, 2, rank)													# use rank to rank the contribution, and dividide the rank by total number of rows (therefore it's a percentile, with the lower values indicating higher contribution)
nutrient_contribution_percentile <- as.data.frame(nutrient_contribution_rank/nrow(nutrient_contribution_rank))			# use rank to rank the contribution, and dividide the rank by total number of rows (therefore it's a percentile, with the lower values indicating higher contribution)
PC_names <- colnames(res_pca$rotation)
nutrient_names <- rownames(res_pca$rotation)
colnames(nutrient_contribution) <- PC_names
rownames(nutrient_contribution) <- nutrient_names
colnames(nutrient_contribution_percentile) <- PC_names
rownames(nutrient_contribution_percentile) <- nutrient_names

# The x-/y-coordinates to plot for each PC
pca_positions <- as.data.frame(res_pca$x)
pca_summary <- summary(res_pca)$importance[2,] * 100

# food relationships
pca_plot_df <- merge(food_annotation, as.data.frame(res_pca$x), by.x="wide_df_colnames_match", by.y=0)

# This block uses a slightly different method to do PCA (note that this library convert all the "PC" into "Dim."
# This is because it's compatible across multiple different PCA packages, and some packages uses "Dim.", dome "PC"
# Contributions for foods
#res_var <- get_pca_var(res_pca)
#res_var$coord          # Coordinates.  Loadings * the component standard deviations
#res_var$contrib        # % Contributions for each nutrient to the PCs (just the strength of contribution, no negative/positive direction). The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
#res_var$cos2           # Quality of representation. var.coord^2
# Contributions for nutrients
#res_ind <- get_pca_ind(res_pca)
#res_ind$coord          # Coordinates.  Loadings * the component standard deviations
#res_ind$contrib        # % Contributions for each sample to the PCs (just the strength of contribution, no negative/positive  direction). The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
#res_ind$cos2           # Quality of representation. var.coord^2



#### Section 2: Plotting every possible combination of dimentions and save the plots
dir.create("pca_output", showWarnings = FALSE)
for(PC_x in 1:ncol(coordinates)){
	for(PC_y in 1:ncol(coordinates)){
		# the if statement is to remove duplicate plots (i.e. it will just plot PC3 vs PC4, and not PC4 vs PC3)
		if(PC_x < PC_y){
			plotting_dimentions <- c(PC_x, PC_y)
			## Find top contributing nutriens for PC_x
			top_contrib_nutrients_x <- data.frame(contribution = nutrient_contribution[,PC_x], 
							nutrient_contribution_percentile = nutrient_contribution_percentile[,PC_x], 
							nutrient=row.names(nutrient_contribution))
			# order by contribution percentile
			top_contrib_nutrients_x <- top_contrib_nutrients_x[order(top_contrib_nutrients_x$nutrient_contribution_percentile),]
			print(paste0("Top contributing nutrients for ", PC_x, ":"))
			print(head(top_contrib_nutrients_x))

			## Find top contributing nutriens for PC_
			top_contrib_nutrients_y <- data.frame(contribution = nutrient_contribution[,PC_y], 
							nutrient_contribution_percentile = nutrient_contribution_percentile[,PC_y], 
							nutrient=row.names(nutrient_contribution))
			# order by contribution percentile
			top_contrib_nutrients_y <- top_contrib_nutrients_y[order(top_contrib_nutrients_y$nutrient_contribution_percentile),]
			print(paste0("Top contributing nutrients for ", PC_y, ":"))
			print(head(top_contrib_nutrients_y))

			## plot PCA
			current_pca_plot_df <- pca_plot_df[,c(paste0("PC", PC_x), paste0("PC", PC_y))]
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




#### Section 3: custom function that plots a PCA with gene expression level
## Function that takes in 
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
	pca_summary <- summary(res_pca)$importance[2,] * 100
	## plot PCA
	# organise a plot dataframe
	plot_df <- merge( in_annotation, as.data.frame(res_pca$x), by.x="wide_df_colnames_match", by.y=0)
	# Plotting using variable names isn't very straightforward in ggplot, so I have taken the columns needed and renamed them
	plot_df <- data.frame(food_id = plot_df[,1],
							PC_x = plot_df[,PC_x], 
							PC_y = plot_df[,PC_y], 
							food_name = plot_df[,"food_name"], 
							query_nutrient = plot_df[,query_nutrient])
	p <- ggplot(data=plot_df, aes(x=PC_x, y=PC_y)) + 
		geom_point(aes(colour=query_nutrient), size=3) +
		geom_text(aes(label=food_name), nudge_y = (max(plot_df$PC_y)-min(plot_df$PC_y))/75 ) +
		xlab(paste0(PC_x, " (", format(pca_summary[PC_x], digit=1), "%)")) +
		ylab(paste0(PC_y, " (", format(pca_summary[PC_y], digit=1), "%)")) +
		scale_colour_viridis(name = query_nutrient) + theme_bw()
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

## Use the function to plot the two PCs where a given nutrient has the strongest contribution to
pca_plot_function(pca_input, food_annotation, "Energy")
pca_plot_function(pca_input, food_annotation, "Protein")
pca_plot_function(pca_input, food_annotation, "Glucose")
pca_plot_function(pca_input, food_annotation, "Calcium..Ca")
pca_plot_function(pca_input, food_annotation, "Lactose")
pca_plot_function(pca_input, food_annotation, "Total.lipid..fat.")

pca_plot_function(pca_input, food_annotation, "Fiber..total.dietary")
pca_plot_function(pca_input, food_annotation, "Sugars..Total")
pca_plot_function(pca_input, food_annotation, "Iron..Fe")
pca_plot_function(pca_input, food_annotation, "Starch")
