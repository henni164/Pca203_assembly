#This script is used to make heatmaps of Pca scoring data for oat differential lines.
#The data should be formatted with the differential line names as the first column, and the isolate names as the first row.
#Make sure that there are no typos, extra characters, or spaces in the differential line or isolate names.
#Also, double check that all scores are in the key, with no ambiguous entries or characters not in the scale (like ?, 1/2?).
#If data is missing, just use NA or leave blank, you may need to modify read.csv line depending on how missing data is handled in your dataset.
#If the data is in excel, save as csv for import to R.
#If you use a different scoring key than the one this script assumes, change the script for your scoring system.
#The ComplexHeatmap package is used to make the heatmap.

#set working directory to where files are
setwd("/Users/evahenningsen/Documents/Research/203_assembly/Heatmap_figure")

#load required packages
#The circlize package contains the colorRamp2 function which I like to use to color the heatmap. Many other packages could be used instead.
#BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize) 

####SOUTH AFRICAN ISOLATES
#check.names is set to FALSE to prevent replace of spaces in column names with a dot (like Differential Line becoming Differential.Line)
#stringsAsFactors is set to FALSE to prevent the scoring data column being converted to a factor
scoring_data <- read.csv("/Users/evahenningsen/Documents/Research/203_assembly/Heatmap_figure/63_isolates_for_supplemental_figure.csv", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)

key <- setNames(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), c("0", "0;", ";", ";C", "1;", "1", "2", "3", "3+", "4"))
scoring_data_numeric <- scoring_data
#The unlist function is required if your dataset has more than one isolate column, but it's fine to use even with just one isolate
scoring_data_numeric[,2:ncol(scoring_data_numeric)] <- key[unlist(scoring_data_numeric[,2:ncol(scoring_data_numeric)])]

#make the heatmap data into a matrix
#leave out the first column since it doesn't belong in the heatmap itself
scoring_matrix <- as.matrix(scoring_data[,c(2:ncol(scoring_data)), drop = FALSE])

#save differential line name, and then add back row names
row_names <- scoring_data$"Differential"
rownames(scoring_matrix) = row_names

#Now, make heatmap

#Different clustering methods with hclust
#The default method for clustering_distance_rows is "euclidean" and the default method for clustering_method_rows is "complete".
scoring_plot <- Heatmap(scoring_matrix,
                        row_names_side = "left",
                        column_names_side = "top",
                        cluster_rows = FALSE,
                        cluster_columns = TRUE, #turn off if you don't want clustering
                        row_names_gp = gpar(fontsize = 10),
                        column_names_gp = gpar(cex= 0.9, fontsize = 10),
                        width = unit(205, "mm"), #adjust the widthe depending on how many isolates you have
                        heatmap_legend_param = list(title = NULL, at = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), 
                                                    labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")), 
                        col = colorRamp2(c(9, 8, 7, 6, 5, 4, 3, 2, 1, 0), 
                                         heat.colors(10)))

# Save the plot to a file
pdf("203_heatmap_supplemental.pdf", width = 9.5, height = 6)
scoring_plot
dev.off()
