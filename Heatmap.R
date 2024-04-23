#Install necessary packages

#install.packages("tidyverse") #package contains ggplot2 package
#install.packages("RcolorBrewer")
#install.packages("vegan")
#install.packages("devtools")

#Load necessary libraries

library(devtools)
library(tidyverse) 
library(RColorBrewer)
library(vegan)

#Read in data
Git_link <- "https://raw.githubusercontent.com/TobeEmma/Bio5202/main/Microbiome_data.csv?token=GHSAT0AAAAAACPATDITRECJFYA66RYWH4GMZQZVLJQ"
Bacteria <- read.csv(Git_link)

#Preview the first four rows and three columns
Bacteria[1:4, 1:3]

# Strip off the sample ids (index) and convert them to row names so that the data matrix contains only sequence count data.

row.names(Bacteria) <- Bacteria$index
Bacteria <- Bacteria[, -1]


#To make a heatmap with microbiome sequencing data, we ought to first transform the raw counts of reads to proportions within a sample

Microbial_Prop <- Bacteria/rowSums(Bacteria)

#Preview the first four rows and three columns
Microbial_Prop[1:4, 1:3]

#Now that the preliminary data transformation has been completed, I'll now begin to make the heatmap

#I'll use the color palette in the RColorBrewer packagehttp://127.0.0.1:44453/graphics/plot_zoom_png?width=1536&height=824

# colorRampPalette is in the RColorBrewer package.  This creates a colour palette that shades from light yellow to red in RGB space with 100 unique colours

scaleyellowmagenta <- colorRampPalette(c("lightyellow", "magenta"), space = "rgb")(100)

#Heatmap

heatmap(as.matrix(Microbial_Prop), Rowv = NA, Colv = NA, col = scaleyellowmagenta, margins = c(14,4))


# determine the maximum relative abundance for each column

max_rel_ab <- apply(Microbial_Prop, 2, max)
head(max_rel_ab)

# remove the genera/taxa with less than 1% as their maximum relative abundance

p1 <- names(which(max_rel_ab < 0.01))
Microbial_Prop_1 <- Microbial_Prop[, -which(names(Microbial_Prop) %in% p1)] #This reduces the microbial population to 8 genera suggesting that most of the taxa sampled occur at very low relative abundances.
#The heatmap with the new data looks like this, with a bit of extra fiddling to get all of the labels displayed.

# The margins command sets the width of the white space around the plot. The first element is the bottom margin and the second is the right margin
heatmap(as.matrix(Microbial_Prop_1), Rowv = NA, Colv = NA, col = scaleyellowmagenta, margins = c(12, 2))

#Adding a dendrogram for the samples using the vegan package as it has more options for distance metrics
#I'll do hierarchical clustering using the full dataset, but only display the more abundant taxa in the heatmap

# calculate the Bray-Curtis dissimilarity matrix on the full dataset:
Microbial_dist <- vegdist(Microbial_Prop, method = "bray")

# Do average linkage hierarchical row clustering. 
row_clus <- hclust(Microbial_dist, "aver")

# make the heatmap with Rowv = as.dendrogram(row.clus)
heatmap(as.matrix(Microbial_Prop_1), Rowv = as.dendrogram(row_clus), Colv = NA, col = scaleyellowmagenta, margins = c(14, 4))


#Adding a column dendrogram to cluster the genera that occur more often together
#I'll have to transpose to get the genera as rows

Microbial_dist_g <- vegdist(t(Microbial_Prop_1), method = "bray")

# Do average linkage hierarchical column clustering.
col_clus <- hclust(Microbial_dist_g, "aver")

#Make heatmap with Rowv - as.dendrogram(row.clus)
heatmap(as.matrix(Microbial_Prop_1), Rowv = as.dendrogram(row_clus), Colv = as.dendrogram(col_clus), col = scaleyellowmagenta, margins = c(14, 3))















