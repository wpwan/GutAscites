# Time: 12-20-2019
# Network analysis for group ABC


## Input library and data
library(tidyverse)
library(qdapTools)
library(readxl)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(nycflights13)
library(igraph)
library(intergraph)
library(sna)
library(ggplot2)
library(ggnetwork)
library(plotly)

dir <- "/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/"
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/TCE_Network")
mat <- read_excel("ALL pathways genes joined together_labels.xlsx")
mat <- as.data.frame(mat)
net <- graph.data.frame(mat[, c(4, 1, 2, 3, 5, 6)], directed = TRUE)
V(net)$degree <- centralization.degree(net)$res
df_net <- ggnetwork(net, layout = "fruchtermanreingold", weights = "log2FC_ti", niter = 5000)
plot <- ggplot(df_net, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(color = factor(Number)), alpha = 0.25) +
    geom_nodes(aes(size = degree, color = "gray")) +
    ggtitle("TCE_significant_pathways_and_genes_plot.pdf") +
    geom_nodetext(aes(color = factor(Number), label = vertex.names),
                  fontface = "plain", size = 2) +
    theme_blank()
plot
dev.print(pdf, file.path(dir, "TCE_groupA_Hallmark_pathways_genes.pdf"), useDingbats = FALSE)


































