# Time: 01-09-2019
# Correlation analysis between protein and RNA
# High RNA vs High protein
# High RNA vs Low protein
# Low RNA vs Low protein
# Low RNA vs High protein


## Input library
library(tidyverse)
library(readxl)
library(qdapTools)
library(ggpubr)


################################## TCE platform #######################################

## High RNA vs High protein
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/Correlations_RNAPro/TCERNA/")
mat <- read_excel("860_HighProteins_HighRNAs_MeanValue.xlsx")
mat <- as.data.frame(mat)
head(mat)
ggscatter(mat, x = "Mean_HighRNAs", y = "Mean_HighProteins",
          palette = "jco", add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          rug = FALSE) +
    stat_cor(method = "spearman", label.x = 8, label.y = 12)
# dev.print(pdf, "TCE_Correlation_plot_860_HighProteins_HighRNAs_MeanValue.pdf", useDingbats = FALSE)


## Low RNA vs Low Protein
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/Correlations_RNAPro/TCERNA/")
mat <- read_excel("744_LowProteins_LowRNAs_MeanValue.xlsx")
mat <- as.data.frame(mat)
head(mat)
ggscatter(mat, x = "LowRNA_mean", y = "LowProtein_mean",
          palette = "jco", add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          rug = FALSE) +
    stat_cor(method = "spearman", label.x = 1, label.y = 1)
# dev.print(pdf, "TCE_Correlation_plot_744_LowProteins_LowRNAs_MeanValue.pdf", useDingbats = FALSE)


## High RNA vs Low protein
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/Correlations_RNAPro/TCERNA/")
mat <- read_excel("442_LowProteins_HighRNAs_MeanValue.xlsx")
mat <- as.data.frame(mat)
head(mat)
ggscatter(mat, x = "HighRNA_mean", y = "LowProtein_mean",
          palette = "jco", add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          rug = FALSE) +
    stat_cor(method = "spearman", label.x = 5, label.y = 1)
# dev.print(pdf, "TCE_Correlation_plot_442_LowProteins_HighRNAs_MeanValue.pdf", useDingbats = FALSE)


## Low RNA vs High protein
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/Correlations_RNAPro/TCERNA/")
mat <- read_excel("487_HighProteins_LowRNAs_MeanValue.xlsx")
mat <- as.data.frame(mat)
head(mat)
ggscatter(mat, x = "LowRNA_Mean", y = "HighProtein_mean",
          palette = "jco", add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          rug = FALSE) +
    stat_cor(method = "spearman", label.x = 1, label.y = 12)
# dev.print(pdf, "TCE_Correlation_plot_487_HighProteins_LowRNAs_MeanValue.pdf", useDingbats = FALSE)



################################# Surface platform #########################################

## High RNA vs High protein
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/Correlations_RNAPro/SurfaceRNA/")
mat <- read_excel("Surface_491_HighRNAs_HighProteins_RowMean.xlsx")
mat <- as.data.frame(mat)
head(mat)
ggscatter(mat, x = "HighRNA_mean", y = "HighProtein_mean",
          palette = "jco", add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          rug = FALSE) +
    stat_cor(method = "spearman", label.x = 8, label.y = 12)
# dev.print(pdf, "Surface_Correlation_plot_491_HighRNAs_HighProteins_RowMean.pdf", useDingbats = FALSE)


## Low RNA vs Low Protein
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/Correlations_RNAPro/SurfaceRNA/")
mat <- read_excel("Surface_368_LowRNAs_LowProteins_RowMean.xlsx")
mat <- as.data.frame(mat)
head(mat)
ggscatter(mat, x = "LowRNA_mean", y = "LowProtein_mean",
          palette = "jco", add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          rug = FALSE) +
    stat_cor(method = "spearman", label.x = 1, label.y = 1)
# dev.print(pdf, "Surface_Correlation_plot_368_LowRNAs_LowProteins_RowMean.pdf", useDingbats = FALSE)


## High RNA vs Low protein
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/Correlations_RNAPro/SurfaceRNA/")
mat <- read_excel("Surface_329_HighRNAs_LowProteins_RowMean.xlsx")
mat <- as.data.frame(mat)
head(mat)
ggscatter(mat, x = "HighRNA_mean", y = "LowProtein_mean",
          palette = "jco", add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          rug = FALSE) +
    stat_cor(method = "spearman", label.x = 6, label.y = 1)
# dev.print(pdf, "Surface_Correlation_plot_329_HighRNAs_LowProteins_RowMean.pdf", useDingbats = FALSE)


## Low RNA vs High protein
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/Correlations_RNAPro/SurfaceRNA/")
mat <- read_excel("Surface_345_LowRNAs_HighProteins_RowMean.xlsx")
mat <- as.data.frame(mat)
head(mat)
ggscatter(mat, x = "LowRNA_mean", y = "HighProtein_mean",
          palette = "jco", add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "blue", fill = "lightgray"),
          rug = FALSE) +
    stat_cor(method = "spearman", label.x = 1, label.y = 10)
# dev.print(pdf, "Surface_Correlation_plot_345_LowRNAs_HighProteins_RowMean.pdf", useDingbats = FALSE)










