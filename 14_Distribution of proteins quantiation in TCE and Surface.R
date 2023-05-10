# Time: 12-12-2019
# Distribution of proteins quantitation in TCE and Surface


## Input library and Surface
library(tidyverse)
library(qdapTools)
library(readxl)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
dir <- "/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-9/"
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Proteins20/")
tceall <- read.csv("TCE_proteins_expressed_all_samples_log2value.csv")
tce20 <- read.csv("TCE_proteins_expressed_at_least_20%_samples_log2value.csv")
sufall <- read.csv("Surface_proteins_expressed_all_samples_log2value.csv")
suf20 <- read.csv("Surface_proteins_expressed_at_least_20%_samples_log2value.csv")
biogroup <- data.frame(patients = c("IPAS0983", "IPAS0993", "IPAS7102", "IPAS0996",
                                    "IPAS0979", "IPAS7106", "IPAS0994", "IPAS0995",
                                    "IPAS0998", "IPAS0997", "IPAS7104",  # Tumor
                                    "IPAS7105", "IPAS0982", "IPAS7100", "IPAS0999",
                                    "IPAS7103", "IPAS0980", "IPAS7117", "IPAS7116",
                                    "IPAS7118", "IPAS7112",  # Mixed
                                    "IPAS7115", "IPAS7107", "IPAS7113", "IPAS7114",
                                    "IPAS0981"),  # Immune
                       group = rep(c("tumorEnriched", "mixed", "immuneEnriched"),
                                   c(11, 10, 5)),
                       stringsAsFactors = FALSE)
exclude <- setdiff(biogroup$patients, colnames(suf20))  # "IPAS7100" "IPAS0979"
biogroup.suf <- biogroup[!(biogroup$patients %in% exclude), ]
dim(biogroup.suf)  # 24*2



############################ TCE platform ###################################

## Plot scatter plot for all the proteins in TCE
df <- tceall %>% tibble::column_to_rownames("X")
df <- df[, biogroup$patients]
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:26]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:26])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1645, 9), 1644))
# df$quat <- with(df, cut(rowMed, breaks = quantile(rowMed, probs = seq(0, 1, by = 0.1), na.rm = TRUE), include.lowest = TRUE))
# write.csv(df, file.path(dir, "TCE_all_proteins_with_quantiles.csv"))

### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of TCE all proteins",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(26, 0),
                    breaks = c(26, 24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of TCE_all proteins.pdf"))


### Histogram of all TCE proteins in 26 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white",
                   binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of TCE all proteins",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of TCE_all proteins.pdf"))


## Plot scatter plot for all the proteins in TCE_tumorEnriched group
df <- tceall %>% tibble::column_to_rownames("X")
df <- df[, biogroup$patients]
df <- df[, 1:11]  # Tumor-enriched group (A)
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:11]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:11])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1645, 9), 1644))
# df$quat <- with(df, cut(rowMed, breaks = quantile(rowMed, probs = seq(0, 1, by = 0.1), na.rm = TRUE), include.lowest = TRUE))
# write.csv(df, file.path(dir, "TCE_all_proteins_with_quantiles in tumorEnriched group.csv"))

### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of TCE all proteins in tumorEnriched group",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(11, 0),
                    breaks = c(10, 8, 6, 4, 2, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of TCE_all proteins in tumorEnriched group.pdf"))

### Histogram of all TCE proteins in 11 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white", binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of TCE all proteins in tumorEnriched group",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of TCE_all proteins in tumorEnriched group.pdf"))


## Plot scatter plot for all the proteins in TCE_Mixed group
df <- tceall %>% tibble::column_to_rownames("X")
df <- df[, biogroup$patients]
df <- df[, 12:21]  # group (B)
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:10]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:10])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1645, 9), 1644))
# df$quat <- with(df, cut(rowMed, breaks = quantile(rowMed, probs = seq(0, 1, by = 0.1), na.rm = TRUE), include.lowest = TRUE))
# write.csv(df, file.path(dir, "TCE_all_proteins_with_quantiles in mixed group.csv"))

### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of TCE all proteins in mixed group",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(10, 0),
                    breaks = c(10, 8, 6, 4, 2, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of TCE_all proteins in mixed group.pdf"))

### Histogram of all TCE proteins in 10 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white", binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of TCE all proteins in mixed group",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of TCE_all proteins in mixed group.pdf"))


## Plot scatter plot for all the proteins in TCE_ImmuneEnriched group
df <- tceall %>% tibble::column_to_rownames("X")
df <- df[, biogroup$patients]
df <- df[, 22:26]  # group (C)
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:5]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:5])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1645, 9), 1644))
# df$quat <- with(df, cut(rowMed, breaks = quantile(rowMed, probs = seq(0, 1, by = 0.1), na.rm = TRUE), include.lowest = TRUE))
# write.csv(df, file.path(dir, "TCE_all_proteins_with_quantiles in immuneEnriched group.csv"))

### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of TCE all proteins in immuneEnriched group",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(5, 0),
                    breaks = c(5, 4, 3, 2, 1, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of TCE_all proteins in immuneEnriched group.pdf"))

### Histogram of all TCE proteins in 5 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white", binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of TCE all proteins in immuneEnriched group",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of TCE_all proteins in immuneEnriched group.pdf"))


## B+C group
## Plot scatter plot
df <- tceall %>% tibble::column_to_rownames("X")
df <- df[, biogroup$patients]
df <- df[, -c(1:11)]  # group B+C
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:15]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:15])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1645, 9), 1644))
# df$quat <- with(df, cut(rowMed, breaks = quantile(rowMed, probs = seq(0, 1, by = 0.1), na.rm = TRUE), include.lowest = TRUE))
# write.csv(df, file.path(dir, "TCE_all_proteins_with_quantiles in group BC.csv"))

### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of TCE all proteins in mixed + immunueEnriched group",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(15, 0),
                    breaks = c(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of TCE_all proteins in group BC.pdf"))

### Histogram of all TCE proteins in 15 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white", binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of TCE all proteins in group B+C",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of TCE_all proteins in group BC.pdf"))



############################ Surface platform ###################################

## Plot scatter plot for all the proteins in Surface
df <- sufall %>% tibble::column_to_rownames("X")
df <- df[, biogroup.suf$patients]
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:24]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:24])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1385, 9), 1381))
# write.csv(df, file.path(dir, "Surface_all_proteins_with_quantiles.csv"))


### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of Surface all proteins",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(24, 0),
                    breaks = c(24, 22, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of Surface_all proteins.pdf"))


### Histogram of all Surface proteins in 24 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white",
                   binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of Surface all proteins",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of Surface_all proteins.pdf"))


## Plot scatter plot for all the proteins in Surface in TumorEnriched group
df <- sufall %>% tibble::column_to_rownames("X")
df <- df[, biogroup.suf$patients]
df <- df[, 1:10]  # TumorEnriched group
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:10]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:10])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1385, 9), 1381))
# write.csv(df, file.path(dir, "Surface_all_proteins_with_quantiles in tumorEnriched group.csv"))


### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of Surface all proteins in tumorEnriched group",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(10, 0),
                    breaks = c(10, 8, 6, 4, 2, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of Surface_all proteins in tumorEnriched group.pdf"))


### Histogram of all Surface proteins in 10 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white",
                   binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of Surface all proteins in tumorEnriched group",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of Surface_all proteins in tumorEnriched group.pdf"))


## Plot scatter plot for all the proteins in Surface in mixed group
df <- sufall %>% tibble::column_to_rownames("X")
df <- df[, biogroup.suf$patients]
df <- df[, 11:19]  # Mixed group
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:9]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:9])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1385, 9), 1381))
# write.csv(df, file.path(dir, "Surface_all_proteins_with_quantiles in mixed group.csv"))


### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of Surface all proteins in mixed group",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(9, 0),
                    breaks = c(9, 8, 7, 6, 5, 4, 3, 2, 1, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of Surface_all proteins in mixed group.pdf"))


### Histogram of all Surface proteins in 9 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white",
                   binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of Surface all proteins in mixed group",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of Surface_all proteins in mixed group.pdf"))


## Plot scatter plot for all the proteins in Surface in immuneEnriched group
df <- sufall %>% tibble::column_to_rownames("X")
df <- df[, biogroup.suf$patients]
df <- df[, 20:24]  # ImmuneEnriched group
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:5]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:5])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1385, 9), 1381))
# write.csv(df, file.path(dir, "Surface_all_proteins_with_quantiles in immuneEnriched group.csv"))


### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of Surface all proteins in immuneEnriched group",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(5, 0),
                    breaks = c(5, 4, 3, 2, 1, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of Surface_all proteins in immuneEnriched group.pdf"))


### Histogram of all Surface proteins in 5 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white",
                   binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of Surface all proteins in immuneEnriched group",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of Surface_all proteins in immuneEnriched group.pdf"))


## Plot scatter plot for all the proteins in Surface in mixed + immuneEnriched group
df <- sufall %>% tibble::column_to_rownames("X")
df <- df[, biogroup.suf$patients]
df <- df[, 11:24]  # Mixed + ImmuneEnriched group
df$psum <- rowSums(df != 0)
df <- df[order(df$psum), ]  # Samples number with protein expression
library(matrixStats)
df$rowMed <- rowMedians(as.matrix(df[, 1:14]))  # Proteins medians
df$rowMean <- rowMeans(df[, 1:14])
df <- df[order(df$rowMean), ]
df$quat <- rep(c(1:9, 10), c(rep(1385, 9), 1381))
# write.csv(df, file.path(dir, "Surface_all_proteins_with_quantiles in mixed + immuneEnriched group.csv"))


### Make a plot
mylabel <- data.frame(name = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                               "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                               "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                               "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      Genes = c("CDH1", "TP53", "FAT1", "KMT2C", "ARID1A", "CREBBP",
                                "CTAGE1", "CTNNA2", "ERBB2", "EGFR", "YAP1",
                                "C10orf54", "ITGAM", "LTF", "LAIR1", "HAVCR2",
                                "FCGR3B", "CD48", "LILRA6", "SIGLEC9"),
                      stringsAsFactors = FALSE)
df <- tibble::rownames_to_column(df)
df <- left_join(df, mylabel, by = c("rowname" = "name"))
df <- df %>% tibble::column_to_rownames("rowname")
library(ggrepel)
df <- df[order(df$psum, decreasing = TRUE), ]
p <- ggplot(df, aes(x = psum, y = rowMean, col = factor(quat))) +
    geom_point(shape = 1, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("violetred4", "red", "orangered", "gold", "yellow",
                                 "snow2", "slategray1", "skyblue", "royalblue", "purple4")) +
    labs(title = "Plot of Surface all proteins in immuneEnriched group",
         x = "Number of observations",
         y = "Mean log2(value)") +
    scale_x_reverse(limits = c(14, 0),
                    breaks = c(14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)) +
    geom_text_repel(data = df, aes(label = Genes), col = "black", box.padding = 1) +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Scatter plot of Surface_all proteins in mixed + immuneEnriched group.pdf"))


### Histogram of all Surface proteins in 14 patients
p <- ggplot(df, aes(x = psum)) +
    geom_histogram(col = "black", fill = "white",
                   binwidth = 1, position = position_dodge()) +
    labs(title = "Plot of Surface all proteins in mixed + immuneEnriched group",
         x = "Number of observations",
         y = "Number of protein groups") +
    scale_x_reverse() +
    theme_classic()
p
# dev.print(pdf, file.path(dir, "Histogram of Surface_all proteins in mixed + immuneEnriched group.pdf"))









































