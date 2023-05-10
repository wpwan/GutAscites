# Time: 11-14-2019
# Sankey plot for 0.25RNA and 0.50Protein


library(tidyverse)
library(qdapTools)
library(readxl)
dir <- "/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-4/"

## Input data with common symbols expressed in 20% samples at least
setwd('/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Proteins20/')
tce20 <- read.csv('TCE_proteins_expressed_at_least_20%_samples_log2value.csv')
tce20 <- as.data.frame(df2matrix(tce20))
dim(tce20)  # 6984*26


## Input Gastric Cancer Ascites RNA data
setwd('/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/RNAsGC/')
rna <- read.csv("All RNAs in the 22 patients same with proteins data.csv")
rna <- as.data.frame(df2matrix(rna))
rna <- log2(rna + 1)
inter <- intersect(colnames(rna), colnames(tce20))
rna <- rna[, inter]
dim(rna)  # 22339*21
rna$rates <- rowSums(rna == 0)/21
rna20 <- rna[rna$rates < 0.8, ]
rna20 <- rna20[, !(colnames(rna20) %in% "rates")]
dim(rna20)  # 18904*21
# write.csv(rna20, file.path(dir, "18904 RNAs expressed at least 0.2 in 21 patients.csv"))


## Common symbols for RNA and protein
com <- intersect(rownames(rna20), rownames(tce20))
length(com)  # 6638

tce20.com <- tce20[rownames(tce20) %in% com, ]
rna20.com <- rna20[rownames(rna20) %in% com, ]
dim(tce20.com)  # 6638*26
dim(rna20.com)  # 6638*21
# write.csv(tce20.com, file.path(dir, "6638 proteins common symbols with RNAs expressed at least 0.2 in 26 patients.csv"))
# write.csv(rna20.com, file.path(dir, "6638 RNAs common symbols with proteins expressed at least 0.2 in 21 patients.csv"))


## Compute quantiles
rna20.com$rowMean <- rowMeans(rna20.com)
rna20.com <- rna20.com[order(rna20.com$rowMean), ]
q1 <- quantile(as.vector(rna20.com$rowMean), 0.25)
q2 <- quantile(as.vector(rna20.com$rowMean), 0.50)
q3 <- quantile(as.vector(rna20.com$rowMean), 0.75)
rna20.q1 <- rna20.com[rna20.com$rowMean <= q1, ]
dim(rna20.q1)  # 1660*22
rna20.q2 <- rna20.com[rna20.com$rowMean <= q2, ]
dim(rna20.q2)  # 3319*22
rna20.q3 <- rna20.com[rna20.com$rowMean <= q3, ]
dim(rna20.q3)  # 4978*22
rna20.com$group <- rep(c("1stQu", "2ndQu", "3rdQu", "4thQu"), c(1660, 1659, 1659, 1660))

tce20.com$rowMean <- rowMeans(tce20.com)
tce20.com <- tce20.com[order(tce20.com$rowMean), ]
q11 <- quantile(as.vector(tce20.com$rowMean), 0.25)
q21 <- quantile(as.vector(tce20.com$rowMean), 0.50)
q31 <- quantile(as.vector(tce20.com$rowMean), 0.75)
tce20.q11 <- tce20.com[tce20.com$rowMean <= q11, ]
dim(tce20.q11)  # 1660*27
tce20.q21 <- tce20.com[tce20.com$rowMean <= q21, ]
dim(tce20.q21)  # 3319*27
tce20.q31 <- tce20.com[tce20.com$rowMean <= q31, ]
dim(tce20.q31)  # 4978*27
tce20.com$group <- rep(c("1stQu", "2ndQu", "3rdQu", "4thQu"), c(1660, 1659, 1659, 1660))

inter1 <- intersect(rownames(rna20.com[rna20.com$group == "1stQu", ]),
                    rownames(tce20.com[tce20.com$group %in% c("3rdQu", "4thQu"), ]))  # 487
inter2 <- intersect(rownames(rna20.com[rna20.com$group %in% c("3rdQu", "4thQu"), ]),
                    rownames(tce20.com[tce20.com$group == "1stQu", ]))  # 442
dat1 <- data.frame(name = inter1, types = rep("HighProtein", 487), stringsAsFactors = FALSE)
dat2 <- data.frame(name = inter2, types = rep("HighRNA", 442), stringsAsFactors = FALSE)
dat <- rbind(dat1, dat2)
dim(dat)  # 929*2


tce20.com2 <- tibble::rownames_to_column(tce20.com)
tce20.com2 <- left_join(tce20.com2, dat, by = c("rowname" = "name"))
rna20.com2 <- tibble::rownames_to_column(rna20.com)
rna20.com2 <- left_join(rna20.com2, dat, by = c("rowname" = "name"))
# write.csv(tce20.com2, file.path(dir, "6638 proteins common symbols with RNAs expressed at least 0.2 in 26 patients_add quantiles.csv"))
# write.csv(rna20.com2, file.path(dir, "6638 RNAs common symbols with proteins expressed at least 0.2 in 21 patients_add quantiles.csv"))


## Sankey plot for 0.25RNA and 0.50Protein
setwd(dir)
df <- read_excel("Protein_RNA_Sanky_plot.xlsx")
df <- as.data.frame(df)
library(ggalluvial)
ggplot(df, aes(x = Survey, y = freq2, stratum = response,
               alluvium = subject, fill = response, label = response)) +
    scale_x_discrete(expand = c(0.1, 0.1)) +
    geom_flow() +
    geom_stratum(alpha = 0.5) +
    geom_text(stat = "stratum", size = 3) +
    theme(legend.position = "none") +
    theme_classic() +
    ggtitle("Proteins_RNAs correlation plot")
# dev.print(pdf, "Protein_RNA_Sanky_plot.pdf")
