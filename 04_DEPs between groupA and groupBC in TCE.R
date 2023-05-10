## Input data
library(tidyverse)
library(qdapTools)
library(readxl)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)
library(survival)
library(survminer)
dir <- "/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-5/"
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Proteins20/")
tceall <- read.csv("TCE_proteins_expressed_all_samples_log2value.csv")
tce20 <- read.csv("TCE_proteins_expressed_at_least_20%_samples_log2value.csv")
sufall <- read.csv("Surface_proteins_expressed_all_samples_log2value.csv")
suf20 <- read.csv("Surface_proteins_expressed_at_least_20%_samples_log2value.csv")
biogroup <- data.frame(patients = c("IPAS7102", "IPAS0996", "IPAS0983", "IPAS0995",
                                    "IPAS0998", "IPAS0994", "IPAS0997", "IPAS7106",
                                    "IPAS0993", "IPAS0979",
                                    "IPAS7104", "IPAS7103", "IPAS0999", "IPAS7105",
                                    "IPAS7112", "IPAS7117", "IPAS7118", "IPAS7116",
                                    "IPAS0980", "IPAS7100", "IPAS0982",
                                    "IPAS7107", "IPAS0981", "IPAS7113", "IPAS7115",
                                    "IPAS7114"),
                       group = rep(c("tumorEnriched", "mixed", "immuneEnriched"),
                                   c(10, 11, 5)),
                       stringsAsFactors = FALSE)
exclude <- setdiff(biogroup$patients, colnames(suf20))  # "IPAS7100" "IPAS0979"
biogroup.suf <- biogroup[!(biogroup$patients %in% exclude), ]
dim(biogroup.suf)  # 24*2


#################### Figure-1 ########################

## Figure-1B

### ANOVA analysis
tce.aov <- as.data.frame(t(df2matrix(tceall)[, biogroup$patients]))
dif <- NULL
for (i in 1:ncol(tce.aov)) {
    res.aov <- aov(tce.aov[, i] ~ biogroup$group)
    pvalue <- summary(res.aov)[[1]][[5]][1]
    pairttest <- pairwise.t.test(tce.aov[, i], biogroup$group, p.adjust.method = "BH")
    adj.pvalue <- cbind(pairttest$p.value[1, 1],
                        pairttest$p.value[2, 1],
                        pairttest$p.value[2, 2])
    log2FC.tm <- mean(tce.aov[1:10, i]) - mean(tce.aov[11:21, i])
    log2FC.ti <- mean(tce.aov[1:10, i]) - mean(tce.aov[22:26, i])
    log2FC.mi <- mean(tce.aov[11:21, i]) - mean(tce.aov[22:26, i])
    dif <- rbind(dif, cbind(pvalue, adj.pvalue, log2FC.tm, log2FC.ti, log2FC.mi))
}
rownames(dif) <- colnames(tce.aov)
colnames(dif) <- c("aov_pvalue", "p_mixed vs p_immune",
                   "p_tumor vs p_immune", "p_tumor vs p_mixed",
                   "log2FC_tm", "log2FC_ti", "log2FC_mi")
dif <- as.data.frame(dif)
dif <- dif[order(dif$aov_pvalue), ]
# write.csv(dif, file.path(dir, "TCE_16449_proteins_ANOVA_result_(IPAS0982, IPAS7104).csv"))
res <- dif
sig.tce <- res[res[, 2] < 0.05 & res[, 3] < 0.05 & res[, 4] < 0.05 & abs(res[, 5]) > 1 &
                  abs(res[, 6]) > 1 & abs(res[, 7]) > 1, ]
dim(sig.tce)  # 510*7
# write.csv(sig.tce, file.path(dir, "TCE_335_SigProteins_ANOVA_result_(IPAS0982, IPAS7104).csv"))


### Plot heatmap
tce.sig <- tceall[tceall$X %in% rownames(sig.tce), ]
dim(tce.sig)  # 335*27
tce.sig2 <- as.data.frame(df2matrix(tce.sig))
tce.sig2 <- tce.sig2[, biogroup$patients]
rfdir <- "/Users/szhao5/Box Sync/UseFiles/Rfunction/"
source(file.path(rfdir, "ComplexHeatmap for TCE.R"))
pdf(file.path(dir, "Heatmap_TCE_335_sigProteins_in_3_groups_(IPAS0982, IPAS7104)-2.pdf"), width = 8, height = 10)
ComHeatmapTCE(ph = tce.sig2,
              cluster_columns = TRUE,
              cluster_rows = FALSE,
              show_column_dend = TRUE,
              show_row_dend = FALSE)
dev.off()
