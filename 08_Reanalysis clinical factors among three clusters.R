# Time: 21:24,12-11-2020
# Compare the tumor purity among three groups
# (for reviewer suggestions)
# (To save life, please complete it within 100 lines)

## Input library
library(tidyverse)
library(qdapTools)
library(readxl)
setwd("/Users/szhao5/Dropbox/MSlandscape/Submit")

## Tumor purity analysis
mat <- read_excel("Supplementary Table-1. Characteristics of 26 ascites samples with gastric adenocarcinoma.xlsx")
mat <- as.data.frame(mat)
mat[, 1:3]
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
matsel <- mat[, 1:3]
colnames(matsel) <- c("samples", "purity", "group")
matsel <- left_join(biogroup, matsel, by = c("patients" = "samples"))
colnames(matsel) <- c("patients", "cluster", "purity", "pathology")
head(matsel)

### violion plot among three clusters
library(ggpubr)
matsel$cluster <- factor(matsel$cluster,
                         levels = c("tumorEnriched", "mixed", "immuneEnriched"))
my.comparisons <- list(c("tumorEnriched", "mixed"),
                       c("mixed", "immuneEnriched"),
                       c("tumorEnriched", "immuneEnriched"))
p <- ggplot(matsel, aes(x = cluster, y = purity, fill = cluster)) +
  geom_violin(aes(color = cluster, trim = FALSE, scale = "width"
                  )) +
  scale_fill_manual(values = c("red", "steelblue", "orange")) +
  scale_color_manual(values = c("red", "steelblue", "orange")) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                             dodge.width = 0.75),
             col = "gray75",
             size = 0.8) +
  labs(title = "Plot of selected genes",
       x = "Groups", y = "log2 expression")
p
setwd("/Users/szhao5/Dropbox/MSlandscape/ReviewersR1")
## dev.print(pdf, "Tumor purity comparisons among three clusters.pdf")

### ANOVA analysis
res.aov <- aov(matsel$purity ~ matsel$cluster)
pvalue <- summary(res.aov)[[1]][[5]][1]
pairttest <- pairwise.t.test(matsel$purity, matsel$cluster, p.adjust.method = "BH")
adj.pvalue <- cbind(pairttest$p.value[1, 1],
                    pairttest$p.value[2, 1],
                    pairttest$p.value[2, 2])
log2FC.tm <- mean(matsel[matsel$cluster == "tumorEnriched", 3]) - mean(matsel[matsel$cluster == "tumorEnriched", 3])
log2FC.ti <- mean(matsel[matsel$cluster == "tumorEnriched", 3]) - mean(matsel[matsel$cluster == "immuneEnriched", 3])
log2FC.mi <- mean(matsel[matsel$cluster == "mixed", 3]) - mean(matsel[matsel$cluster == "immuneEnriched", 3])
dif <- cbind(pvalue, adj.pvalue, log2FC.tm, log2FC.ti, log2FC.mi)
colnames(dif) <- c("aov_pvalue", "p_tumor vs p_mixed", "p_tumor vs p_immune",
                   "p_mixed vs p_immune", "log2FC_tm", "log2FC_ti", "log2FC_mi")
dif <- as.data.frame(dif)
## write.csv(dif, "Tumor purity comparisons among three clusters.csv")

### Fisher exact test
temp <- data.frame(A = c(1, 0, 10),
                   B = c(1, 5, 4),
                   C = c(0, 4, 1),
                   stringsAsFactors = FALSE)
rownames(temp) <- c("Low", "Mid", "High")
fisher.test(temp)  # p = 0.004

### Heatmap for tumor purity
ph <- cbind(matsel[, c(1, 3)], matsel[, 3], matsel[, 3])
ph <- df2matrix(ph)
ph <- t(ph)
library(pheatmap)
ph.scale <- t(scale(t(ph)))
ph.scale[is.na(ph.scale)] <- 0
ph.scale.new <- ifelse(abs(ph.scale) > 1, sign(ph.scale)*1, ph.scale)
pheatmap(ph.scale.new, cluster_cols = FALSE)
## dev.print(pdf, "Tumor purity comparisons among three clusters_heatmap_scale.pdf")

pheatmap(ph, cluster_cols = FALSE)
## dev.print(pdf, "Tumor purity comparisons among three clusters_heatmap.pdf")








