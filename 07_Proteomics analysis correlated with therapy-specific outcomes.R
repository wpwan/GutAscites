#Time： 21：22, Nov25-2022
#（What to do） Perform the proteomic analysis correlated with therapy
# Complete it within 200 lines


# Load libraries and data
library(qdapTools)
library(tidyverse)
library(readxl)
setwd("/Volumes/Robert5T/01_Research/21_MDACC/1_Gastric_cancer/AscitesPatient/MS_data/Samples20%")
tce20 <- read.csv("TCE_proteins_expressed_at_least_20%_samples_log2value.csv")
tce20 <- as.data.frame(df2matrix(tce20))  # dim = 6,984*26
## Chemotherapy after peritoneal metastasis before ascites collection_groupY(Yes) and groupN(No)
chemogroup <- data.frame(patients = c("IPAS0993", "IPAS7116", "IPAS0981",
                                      "IPAS0998", "IPAS7104", "IPAS7118",
                                      "IPAS7115", "IPAS0983", "IPAS7117",
                                      "IPAS7105", "IPAS0996", "IPAS0995",
                                      "IPAS7112", "IPAS0999", "IPAS0982",
                                      "IPAS0979", # groupY
                                      "IPAS7107", "IPAS7106", "IPAS7114",
                                      "IPAS7113", "IPAS0994", "IPAS0997",
                                      "IPAS0980", "IPAS7103", "IPAS7100",
                                      "IPAS7102"), # groupN
                         group = rep(c("Chemotherapy", "Unchemotherapy"), c(16, 10)),
                         stringsAsFactors = FALSE)
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
comgroup <- left_join(chemogroup, biogroup, by = "patients")
head(comgroup)
## Statistics of difference in therapy-specific group
mat <- table(comgroup$group.x, comgroup$group.y)
fisher.test(mat)  # p-value = 0.5173
chisq.test(mat)  # p-value = 0.5213


# DEPs analysis between chemotherapy and unchemotherapy group
tce20 <- tce20[, chemogroup$patients]
dim(tce20)  # dim = 6,984*26
head(tce20[, 1:6])
data <- t(tce20)
dif <- NULL
for (i in 1:ncol(data)) {
  ttest <- t.test(data[, i] ~ chemogroup$group)
  log2FC <- mean(data[1:16, i]) - mean(data[17:26, i])  # Chemotherapy - Unchemotherapy
  dif <- rbind(dif, cbind(ttest$p.value, log2FC))
}
rownames(dif) <- colnames(data)
colnames(dif) <- c("pvalue", "log2FC")
dif <- as.data.frame(dif)
dif$adj.pval <- p.adjust(dif$pvalue, method = "BH")
dif <- dif[order(dif$adj.pval), ]
setwd("/Users/robert/Dropbox/Project/Ascites_Revision/Data&Results")
## write.csv(dif, "T.test result between chemotherapy and unchemotherapy in 6984 proteins(0.2) in 26samples from TCE.csv")


### Plot volcano plot
vpd <- read_excel("T.test result between chemotherapy and unchemotherapy in 6984 proteins(0.2) in 26samples from TCE_Sig237_Join.xlsx")
vpd <- as.data.frame(vpd)
library(ggrepel)
diff <- ifelse(vpd$pvalue < 0.05 & vpd$log2FC > 1, "p < 0.05 & log2FC > 1",
        ifelse(vpd$pvalue < 0.05 & vpd$log2FC < -1, "p < 0.05 & log2FC < -1",
               "p >= 0.05 & abs(log2FC) <= 1"))
p <- ggplot(vpd, aes(log2FC, -log10(pvalue))) +
  geom_point(aes(col = diff), alpha = 0.8, size = 6) +
  scale_color_manual(values = c("steelblue", "red", "gray60"))
p + geom_text_repel(data = vpd[1:24,], aes(label = KeyGene), box.padding = 0.8) +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
  theme_light()
## dev.print(pdf, "Volcano plot of DEPs between chemotherapy and unchemotherapy in 26samples in TCE.pdf", useDingbats = FALSE)


## Develop bar plot function
ggbarplot <- function(data) {
  data <- within(data, pathways <- factor(pathways,
                                          levels = data[order(data$qvalue), "pathways"]))
  p <- ggplot(data, aes(x = pathways, y = -log10(qvalue))) +
    geom_bar(stat = "identity", width = 0.6, col = "gray", fill = "gray")
  p + coord_flip() +
    theme_classic() +
    scale_x_discrete(limits = rev(levels(data$pathways)))  # Reverse order
}


## Bar plot for pathways enriched in chemogroup
setwd("/Users/robert/Dropbox/Project/Ascites_Revision/Data&Results")
df <- read_excel("Hallmark_KEGG_Upregulated in chemotherapy group_Pathways_Sel.xlsx")
df <- as.data.frame(df)
df$qvalue <- as.numeric(df$qvalue)
ggbarplot(df)
## dev.print(pdf, "Hallmark_KEGG_Upregulated in chemotherapy group_Pathways_Sel.pdf")


## Bar plot for pathways enriched in unchemotherapy group
df <- read_excel("Hallmark_KEGG_upregulated_in_unchemotherapy_group_Pathways_Sel.xlsx")
df <- as.data.frame(df)
df$qvalue <- as.numeric(df$qvalue)
ggbarplot(df)
## dev.print(pdf, "Hallmark_KEGG_upregulated_in_unchemotherapy_group_Pathways_Sel.pdf")





