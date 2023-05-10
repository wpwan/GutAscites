#Time: 22:31, Dec24-2022
#(What to do) Perform the proteomic analysis between NOS and SRC
#Complete it within 300 lines


# Load libraries and data
library(qdapTools)
library(tidyverse)
library(readxl)
setwd("/Users/robert/Dropbox/Project/Ascites_Revision/mybk/AscitesPatient/MS_data/Samples20%")
tce20 <- read.csv("TCE_proteins_expressed_at_least_20%_samples_log2value.csv")
tce20 <- as.data.frame(df2matrix(tce20))  # dim = 6,984*26
dim(tce20)
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
setwd("/Users/robert/Dropbox/Project/Ascites_Revision/Data&Results")
cld <- read_excel("Supplementary Table S1 for NOS and SRC analysis.xlsx", col_names = TRUE)
cld <- as.data.frame(cld)  # dim = 26*45
colnames(cld)[3] <- "Histology"
colnames(cld)[5] <- "Status"
colnames(cld)[6] <- "Ascites_collection_months"
colnames(cld)[7] <- "Pertoneal_metastasis_months"
biocld <- left_join(biogroup, cld, by = c("patients" = "Sample"))
head(biocld[, 1:8])


# Survival analysis between NOS and SRC among three clusters
library(survival)
library(survminer)
table(biocld$Histology)
biocld.his <- biocld[biocld$Histology != "N/A", ]
dim(biocld.his)  # dim = 24*46
biocld.his$Ascites_collection_months <- as.numeric(biocld.his$Ascites_collection_months)
biocld.his$Pertoneal_metastasis_months <- as.numeric(biocld.his$Pertoneal_metastasis_months)
biocld.his$Status <- ifelse(biocld.his$Status == "Dead", 1, 0)
survdiff(Surv(Pertoneal_metastasis_months, Status) ~ Histology, data = biocld.his)
sfit <- survfit(Surv(Pertoneal_metastasis_months, Status) ~ Histology, data = biocld.his)
mycol <- c("red", "black")
p <- ggsurvplot(sfit,
                data = biocld.his,
                size = 1,                 # change line size
                conf.int = FALSE,         # add confidence interval
                pval = TRUE,              # add p-value
                pval.method = TRUE,
                test.for.trend = FALSE,   # group > 2 to perform
                surv.median.line = "hv",
                palette = mycol,
                risk.table = TRUE,         # add risk table 
                risk.table.col = "strata", # risk table color by groups
                xlab = "Months",
                risk.table.height = 0.25,  # useful when multiple groups
                ggtheme = theme_classic())      # change ggplot2 theme
p
## dev.print(pdf, "Survival plot of peritoneal metastasis between NOS and SRC group in TCE 22 patients.pdf", useDingbats = FALSE)


survdiff(Surv(Ascites_collection_months, Status) ~ Histology, data = biocld.his)
sfit <- survfit(Surv(Ascites_collection_months, Status) ~ Histology, data = biocld.his)
mycol <- c("red", "black")
p <- ggsurvplot(sfit,
                data = biocld.his,
                size = 1,                 # change line size
                conf.int = FALSE,         # add confidence interval
                pval = TRUE,              # add p-value
                pval.method = TRUE,
                test.for.trend = FALSE,   # group > 2 to perform
                surv.median.line = "hv",
                palette = mycol,
                risk.table = TRUE,         # add risk table 
                risk.table.col = "strata", # risk table color by groups
                xlab = "Months",
                risk.table.height = 0.25,  # useful when multiple groups
                ggtheme = theme_classic())      # change ggplot2 theme
p
## dev.print(pdf, "Survival plot of ascites collection between NOS and SRC group in TCE 24 patients.pdf", useDingbats = FALSE)


# DEPs analysis between NOS and SRC group in 24 TCE patients
biocld.his.order <- biocld.his[order(biocld.his$Histology), ]
head(biocld.his.order[, 1:8])
tce20.sel <- tce20[, biocld.his.order$patients]
dim(tce20.sel)  # dim = 6984*24
data <- t(tce20.sel)
dif <- NULL
for (i in 1:ncol(data)) {
  ttest <- t.test(data[, i] ~ biocld.his.order$Histology)
  log2FC <- mean(data[14:24, i]) - mean(data[1:13, i])  # SRC - NOS
  dif <- rbind(dif, cbind(ttest$p.value, log2FC))
}
rownames(dif) <- colnames(data)
colnames(dif) <- c("pvalue", "log2FC")
dif <- as.data.frame(dif)
dif$adj.pval <- p.adjust(dif$pvalue, method = "BH")
dif <- dif[order(dif$adj.pval), ]
## write.csv(dif, "T.test result between SRC and NOS in 6984 proteins(0.2) in 24 samples from TCE.csv")


# Plot volcano plot
mydif <- read_excel("T.test result between SRC and NOS in 6984 proteins(0.2) in 24 samples from TCE_add keygene.xlsx")
vpd <- as.data.frame(mydif)
library(ggrepel)
diff <- ifelse(vpd$pvalue < 0.05 & vpd$log2FC > 1, "p < 0.05 & log2FC > 1",
        ifelse(vpd$pvalue < 0.05 & vpd$log2FC < -1, "p < 0.05 & log2FC < -1",
               "p >= 0.05 & abs(log2FC) <= 1"))
p <- ggplot(vpd, aes(log2FC, -log10(pvalue))) +
  geom_point(aes(col = diff), alpha = 0.8, size = 6) +
  scale_color_manual(values = c("steelblue", "red", "gray60"))
p + geom_text_repel(data = vpd[1:38,], aes(label = keygene), box.padding = 0.8) +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
  theme_light()
## dev.print(pdf, "Volcano plot between SRC and NOS group in 24 TCE patients.pdf", useDingbats = FALSE)
## dev.print(pdf, "Volcano plot between SRC and NOS group in 24 TCE patients2.pdf", useDingbats = FALSE)



















