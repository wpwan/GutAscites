library(tidyverse)
library(qdapTools)
library(readxl)
library(survival)
library(survminer)
dir <- "/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Reanalysis-5/"
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
dir2 <- "/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Submit"
df <- read_excel(file.path(dir2, "Supplementary Table-1. Characteristics of 26 ascites samples with gastric adenocarcinoma.xlsx"))
surdf <- df[, c(1, 5, 7)]
colnames(surdf) <- c("patients", "status", "months")
surdf <- as.data.frame(surdf)
surdf <- left_join(biogroup, surdf, by = "patients")
surdf$group2 <- surdf$group
surdf$group2 <- ifelse(surdf$group2 == "immuneEnriched", "mixed", surdf$group2)
surdf$months <- as.numeric(surdf$months)
surdf$status <- ifelse(surdf$status == "Dead", 1, 0)
survdiff(Surv(months, status) ~ group2, data = surdf)
sfit <- survfit(Surv(months, status) ~ group2, data = surdf)
mycol <- c("steelblue", "red")
p <- ggsurvplot(sfit,
                data = surdf,
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
write.csv(surdf, "My survival data_plot version.csv")
dev.print(pdf, "TCE_Survival plot of peritoneal metastasis in tumorEnriched and mixed group_final version_09-21-2020.pdf", useDingbats = FALSE)

