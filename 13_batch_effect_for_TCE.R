## Time: 08-01-2020
## Compute the tumor burdern
## (Supplementary Figure-1 and methods)
## (To save life, please complete it within 200 lines)


### Tumor burdern analysis_gain
library(tidyverse)
library(qdapTools)
library(readxl)
library(reshape2)
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/DNAsGC/Selected_Exomecn")
tum <- read.table("22_patients_exomecn_tumorEnriched.seg.txt", header = TRUE)
mix <- read.table("22_patients_exomecn_mixed.seg.txt", header = TRUE, sep = ",")
imm <- read.table("22_patients_exomecn_immuneEnriched.seg.txt", header = TRUE, sep = ",")
n1 <- length(levels(as.factor(tum$Patients))); n1  # n=10
n2 <- length(levels(as.factor(mix$Patients))); n2  # n=8
n3 <- length(levels(as.factor(imm$Patients))); n3  # n=4


### Gain in three clusters_samples means
imm.gain <- imm[imm$seg.mean > 0, ]  # 538
summary(imm.gain$num.mark)
quantile(imm.gain$num.mark, c(0.01, 0.05, 0.85, 0.95, 0.99))
imm.gain$num.mark <- ifelse(imm.gain$num.mark > 400, 400, imm.gain$num.mark)
imm.mean <- NULL
for (i in levels(as.factor(imm$Patients))) {
    value <- mean(imm.gain[imm.gain == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    imm.mean <- rbind(imm.mean, temp)
}

mix.gain <- mix[mix$seg.mean > 0, ]  # 1344
summary(mix.gain$num.mark)
quantile(mix.gain$num.mark, c(0.01, 0.05, 0.85, 0.95, 0.99))
mix.gain$num.mark <- ifelse(mix.gain$num.mark > 2500, 2500,
                     ifelse(mix.gain$num.mark < 3, 3, mix.gain$num.mark))
mix.mean <- NULL
for (i in levels(as.factor(mix$Patients))) {
    value <- mean(mix.gain[mix.gain == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    mix.mean <- rbind(mix.mean, temp)
}

tum.gain <- tum[tum$seg.mean > 0, ]  # 989
tum.mean <- NULL
for (i in levels(as.factor(tum$Patients))) {
    value <- mean(tum.gain[tum.gain == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    tum.mean <- rbind(tum.mean, temp)
}
dim(tum.mean)
head(tum.gain)
tm <- rbind(tum.mean, mix.mean)
tmi <- rbind(tm, imm.mean)
library(ggpubr)
tmi$group <- rep(c("tumor", "mixed", "immune"), c(n1, n2, n3))
my.comparisons <- list(c("tumor", "mixed"),
                       c("tumor",  "immune"),
                       c("mixed", "immune"))
tmi <- within(tmi, group <- factor(group, levels = c("tumor", "mixed", "immune")))
p <- ggplot(tmi, aes(x = group, y = value, fill = group)) +
        geom_violin(aes(color = group),
                trim = FALSE, scale = "width",
                position = position_dodge(0.9)) +
        geom_boxplot(# aes(color = "gray"),
            width = 0.1,
            fill = "white",
            position = position_dodge(0.9),
            outlier.shape = NA) +
        scale_fill_manual(values = c("red", "steelblue", "orange")) +
        scale_color_manual(values = c("red", "steelblue", "orange")) +
        geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                                   dodge.width = 0.75),
                   col = "gray75",
                   size = 0.8) +
        labs(title = "Plot of selected genes",
             x = "Groups", y = "log2 expression") +
        theme(plot.title = element_text(color = "black",
                                        size = 14,
                                        face = "bold.italic"),
              axis.title.x = element_text(color = "black",
                                          size = 10,
                                          face = "bold"),
              axis.title.y = element_text(color = "black",
                                          size = 10,
                                          face = "bold")) +
#       geom_text(x = 4, y = 14, label = "p < 0.001", color = "black") +
        ## facet_wrap(~ variable, scales = "free") + # fixed
        stat_compare_means(comparisons = my.comparisons,
                           ## method = "anova",
                           hide.ns = FALSE) +
        stat_compare_means()  # Add global value
    p + theme_classic()
dev.print(pdf, "Violin plot of gain in three response groups_22samples_means.pdf", useDingbats = FALSE)


### Loss in three clusters_samples means
imm.gain <- imm[imm$seg.mean < 0, ]  # 538
summary(imm.gain$num.mark)
quantile(imm.gain$num.mark, c(0.01, 0.05, 0.85, 0.95, 0.99))
imm.gain$num.mark <- ifelse(imm.gain$num.mark > 1000, 1000, imm.gain$num.mark)
imm.mean <- NULL
for (i in levels(as.factor(imm$Patients))) {
    value <- mean(imm.gain[imm.gain == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    imm.mean <- rbind(imm.mean, temp)
}

mix.gain <- mix[mix$seg.mean < 0, ]  # 1344
summary(mix.gain$num.mark)
quantile(mix.gain$num.mark, c(0.01, 0.05, 0.85, 0.95, 0.99))
mix.gain$num.mark <- ifelse(mix.gain$num.mark > 2500, 2500,
                     ifelse(mix.gain$num.mark < 3, 3, mix.gain$num.mark))
mix.mean <- NULL
for (i in levels(as.factor(mix$Patients))) {
    value <- mean(mix.gain[mix.gain == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    mix.mean <- rbind(mix.mean, temp)
}

tum.gain <- tum[tum$seg.mean < 0, ]  # 989
tum.mean <- NULL
for (i in levels(as.factor(tum$Patients))) {
    value <- mean(tum.gain[tum.gain == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    tum.mean <- rbind(tum.mean, temp)
}
dim(tum.mean)
head(tum.gain)
tm <- rbind(tum.mean, mix.mean)
tmi <- rbind(tm, imm.mean)
library(ggpubr)
tmi$group <- rep(c("tumor", "mixed", "immune"), c(n1, n2, n3))
my.comparisons <- list(c("tumor", "mixed"),
                       c("tumor",  "immune"),
                       c("mixed", "immune"))
tmi <- within(tmi, group <- factor(group, levels = c("tumor", "mixed", "immune")))
p <- ggplot(tmi, aes(x = group, y = value, fill = group)) +
        geom_violin(aes(color = group),
                trim = FALSE, scale = "width",
                position = position_dodge(0.9)) +
        geom_boxplot(# aes(color = "gray"),
            width = 0.1,
            fill = "white",
            position = position_dodge(0.9),
            outlier.shape = NA) +
        scale_fill_manual(values = c("red", "steelblue", "orange")) +
        scale_color_manual(values = c("red", "steelblue", "orange")) +
        geom_point(position = position_jitterdodge(jitter.width = 0.2,
                                                   dodge.width = 0.75),
                   col = "gray75",
                   size = 0.8) +
        labs(title = "Plot of selected genes",
             x = "Groups", y = "log2 expression") +
        theme(plot.title = element_text(color = "black",
                                        size = 14,
                                        face = "bold.italic"),
              axis.title.x = element_text(color = "black",
                                          size = 10,
                                          face = "bold"),
              axis.title.y = element_text(color = "black",
                                          size = 10,
                                          face = "bold")) +
#       geom_text(x = 4, y = 14, label = "p < 0.001", color = "black") +
        ## facet_wrap(~ variable, scales = "free") + # fixed
        stat_compare_means(comparisons = my.comparisons,
                           ## method = "anova",
                           hide.ns = FALSE) +
        stat_compare_means()  # Add global value
    p + theme_classic()
dev.print(pdf, "Violin plot of loss in three response groups_22samples_means.pdf", useDingbats = FALSE)










