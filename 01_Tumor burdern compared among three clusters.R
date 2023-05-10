### Tumor burdern analysis_gain
library(tidyverse)
library(qdapTools)
library(readxl)
library(reshape2)
setwd("/Users/robert/Dropbox/Revised_MS/ReviewersR2/ResultsR2/Figure")
tum <- read.table("patients_exomecn_tumorEnriched.seg.txt", header = TRUE)
mix <- read.table("patients_exomecn_mixed.txt", header = TRUE)
imm <- read.table("patients_exomecn_immuneEnriched.seg.txt", header = TRUE, sep = ",")
n1 <- length(levels(as.factor(tum$Patients))); n1
n2 <- length(levels(as.factor(mix$Patients))); n2
n3 <- length(levels(as.factor(imm$Patients))); n3


### Gain in three clusters_samples means
imm.gain <- imm[imm$seg.mean > 0, ]
summary(imm.gain$num.mark)
quantile(imm.gain$num.mark, c(0.01, 0.05, 0.85, 0.95, 0.99))
imm.gain$num.mark <- ifelse(imm.gain$num.mark > 400, 400, imm.gain$num.mark)
imm.mean <- NULL
for (i in levels(as.factor(imm$Patients))) {
    value <- mean(imm.gain[imm.gain == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    imm.mean <- rbind(imm.mean, temp)
}

mix.gain <- mix[mix$seg.mean > 0, ]
summary(mix.gain$num.mark)
quantile(mix.gain$num.mark, c(0.01, 0.05, 0.85, 0.95, 0.99))
mix.gain$num.mark <- ifelse(mix.gain$num.mark > 2000, 2000,
                     ifelse(mix.gain$num.mark < 3, 3, mix.gain$num.mark))
mix.mean <- NULL
for (i in levels(as.factor(mix$Patients))) {
    value <- mean(mix.gain[mix.gain == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    mix.mean <- rbind(mix.mean, temp)
}

tum.gain <- tum[tum$seg.mean > 0, ]
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
dev.print(pdf, "Violin plot of gain in three response groups.pdf", useDingbats = FALSE)


### Loss in three clusters_samples means
imm.loss <- imm[imm$seg.mean < 0, ]
summary(imm.loss$num.mark)
quantile(imm.loss$num.mark, c(0.01, 0.05, 0.85, 0.95, 0.99))
imm.loss$num.mark <- ifelse(imm.loss$num.mark > 1000, 1000, imm.loss$num.mark)
imm.mean <- NULL
for (i in levels(as.factor(imm$Patients))) {
    value <- mean(imm.loss[imm.loss == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    imm.mean <- rbind(imm.mean, temp)
}

mix.loss <- mix[mix$seg.mean < 0, ]
summary(mix.loss$num.mark)
quantile(mix.loss$num.mark, c(0.01, 0.05, 0.85, 0.95, 0.99))
mix.loss$num.mark <- ifelse(mix.loss$num.mark > 2500, 2500,
                     ifelse(mix.loss$num.mark < 3, 3, mix.loss$num.mark))
mix.mean <- NULL
for (i in levels(as.factor(mix$Patients))) {
    value <- mean(mix.loss[mix.loss == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    mix.mean <- rbind(mix.mean, temp)
}

tum.loss <- tum[tum$seg.mean < 0, ]
tum.mean <- NULL
for (i in levels(as.factor(tum$Patients))) {
    value <- mean(tum.loss[tum.loss == i, "num.mark"])
    temp <- data.frame(Patients = i, value = value)
    tum.mean <- rbind(tum.mean, temp)
}
dim(tum.mean)
head(tum.loss)
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
dev.print(pdf, "Violin plot of loss in three response groups.pdf", useDingbats = FALSE)




















