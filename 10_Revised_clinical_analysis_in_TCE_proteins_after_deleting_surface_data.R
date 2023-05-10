()## Time: 08-19-2020 04:11
## Clinical analysis for completing data
## (For figure-1, 2, 3, 4 and methods)
## (To save life, please complete it as short as possible)


#########################
### Figure-1 data
#########################

library(tidyverse)
library(qdapTools)
library(readxl)
dir <- "/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Proteins_CombinedOC/clinical_analysis_after_deleting_surface"
setwd("/Users/szhao5/Box Sync/Wdir/1_GastricCancer/Ascites/Proteins20/")
tceall <- read.csv("TCE_proteins_expressed_all_samples_log2value.csv")
head(tceall[, 1:4])

tceall <- as.data.frame(df2matrix(tceall))
head(tceall[, 1:4])
setwd(dir)
write.csv(tceall, "TCE_proteins_expressed_all_samples_log2value_add_proteins_sum.csv")


### Intetinal and diffuse
intestinal <- read_excel("Supplementary Table-1. Characte_Intestinal.xlsx")
diffuse <- read_excel("Supplementary Table-1. Characte 2_Diffuse.xlsx")
intestinal.tce <- tceall[, intestinal$Sample]
diffuse.tce <- tceall[, diffuse$Sample]

intestinal.tce$sum <- rowSums(intestinal.tce != 0)
diffuse.tce$sum <- rowSums(diffuse.tce != 0)
write.csv(intestinal.tce, "TCE_proteins_expressed_all_samples_log2value_add_proteins_sum_intestinal.csv")
write.csv(diffuse.tce, "TCE_proteins_expressed_all_samples_log2value_add_proteins_sum_diffuse.csv")

diff.intes <- read_excel("TCE_diffuse_intestinal_combined.xlsx")
diff.intes <- as.data.frame(diff.intes)
library(ggrepel)
p <- ggplot(data = diff.intes, aes(x = sum, y = NRows, col = Histology)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_reverse(limits = c(20, 0),
                  breaks = c(20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0)) +
  geom_text_repel(data = diff.intes, aes(label = NRows), col = "black", box.padding = 1) +
  theme_classic()
p
dev.print(pdf, "TCE_diffuse_intestinal_point_line_plot.pdf", useDingbats = FALSE)

# T test analysis between intestinal and diffuse
t.test(diff.intes$NRows ~ diff.intes$Histology)  # p-value = 0.08793
shapiro.test(diff.intes$NRows[diff.intes$Histology == "diffuse"])  # p-value = 2.839e-06
shapiro.test(diff.intes$NRows[diff.intes$Histology == "intestinal"])  # p-value = 0.1148
wilcox.test(diff.intes$NRows ~ diff.intes$Histology)  # p-value = 0.004

# Pipe plot for all intestinal and diffuse subtype proteins
fam <- read_excel("pipe_plot_data_for_diffuse.xlsx")
fam <- as.data.frame(fam)
fam$pop <- fam$Nrows / sum(fam$Nrows)
fam <- fam[order(fam$Nrows), ]
fam
fam <- within(fam, Family <- factor(Family, levels = fam[order(fam$Nrows), "Family"]))
mycols <- c("darkorchid", "yellow", "chartreuse3", "cyan", "dodgerblue4", "red", "gray")
ggplot(fam, aes(x = "", y = pop, fill = Family)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = pop, label = Nrows), color = "white") +
  scale_fill_manual(values = mycols) +
  theme_void()
# dev.print(pdf, file.path(dir, "Pipe chart for 13279 proteins in TCE_Diffuse_sum>0.pdf"))


fam <- read_excel("pipe_plot_data_for_intestinal.xlsx")
fam <- as.data.frame(fam)
fam$pop <- fam$Nrows / sum(fam$Nrows)
fam <- fam[order(fam$Nrows), ]
fam
fam <- within(fam, Family <- factor(Family, levels = fam[order(fam$Nrows), "Family"]))
mycols <- c("darkorchid", "yellow", "chartreuse3", "cyan", "dodgerblue4", "red", "gray")
ggplot(fam, aes(x = "", y = pop, fill = Family)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = pop, label = Nrows), color = "white") +
  scale_fill_manual(values = mycols) +
  theme_void()
# dev.print(pdf, file.path(dir, "Pipe chart for 8951 proteins in TCE_Intestinal_sum>0.pdf"))


### NOS and SRC
nos <- read_excel("Supplementary Table-1. Characte_NOS.xlsx")
src <- read_excel('Supplementary Table-1. Characte_SRC.xlsx')
nos.tce <- tceall[, nos$Sample]
src.tce <- tceall[, src$Sample]

nos.tce$sum <- rowSums(nos.tce != 0)
src.tce$sum <- rowSums(src.tce != 0)
write.csv(nos.tce, "TCE_proteins_expressed_all_samples_log2value_add_proteins_sum_NOS.csv")
write.csv(src.tce, "TCE_proteins_expressed_all_samples_log2value_add_proteins_sum_SRC.csv")

nos.src <- read_excel("TCE_NOS_SRC_point_line_plot.xlsx")
nos.src <- as.data.frame(nos.src)
library(ggrepel)
p <- ggplot(data = nos.src, aes(x = sum, y = NRows, col = group)) +
  geom_point(size = 3) +
  geom_line() +
  scale_x_reverse(limits = c(13, 0),
                  breaks = c(12, 10, 8, 6, 4, 2, 0)) +
  geom_text_repel(data = nos.src, aes(label = NRows), col = "black", box.padding = 1) +
  theme_classic()
p
dev.print(pdf, "TCE_NOS_SRC_point_line_plot.pdf", useDingbats = FALSE)

## difference analysis between NOS and SRC
shapiro.test(nos.src$NRows[nos.src$group == "NOS"])  # p-value = 0.0002115
shapiro.test(nos.src$NRows[nos.src$group == "SRC"])  # p-value = 0.0005613
wilcox.test(nos.src$NRows ~ nos.src$group)  # p-value = 0.5267


# Pipe plot for all NOS and SRC subtype proteins
fam <- read_excel("Pipe chart for 12085 proteins in TCE_NOS_sum>0.xlsx")
fam <- as.data.frame(fam)
fam$pop <- fam$Nrows / sum(fam$Nrows)
fam <- fam[order(fam$Nrows), ]
fam
fam <- within(fam, Family <- factor(Family, levels = fam[order(fam$Nrows), "Family"]))
mycols <- c("darkorchid", "yellow", "chartreuse3", "cyan", "dodgerblue4", "red", "gray")
ggplot(fam, aes(x = "", y = pop, fill = Family)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = pop, label = Nrows), color = "white") +
  scale_fill_manual(values = mycols) +
  theme_void()
# dev.print(pdf, file.path(dir, "Pipe chart for 12085 proteins in TCE_NOS_sum>0.pdf"))


fam <- read_excel("Pipe chart for 11399 proteins in TCE_SRC_sum>0.xlsx")
fam <- as.data.frame(fam)
fam$pop <- fam$Nrows / sum(fam$Nrows)
fam <- fam[order(fam$Nrows), ]
fam
fam <- within(fam, Family <- factor(Family, levels = fam[order(fam$Nrows), "Family"]))
mycols <- c("darkorchid", "yellow", "chartreuse3", "cyan", "dodgerblue4", "red", "gray")
ggplot(fam, aes(x = "", y = pop, fill = Family)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = pop, label = Nrows), color = "white") +
  scale_fill_manual(values = mycols) +
  theme_void()
# dev.print(pdf, file.path(dir, "Pipe chart for 11399 proteins in TCE_SRC_sum>0.pdf"))


