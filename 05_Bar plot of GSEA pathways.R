

### Input library and data
setwd("/Users/robert/Dropbox/Revised_MS/ReviewersR2/ResultsR2/Figure-3")
library(tidyverse)
library(qdapTools)
library(readxl)


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


### Bar plot for pathways enriched in groupA
df <- read_excel("GSEA_SigHigh_in_groupA_compared_with_groupBC_in_TCE_23_patients_Hallmark_Sel.xlsx")
df <- as.data.frame(df)
df$qvalue <- as.numeric(df$qvalue)
ggbarplot(df)



df <- read_excel("GSEA_SigHigh_in_groupA_compared_with_groupBC_in_TCE_23_patients_KEGG_Sel.xlsx")
df <- as.data.frame(df)
df$qvalue <- as.numeric(df$qvalue)
ggbarplot(df)



### Bar plot for pathways enriched in groupBC
df <- read_excel("GSEA_SigHigh_in_groupBC_compared_with_groupA_in_TCE_23_patients_GOBP_Sel.xlsx")
df <- as.data.frame(df)
df$qvalue <- as.numeric(df$qvalue)
ggbarplot(df)



df <- read_excel("GSEA_SigHigh_in_groupBC_compared_with_groupA_in_TCE_23_patients_Hallmark_KEGG_Sel.xlsx")
df <- as.data.frame(df)
df$qvalue <- as.numeric(df$qvalue)
ggbarplot(df)





