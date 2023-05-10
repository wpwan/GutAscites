## time:4-3-2019
## Heatmap with oncoPrint (mutation genes n = 24)
## Select patients n = 22
## Heatmap with genes (n=24) in protein and mRNA levels between pek2 clusters

## Input library
library(tidyverse)
library(qdapTools)
library(readxl)
library(ComplexHeatmap)
library(magrittr)
library(qwraps2)

## Input the patients group information
setwd('/Users/szhao5/Box Sync/Workingdirectory/RNA Gastric Cancer/')
group <- read.csv('RNA and Proteins group information.csv')
setwd('/Users/szhao5/Box Sync/Workingdirectory/DNA Gastric Cancer/Patients in Proteins RNAs')
prpatients <- read.csv('22 overlap patients between all Proteins and RNAs.csv')
group$Patients <- as.character(group$Patients)
prpatients$Mass.spect.number <- as.character(prpatients$Mass.spect.number)
gpr <- left_join(group, prpatients, by = c('Patients' = 'Mass.spect.number'))
gpr2 <- gpr[, c(1:4, 9)]

## Deal with mut data
setwd('/Users/szhao5/Box Sync/Workingdirectory/DNA Gastric Cancer/')
mut <- read_excel('for-clustering.xlsx')
mut <- as.data.frame(mut)
mut <- mut %>% tibble::column_to_rownames('gene.knowngene')
colnames(mut) <- apply(as.data.frame(colnames(mut)), 1, function(x)
                       substring(x[1], gregexpr(pattern = 'I', x[1])[[1]][1]))
colnames(mut) <- gsub('\\IP-58-2', 'IP-058-2', colnames(mut))
colnames(mut) <- gsub('\\IP-63-4', 'IP-063-4', colnames(mut))
colnames(mut) <- gsub('\\IP-0107-1', 'IP-107-1', colnames(mut))
colnames(mut) <- gsub('\\IP-0107-2', 'IP-107-2', colnames(mut))
colnames(mut) <- apply(as.data.frame(colnames(mut)), 1, function(x)
                       substr(x[1], 1, gregexpr(pattern = '-', x[1])[[1]][1]+4))
colnames(mut) <- gsub('-', '', colnames(mut))
colnames(mut) <- gsub('^(.{2})(.*)$', '\\1-\\2', colnames(mut)) #dim=24*40
write.csv(mut, 'for-clustering_colnames-adjusted.csv')
# Look for duplicated patients number ("IP-107" "IP-067")
colnames(mut)[duplicated(colnames(mut))]
# Delete the duplicated columns and re-input
mut <- read.csv('for-clustering_colnames-adjusted.csv')
mut <- mut %>% tibble::column_to_rownames('X')
colnames(mut) <- gsub('\\.', '-', colnames(mut))
t_mut <- matrix2df(t(mut), col1 = 'Patients')

## Mutation data join with group information
gpr2$PTID.x <- as.character(gpr2$PTID.x)
gmut <- left_join(gpr2, t_mut, by = c('PTID.x' = 'Patients'))
gmut <- gmut %>% tibble::column_to_rownames('Patients')
#write.csv(gmut, 'for-clustering_colnames-adjusted_group.csv')
# Delete rows ('G_IPAS7116', 'G_IPAS7117', 'G_IPAS7115')
gmut2 <- gmut[-c(15, 16, 20), ]
#write.csv(gmut2, 'for-clustering_colnames-adjusted_group_18Patients.csv')
gmut2[is.na(gmut2)] <- 0

## Deal with pindel data
pid <- read.csv('New Pindel combined in study1 and study2 of gastric cancer.csv')
mut1 <- read.csv('for-clustering_colnames-adjusted.csv')
pid$gene <- as.character(pid$gene)
mut1$X <- as.character(mut1$X)
pid_select <- inner_join(mut1[, 1:2], pid, by = c('X' = 'gene')) #TAF1, NOTCH3, FANCM are NA
pid24 <- pid_select[, c(1:9, 131)]
pid24$sample_name <- apply(as.data.frame(pid24$sample_name), 1,
                           function(x) substring(x[1], gregexpr(pattern = 'I', x[1])[[1]][1]))
pid24$sample_name <- gsub('\\IP-0107-1', 'IP-107-1', pid24$sample_name)
pid24$sample_name <- gsub('\\IP-0107-2', 'IP-107-2', pid24$sample_name)
pid24$sample_name <- gsub('\\IP-58-2', 'IP-058-2', pid24$sample_name)
pid24$sample_name <- gsub('\\IP-63-4', 'IP-063-4', pid24$sample_name)
pid24$sample_name <- apply(as.data.frame(pid24$sample_name), 1,
                           function(x) substr(x[1], 1, gregexpr(pattern = '-', x[1])[[1]][1] + 3))
df <- pid24[, c(1, 9, 10)]
for(i in levels(df$type)){
df <- pid24[, c(1, 9, 10)]
df <- df %>% filter(type == i) %>% group_by(X, sample_name) %>% summarise(n = n())
df <- spread(df, X, n)
#rownames(df) <- NULL
#df <- df %>% tibble::column_to_rownames('sample_name')
write.table(df, paste0(i, '_table.csv'), sep=',', quote = FALSE)
}
# TD is None

## Group for pindel data selected
insert <- read.csv('I_table.csv')
deletion <- read.csv('D_table.csv')
insert[is.na(insert)] <- 0
deletion[is.na(deletion)] <- 0

## Join with gmut2
insert$sample_name <- as.character(insert$sample_name)
gmut2 <- gmut2 %>% tibble::rownames_to_column()
gmut2_insert <- left_join(gmut2[, c(1, 2, 5)], insert, by = c('PTID.x' = 'sample_name'))
deletion$sample_name <- as.character(deletion$sample_name)
gmut2_deletion <- left_join(gmut2[, c(1, 2, 5)], deletion, by = c('PTID.x' = 'sample_name'))
t_insert <- matrix2df(t(df2matrix(gmut2_insert[, -c(2:3)])), col1 = 'Gene')
t_deletion <- matrix2df(t(df2matrix(gmut2_deletion[, -c(2:3)])), col1 = 'Gene')
t_gmut2 <- matrix2df(t(df2matrix(gmut2[, -c(2:5)])), col1 = 'Gene')
t_gmut2_insert <- left_join(t_gmut2, t_insert, by = 'Gene')
t_gmut2_insert_deleltion <- left_join(t_gmut2_insert, t_deletion, by = 'Gene')
# Insert data with pek2 cluster order
FinalInsert <- t_gmut2_insert_deleltion[, c(1, 20:37)]
FinalDeletion <- t_gmut2_insert_deleltion[, c(1, 38:55)]
#write.csv(FinalInsert, 'FinalInsertData.csv')
#write.csv(FinalDeletion, 'FinalDeletionData.csv')

## Mutation count with R
mutectDataM4 <- read_excel('Mutation_cluster.xlsx', col_names = TRUE)
mutectDataM4$Deleterious <- 'NotDele'
row.names(mutectDataM4) <- paste(mutectDataM4$sample_name, mutectDataM4$chr,
                                 mutectDataM4$start, mutectDataM4$ref_allele, sep = '_')
deleName <- c('mutationassessor_pred', 'cadd_phred', 'sift_pred', 'polyphen2_hvar_pred',
              'lrt_pred', 'mutationtaster_pred', 'fathmm_pred', 'provean_pred', 'metalr_pred')
deleCountM <- matrix(0, nrow = dim(mutectDataM4)[1], ncol = 9)
row.names(deleCountM) <- row.names(mutectDataM4)
deleCountM[mutectDataM4[, deleName] == "D"] <- 1
deleCountM[mutectDataM4$mutationassessor_pred == "H", 1] <- 1
deleCountM[as.numeric(as.character(mutectDataM4$cadd_phred)) > 20, 2] <- 1
deleCountM[is.na(deleCountM)] <- 0
mutectDataM4$Deleterious[rowSums(deleCountM) > 4] <- 'Deleterious'
colnames(mutectDataM4)
mutectDataM4$sample_name2 <- apply(as.data.frame(mutectDataM4$sample_name), 1, function(x)
                             substring(x[1], gregexpr(pattern = 'I', x[1])[[1]][1]))
mutectDataM4$sample_name2 <- gsub('\\IP-0107-2', 'IP-107-2', mutectDataM4$sample_name2)
mutectDataM4$sample_name2 <- gsub('\\IP-58-2', 'IP-058-2', mutectDataM4$sample_name2)
mutectDataM4$sample_name2 <- gsub('\\IP-63-4', 'IP-063-4', mutectDataM4$sample_name2)
mutectDataM4$sample_name2 <- gsub('\\IP-0107-1', 'IP-107-1', mutectDataM4$sample_name2)
#write.csv(mutectDataM4, 'Mutation count for Deleterious missense.csv')


##################### Plot with oncoPrint #############################

## Step1 Form data
#dpek1 <- read.csv('oncoprintDataGastricCancerAscites_PEK2_1.csv')
#dpek2 <- read.csv('oncoprintDataGastricCancerAscites_PEK2_2.csv')
dpek1 <- read.csv('oncoprintDataGastricCancerAscites_PEK2_1_eleted gain and shallow deletion.csv')
dpek2 <- read.csv('oncoprintDataGastricCancerAscites_PEK2_2_deleted gain and shallow deletion.csv')

### Plot with dpek1
dpek1$tumor_name <- gsub('\\G_', '', dpek1$tumor_name)
a <- table(dpek1$gene)
#sampleinfoData <- read_excel('oncoprintDataGastricCancerAscites.xlsx')
sample <- unique(dpek1$tumor_name)
sample <- gsub('\\G_', '', sample)
gene <- unique(dpek1$gene)
mat.u <- matrix(NA, ncol = length(sample), nrow = length(gene))
for(i in 1:length(gene)){
    print(i)
    for(j in 1:length(sample)){
        data2 <- dpek1[dpek1$gene == gene[i] & dpek1$tumor_name == sample[j], ]
        input <- unlist(data2$INPUT)
        if(length(input) > 0){
            INPUT <- paste(input[1], ';', sep = '')
            if(length(input) > 1){
                for(k in 2:length(input)){
                    INPUT <- paste(INPUT, paste(input[k], ';', sep = ''), sep = '')
                }
            }
            mat.u[i, j] <- INPUT
        }
    }
}

mat.u[is.na(mat.u)] <- ''
rownames(mat.u) <- gene
colnames(mat.u) <- sample
alter_fun <- list(background = function(x, y, w, h){
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = '#F0F0F0', col = NA))
    },
    HOMDEL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#021AFF", col = NA))
    },
    HETLOSS = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#73FEFF", col = NA))
    },
    AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#C51A8A", col = NA))
    },
    GAIN = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FB9FB5", col = NA))
    },
    MISSENSE = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#34A02C", col = NA))
    },
    TRUNC = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = 'black', col = NA))
    },
    DelMISSENSE = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#89419D", col = NA))
    }
   )
col = c('HOMDEL' = "#021AFF",
        #'HETLOSS' = "#73FEFF",
        'AMP' = "#C51A8A",
        #'GAIN' = "#FB9FB5",
        'MISSENSE' = "#34A02C",
        'TRUNC' = 'black',
        'DelMISSENSE' = "#89419D")
oncoPrint(mat.u,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun,
          col = col,
          row_names_gp = gpar(fontsize = 8),
          show_column_names = TRUE,
          column_title = "Landscape of Genomic Alteration in Gastric Ascites samples in PEK2_1 group",
          heatmap_legend_param = list(title = "Alterations",
                                      at = c("AMP",
                                             #'GAIN',
                                             "HOMDEL",
                                             #'HETLOSS',
                                             'TRUNC',
                                             'DelMISSENSE',
                                             'MISSENSE'),
                                      labels = c('Amplification',
                                                 #'Gain',
                                                 'Deep deletion',
                                                 #'Shallow deletion',
                                                 'Truncating',
                                                 'Deleterious missense',
                                                 'Missense')))


### Plot with dpek2
dpek2$tumor_name <- gsub('\\G_', '', dpek2$tumor_name)
a <- table(dpek2$gene)
#sampleinfoData <- read_excel('oncoprintDataGastricCancerAscites.xlsx')
sample <- unique(dpek2$tumor_name)
sample <- gsub('\\G_', '', sample)
gene <- unique(dpek2$gene)
mat.u <- matrix(NA, ncol = length(sample), nrow = length(gene))
for(i in 1:length(gene)){
    print(i)
    for(j in 1:length(sample)){
        data2 <- dpek2[dpek2$gene == gene[i] & dpek2$tumor_name == sample[j], ]
        input <- unlist(data2$INPUT)
        if(length(input) > 0){
            INPUT <- paste(input[1], ';', sep = '')
            if(length(input) > 1){
                for(k in 2:length(input)){
                    INPUT <- paste(INPUT, paste(input[k], ';', sep = ''), sep = '')
                }
            }
            mat.u[i, j] <- INPUT
        }
    }
}

mat.u[is.na(mat.u)] <- ''
rownames(mat.u) <- gene
colnames(mat.u) <- sample
alter_fun <- list(background = function(x, y, w, h){
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = '#F0F0F0', col = NA))
    },
    HOMDEL = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#021AFF", col = NA))
    },
    HETLOSS = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#73FEFF", col = NA))
    },
    AMP = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#C51A8A", col = NA))
    },
    GAIN = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#FB9FB5", col = NA))
    },
    MISSENSE = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#34A02C", col = NA))
    },
    TRUNC = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = 'black', col = NA))
    },
    DelMISSENSE = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#89419D", col = NA))
    }
   )
col = c('HOMDEL' = "#021AFF",
        #'HETLOSS' = "#73FEFF",
        'AMP' = "#C51A8A",
        #'GAIN' = "#FB9FB5",
        'MISSENSE' = "#34A02C",
        'TRUNC' = 'black',
        'DelMISSENSE' = "#89419D")
oncoPrint(mat.u,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun,
          col = col,
          row_names_gp = gpar(fontsize = 8),
          show_column_names = TRUE,
          column_title = "Landscape of Genomic Alteration in Gastric Ascites samples in PEK2_2 group",
          heatmap_legend_param = list(title = "Alterations",
                                      at = c("AMP",
                                             #'GAIN',
                                             "HOMDEL",
                                             #'HETLOSS',
                                             'TRUNC',
                                             'DelMISSENSE',
                                             'MISSENSE'),
                                      labels = c('Amplification',
                                                 #'Gain',
                                                 'Deep deletion',
                                                 #'Shallow deletion',
                                                 'Truncating',
                                                 'Deleterious missense',
                                                 'Missense')))


#### Mutation signature plot
##### Make order of samples which is agree with heatmap column names
#oncop <- oncoPrint(mat.u, get_type = function(x) strsplit(x, ";")[[1]],
#                   alter_fun = alter_fun, col = col,
#                   row_names_gp = gpar(fontsize = 8),
#                   column_title = "Landscape of Genomic Alteration in Gastric Ascites samples",
#                   heatmap_legend_param = list(
#                       title = "Alterations",
#                       at = c("AMP",'GAIN', "HOMDEL", 'HETLOSS','TRUNC','DelMISSENSE','MISSENSE'),
#                       labels = c('Amplification', 'Gain', 'Deep deletion', 'Shallow deletion',
#                                  'Truncating', 'Deleterious missense', 'Missense')))
#ht_list <- oncop@ht_list
#sample.order <- (ht_list$matrix_4)@column_order

##### Construct mutation signature
mut <- read_excel('16510JoinedAncotator-filter-v1.xlsx')
mut$ref <- mut$ref_allele
mut$alt <- mut$alt_allele

mut[mut$ref=="A" ,"ref"] <- "T"
mut[mut$ref_allele=="A" & mut$alt_allele=="G", "alt"] <- "C"
mut[mut$ref_allele=="A" & mut$alt_allele=="T", "alt"] <- "A"
mut[mut$ref_allele=="A" & mut$alt_allele=="C", "alt"] <- "G"

mut[mut$ref=="G" ,"ref"] <- "C"
mut[mut$ref_allele=="G" & mut$alt_allele=="A", "alt"] <- "T"
mut[mut$ref_allele=="G" & mut$alt_allele=="T", "alt"] <- "A"
mut[mut$ref_allele=="G" & mut$alt_allele=="C", "alt"] <- "G"


###### Patients names to change
mut$tumor_name.v2 <- gsub('.*I', 'I', mut$tumor_name)
mut$tumor_name.v2 <- gsub('\\IP-0107-2', 'IP-107-2', mut$tumor_name.v2)
mut$tumor_name.v2 <- gsub('\\IP-58-2', 'IP-058-2', mut$tumor_name.v2)
mut$tumor_name.v2 <- gsub('\\IP-63-4', 'IP-063-4',mut$tumor_name.v2)
mut$tumor_name.v2 <- gsub('\\IP-0107-1', 'IP-107-1', mut$tumor_name.v2)
mut$tumor_name.v2 <- apply(as.data.frame(mut$tumor_name.v2), 1,
                           function(x) substr(x[1], 1, gregexpr(pattern = '-', x[1])[[1]][1] + 3))
pek2order <- gmut2[, c(1, 5)]
colnames(pek2order) <- c('Patients', 'tumor_name.v2')
mut <- left_join(pek2order, mut, by = 'tumor_name.v2')
mut.types <- paste(mut$ref, '>', mut$alt, sep = '')
mut.u <- data.frame(mut[, c('Patients', 'ref', 'alt')], mut.types)
mut.u2 <- aggregate(mut.u$Patients, mut.u, length)
sample.order <- c('IPAS0996', 'IPAS0997', 'IPAS0980', 'IPAS0994', 'IPAS0982', 'IPAS0998',
                  'IPAS0995', 'IPAS0979', 'IPAS7102', 'IPAS7100', 'IPAS0999', 'IPAS7104',
                  'IPAS7103', 'IPAS7106', 'IPAS7107', 'IPAS7113', 'IPAS0981', 'IPAS7114')
mut.u2$Patients <- gsub('\\G_', '', mut.u2$Patients)
mut.u2$Patients <- factor(mut.u2$Patients, levels = sample.order)

mut.u2$color <- 'black'
mut.u2[mut.u2$mut.types == 'C>A', 'color'] <- "#1EBFF0"
mut.u2[mut.u2$mut.types == 'C>T', 'color'] <- "#E62725"
mut.u2[mut.u2$mut.types == 'T>A', 'color'] <- "#CBCACB"
mut.u2[mut.u2$mut.types == 'T>C', 'color'] <- "#A1CF64"
mut.u2[mut.u2$mut.types == 'T>G', 'color'] <- "#EDC8C5"
color <- c("#1EBFF0", 'black', "#E62725", "#CBCACB", "#A1CF64", "#EDC8C5")
p2 <- ggplot(mut.u2, aes(x = Patients, y = x, fill = mut.types)) +
    geom_bar(stat = "identity", position = 'fill') +
    scale_fill_manual(values =  color) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="grey20", size=8, angle=90, hjust=.5, vjust=.5))
p2













