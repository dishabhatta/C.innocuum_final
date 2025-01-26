# Paper: Cinnocuum and its diversity
# Author: DB
# Figure 2, Figure S2: pangenome analysis

# Figure 2A : Number of genes vs total no of genomes for all genomes
  
# Figure 2B : No of unique genes vs total no of genomes for all genomes
  
# Figure S2A : No of genes vs total no of genomes for clade 1 & 2
  
# Figure S2B : No of genes vs total no of genomes for clade 1 & 2
  
# Input files :  number_of_genes_in_pan_genome.Rtab, gene_presence_absence.Rtab, number_of_new_genes.Rtab, gene_presence_absence.csv, clade_12_number_of_genes_in_pan_genome.Rtab, clade_12_gene_presence_absence.Rtab, clade_12_gene_presence_absence.csv
  
# Output files: plot figures, no_of_new_genes_vs_genomes.txt, no_of_new_genes_vs_genomes_clade12.txt

# Libraries used in the entire script

library(stringr)
library(stringi)
library(readxl)
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(plyr)
library(vegan)
library(micropan)


### fig 2A
mydata = read.table("/Data/fig_2_S2/number_of_genes_in_pan_genome.Rtab")

data1 <- mydata
colnames(data1) <- 1:232
data2 <- pivot_longer(data1, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
data2$no_of_genomes <- as.numeric(data2$no_of_genomes)

### calculating 
model <- lm(log(data2$no_of_genes) ~ log(data2$no_of_genomes))
summary(model)

## Result:
##Call:
##lm(formula = log(data2$no_of_genes) ~ log(data2$no_of_genomes))

##Residuals:
##  Min       1Q   Median       3Q      Max 
##-0.48256 -0.03533 -0.00524  0.03378  0.23768 

##Coefficients:
##  Estimate Std. Error t value Pr(>|t|)    
##(Intercept)              8.691715   0.007460  1165.1   <2e-16 ***
##  log(data2$no_of_genomes) 0.332236   0.001635   203.3   <2e-16 ***
##  ---
##  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.07533 on 2318 degrees of freedom
##Multiple R-squared:  0.9469,	Adjusted R-squared:  0.9468 
##F-statistic: 4.131e+04 on 1 and 2318 DF,  p-value: < 2.2e-16

#library(micropan)
gene1 <- read.table("/Data/fig_2_S2/gene_presence_absence.Rtab")

gene2 <- gene1
colnames(gene2) <- gene1[1, ]
rownames(gene2) <- gene1$V1
gene3 <- gene2[-1, ]
gene4 <- gene3[ , -1]
gene5 <- t(gene4)
set.seed(1234)
h.est <- heaps(gene5, n.perm = 500) 

## Result of hest
## Named num [1:2] 3254.965 0.794

### plotting Figure 2A; contains all strains

ggplot(data2, aes(x=factor(no_of_genomes), y=no_of_genes)) + 
  #geom_point() +
  geom_boxplot(outlier.shape = NA, fill= "#E8A419", lwd =0.2) +
  geom_smooth(method = 'nls', formula = 'y~a*x^b') +
  # geom_text(x = 600, y = 1, label = power_eqn(DD), parse = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=1), text = element_text(size = 4)) +
  scale_x_discrete(breaks = seq(1, length(data2$no_of_genomes), by =2)) +
  geom_text(x = 40, y= 34000, label= "N = 5953.3835n^0.3322", size = 2) +
  geom_text(x= 50, y = 32000, label = "alpha = 0.794 (Heap's law), p-value < 2e-16", size = 2) +
  ggtitle("Total number of genes vs number of genomes")


data3 <- data2 %>% group_by(no_of_genomes) %>% summarise_at(vars(no_of_genes), list(mean = mean, sd = sd)) %>% as.data.frame()
data3$mino3 <- data3$mean - data3$sd
data3$maxo3 <- data3$mean + data3$sd

ggplot(data3, aes(x=(no_of_genomes), y=mean)) + 
  geom_ribbon(aes(ymin = (mean - sd), ymax = (mean + sd)), alpha = .5, fill = "darkseagreen3", color = "transparent") +
  geom_line(size = 0.5, color = "aquamarine4") +
  #geom_boxplot(outlier.shape = NA, fill= "#E8A419") +
  #geom_point() +
  # geom_text(x = 600, y = 1, label = power_eqn(DD), parse = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=0.3), text = element_text(size = 3)) +
  scale_x_continuous(breaks = seq(1, length(data3$no_of_genomes), by =2)) +
  geom_text(x = 40, y= 34000, label= "N = 5953.3835n^0.3322", size = 2) +
  geom_text(x= 50, y = 32000, label = "alpha = 0.794 (Heap's law), p-value < 2e-16", size = 2) +
  ggtitle("Total number of genes vs number of genomes")



new_genes = read.table(file = "/Data/fig_2_S2/number_of_new_genes.Rtab")

data_new <- new_genes
colnames(data_new) <- 1:232
data_new3 <- as.data.frame(colMeans(data_new))
data_new3$no_of_genomes <- 1:232
colnames(data_new3) <- c('avg_no_of_new_genes', 'no_of_genomes')

write_tsv(data_new3, file = "/Data/fig_2_S2/no_of_new_genes_vs_genomes.txt")
#data_new2 <- pivot_longer(data_new, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
#data_new2$no_of_genomes <- as.numeric(data_new2$no_of_genomes)

#ggplot(data_new3, aes(x=factor(no_of_genomes), y=avg_no_of_new_genes)) + 
#geom_point() +
#  geom_boxplot(outlier.shape = NA, fill= "#E8A419")

### fig 2B
gene_presence_absence <- read_csv("/Data/fig_2_S2/gene_presence_absence.csv")

colnames(gene_presence_absence)[colnames(gene_presence_absence) == "No. isolates"] <- "no_isolates" 
gene_presence_absence$`Avg sequences per isolate`[gene_presence_absence$`Avg sequences per isolate` == 1.00] <- 1
gene_presence_absence1 <- subset(gene_presence_absence, gene_presence_absence$`Avg sequences per isolate` < 1.02)

gene_count <- gene_presence_absence1 %>% group_by(no_isolates) %>% tally()
ggplot(gene_count, aes(x=no_isolates, y=n)) + geom_col(fill = "#E8A419") + theme_classic() +
  scale_x_continuous(breaks = seq(1, 232)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=1), text = element_text(size = 3)) +
  ggtitle("Gene frequency vs number of genomes") +
  xlab("No. of genomes") +
  ylab("No. of genes") +
  geom_text(x = 231, y= 4000, label = "core", size =3)


## Fig S2A

clade12 = read.table("/Data/fig_2_S2/clade12_number_of_genes_in_pan_genome.Rtab")

clade12_2 <- clade12
colnames(clade12_2) <- 1:193
data_clade12 <- pivot_longer(clade12_2, cols = everything(), names_to = "no_of_genomes", values_to = "no_of_genes")
data_clade12$no_of_genomes <- as.numeric(data_clade12$no_of_genomes)


## calculating equation

model_clade12 <- lm(log(data_clade12$no_of_genes) ~ log(data_clade12$no_of_genomes))
summary(model_clade12)

## Result:
##Call:
##  lm(formula = log(data_clade12$no_of_genes) ~ log(data_clade12$no_of_genomes))

##Residuals:
## Min        1Q    Median        3Q       Max 
## -0.209351 -0.014232 -0.002964  0.016947  0.099188 

##Coefficients:
##  Estimate Std. Error t value Pr(>|t|)    
##(Intercept)                     8.5088855  0.0032678  2603.9   <2e-16 ***
##  log(data_clade12$no_of_genomes) 0.3031731  0.0007451   406.9   <2e-16 ***
##  ---
##  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Residual standard error: 0.03113 on 1928 degrees of freedom
##Multiple R-squared:  0.9885,	Adjusted R-squared:  0.9885 
##F-statistic: 1.655e+05 on 1 and 1928 DF,  p-value: < 2.2e-16


## calculating hest
gene12 <- read.table("/Data/fig_2_S2/clade12_gene_presence_absence.Rtab")

gene12_2 <- gene12
colnames(gene12_2) <- gene12[1, ]
rownames(gene12_2) <- gene12$V1
gene12_3 <- gene12_2[-1, ]
gene12_4 <- gene12_3[ , -1]
gene12_5 <- t(gene12_4)
set.seed(1234)
h.est <- heaps(gene12_5, n.perm = 500)

## Result:
## Named num [1:2] 1082.842 0.571

## plot FigS2A

ggplot(data_clade12, aes(x=factor(no_of_genomes), y=no_of_genes)) + 
  #geom_point() +
  geom_boxplot(outlier.shape = NA, fill= "#E8A419", lwd = 0.2) +
  geom_smooth(method = 'nls', formula = 'y~a*x^b') +
  # geom_text(x = 600, y = 1, label = power_eqn(DD), parse = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=0.5), text = element_text(size = 4)) +
  scale_x_discrete(breaks = seq(1, length(data_clade12$no_of_genomes), by =2)) +
  geom_text(x = 40, y= 24000, label= "N = 4958.6336n^0.3031", size = 2) +
  geom_text(x= 50, y = 22000, label = "alpha = 0.571 (Heap's law), p-value < 2e-16", size =2) +
  ggtitle("Total number of genes vs number of genomes")


data2_clade12 <- data_clade12 %>% group_by(no_of_genomes) %>% summarise_at(vars(no_of_genes), list(mean = mean, sd = sd)) %>% as.data.frame()
#data3$mino3 <- data3$mean - data3$sd
#data3$maxo3 <- data3$mean + data3$sd

ggplot(data2_clade12, aes(x=(no_of_genomes), y=mean)) + 
  geom_ribbon(aes(ymin = (mean - sd), ymax = (mean + sd)), alpha = .5, fill = "darkseagreen3", color = "transparent") +
  geom_line(size = 0.5, color = "aquamarine4") +
  #geom_boxplot(outlier.shape = NA, fill= "#E8A419") +
  #geom_point() +
  # geom_text(x = 600, y = 1, label = power_eqn(DD), parse = TRUE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=0.3), text = element_text(size = 3)) +
  scale_x_continuous(breaks = seq(1, 193, by =2)) +
  geom_text(x = 40, y= 24000, label= "N = N = 4958.6336n^0.3031", size = 2) +
  geom_text(x= 50, y = 22000, label = "alpha = 0.571 (Heap's law), p-value < 2e-16", size = 2) +
  ggtitle("Total number of genes vs number of genomes")



# fig S2B

new_genes12 = read.table(file = "/Data/fig_2_S2/clade12_number_of_new_genes.Rtab")

data_new12 <- new_genes12
colnames(data_new12) <- 1:193
data_new12_3 <- as.data.frame(colMeans(data_new12))
data_new12_3$no_of_genomes <- 1:193
colnames(data_new12_3) <- c('avg_no_of_new_genes', 'no_of_genomes')
write_tsv(data_new12_3, file = "/Data/fig_2_S2/no_of_new_genes_vs_genomes_clade12.txt")


gene_presence_absence_12 <- read_csv("/Data/fig_2_S2/gene_presence_absence.csv")
colnames(gene_presence_absence_12)[colnames(gene_presence_absence_12) == "No. isolates"] <- "no_isolates" 

gene_count_12 <- subset(gene_presence_absence_12, gene_presence_absence_12$`Avg sequences per isolate` == 1)
gene_count12 <- gene_count_12 %>% group_by(no_isolates) %>% tally()
ggplot(gene_count12, aes(x=no_isolates, y=n)) + geom_col(fill = "#E8A419") + theme_classic() +
  scale_x_continuous(breaks = seq(1, length(gene_count12$no_isolates), by =2)) +
  theme(axis.text.x = element_text(angle = 90), panel.background = element_rect(colour = "black", size=1), text = element_text(size = 3)) +
  ggtitle("Gene frequency vs number of genomes") +
  xlab("No. of genomes") +
  ylab("No. of unique genes") +
  geom_text(x = 190, y= 4000, label = "core", size = 4)



