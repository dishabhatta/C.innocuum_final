# Paper: Cinnocuum and its diversity
# Author: DB
# Figure 3: core genome tree, PCA plot , cog % in pangenome, functional enrichment
# Figure S3: module completeness


# Fig 3A: core genome tree- from graphlan; 

# Fig 3B: PCoA plot; bray-curtis

# Fig 3C: percentage of COG assignments in core, softshell, accessory-singletons

# Fig 3D: functional enrichment

# Fig S3ABC:  module completeness

# Input files:
## Fig 3B: cinn_latest.txt; all_cinn_pcoa.tsv; cog-20.def.txt
## Fig 3C: NEWCINNWGS_gene_clusters_summary.txt; gene_cluster_summary_lite.txt
## Fig 3D: enriched_clade_module,
## Fig S3: full_modules.txt


# Output files: plot files, all_cinn_pcoa.tsv, gene_cluster_summary_lite.txt, full_summary.txt


# Libraries

library(stringr)
library(stringi)
library(readxl)
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(philentropy)
library(vegan)
library(ape)



### Fig 3B: PCoA plot of gene presence and absence with bray curtis

#file_list_cinn <- read_tsv(file = "./cinn_latest")
#file_list_cinn <- as.vector(file_list_cinn$file_list)
### create an empty list which will be populated with all the overview.txt
#list.data.pcoa.cinn<-list()
### file path to overview.txt from run_dbcan
#filepath_pcoa_cinn = "./prokka_tsv/"

#for(i in (1:length(file_list_cinn))){
  #print(i)
#  list.data.pcoa.cinn[[i]] <- read_tsv(paste0(filepath_pcoa_cinn, file_list_cinn[i], '.tsv'))
#}

### combining all the overview.txt from a list to one large data frame and simultaneously adding their species name

#all_cinn_pcoa <- list.data.pcoa.cinn[[1]]
#all_cinn_pcoa$genome_name <- "CM1_52_S208"

#for (i in 2:232) {
#  all_cinn_pcoa <- rbind.fill(all_cinn_pcoa, list.data.pcoa.cinn[[i]])
#  all_cinn_pcoa <- mutate(all_cinn_pcoa, genome_name = replace(genome_name, is.na(genome_name), file_list_cinn[[i]]))
#}

#write_tsv(all_cinn_pcoa, file = "/Data/fig_3_S3/all_cinn_pcoa.tsv")

all_cinn_pcoa <- read_tsv(file = "/Data/fig_3_S3/all_cinn_pcoa.tsv")

prokka_cog_cinn <- subset(all_cinn_pcoa, !is.na(all_cinn_pcoa$COG))
pcoa_cog_cinn <- select(prokka_cog_cinn, COG, genome_name) %>% distinct()
pcoa_cog_cinn$prab <- 1
pcoa_cog_cinn2 <- pcoa_cog_cinn %>% pivot_wider(names_from = COG, values_from = prab, values_fill = 0)

pcoa_cog_cinn3 <- pcoa_cog_cinn2 

pcoa_cog_cinn3 <- pcoa_cog_cinn3[,-1]
rownames(pcoa_cog_cinn3) <- pcoa_cog_cinn2$genome_name

pcoa_cog_cinn4 <- as.matrix(pcoa_cog_cinn3)

#library(vegan)
bray_dist_pcoa_cinn <- (vegdist(pcoa_cog_cinn4, method = "bray", binary = TRUE))

#library(ape)
pcoa_vals_cinn <- pcoa(bray_dist_pcoa_cinn)
pcoa_axes_cinn <- pcoa_vals_cinn$vectors[,1:2]
pcoa_axes_cinn2 <- as.data.frame(pcoa_axes_cinn)
pcoa_axes_cinn2$genome_name <- (rownames(pcoa_axes_cinn))

pcoa_axes_cinn2$newgenome_name <- sub("^([^_]*_[^_]*)_.*$", "\\1", pcoa_axes_cinn2$genome_name)
pcoa_axes_cinn2$newgenome_name <- gsub("\\..*","",pcoa_axes_cinn2$newgenome_name)
pcoa_axes_cinn2$newgenome_name[pcoa_axes_cinn2$genome_name == "Ref_Cinnocuum_14501"] <- "Cinnocuum_14501"
pcoa_axes_cinn2$newgenome_name[pcoa_axes_cinn2$genome_name == "Ref_Cinnocuum_LCLUMC"] <- "Cinnocuum_LCLUMC"
pcoa_axes_cinn2$newgenome_name[pcoa_axes_cinn2$genome_name == "Ref_Cinnocuum_2959"] <- "Cinnocuum_2959"
pcoa_axes_cinn2$newgenome_name[pcoa_axes_cinn2$genome_name == "Ref_Cinnocuum_I46"] <- "Cinnocuum_I46"

cluster_group <- read_tsv(file = "/Data/fig_3_S3/cluster_groups.txt")
cluster_group1 <- cluster_group %>%  arrange(clade)

cluster_group1$color[cluster_group1$clade == "clade_1"] <- "#E8A419"
cluster_group1$color[cluster_group1$clade == "clade_2"] <- "#9FC095"
cluster_group1$color[cluster_group1$clade == "clade_3"] <- "#3B99B1"
cluster_group1$color[cluster_group1$clade == "clade_4"] <- "#F5191C"

pcoa_axes_meta_cinn <- inner_join(pcoa_axes_cinn2, cluster_group1, by = c('genome_name' = 'isolate'))

pcoa_eigen_cinn <- as.data.frame(pcoa_vals_cinn$values$Eigenvalues)
prop_var_cinn <- pcoa_eigen_cinn$`pcoa_vals_cinn$values$Eigenvalues`/sum(pcoa_eigen_cinn$`pcoa_vals_cinn$values$Eigenvalues`) * 100

pcoa_axes_meta_cinn$isolate_id <- "0"
pcoa_axes_meta_cinn$isolate_id[grep("^CM", pcoa_axes_meta_cinn$newgenome_name)] <- "this_study"
pcoa_axes_meta_cinn$isolate_id[grep("^d22", pcoa_axes_meta_cinn$newgenome_name)] <- "this_study"
pcoa_axes_meta_cinn$isolate_id[pcoa_axes_meta_cinn$isolate_id == "0"] <- "other"

ggplot(pcoa_axes_meta_cinn, aes(x=Axis.1, y=Axis.2, fill=clade, color=isolate_id)) + geom_point(shape = 21) + 
  theme_classic() +
  scale_fill_manual(breaks = pcoa_axes_meta_cinn$clade, values =  pcoa_axes_meta_cinn$color) +
  scale_color_manual(breaks = c('this_study', 'other'), values = c("black", rgb(0, 0, 0, alpha=0))) +
  #scale_color_manual(values =) +
  theme(panel.grid = element_blank()) + xlab("Axis.1 (60.42%)") + ylab("Axis.2 (11.2%)") +
  theme(text = element_text(size = 5), panel.background = element_rect(colour = "black", size=0.3), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  coord_fixed() +
  scale_x_continuous(limits = c(-0.04,0.08)) +
  scale_y_continuous(limits = c(-0.04,0.08)) 

### Figure 3C

## soft-core = soft-shell core

#anvio_pan_summary <- read_tsv(file = "./anvio/CINNLATESTWGS_gene_clusters_summary.txt")

#anvio_pan_summary1 <-  anvio_pan_summary %>% select(-c(aa_sequence))
#anvio_pan_summary1$bin_name[is.na(anvio_pan_summary1$bin_name)] <- "Core"
#write_tsv(anvio_pan_summary1, file = "./anvio/gene_cluster_summary_lite.txt")

pan_sum <- read_tsv(file = "/Data/fig_3_S3/gene_cluster_summary_lite.txt") 
pan_sum$COG20_CATEGORY_ACC <- gsub("\\|.*", "", pan_sum$COG20_CATEGORY_ACC)
pan_sum$COG20_CATEGORY <- gsub("\\|.*", "", pan_sum$COG20_CATEGORY)

summ_anvio <- pan_sum %>% group_by(bin_name, COG20_CATEGORY_ACC, COG20_CATEGORY) %>% tally()

summ_anvio1 <- summ_anvio
summ_anvio1$COG20_CATEGORY_ACC[is.na(summ_anvio1$COG20_CATEGORY_ACC)] <- "Not_Assigned"
summ_anvio1$COG20_CATEGORY[is.na(summ_anvio1$COG20_CATEGORY)] <- "Not_Assigned"

summ_anvio1$newcog <- summ_anvio1$COG20_CATEGORY_ACC

core_summary <- subset(summ_anvio1, summ_anvio1$bin_name == "Core")
#core_summary <- core_summary %>% group_by(bin_name, COG20_CATEGORY, newcog) %>% tally()
core_summary$perc <- round((core_summary$n/sum(core_summary$n)) *100, digits = 6)

accsin_summary <- subset(summ_anvio1, summ_anvio1$bin_name == "Accessory_singleton")
#accsin_summary <- accsin_summary %>% group_by(bin_name, newcog) %>% tally()
accsin_summary$perc <- round((accsin_summary$n/sum(accsin_summary$n)) *100, digits = 6)

soft_summary <- subset(summ_anvio1, summ_anvio1$bin_name == "Soft_shell")
#soft_summary <- soft_summary %>% group_by(bin_name, newcog) %>% tally()
soft_summary$perc <- round((soft_summary$n/sum(soft_summary$n)) *100, digits = 6)

final_summary <- rbind(core_summary, accsin_summary, soft_summary)

write_tsv(final_summary, "/Data/fig_3_S3/full_summary.txt")

## making the color palette and sorting the data

color_pal <- as.data.frame(unique(final_summary$newcog))
colnames(color_pal) <- "newcog"
color_pal$color[color_pal$newcog == "R" | color_pal$newcog == "S" | color_pal$newcog == "W" | color_pal$newcog == "Not_Assigned"] <- "grey67"
color_pal$color[color_pal$newcog == "J" | color_pal$newcog == "K" | color_pal$newcog == "L"] <- "#A71B4B"
color_pal$color[color_pal$newcog == "D" | color_pal$newcog == "V" | color_pal$newcog == "T" | color_pal$newcog == "M" | color_pal$newcog == "N" | color_pal$newcog == "O" | color_pal$newcog == "U"] <- "#D04939"
color_pal$color[color_pal$newcog == "C"] <- "#EB7803"
color_pal$color[color_pal$newcog == "G"] <- "#F9BC53"
color_pal$color[color_pal$newcog == "E"] <- "#FEF1A6"
color_pal$color[color_pal$newcog == "F"] <- "#E2F8B5"
color_pal$color[color_pal$newcog == "H"] <- "#9CE5AD"
color_pal$color[color_pal$newcog == "I"] <- "#43CBB1"
color_pal$color[color_pal$newcog == "P"] <- "#00AAB6"
color_pal$color[color_pal$newcog == "Q"] <- "#0080B2"
color_pal$color[color_pal$newcog == "X"] <- "#584B9F"

fin_sum <- inner_join(final_summary, color_pal, by = "newcog")
ord <- c("Not_Assigned", "R", "S", "W", "J", "K", "L", "D", "V", "T", "M", "N", "O", "U", "C", "G", "E", "F", "H", "I", "P", "Q", "X")
fin_sum1 <- fin_sum %>% arrange(factor(newcog, levels = ord))

ggplot(fin_sum1, aes(x=factor(bin_name, levels = c('Core', 'Soft_shell', 'Accessory_singleton')), y=perc, fill =factor(newcog, levels = ord))) + geom_col(width = 0.4) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size=0.5), text = element_text(size = 5), legend.key.size = unit(5, "mm")) +
  scale_fill_manual(breaks = fin_sum1$newcog, values = fin_sum1$color)

### Figure 3D

## read in data files for functional enrichment

enriched_KEGG_Modules <- read_tsv(file = "/Data/fig_3_S3/enriched_clade_module.txt")

## data modification of kegg_modules

enr_mods <- enriched_KEGG_Modules %>% subset(!(is.na(associated_groups)))
enr_mods1 <- enr_mods[ ,1:11]
enr_mods2 <- pivot_longer(enr_mods1, cols = !c("KEGG_Module", "enrichment_score", "unadjusted_p_value", "adjusted_q_value", "associated_groups", "accession", "gene_clusters_ids"), names_to = "clade", values_to = "fraction")
enr_mods3 <- separate_rows(enr_mods2, KEGG_Module, sep = "!!!")
enr_mods4 <- subset(enr_mods3, enr_mods3$adjusted_q_value < 0.05)

## ggplot

ggplot(enr_mods4, aes(x=clade, y=KEGG_Module, fill=fraction)) + geom_tile(color = "black") + theme_bw() +
  scale_fill_gradient2(high = "navyblue", low = "white") + coord_fixed(ratio = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 6), legend.key.size = unit(4, 'mm')) +
  scale_y_discrete(limits = unique(enr_mods4$KEGG_Module))


### Figure S3

# Data input

est_met <- read_tsv(file = "/Data/fig_3_S3/full_modules.txt")

# Data modification

est_met1 <- est_met
est_met1$isolate <- sub("^([^_]*_[^_]*)_.*$", "\\1", est_met1$genome_name)
est_met1$isolate <- gsub("\\..*","",est_met1$isolate)
est_met1$isolate[est_met1$genome_name == "Ref_Cinnocuum_14501"] <- "Cinnocuum_14501"
est_met1$isolate[est_met1$genome_name == "Ref_Cinnocuum_LCLUMC"] <- "Cinnocuum_LCLUMC"
est_met1$isolate[est_met1$genome_name == "Ref_Cinnocuum_2959"] <- "Cinnocuum_2959"
est_met1$isolate[est_met1$genome_name == "Ref_Cinnocuum_I46"] <- "Cinnocuum_I46"

cluster_group <- read_tsv(file = "/Data/fig_3_S3/cluster_groups.txt")
cluster_group1 <- cluster_group %>%  arrange(clade)

cluster_group1$color[cluster_group1$clade == "clade_1"] <- "#E8A419"
cluster_group1$color[cluster_group1$clade == "clade_2"] <- "#9FC095"
cluster_group1$color[cluster_group1$clade == "clade_3"] <- "#3B99B1"
cluster_group1$color[cluster_group1$clade == "clade_4"] <- "#F5191C"

# Data plotting

merge_met <- inner_join(cluster_group1, est_met1, by= c('isolate'='genome_name'))

merge_met_carb <- subset(merge_met, merge_met$module_category == "Carbohydrate metabolism")
avg_carb <- merge_met_carb %>% select(clade, module_name, module_completeness)
avg_carb_perc <- aggregate(avg_carb$module_completeness, list(avg_carb$clade, avg_carb$module_name), FUN=mean)
avg_carb_perc$module_present[avg_carb_perc$x >= 0.75] <- "TRUE"
avg_carb_perc$module_present[is.na(avg_carb_perc$module_present)] <- "FALSE"
ggplot(avg_carb_perc, aes(x= Group.1, y = Group.2, fill= module_present)) +geom_tile(color = "black") + theme_bw() +
  scale_fill_manual(values = c("gold", "navy")) + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 4)) +
  ggtitle("Carbohydrate")

merge_met_aa <- subset(merge_met, merge_met$module_category == "Amino acid metabolism")
avg_aa <- merge_met_aa %>% select(clade, module_name, module_completeness)
avg_aa_perc <- aggregate(avg_aa$module_completeness, list(avg_aa$clade, avg_aa$module_name), FUN=mean)
avg_aa_perc$module_present[avg_aa_perc$x >= 0.75] <- "TRUE"
avg_aa_perc$module_present[is.na(avg_aa_perc$module_present)] <- "FALSE"
ggplot(avg_aa_perc, aes(x= Group.1, y = Group.2, fill= module_present)) +geom_tile(color = "black") + theme_bw() +
  scale_fill_manual(values = c("gold","navy")) + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 4)) +
  ggtitle("Amino acid")


merge_met_lipid <- subset(merge_met, merge_met$module_category == "Lipid metabolism")
avg_lipid <- merge_met_lipid %>% select(clade, module_name, module_completeness)
avg_lipid_perc <- aggregate(avg_lipid$module_completeness, list(avg_lipid$clade, avg_lipid$module_name), FUN=mean)
avg_lipid_perc$module_present[avg_lipid_perc$x >= 0.75] <- "TRUE"
avg_lipid_perc$module_present[is.na(avg_lipid_perc$module_present)] <- "FALSE"
ggplot(avg_lipid_perc, aes(x= Group.1, y = Group.2, fill= module_present)) +geom_tile(color = "black") + theme_bw() +
  scale_fill_manual(values = c("gold","navy")) + coord_fixed() + theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 4)) +
  ggtitle("Lipid")

