# Paper : C. innocum and its diversity
# Figure 6: toxins, virulence factors + blast output
# Author: DB

## Input files

### Figure 6A: Folder: figure_6
#    - pathofact_pred.tsv
#    - tox_lib.tsv
#    - cluster_groups.txt

### Figure 6B: Folder: figure_6
#    - vir_facs_orf_id_14501.tsv
#    - tox_orf_id_14501.tsv
#    - blast_Ref_Cinnocuum_14501.txt; blast_amr_mge_Ref_Cinnocuum_14501.txt; 
#    - Toxin_gene_library_Ref_Cinnocuum_14501_dvf_report.tsv
#    - Ref_Cinnocuum_14501.gtf

## Output files : vir_facs_orf_id_14501.tsv; tox_orf_id_14501.tsv; amr_pathofact_orf_id.tsv; path_vir2.txt; path_amr.txt; path_tox.txt; gi.csv; Cinn14501.xls


rm(list = ls())
dev.off()

# Libraries

library(stringr)
library(stringi)
library(readxl)
library(readr)
library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(plyr)
library(RColorBrewer)
library(dplyr)

# Toxins: use toxin data from Pathofact
# data input, using for loop to catch all the files 
# with all predictions first

#file_list_cinn <- read_tsv(file = "/Data/fig_6_S5/cinn_latest")
#file_list_cinn <- as.vector(file_list_cinn$file_list)
#list.data.tox.cinn<-list()
### file path to overview.txt from run_dbcan
#filepath_tox = "./pathofact_new/Pathofact_"

#for(i in (1:length(file_list_cinn))){
  #print(i)
#  list.data.tox.cinn[[i]] <- read_tsv(paste0(filepath_tox, file_list_cinn[i], '_predictions.tsv'))
#}

### combining all the overview.txt from a list to one large data frame and simultaneously adding their species name

#pathofact_pred <- list.data.tox.cinn[[1]]
#pathofact_pred$genome_name <- "CM1_52_S208"

#for (i in 2:232) {
#  pathofact_pred <- rbind.fill(pathofact_pred, list.data.tox.cinn[[i]])
#  pathofact_pred <- mutate(pathofact_pred, genome_name = replace(genome_name, is.na(genome_name), file_list_cinn[[i]]))
#}
#write_tsv(pathofact_pred, file = "~/Box Sync/SeekatzLab/Presentations/Manuscripts/Ci_Paper/Revision/pathofact_pred.tsv")

#list.data.tox.lib <- list()
#filepath_tox.lib <- "./pathofact_new/Toxin_gene_library_"

#for(i in (1:length(file_list_cinn))){
  #print(i)
#  list.data.tox.lib[[i]] <- read_tsv(paste0(filepath_tox.lib, file_list_cinn[i], '_report.tsv'))
#}

#tox.lib <- list.data.tox.lib[[1]]
#tox.lib$genome_name <- "CM1_52_S208"

#for (i in 2:232) {
#  tox.lib <- rbind.fill(tox.lib, list.data.tox.lib[[i]])
#  tox.lib <- mutate(tox.lib, genome_name = replace(genome_name, is.na(genome_name), file_list_cinn[[i]]))
#}
#write_tsv(tox.lib, file = "~/Box Sync/SeekatzLab/Presentations/Manuscripts/Ci_Paper/Revision/tox_library.tsv")

pathofact_pred <- read_csv(file = "/Data/fig_6_S5/pathofact_pred.tsv")
pathofact_pred$newgenome_name <- sub("^([^_]*_[^_]*)_.*$", "\\1", pathofact_pred$genome_name)
pathofact_pred$newgenome_name <- gsub("\\..*","",pathofact_pred$newgenome_name)
pathofact_pred$newgenome_name[pathofact_pred$genome_name == "Ref_Cinnocuum_14501"] <- "Cinnocuum_14501"
pathofact_pred$newgenome_name[pathofact_pred$genome_name == "Ref_Cinnocuum_LCLUMC"] <- "Cinnocuum_LCLUMC"
pathofact_pred$newgenome_name[pathofact_pred$genome_name == "Ref_Cinnocuum_2959"] <- "Cinnocuum_2959"
pathofact_pred$newgenome_name[pathofact_pred$genome_name == "Ref_Cinnocuum_I46"] <- "Cinnocuum_I46"

secreted_toxins <- subset(pathofact_pred, pathofact_pred$Toxin_confidence_level == "1: Secreted Toxin")
non_secreted_toxins <- subset(pathofact_pred, pathofact_pred$Toxin_confidence_level == "2: Non-secreted Toxin")

secreted_vf <- subset(pathofact_pred, pathofact_pred$Virulence_confidence_level == "1: Secreted Virulence factor")
non_secreted_vf <- subset(pathofact_pred, pathofact_pred$Virulence_confidence_level == "2: Non-secreted Virulence factor")
potential_vf <- subset(pathofact_pred, pathofact_pred$Virulence_confidence_level == "3: Potential Secreted Virulence factor" | pathofact_pred$Virulence_confidence_level == "4: Potential Non-secreted Virulence factor")

count_sec_tox <- secreted_toxins %>% group_by(newgenome_name) %>% tally() ##sum = 1607
count_nonsec_tox <- non_secreted_toxins %>% group_by(newgenome_name) %>% tally() ##sum = 11390
count_sec_vf <- secreted_vf %>% group_by(newgenome_name) %>% tally() ##sum = 32442
count_nonsec_vf <- non_secreted_vf %>% group_by(newgenome_name) %>% tally() ##sum = 32442
count_pot_vf <- potential_vf %>% group_by(newgenome_name) %>% tally() ##sum = 480202

tox.lib <- read_tsv(file = "/Data/fig_6_S5/tox.lib.tsv")
secreted_tox_lib <- inner_join(secreted_toxins, tox.lib[,-2], by = c("genome_name" = "genome_name", "ORF_ID" = "ORF_ID")) %>% distinct()
clust_secreted_tox_lib <- inner_join(secreted_tox_lib, cluster_group1, by = c('genome_name'='isolate'))
cut_secreted_tox_lib <- subset(clust_secreted_tox_lib, clust_secreted_tox_lib$Score > 50)

non_secreted_tox_lib <- inner_join(non_secreted_toxins, tox.lib[,-2], by = c("genome_name" = "genome_name", "ORF_ID" = "ORF_ID")) %>% distinct()
clust_non_secreted_tox_lib <- inner_join(non_secreted_tox_lib, cluster_group1, by = c('genome_name'='isolate'))
cut_non_secreted_tox_lib <- subset(clust_non_secreted_tox_lib, clust_non_secreted_tox_lib$Score > 50) %>% arrange(clade)

full_tox <- rbind(cut_secreted_tox_lib, cut_non_secreted_tox_lib) %>% arrange(clade)

full_sum <- aggregate(full_tox$Score, list(full_tox$clade, full_tox$NAME), FUN=mean)

ggplot(full_sum, aes(x=Group.1, y=Group.2, fill=x)) + geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90), text = element_text(size = 3)) +
  #scale_x_discrete(limits = unique(full_tox$newgenome_name)) +
  scale_fill_gradient(high = "#000080", low = "#E78100") +
  #scale_fill_gradient2(low = "yellow", high= "navy") +
  coord_fixed() + 
  ggtitle("Secreted + Non Secreted toxins; Score > 50")


# Figure 6B: Circos plot for C. innocuum 14501


## subset pathofact_pred to Ref Cinnocuum 14501 to get just the 14501
## make orf_id list to blast against the genome to get coordinates on the genome

### virulence factors
cinn14501_vir <- subset(pathofact_pred, pathofact_pred$newgenome_name == "Cinnocuum_14501")
cinn14501_vir2 <- subset(cinn14501_vir, cinn14501_vir$Virulence_confidence_level == "1: Secreted Virulence factor" | cinn14501_vir$Virulence_confidence_level == "2: Non-secreted Virulence factor")

vir_df <- as.data.frame(cinn14501_vir2$ORF_ID) %>% distinct()
colnames(vir_df) <- "orf_id"

write_tsv(vir_df, file = "/Data/fig_6_S5/vir_facs_orf_id_14501.tsv")

### toxin

cinn14501_tox2 <- subset(cinn14501_vir, cinn14501_vir$Toxin_confidence_level == "1: Secreted Toxin" | cinn14501_vir$Toxin_confidence_level == "2: Non-secreted Toxin")

tox_df <- as.data.frame(cinn14501_tox2$ORF_ID) %>% distinct()
colnames(tox_df) <- "orf_id"

write_tsv(tox_df, file = "/Data/fig_6_S5/tox_orf_id_14501.tsv") 


## generating the blast file from Ref_cinnocuum_14501_dvf_sl.faa to find the genes pathofact in the .faa is talking about; done on Palmetto
amr_mge2 <- read_tsv(file = "/Data/fig_6_S5/antibiotic_resistance.tsv")
amr_mge2 <- subset(amr_mge2, amr_mge2$newgenome_name == "Cinnocuum_14501")
amr_orf_id <- amr_mge2 %>% select(ORF_ID) %>% distinct()

write_tsv(amr_orf_id, file = "/Data/fig_6_S5/amr_pathofact_orf_id.tsv")
#### I need this file for blast against the Ref_Cinnocuum_14501.faa


## compare it to coordinates and produce file for circos to use

## for virulence factors after blast
vir_14501 <- read_tsv(file = "/Data/fig_6_S5/blast_Ref_Cinnocuum_14501.txt", col_names = c("Query_ID", "gtf_id", "Percentage_of_identical_matches", "Alignment_length", "Number_of_mismatches", "Numberofgapopenings", "Startofalignmentinquery", "Endofalignmentinquery", "Startofalignmentinsubject", "Endofalignmentinsubject", "Expected_value", "Bit_score"))
final_vir <- subset(vir_14501, vir_14501$Percentage_of_identical_matches >= 99.00)
final_vir2 <- final_vir %>% select(Query_ID, gtf_id)

### for coordinates of virulence factors from pathofact
gtf_14501 <- read_tsv(file = "/Data/fig_6_S5/Ref_Cinnocuum_14501.gtf", col_names = c("name", "software", "gene_type", "start", "end", "blah1", "orientation", "blah2", "gene_id"))
gtf_14501<- separate(gtf_14501, col = gene_id, into = c("blah3", "gtf_id"), sep = " ")

final_vir3 <- inner_join(final_vir2, gtf_14501, by = "gtf_id")

write_tsv(final_vir3, file = "/Data/fig_6_S5/path_vir2.txt") # <------ Goes into circos

## after blast antibiotic resistance

amr_14501 <- read_tsv(file = "/Data/fig_6_S5/blast_amr_mge_Ref_Cinnocuum_14501.txt", col_names = c("Query_ID", "gtf_id", "Percentage_of_identical_matches", "Alignment_length", "Number_of_mismatches", "Numberofgapopenings", "Startofalignmentinquery", "Endofalignmentinquery", "Startofalignmentinsubject", "Endofalignmentinsubject", "Expected_value", "Bit_score"))
final_amr <- subset(amr_14501, amr_14501$Percentage_of_identical_matches >= 99.00)
final_amr2 <- final_amr %>% select(Query_ID, gtf_id)


final_amr3 <- inner_join(final_amr2, gtf_14501, by = "gtf_id")
# this goes into circos
write_tsv(final_amr3, file = "/Data/fig_6_S5/path_amr.txt")

## after blast toxins

tox_14501 <- read_tsv(file = "/Data/fig_6_S5/blast_tox_Ref_Cinnocuum_14501.txt", col_names = c("Query_ID", "gtf_id", "Percentage_of_identical_matches", "Alignment_length", "Number_of_mismatches", "Numberofgapopenings", "Startofalignmentinquery", "Endofalignmentinquery", "Startofalignmentinsubject", "Endofalignmentinsubject", "Expected_value", "Bit_score"))
final_tox <- subset(tox_14501, tox_14501$Percentage_of_identical_matches >= 99.00)
final_tox2 <- final_tox %>% select(Query_ID, gtf_id)
final_tox3 <- inner_join(final_tox2, gtf_14501, by = "gtf_id")

# this goes into circos
write_tsv(final_tox3, file = "/Data/fig_6_S5/path_tox.txt") 


## for genomic islands
cinn14501_gi <- read_excel(path = "/Data/fig_6_S5/Cinn14501.xls")

cinn14501_gi2 <- cinn14501_gi %>% select(`Island start`, `Island end`) %>% distinct()

#this goes into circos after modifications

write_csv(cinn14501_gi2, file = "/Data/fig_6_S5/gi.csv")

## mods to output files before circos
#### remove the extra columns in Excel, put in chr1 on the first column and fill_color for the last column forr all output files before running through circos
### with given output files + circos config set to user preference visualization






### Fig S5

#filepath_amr = "./AMR_MGE_prediction_"

#list.data.amr <- list()
#for(i in (1:length(file_list_cinn))){
  #  #print(i)
#  list.data.amr[[i]] <- read_tsv(paste0(filepath_amr, file_list_cinn[i], '_report.tsv'))
#}

# combining all the predictions.tsv from a list to one large data frame and simultaneously adding their species name

#all.amr.mge <- list.data.amr[[1]]
#all.amr.mge$genome_name <- "CM1_52_S208"


#for (i in 2:232) {
#  all.amr.mge <- rbind.fill(all.amr.mge, list.data.amr[[i]])
#  all.amr.mge <- mutate(all.amr.mge, genome_name = replace(genome_name, is.na(genome_name), file_list_cinn[[i]]))
#}


#all.amr.mge$newgenome_name <- sub("^([^_]*_[^_]*)_.*$", "\\1", all.amr.mge$genome_name)
#all.amr.mge$newgenome_name <- gsub("\\..*","",all.amr.mge$newgenome_name)
#all.amr.mge$newgenome_name[all.amr.mge$genome_name == "Ref_Cinnocuum_14501"] <- "Cinnocuum_14501"
#all.amr.mge$newgenome_name[all.amr.mge$genome_name == "Ref_Cinnocuum_LCLUMC"] <- "Cinnocuum_LCLUMC"
#all.amr.mge$newgenome_name[all.amr.mge$genome_name == "Ref_Cinnocuum_2959"] <- "Cinnocuum_2959"
#all.amr.mge$newgenome_name[all.amr.mge$genome_name == "Ref_Cinnocuum_I46"] <- "Cinnocuum_I46"

#write_tsv(all.amr.mge, file = "~/Box Sync/SeekatzLab/Presentations/Manuscripts/Ci_Paper/Revision/antibiotic_resistance.tsv")

all.amr.mge <- read_tsv(file = "/Data/fig_6_S5/all.amr.mge.tsv")
m.amr_mge <- inner_join(all.amr.mge, cluster_group1, by = c('genome_name' = 'isolate'))
m.amr_mge2 <- subset(m.amr_mge, !(m.amr_mge$AMR_category == "-"))
#plot
ggplot(m.amr_mge2, aes(x= clade, y = AMR_category)) + geom_tile() + coord_fixed() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


