path <- '~/Documents/Medical School/Research/TCGA-BRCA/'
filtered_file <- ''
TR_file <-'SearchTable-2024-10-09 18_35_45.642.tsv'
TR <- 'TRB'
chimeric <- T
library(tidyverse)

## read in the matching TR reads and pick receptor type
Konrad_unfiltered <- read_csv(paste(path, filtered_file, sep=''))
Konrad_unfiltered <- Konrad_unfiltered %>% filter(receptor == 'TRA' | receptor == 'TRB')
if(!chimeric){
  Konrad_unfiltered <- Konrad_unfiltered  %>% filter(!grepl('chimeric', Filename))
}
total_pt <- unique(Konrad_unfiltered$Filename)
Konrad_filtered <- Konrad_unfiltered
Konrad_filtered <- distinct(Konrad_unfiltered, CDR3, .keep_all = T)
TR_files <- read_tsv(paste(path, TR_file, sep=''))
TR_files <- TR_files %>% filter((Gene == 'TRA' | Gene == 'TRB'))
TR_files <- distinct(TR_files, CDR3, .keep_all = T)

##find read numbers
reads_meta <- Konrad_filtered %>% filter(`Sample Type` == 'Metastatic')
reads_prim <- Konrad_filtered %>% filter(`Sample Type` == 'Primary Tumor')

##find patients numbers
pt_CDR3_remove <- unique(Konrad_filtered$Filename)
pt_meta <- unique((Konrad_filtered %>% filter(`Sample Type` == 'Metastatic'))$Filename)
pt_prim <- unique((Konrad_filtered %>% filter(`Sample Type` == 'Primary Tumor'))$Filename)

## isolate the important columns
TR_files <- TR_files %>% select(CDR3, Epitope,`Epitope gene`, `Epitope species`)
Konrad_filtered <- Konrad_filtered[, apply(Konrad_filtered, 2, function(x) !all(is.na(x)))]

##join together based on CDR3
final <- merge(Konrad_filtered, TR_files, by = 'CDR3')
write.csv(final, paste(path, 'TCGA_matched.csv', sep = ''))

##find the metastatic and primary tumor reads
final_metastatic <- final %>% filter(`Sample Type` == "Metastatic")
final_primary <- final %>% filter(`Sample Type` == "Primary Tumor")

#identify the epitopes associated with each read and with each patient & remove chimeric
final_metastatic %>% count(`Epitope species`)
pt_epitopes_meta <- final_metastatic %>% group_by(Filename) %>% count(`Epitope species`)

final_primary %>% count(`Epitope species`)
pt_epitopes_prim <- final_primary %>% group_by(Filename) %>% count(`Epitope species`)

##Find number of patients with matches
unique(pt_epitopes_prim$Filename)
unique(pt_epitopes_meta$Filename)

## To find unique reads per patient
pt_epitopes_prim %>% group_by(Filename, `Epitope species`) %>% summarise(total_count=n(),.groups = 'drop') %>% count(`Epitope species`)
pt_epitopes_meta %>% group_by(Filename, `Epitope species`) %>% summarise(total_count=n(),.groups = 'drop') %>% count(`Epitope species`)



