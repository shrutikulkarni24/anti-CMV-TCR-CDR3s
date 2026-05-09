path <- '~/Documents/Medical School/Research/'
filtered_file <- ''
TR_file <- 'anti viral CDR3s from vdjdb Oct 9/SearchTable-2024-10-09 18_35_45.642.tsv'
TRA <- 'TRA'
TR <- 'TRB'
chimeric <- F
wxs <- F
library(tidyverse)

## read in the matching TR reads and pick receptor type
Konrad_unfiltered <- read_csv(paste(path, filtered_file, sep=''))
Konrad_unfiltered <- Konrad_unfiltered %>% filter(Receptor == TR | Receptor == paste(TR, '_UM', sep = ''))
if(!chimeric){
  Konrad_unfiltered <- Konrad_unfiltered  %>% filter(!grepl('chimeric', Filename))
}
if(wxs){
  Konrad_unfiltered <- Konrad_unfiltered %>% filter(`Assay Type` == 'WXS')
} else{
  Konrad_unfiltered <- Konrad_unfiltered %>% filter(`Assay Type` == 'RNA-Seq')
}
total_pt <- unique(Konrad_unfiltered$submitted_subject_id)
Konrad_filtered <- Konrad_unfiltered
Konrad_filtered <- distinct(Konrad_unfiltered, CDR3, .keep_all = T)
TR_files <- read_tsv(paste(path, TR_file, sep=''))
TR_files <- distinct(TR_files %>% filter(Gene == TR| Gene == TRA), CDR3, .keep_all = T)

##find read numbers
reads_meta <- Konrad_filtered %>% filter(`Sample Type` == 'Metastatic')
reads_meta <- rbind(reads_meta, Konrad_filtered %>% filter(is.na(samples.sample_type)) %>% filter(body_site != 'Breast'))
reads_prim <- Konrad_filtered %>% filter(`Sample Type` == 'Primary Tumor')
reads_prim <- rbind(reads_prim, Konrad_filtered %>% filter(is.na(samples.sample_type)) %>% filter(body_site == 'Breast'))

##find patients numbers
pt_CDR3_remove <- unique(Konrad_filtered$submitted_subject_id)
pt_meta <- unique((reads_meta$submitted_subject_id))
pt_prim <- unique(reads_prim$submitted_subject_id)

## isolate the important columns
TR_files <- TR_files %>% select(CDR3, Epitope,`Epitope gene`, `Epitope species`)
Konrad_filtered <- Konrad_filtered[, apply(Konrad_filtered, 2, function(x) !all(is.na(x)))]

##join together based on CDR3
final <- merge(Konrad_filtered, TR_files, by = 'CDR3')
write.csv(final, paste(path, 'TCGA_matched_TRA.csv', sep = ''))

##find the metastatic and primary tumor reads
final_metastatic <- final %>% filter(`Sample Type`== 'Metastatic')
final_metastatic <- rbind(final_metastatic, final %>% filter(is.na(samples.sample_type)) %>% filter(body_site != 'Breast'))
final_primary <- final %>%  filter(`Sample Type` == 'Primary Tumor')
final_primary <- rbind(final_primary, final %>% filter(is.na(samples.sample_type)) %>% filter(body_site == 'Breast'))

#identify the epitopes associated with each read and with each patient & remove chimeric
final_metastatic %>% count(`Epitope species`)
pt_epitopes_meta <- final_metastatic %>% group_by(submitted_subject_id) %>% count(`Epitope species`)

final_primary %>% count(`Epitope species`)
pt_epitopes_prim <- final_primary %>% group_by(submitted_subject_id) %>% count(`Epitope species`)

write_csv(pt_epitopes_meta, paste(path, 'PT_meta_wxs_TRB_matches_05-27-25.csv', sep = ''))
write_csv(pt_epitopes_prim, paste(path, 'PT_prim_wxs_TRB_matches_05-27-25.csv', sep = ''))

##Find number of patients with matches
unique(pt_epitopes_prim$submitted_subject_id)
unique(pt_epitopes_meta$submitted_subject_id)

## To find unique reads per patient
pt_epitopes_prim %>% group_by(submitted_subject_id, `Epitope species`) %>% summarise(total_count=n(),.groups = 'drop') %>% count(`Epitope species`)
pt_epitopes_meta %>% group_by(submitted_subject_id, `Epitope species`) %>% summarise(total_count=n(),.groups = 'drop') %>% count(`Epitope species`)