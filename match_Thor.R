path <- '~/Documents/Medical School/Research/'
Thor_file <- 'Thorsson/TCGA_mitcr_cdr3_result_BRCA.csv'
TR_file <-'anti viral CDR3s from vdjdb Oct 9/SearchTable-2024-10-09 18_35_45.642.tsv'
TR <- 'TRB'
TR_a <- 'beta'
library(tidyverse)
library(tidyverse)
library(readxl)
library(stringr)

##filter and clean sheet data
Thor_full <- read_csv(paste(path, Thor_file, sep = ''))
Thor_BRCA <- Thor_full %>% filter(Study == 'BRCA')
Thor_BRCA <- Thor_BRCA %>% filter(chain == TR_a)
Thor_BRCA$aaSeq <- stringr::str_replace(Thor_BRCA$aaSeq, '\\*', '')
Thor_BRCA$aaSeq <- stringr::str_replace(Thor_BRCA$aaSeq, '~', '')

##remove duplicates
total_pt <- unique(Thor_BRCA$ParticipantBarcode)
Thor_filtered <- distinct(Thor_BRCA, aaSeq, .keep_all = T)
TR_files <- read_tsv(paste(path, TR_file, sep=''))
TR_files <- distinct(TR_files %>% filter(Gene == TR), CDR3, .keep_all = T)
names(TR_files)[names(TR_files) == 'CDR3'] <- 'aaSeq'

##find read numbers
reads_meta <- Thor_filtered %>% filter(SampleTypeLetterCode == 'TM')
reads_prim <- Thor_filtered %>% filter(SampleTypeLetterCode == 'TP')

##find patients numbers
pt_CDR3_remove <- unique(Thor_filtered$ParticipantBarcode)
pt_meta <- unique((Thor_filtered %>% filter(SampleTypeLetterCode == 'TM'))$ParticipantBarcode)
pt_prim <- unique((Thor_filtered %>% filter(SampleTypeLetterCode == 'TP'))$ParticipantBarcode)

## isolate the important columns
TR_files <- TR_files %>% select(aaSeq, Epitope,`Epitope gene`, `Epitope species`)
Thor_filtered <-Thor_filtered[, apply(Thor_filtered, 2, function(x) !all(is.na(x)))]

##join together based on CDR3
final <- merge(Thor_filtered, TR_files, by = 'aaSeq')

final_metastatic <- final %>% filter(SampleTypeLetterCode == 'TM')
final_primary <- final %>% filter(SampleTypeLetterCode == 'TP')

final_metastatic %>% count(`Epitope species`)
pt_epitopes_meta <- final_metastatic %>% group_by(ParticipantBarcode) %>% count(`Epitope species`)

final_primary %>% count(`Epitope species`)
pt_epitopes_prim <- final_primary %>% group_by(ParticipantBarcode) %>% count(`Epitope species`)

##Find number of patients with matches
pt_final_p <- unique(pt_epitopes_prim$ParticipantBarcode)
pt_final_m <- unique(pt_epitopes_meta$ParticipantBarcode)

## To find unique reads per patient
pt_epitopes_prim %>% group_by(ParticipantBarcode, `Epitope species`) %>% summarise(total_count=n(),.groups = 'drop') %>% count(`Epitope species`)
pt_epitopes_meta %>% group_by(ParticipantBarcode, `Epitope species`) %>% summarise(total_count=n(),.groups = 'drop') %>% count(`Epitope species`)

