path <- '~/Documents/Medical School/Research/CMI-MBC/'
Konrad_file <- 'CMI-MBC VDJ_Recoveries Konrad.csv'
gdc_sample <-'gdc_sample_sheet-2021-05-13.tsv.xlsx'
library(tidyverse)
library(readxl)

Konrad_full <- read_csv(paste(path, Konrad_file, sep = ''))
Konrad_nochim <- Konrad_full %>% filter(!grepl('chimeric', Filename))
Konrad_wxs <- Konrad_full %>% filter(grepl('wxs', Filename))
Konrad_nochim$Filename <- gsub('.tsv', '', as.character(Konrad_nochim$Filename))
Konrad_wxs$Filename <- gsub('.tsv', '', as.character(Konrad_wxs$Filename))

receptors <- c('TRB_UM','TRB')
Konrad_nochim_recep <- Konrad_nochim %>% filter(Receptor %in% receptors)
Konrad_wxs_recep <- Konrad_wxs %>% filter(Receptor %in% receptors)

gdc_full <- read_excel(paste(path, gdc_sample, sep = ''))
gdc_full <- rename(gdc_full, Filename = `File Name`)

Konrad_nochim_filenames <- Konrad_nochim_recep %>% select(Filename)
Konrad_wxs_filenames <- Konrad_wxs_recep %>% select(Filename)
gdc_filenames <- gdc_full %>% select(`Filename`)

matching_nochim_Konrad <- Konrad_nochim_recep[Konrad_nochim_recep$Filename %in% gdc_full$Filename, ]
matching_wxs_Konrad <- Konrad_wxs_recep[Konrad_wxs_recep$Filename %in% gdc_full$Filename, ]
nt_gdc <- gdc_full %>% select(c(Filename, `Sample Type`))

final_nochim <- left_join(matching_nochim_Konrad, nt_gdc, by = 'Filename')
final_wxs <-left_join(matching_wxs_Konrad, nt_gdc, by = 'Filename')
final_nochim <- final_nochim %>% select(-c(`Unnamed: 0`, `Unnamed: 0.1`))
final_wxs <- final_wxs %>% select(-c(`Unnamed: 0`, `Unnamed: 0.1`))
write.csv(final_wxs, paste(path, 'filtered_Konrad_TRB.csv', sep = ''))


