path <- '~/Documents/Medical School/Research/Boris LGG WXS recombination reads/'
filtered_file <- 'VDJ_TRA+TRB_filtered.csv'
TR_file <-'VDJ_matched_TRA+TRB.csv'
library(tidyverse)

## read in the matching TR reads
VDJ_filtered <- read_csv(paste(path, filtered_file, sep=''))
TR_files <- read_csv(paste(path, TR_file, sep=''))

## isolate the important columns
TR_files <- TR_files %>% select(CDR3, Epitope,`Epitope gene`, `Epitope species`)
VDJ_filtered <- VDJ_filtered[, apply(VDJ_filtered, 2, function(x) !all(is.na(x)))]

##join together based on CDR3
final <- merge(VDJ_filtered, TR_files, by = 'CDR3')
final_CMV <- final %>% filter(`Epitope species` == 'CMV')
write_csv(final_CMV, paste(path, 'VDJ_matched_CMV.csv', sep = ''))
