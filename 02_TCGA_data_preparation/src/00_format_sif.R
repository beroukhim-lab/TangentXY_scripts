library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

sif <- sif %>%
  rename(Gender=gender) %>%
  rename(Age=age) %>%
  rename(Ethnicity=ethnicity) %>%
  rename(Race=race) %>%
  mutate(Gender=case_when(Gender=='na' ~ NA, TRUE ~ str_to_title(Gender))) %>%
  mutate(Age=as.numeric(Age)) %>%
  mutate(tss=case_when(is.na(tss) ~ 'NA', TRUE ~ tss))

saveRDS(sif, file=here('02_TCGA_data_preparation/output/00_format_sif', 'sif.rds'), compress=FALSE)
