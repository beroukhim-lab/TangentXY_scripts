library(tidyverse)
library(here)

library(rvest)

# sif <- read_csv(here('05_CCLE_data_preparation/data', 'internal-22q2_v15-sample-info.csv')) # Downloaded from Skyros (https://cds.team/depmap/)
sif <- read_csv(here('05_CCLE_data_preparation/data', 'Model.csv')) # Downloaded from DepMap (https://depmap.org/portal/download/all/)

## Firstly, retrieve ATCC IDs and DSMZ IDs from cellosaurus
for (i in 1:nrow(sif)) {
  rrid <- sif[i, 'RRID']
  url <- paste0('https://web.expasy.org/cellosaurus/', rrid)
  html <- read_html(url)
  td <- html %>%
    html_nodes('td') %>%
    html_text()

  crl <- td[grepl(' CRL-', td)]
  acc <- td[grepl(' ACC-', td)]

  atcc <- crl[grepl('ATCC;', crl)]
  dsmz <- acc[grepl('DSMZ;', acc)]

  atcc.id <- str_extract(atcc, 'CRL-[0-9]*')
  dsmz.id <- str_extract(dsmz, 'ACC-[0-9]*')

  if (length(atcc.id)!=0) {
    sif[i, 'ATCC.ID'] <- atcc.id[1]
  }
  if (length(dsmz.id)!=0) {
    sif[i, 'DSMZ.ID'] <- dsmz.id[1]
  }

  print(paste0(i, ' ATCC:', atcc.id, ', DSMZ:', dsmz.id))
}

## Secondly, retrieve karyotype information from ATCC and DSMZ
for (i in 1:nrow(sif)) {
  ## Karyotype from ATCC
  atcc.id <- sif[i, 'ATCC.ID'] %>% tolower()
  if (!is.na(atcc.id)) {
    url <- paste0('https://www.atcc.org/products/', atcc.id)

    html <- try(read_html(url), silent=TRUE)
    if (class(html)=='try-error') {
      next
    } else {
      title <- html %>%
        html_nodes('[class="product-information__title"]') %>%
        html_text()
      data <- html %>%
        html_nodes('[class="product-information__data"]') %>%
        html_text()
      if (length(title)==length(data)) {
        df <- data.frame(title=title, data=data)
        karyotype <- df %>%
          filter(title=='Karyotype') %>%
          pull(data) %>%
          gsub('\r\n *', '', .)
        if (length(karyotype)!=0) {
          sif[i, 'karyotype.atcc'] <- karyotype
          print(paste0(i, ' ATCC: ', karyotype))
        } else {
          next
        }
      }
    }
  }

  ## Karyotype from DSMZ
  dsmz.id <- sif[i, 'DSMZ.ID']
  if(!is.na(dsmz.id)) {
    url <- paste0('https://www.dsmz.de/collection/catalogue/details/culture/', dsmz.id)

    html <- try(read_html(url), silent=TRUE)
    if (class(html)=='try-error') {
      next
    } else {
      class <- html %>%
        html_nodes('[class="label"]') %>%
        html_text()
      value <- html %>%
        html_nodes('[class="value"]') %>%
        html_text()
      if (length(class)==length(value)) {
        df <- data.frame(class=class, value=value)
        karyotype <- df %>%
          filter(class=='Cytogenetics: ') %>%
          pull(value)
        if (length(karyotype)!=0) {
          sif[i, 'karyotype.dsmz'] <- karyotype
          print(paste0(i, ' DSMZ: ', karyotype))
        } else {
          next
        }
      }
    }
  }

  if (is.na(sif[i, 'karyotype.atcc']) & is.na(sif[i, 'karyotype.dsmz'])) {
    sif[i, 'karyotype.info'] <- 'N'
  } else {
    sif[i, 'karyotype.info'] <- 'Y'
  }
}

data.list.file <- here('05_CCLE_data_preparation/data', 'sample.tsv') # Downloaded from Terra (https://app.terra.bio/#workspaces/fccredits-silver-tan-7621/CCLE_v2/data)
data.list <- read.delim(data.list.file)

wes.bam.list <- data.list %>%
  select(entity.sample_id, stripped_cell_line_name, hg38_wes_bam, hg38_wes_bai)
  filter(hg38_wes_bam!='')

ids.with.bam <- wes.bam.list %>%
  filter(hg38_wes_bam!='') %>%
  pull(entity.sample_id)

sif <- sif %>%
  mutate(bam.exists=case_when(DepMap_ID %in% ids.with.bam ~ 'Y', TRUE ~ 'N'))

write.csv(sif, file=here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na='', row.names=FALSE)
