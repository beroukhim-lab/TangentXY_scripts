library(tidyverse)
library(here)

annot.file <- here('11_chrX_allele_specific_expression/data', 'Ensembl_TableBrowser_Exon_hg38.bed')
annot <- read.delim(annot.file, header=FALSE) %>%
  rename(chr=V1, start=V2, stop=V3) %>%
  mutate(CHROM=sub('^chr', '', chr)) %>%
  mutate(CHROM=case_when(CHROM=='M' ~ 'MT', TRUE ~ CHROM)) %>%
  filter(CHROM %in% c(1:22, 'X', 'Y', 'MT'))

propagate.region <- function(df) {
  chr <- df$CHROM
  new.df <- data.frame(CHROM=chr, POS=seq(df$start, df$stop, by=1))
  return(new.df)
}

annot.split <- split(annot, seq(nrow(annot)))
exon.positions.list <- parallel::mclapply(annot.split, propagate.region, mc.cores=parallel::detectCores() - 2)
exon.positions.df <- exon.positions.list %>%
  bind_rows() %>%
  distinct()
saveRDS(exon.positions.df, file=here('11_chrX_allele_specific_expression/output/07_CCLE_annotate_SNPs', 'exon.positions.df.rds'), compress=FALSE)
