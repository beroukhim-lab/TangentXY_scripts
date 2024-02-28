library(tidyverse)
library(here)

vcf.file <- here('11_chrX_allele_specific_expression/output/01_TCGA_prepare_vcf_for_ASEReadCounter', 'TCGA-FU-A3WB-10A-01D-A22X-09.vcf.gz')
snp <- vcfR::read.vcfR(vcf.file)
snp.fix <- snp@fix %>%
  as.data.frame() %>%
  mutate(POS=as.numeric(POS))

snp.probes <- snp.fix %>%
  select(c('CHROM', 'POS')) %>%
  distinct()

annot.file <- here('11_chrX_allele_specific_expression/data', 'Ensembl_TableBrowser_ExonPlus2bases.bed')
annot <- read.delim(annot.file, skip=1, header=FALSE) %>%
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
  bind_rows()

exon.positions.df.distinct <- exon.positions.df %>%
  distinct() %>%
  mutate(region='exon')

snp.annot.df <- snp.probes %>%
  left_join(exon.positions.df.distinct, by=c('CHROM', 'POS')) %>%
  mutate(region=case_when(is.na(region) ~ 'intron', TRUE ~ region))
saveRDS(snp.annot.df, file=here('11_chrX_allele_specific_expression/output/03_TCGA_annotate_SNPs', 'snp.annot.df.rds'), compress=FALSE)s
