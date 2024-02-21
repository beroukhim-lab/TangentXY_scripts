library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

load(file=here('04_TCGA_Tangent_vs_TangentXY/output/01_sample_selection', 'samples.list.RData'))

hg19 <- rCGH::hg19 %>%
  mutate(chrom=case_when(chrom==23 ~ 'X', chrom==24 ~ 'Y', TRUE ~ as.character(chrom)))

probes <- doc.t %>%
  as.data.frame() %>%
  mutate(locus=rownames(.)) %>%
  select(locus) %>%
  separate(col=locus, into=c('chr', 'pos'), sep=':', remove=FALSE) %>%
  separate(col=pos, into=c('start', 'end'), sep='-') %>%
  mutate(start=as.numeric(start), end=as.numeric(end)) %>%
  mutate(probe=1:n()) %>%
  mutate(chrom=paste0('chr', chr)) %>%
  left_join(hg19 %>% select(chrom, centromerStart, centromerEnd), by=c('chr'='chrom')) %>%
  mutate(arm=case_when(end < centromerStart ~ 'p',
                        start > centromerEnd ~ 'q'))

## Signal and noise
## Make function for calculating signal and noise
calc.signal.noise <- function(list) {
  data <- list %>%
    as.data.frame() %>%
    setNames('signal') %>%
    bind_cols(probes)

  signal.auto <- data %>%
    filter(!chr %in% c('X', 'Y')) %>%
    group_by(chr, arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sds(arm.median)) %>%
    pull(signal)

  noise.auto <- data %>% 
    filter(!chr %in% c('X', 'Y')) %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  result.df <- data.frame(signal.auto=signal.auto,
                          noise.auto=noise.auto)

  return(result.df)
}

n.num <- c(10, 50, 100, 200, 500)
for (i in seq_along(n.num)) {
  n.num.i <- n.num[i]
  print(paste(i, n.num.i))

  ## Prototype Tangent
  Tn.pt <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/10_PrototypeTangent', paste0('Tn_PrototypeTangent_', n.num.i, '.rds'))) %>%
    as.data.frame() %>%
    rownames_to_column('locus') %>%
    right_join(probes %>% select(locus), by='locus') %>%
    column_to_rownames('locus')

  Tn.pt.sn.list.i <- parallel::mclapply(Tn.pt, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  Tn.pt.sn.df.i <- Tn.pt.sn.list.i %>%
    bind_rows() %>%
    mutate(SampleID=names(Tn.pt.sn.list.i)) %>%
    mutate(method='Tangent')

  ## TangentXY
  Tn.xy <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/04_TangentXY/TangentXY', paste0('Tn_TangentXY_', n.num.i, '.rds'))) %>%
    as.data.frame() %>%
    rownames_to_column('locus') %>%
    right_join(probes %>% select(locus), by='locus') %>%
    column_to_rownames('locus')

  Tn.xy.sn.list.i <- parallel::mclapply(Tn.xy, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  Tn.xy.sn.df.i <- Tn.xy.sn.list.i %>%
    bind_rows() %>%
    mutate(SampleID=names(Tn.xy.sn.list.i)) %>%
    mutate(method='TangentXY')

  sn.df.i <- Tn.pt.sn.df.i %>%
    bind_rows(Tn.xy.sn.df.i) %>%
    left_join(sif, by='SampleID') %>%
    mutate(n.num=n.num.i)

  if (i==1) {
    sn.df <- sn.df.i
  } else {
    sn.df <- sn.df %>% bind_rows(sn.df.i)
  }
}
saveRDS(sn.df, file=here('04_TCGA_Tangent_vs_TangentXY/output/11_PrototypeTangent_vs_TangentXY', 'sn.df.rds'), compress=FALSE)


signal.noise <- sn.df %>%
  mutate(n.num=factor(.$n.num, levels=.$n.num %>% unique())) %>%
  mutate(sn.auto=signal.auto/noise.auto) %>%
  mutate(n.num2=case_when(n.num=='10' ~ '5F + 5M',
                          n.num=='50' ~ '25F + 25M',
                          n.num=='100' ~ '50F + 50M',
                          n.num=='200' ~ '100F + 100M',
                          n.num=='500' ~ '250F + 250M')) %>%
  mutate(n.num2=factor(.$n.num2, levels=c('Pre-normalization', '5F + 5M', '25F + 25M', '50F + 50M', '100F + 100M', '250F + 250M')))


## Noise
noise.auto.wider <- signal.noise %>%
  select(SampleID, method, n.num, n.num2, noise.auto) %>%
  pivot_wider(names_from=method, values_from=noise.auto)

g <- ggplot(noise.auto.wider, aes(x=Tangent, y=TangentXY)) +
  geom_abline(intercept=0, slope=1, linetype='dashed', col='gray') +
  geom_point() +
  lemon::facet_rep_wrap(~ n.num2, nrow=1, repeat.tick.labels=TRUE) +
  coord_fixed() +
  labs(title='Noise (Autosomes)') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/11_PrototypeTangent_vs_TangentXY', 'FigS3a.png'), dpi=100, width=14, height=4)
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/11_PrototypeTangent_vs_TangentXY', 'FigS3a.pdf'), width=14, height=4, useDingbats=TRUE)

## SN
sn.auto.wider <- signal.noise %>%
  select(SampleID, method, n.num, n.num2, sn.auto) %>%
  pivot_wider(names_from=method, values_from=sn.auto)

g <- ggplot(sn.auto.wider, aes(x=Tangent, y=TangentXY)) +
  geom_abline(intercept=0, slope=1, linetype='dashed', col='gray') +
  geom_point() +
  lemon::facet_rep_wrap(~ n.num2, nrow=1, repeat.tick.labels=TRUE) +
  coord_fixed() +
  labs(title='Signal-to-noise ratio (Autosomes)') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/11_PrototypeTangent_vs_TangentXY', 'FigS3b.png'), dpi=100, width=14, height=4)
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/11_PrototypeTangent_vs_TangentXY', 'FigS3b.pdf'), width=14, height=4, useDingbats=TRUE)




