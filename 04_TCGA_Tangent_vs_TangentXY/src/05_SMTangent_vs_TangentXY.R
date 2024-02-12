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

  signal.x <- data %>%
    filter(chr=='X') %>%
    group_by(arm) %>%
    summarize(arm.median=median(signal)) %>%
    ungroup() %>%
    summarize(signal=sds(arm.median)) %>%
    pull(signal)

  noise.x <- data %>% 
    filter(chr=='X') %>%
    pull(signal) %>%
    diff() %>%
    abs() %>%
    median()

  result.df <- data.frame(signal.x=signal.x,
                          noise.x=noise.x)

  return(result.df)
}

n.num <- c(10, 50, 100, 200, 500)
for (i in seq_along(n.num)) {
  n.num.i <- n.num[i]
  print(paste(i, n.num.i))

  ## Sex matched Tangent - Female
  Tn.f <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/02_SexMatchedTangent_female', paste0('Tn_SMTangent_female_', n.num.i/2, '.rds'))) %>%
    as.data.frame() %>%
    rownames_to_column('locus') %>%
    right_join(probes %>% select(locus), by='locus') %>%
    column_to_rownames('locus')

  Tn.f.sn.list.i <- parallel::mclapply(Tn.f, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  Tn.f.sn.df.i <- Tn.f.sn.list.i %>%
    bind_rows() %>%
    mutate(SampleID=names(Tn.f.sn.list.i)) %>%
    mutate(method='Sex-matched Tangent')

  ## Sex matched Tangent - Male
  Tn.m <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/03_SexMatchedTangent_male', paste0('Tn_SMTangent_male_', n.num.i/2, '.rds'))) %>%
    as.data.frame() %>%
    rownames_to_column('locus') %>%
    right_join(probes %>% select(locus), by='locus') %>%
    column_to_rownames('locus')

  Tn.m.sn.list.i <- parallel::mclapply(Tn.m, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  Tn.m.sn.df.i <- Tn.m.sn.list.i %>%
    bind_rows() %>%
    mutate(SampleID=names(Tn.m.sn.list.i)) %>%
    mutate(method='Sex-matched Tangent')

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

  sn.df.i <- Tn.f.sn.df.i %>%
    bind_rows(Tn.m.sn.df.i) %>%
    bind_rows(Tn.xy.sn.df.i) %>%
    left_join(sif, by='SampleID') %>%
    mutate(n.num=n.num.i)

  if (i==1) {
    sn.df <- sn.df.i
  } else {
    sn.df <- sn.df %>% bind_rows(sn.df.i)
  }
}
saveRDS(sn.df, file=here('04_TCGA_Tangent_vs_TangentXY/output/05_SMTangent_vs_TangentXY', 'sn.df.rds'), compress=FALSE)

T.sn.list <- parallel::mclapply(doc.t %>% select(c(tumors.female, tumors.male)), calc.signal.noise, mc.cores=parallel::detectCores()-2)
T.sn.df <- T.sn.list %>%
  bind_rows() %>%
  mutate(SampleID=names(T.sn.list)) %>%
  mutate(method='-') %>%
  mutate(n.num='Pre-normalization') %>%
  left_join(sif, by='SampleID')

## Visualization of signal and noise with violin plot
signal.noise <- T.sn.df %>%
  bind_rows(sn.df %>% mutate(n.num=as.character(n.num))) %>%
  mutate(n.num=factor(.$n.num, levels=.$n.num %>% unique())) %>%
  mutate(n.num2=case_when(n.num=='Pre-normalization' ~ 'Pre-normalization',
                          n.num=='10' ~ '5F + 5M',
                          n.num=='50' ~ '25F + 25M',
                          n.num=='100' ~ '50F + 50M',
                          n.num=='200' ~ '100F + 100M',
                          n.num=='500' ~ '250F + 250M')) %>%
  mutate(n.num2=factor(.$n.num2, levels=c('Pre-normalization', '5F + 5M', '25F + 25M', '50F + 50M', '100F + 100M', '250F + 250M'))) %>%
  mutate(method2=case_when(method=='Sex-matched Tangent' & Gender=='Female' ~ 'Sex-matched Tangent_Female',
                            method=='Sex-matched Tangent' & Gender=='Male' ~ 'Sex-matched Tangent_Male',
                            TRUE ~ method))
saveRDS(signal.noise, file=here('04_TCGA_Tangent_vs_TangentXY/output/05_SMTangent_vs_TangentXY', 'signal.noise.rds'), compress=FALSE)


## Noise (chrX, female & male merged)
stat.test.noise.chrX.FMseparated <- signal.noise %>%
  filter(method!='-') %>%
  mutate(n.num=as.character(n.num)) %>%
  filter(n.num %in% c('10', '50', '100', '200', '500')) %>%
  mutate(n.num=factor(.$n.num, levels=c('Pre-normalization', '10', '50', '100', '200', '500'))) %>%
  group_by(Gender, n.num) %>%
  rstatix::wilcox_test(noise.x ~ method, paired=TRUE) %>%
  rstatix::adjust_pvalue(method = "bonferroni") %>%
  rstatix::add_significance()

g <- ggplot(signal.noise, aes(x=n.num2, y=noise.x)) +
  geom_boxplot(aes(group=interaction(Gender, method, n.num2), fill=method2)) +
  scale_fill_manual(values=list('-'='gray', 'Sex-matched Tangent_Female'='red', 'Sex-matched Tangent_Male'='royalblue', 'TangentXY'='green')) +
  lemon::facet_rep_wrap(~Gender, nrow=1, repeat.tick.labels=TRUE) +
  ggpubr::stat_pvalue_manual(stat.test.noise.chrX.FMseparated %>% rstatix::add_xy_position(x='n.num'), label='p.adj.signif', col='black', size=5) +
  ylim(0, NA) +
  labs(title='Noise (ChrX)', x='# normal samples in reference plane', y='Noise', fill='Method') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/05_SMTangent_vs_TangentXY', 'FigS2a.png'), dpi=100, width=14, height=6)
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/05_SMTangent_vs_TangentXY', 'FigS2a.pdf'), width=14, height=6)
