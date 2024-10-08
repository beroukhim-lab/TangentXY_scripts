library(tidyverse)
library(here)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

doc.t <- readRDS(file=here('02_TCGA_data_preparation/output/02_6_DOC_Preprocessing_removeCommonGermlineCNVs', 'TCGA_WES_hg19_T_QCed_commonCNVremoved.rds'))

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

fm.table <- data.frame(case=1:8, female=c(5, 25, 50, 150, 350, 450, 475, 495), male=c(495, 475, 450, 350, 150, 50, 25, 5))
for (i in 1:nrow(fm.table)) {
  n.num.f <- fm.table[i, 'female']
  n.num.m <- fm.table[i, 'male']
  print(paste0(i, ' F:', n.num.f, ', M:', n.num.m))

  fm <- paste0('F', n.num.f, 'M', n.num.m)

  ## Sex matched Tangent - Female
  Tn.f <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/06_SexMatchedTangent_female_FMratio', paste0('Tn_SMTangent_female_', n.num.f, '.rds'))) %>%
    as.data.frame() %>%
    rownames_to_column('locus') %>%
    right_join(probes %>% select(locus), by='locus') %>%
    column_to_rownames('locus')

  Tn.f.sn.list.i <- parallel::mclapply(Tn.f, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  Tn.f.sn.df.i <- Tn.f.sn.list.i %>%
    bind_rows() %>%
    mutate(SampleID=names(Tn.f.sn.list.i)) %>%
    mutate(method='Sex-matched Tangent_Female')


  ## Sex matched Tangent - Male
  Tn.m <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/07_SexMatchedTangent_male_FMratio', paste0('Tn_SMTangent_male_', n.num.m, '.rds'))) %>%
    as.data.frame() %>%
    rownames_to_column('locus') %>%
    right_join(probes %>% select(locus), by='locus') %>%
    column_to_rownames('locus')

  Tn.m.sn.list.i <- parallel::mclapply(Tn.m, calc.signal.noise, mc.cores=parallel::detectCores()-2)
  Tn.m.sn.df.i <- Tn.m.sn.list.i %>%
    bind_rows() %>%
    mutate(SampleID=names(Tn.m.sn.list.i)) %>%
    mutate(method='Sex-matched Tangent_Male')


  ## TangentXY
  Tn.xy <- readRDS(file=here('04_TCGA_Tangent_vs_TangentXY/output/08_TangentXY_FMratio/TangentXY', paste0('Tn_TangentXY_', fm, '.rds'))) %>%
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
    mutate(fm=fm)

  if (i==1) {
    sn.df <- sn.df.i
  } else {
    sn.df <- sn.df %>% bind_rows(sn.df.i)
  }
}
saveRDS(sn.df, file=here('04_TCGA_Tangent_vs_TangentXY/output/09_SMTangent_vs_TangentXY_FMratio', 'sn.df.rds'), compress=FALSE)

## Visualization of signal and noise with violin plot
signal.noise <- sn.df %>%
  mutate(fm=factor(.$fm, levels=.$fm %>% unique())) %>%
  mutate(sn.x=signal.x/noise.x) %>%
  mutate(fm2=case_when(fm=='F5M495' ~ '5F + 495M',
                        fm=='F25M475' ~ '25F + 475M',
                        fm=='F50M450' ~ '50F + 450M',
                        fm=='F150M350' ~ '150F + 350M',
                        fm=='F350M150' ~ '350F + 150M',
                        fm=='F450M50' ~ '450F + 50M',
                        fm=='F475M25' ~ '475F + 25M',
                        fm=='F495M5' ~ '495F + 5M')) %>%
  mutate(fm2=factor(.$fm2, levels=c('5F + 495M', '25F + 475M', '50F + 450M', '150F + 350M', '350F + 150M', '450F + 50M', '475F + 25M', '495F + 5M')))

## Noise (ChrX)
stat.test.noise.auto.FMmerged.female <- signal.noise %>%
  filter(Gender=='Female') %>%
  filter(method!='-') %>%
  group_by(fm2) %>%
  rstatix::wilcox_test(noise.x ~ method, paired=TRUE) %>%
  mutate(Gender=case_when(group1=='Sex-matched Tangent_Female' ~ 'Female', group1=='Sex-matched Tangent_Male' ~ 'Male')) %>%
  rstatix::add_xy_position(x='fm2')

stat.test.noise.auto.FMmerged.male <- signal.noise %>%
  filter(Gender=='Male') %>%
  filter(method!='-') %>%
  group_by(fm2) %>%
  rstatix::wilcox_test(noise.x ~ method, paired=TRUE) %>%
  mutate(Gender=case_when(group1=='Sex-matched Tangent_Female' ~ 'Female', group1=='Sex-matched Tangent_Male' ~ 'Male')) %>%
  rstatix::add_xy_position(x='fm2')

stat.test.noise.auto.FMmerged <- stat.test.noise.auto.FMmerged.female %>%
  bind_rows(stat.test.noise.auto.FMmerged.male) %>%
  rstatix::adjust_pvalue(method = "bonferroni") %>%
  rstatix::add_significance() %>%
  mutate(y.position=case_when(Gender=='Female' ~ y.position+0.03, TRUE ~ y.position))

g <- ggplot(signal.noise, aes(x=fm2, y=noise.x)) +
  geom_boxplot(aes(group=interaction(fm2, method), fill=method)) +
  ggpubr::stat_pvalue_manual(stat.test.noise.auto.FMmerged, label='p.adj.signif', col='black', size=5) +
  scale_fill_manual(values=list('Sex-matched Tangent_Female'='red', 'Sex-matched Tangent_Male'='royalblue', 'TangentXY'='green')) +
  ylim(0, NA) +
  lemon::facet_rep_wrap(~Gender, nrow=1, repeat.tick.labels=TRUE) +
  labs(title='Noise (ChrX)', x='# normal samples in reference plane', y='Noise', fill='Method') +
  theme_classic(base_size=20) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/09_SMTangent_vs_TangentXY_FMratio', 'FigS2b.png'), dpi=100, width=14, height=6)
ggsave(g, file=here('04_TCGA_Tangent_vs_TangentXY/output/09_SMTangent_vs_TangentXY_FMratio', 'FigS2b.pdf'), width=14, height=6)
