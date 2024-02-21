library(tidyverse)
library(here)

library(survival)
library(ggsurvfit)

sif <- read.delim(file=here('02_TCGA_data_preparation/data', 'sif.txt'))

sample.amp.del <- readRDS(file=here('07_SCNAs_in_chrX_and_chrY/output/01_TCGA_SCNA_classification', 'sample.amp.del.rds'))

cli1 <- readxl::read_xlsx(path=here('09_TCGA_SCNA_vs_prognosis/data', 'TCGA-CDR-SupplementalTableS1.xlsx'), sheet='TCGA-CDR') %>%
  select(-1) %>%
  mutate(TCGA.ID=gsub('-', '.', bcr_patient_barcode)) %>%
  as.data.frame()

cli2 <- readxl::read_xlsx(path=here('09_TCGA_SCNA_vs_prognosis/data', 'TCGA-CDR-SupplementalTableS1.xlsx'), sheet='ExtraEndpoints') %>%
  select(-1) %>%
  mutate(TCGA.ID=gsub('-', '.', bcr_patient_barcode)) %>%
  as.data.frame()

cli <- cli1 %>%
  left_join(cli2 %>% select(TCGA.ID, PFS, PFS.time), by='TCGA.ID')

chrXY.karyo <- sample.amp.del %>%
  filter(chr %in% c('X', 'Y')) %>%
  select(SampleID, chr, karyo.class, TCGA.ID, Gender, project)

cli.karyo <- cli %>%
  select(-c('gender', 'type')) %>%
  inner_join(chrXY.karyo, by='TCGA.ID') %>%
  mutate(if.amp=case_when(karyo.class %in% c('No.Alt', 'Whole.Amp') ~ karyo.class)) %>%
  mutate(if.del=case_when(karyo.class %in% c('No.Alt', 'Whole.Del') ~ karyo.class))

## Pan-cancer
## OS ~ Amp
surv.os.amp <- ggsurvfit::survfit2(Surv(OS.time, OS) ~ if.amp, data=cli.karyo %>% filter(!is.na(if.amp)))
survdiff.os.amp <- survdiff(Surv(OS.time, OS) ~ if.amp, data=cli.karyo %>% filter(!is.na(if.amp)))
g <- ggsurvfit(surv.os.amp) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.os.amp$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Amp'='red')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Amp'='red')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='OS ~ ChrX Amp (Female + Male)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/OS', 'OS_Amp_all.png'), dpi=100, width=10, height=8)

## PFS ~ Amp
surv.pfs.amp <- ggsurvfit::survfit2(Surv(PFS.time, PFS) ~ if.amp, data=cli.karyo %>% filter(!is.na(if.amp)))
survdiff.pfs.amp <- survdiff(Surv(PFS.time, PFS) ~ if.amp, data=cli.karyo %>% filter(!is.na(if.amp)))
g <- ggsurvfit(surv.pfs.amp) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfs.amp$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Amp'='red')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Amp'='red')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='PFS ~ ChrX Amp (Female + Male)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/PFS', 'PFS_Amp_all.png'), dpi=100, width=10, height=8)

## PFI ~ Amp
surv.pfi.amp <- ggsurvfit::survfit2(Surv(PFI.time, PFI) ~ if.amp, data=cli.karyo %>% filter(!is.na(if.amp)))
survdiff.pfi.amp <- survdiff(Surv(PFI.time, PFI) ~ if.amp, data=cli.karyo %>% filter(!is.na(if.amp)))
g <- ggsurvfit(surv.pfi.amp) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfi.amp$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Amp'='red')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Amp'='red')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='PFI ~ ChrX Amp (Female + Male)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/PFI', 'PFI_Amp_all.png'), dpi=100, width=10, height=8)

## OS ~ Del
surv.os.del <- ggsurvfit::survfit2(Surv(OS.time, OS) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)))
survdiff.os.del <- survdiff(Surv(OS.time, OS) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)))
g <- ggsurvfit(surv.os.del) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.os.del$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='OS ~ ChrX Del (Female + Male)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/OS', 'OS_Del_all.png'), dpi=100, width=10, height=8)

## PFS ~ Del
surv.pfs.del <- ggsurvfit::survfit2(Surv(PFS.time, PFS) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)))
survdiff.pfs.del <- survdiff(Surv(PFS.time, PFS) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)))
g <- ggsurvfit(surv.pfs.del) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfs.del$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='PFS ~ ChrX Del (Female + Male)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/PFS', 'PFS_Del_all.png'), dpi=100, width=10, height=8)

## PFI ~ Del
surv.pfi.del <- ggsurvfit::survfit2(Surv(PFI.time, PFI) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)))
survdiff.pfi.del <- survdiff(Surv(PFI.time, PFI) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)))
g <- ggsurvfit(surv.pfi.del) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfi.del$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='PFI ~ ChrX Del (Female + Male)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/PFI', 'PFI_Del_all.png'), dpi=100, width=10, height=8)

## OS ~ Del (Female only)
surv.os.del <- ggsurvfit::survfit2(Surv(OS.time, OS) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)) %>% filter(Gender=='Female'))
survdiff.os.del <- survdiff(Surv(OS.time, OS) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)) %>% filter(Gender=='Female'))
g <- ggsurvfit(surv.os.del) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.os.del$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='OS ~ ChrX Del (Female only)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/OS', 'OS_Del_female.png'), dpi=100, width=10, height=8)

## PFS ~ Del (Female only)
surv.pfs.del <- ggsurvfit::survfit2(Surv(PFS.time, PFS) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)) %>% filter(Gender=='Female'))
survdiff.pfs.del <- survdiff(Surv(PFS.time, PFS) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)) %>% filter(Gender=='Female'))
g <- ggsurvfit(surv.pfs.del) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfs.del$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='PFS ~ ChrX Del (Female only)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/PFS', 'PFS_Del_female.png'), dpi=100, width=10, height=8)

## PFI ~ Del (Female only)
surv.pfi.del <- ggsurvfit::survfit2(Surv(PFI.time, PFI) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)) %>% filter(Gender=='Female'))
survdiff.pfi.del <- survdiff(Surv(PFI.time, PFI) ~ if.del, data=cli.karyo %>% filter(!is.na(if.del)) %>% filter(Gender=='Female'))
g <- ggsurvfit(surv.pfi.del) +
  add_confidence_interval() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfi.del$pvalue, digit=3))), size=10, hjust=1) +
  add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
  scale_color_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_fill_manual(values=c('WT'='black', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='PFI ~ ChrX Del (Female only)', x='Days') +
  theme_classic(base_size=20) +
  theme(legend.position = c(0.8, 0.7))
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/PFI', 'PFI_Del_female.png'), dpi=100, width=10, height=8)



## OS, PFS, and PFI ~ Amp/Del for each tumor type (Females and males combined)
projects <- cli.karyo$project %>% unique()
for (i in 1:length(projects)) {
  project.i <- projects[i]

  for (chr.i in c('X', 'Y')) {
    for (alt.i in c('Amp', 'Del')) {
      cli.karyo.i <- cli.karyo %>%
        filter(project==project.i) %>%
        filter(chr==chr.i) %>%
        filter(if (alt.i=='Amp') !is.na(if.amp) else if (alt.i=='Del') !is.na(if.del))

      if (alt.i=='Amp') {
        cli.karyo.i <- cli.karyo.i %>%
          rename(if.alt='if.amp')
      } else if (alt.i=='Del') {
        cli.karyo.i <- cli.karyo.i %>%
          rename(if.alt='if.del')
      }

      os.data <- cli.karyo.i %>%
        filter(!is.na(OS))
      pfs.data <- cli.karyo.i %>%
        filter(!is.na(PFS))
      pfi.data <- cli.karyo.i %>%
        filter(!is.na(PFI))

      os.df.i <- data.frame(tumor.type=project.i,
                            chr=chr.i,
                            alt.type=alt.i,
                            wt=os.data %>% filter(if.alt=='No.Alt') %>% nrow(),
                            alt=os.data %>% filter(if (alt.i=='Amp') if.alt=='Whole.Amp' else if (alt.i=='Del') if.alt=='Whole.Del') %>% nrow())

      pfs.df.i <- data.frame(tumor.type=project.i,
                            chr=chr.i,
                            alt.type=alt.i,
                            wt=os.data %>% filter(if.alt=='No.Alt') %>% nrow(),
                            alt=os.data %>% filter(if (alt.i=='Amp') if.alt=='Whole.Amp' else if (alt.i=='Del') if.alt=='Whole.Del') %>% nrow())

      pfi.df.i <- data.frame(tumor.type=project.i,
                            chr=chr.i,
                            alt.type=alt.i,
                            wt=os.data %>% filter(if.alt=='No.Alt') %>% nrow(),
                            alt=os.data %>% filter(if (alt.i=='Amp') if.alt=='Whole.Amp' else if (alt.i=='Del') if.alt=='Whole.Del') %>% nrow())

      print(paste(i, project.i, chr.i, alt.i))

      if (os.df.i$wt >=5 & os.df.i$alt >= 5) {
        surv.os.i <- ggsurvfit::survfit2(Surv(OS.time, OS) ~ if.alt, data=os.data)
        survdiff.os.i <- survdiff(Surv(OS.time, OS) ~ if.alt, data=os.data)
        g <- ggsurvfit(surv.os.i) +
          add_confidence_interval() +
          add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.os.i$pvalue, digit=3))), size=10, hjust=1) +
          add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
          scale_color_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
          scale_fill_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
          scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
          labs(title=paste0('OS ~ Chr', chr.i, ' ', alt.i, ' (', project.i, ')'), x='Days') +
          theme_classic(base_size=30) +
          theme(axis.line.x=element_line(linewidth=0.5)) +
          theme(axis.line.y=element_line(linewidth=0.5))
          theme(legend.position = c(0.8, 0.7))
        ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/OS', paste0('OS_', project.i, '_', chr.i, '_', alt.i, '.png')), dpi=100, width=10, height=10)

        os.df.i <- os.df.i %>%
            mutate(pval=survdiff.os.i$pvalue)
      }

      if (pfs.df.i$wt >=5 & pfs.df.i$alt >= 5) {
        surv.pfs.i <- ggsurvfit::survfit2(Surv(PFS.time, PFS) ~ if.alt, data=pfs.data)
        survdiff.pfs.i <- survdiff(Surv(PFS.time, PFS) ~ if.alt, data=pfs.data)
        g <- ggsurvfit(surv.pfs.i) +
          add_confidence_interval() +
          add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfs.i$pvalue, digit=3))), size=10, hjust=1) +
          add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
          scale_color_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
          scale_fill_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
          scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
          labs(title=paste0('PFS ~ Chr', chr.i, ' ', alt.i, ' (', project.i, ')'), x='Days') +
          theme_bw(base_size=30) +
          theme(panel.grid.major=element_blank()) +
          theme(panel.grid.minor=element_blank()) +
          theme(legend.position = c(0.8, 0.7))
        ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/PFS', paste0('PFS_', project.i, '_', chr.i, '_', alt.i, '.png')), dpi=100, width=10, height=10)

        pfs.df.i <- pfs.df.i %>%
            mutate(pval=survdiff.pfs.i$pvalue)
      }

      if (pfi.df.i$wt >=5 & pfi.df.i$alt >= 5) {
        surv.pfi.i <- ggsurvfit::survfit2(Surv(PFI.time, PFI) ~ if.alt, data=pfi.data)
        survdiff.pfi.i <- survdiff(Surv(PFI.time, PFI) ~ if.alt, data=pfi.data)
        g <- ggsurvfit(surv.pfi.i) +
          add_confidence_interval() +
          add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfi.i$pvalue, digit=3))), size=10, hjust=1) +
          add_risktable(risktable_stats='n.risk', size=7, theme=theme_risktable_default(plot.title.size=10, axis.text.y.size=10)) +
          scale_color_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
          scale_fill_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
          scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
          labs(title=paste0('PFI ~ Chr', chr.i, ' ', alt.i, ' (', project.i, ')'), x='Days') +
          theme_bw(base_size=30) +
          theme(panel.grid.major=element_blank()) +
          theme(panel.grid.minor=element_blank()) +
          theme(legend.position = c(0.8, 0.7))
        ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/PFI', paste0('PFI_', project.i, '_', chr.i, '_', alt.i, '.png')), dpi=100, width=10, height=10)

        pfi.df.i <- pfi.df.i %>%
            mutate(pval=survdiff.pfi.i$pvalue)
      }

      if (i==1 & chr.i=='X' & alt.i=='Amp') {
        os.df <- os.df.i
        pfs.df <- pfs.df.i
        pfi.df <- pfi.df.i
      } else {
        os.df <- os.df %>% bind_rows(os.df.i)
        pfs.df <- pfs.df %>% bind_rows(pfs.df.i)
        pfi.df <- pfi.df %>% bind_rows(pfi.df.i)
      }
    }
  }
}

os.df <- os.df %>%
  mutate(adj.pval=p.adjust(pval, method='bonferroni'))
pfs.df <- pfs.df %>%
  mutate(adj.pval=p.adjust(pval, method='bonferroni'))
pfi.df <- pfi.df %>%
  mutate(adj.pval=p.adjust(pval, method='bonferroni'))

saveRDS(os.df, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis', 'os.df.rds'), compress=FALSE)
saveRDS(pfs.df, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis', 'pfs.df.rds'), compress=FALSE)
saveRDS(pfi.df, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis', 'pfi.df.rds'), compress=FALSE)

pfs.p.adj.table <- pfs.df %>%
  unite(col=chr.alt, c('chr', 'alt.type'), sep='_') %>%
  select(tumor.type, chr.alt, adj.pval) %>%
  pivot_wider(names_from=chr.alt, values_from=adj.pval)

## Females and males separately
for (i in 1:length(projects)) {
  project.i <- projects[i]
  print(project.i)

  data.female <- cli.karyo %>%
    filter(project==project.i, Gender=='Female')

  data.male <- cli.karyo %>%
    filter(project==project.i, Gender=='Male')

  if (nrow(data.female)!=0) {
    surv.female <- ggsurvfit::survfit2(Surv(OS.time, OS) ~ karyo.class, data=data.female)
    survdiff.female <- survdiff(Surv(OS.time, OS) ~ karyo.class, data=data.female)

    g.female <- ggsurvfit(surv.female) +
      add_confidence_interval() +
      add_risktable() +
      scale_color_manual(values=c('WT'='black', 'Whole.Amp'='red', 'Arm.Amp'='pink', 'Whole.Del'='blue', 'Arm.Del'='cyan', 'Amp.Del'='purple')) +
      scale_fill_manual(values=c('WT'='black', 'Whole.Amp'='red', 'Arm.Amp'='pink', 'Whole.Del'='blue', 'Arm.Del'='cyan', 'Amp.Del'='purple')) +
      labs(title=paste0(project.i, ', Female')) +
      theme_bw(base_size=30)

    ggsave(g.female, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/OS', paste0(project.i, '_female.png')), dpi=100, width=12, height=8)
  }

  if (nrow(data.male)!=0) {
    surv.male <- ggsurvfit::survfit2(Surv(OS.time, OS) ~ karyo.class, data=data.male)

    g.male <- ggsurvfit(surv.male) +
      add_confidence_interval() +
      add_risktable() +
      scale_color_manual(values=c('WT'='black', 'Whole.Amp'='red', 'Arm.Amp'='pink', 'Whole.Del'='blue', 'Arm.Del'='cyan', 'Amp.Del'='purple')) +
      scale_fill_manual(values=c('WT'='black', 'Whole.Amp'='red', 'Arm.Amp'='pink', 'Whole.Del'='blue', 'Arm.Del'='cyan', 'Amp.Del'='purple')) +
      labs(title=paste0(project.i, ', Male')) +
      theme_bw(base_size=30)

    ggsave(g.male, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis/KMplot/OS', paste0(project.i, '_male.png')), dpi=100, width=12, height=8)
  }
}


## Make figures for paper
## PFS ~ ChrX Del in OV
cli.karyo.ov.chrX.del <- cli.karyo %>%
  filter(project=='OV') %>%
  filter(chr=='X') %>%
  filter(!is.na(if.del)) %>%
  filter(!is.na(PFS)) %>%
  rename(if.alt=if.del)

surv.pfs.ov.chrX.del <- ggsurvfit::survfit2(Surv(PFS.time, PFS) ~ if.alt, data=cli.karyo.ov.chrX.del)
survdiff.pfs.ov.chrX.del <- survdiff(Surv(PFS.time, PFS) ~ if.alt, data=cli.karyo.ov.chrX.del)
g <- ggsurvfit(surv.pfs.ov.chrX.del) +
  add_censor_mark() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfs.ov.chrX.del$pvalue, digit=3))), size=5, hjust=1, y=1) +
  scale_color_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
  scale_fill_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='OV ~ ChrX deletion', x='Progression-free survival (Days)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(legend.position = c(0.8, 0.85)) +
  theme(legend.background=element_blank())
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis', 'Fig4a.png'), dpi=100, width=6, height=6)
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis', 'Fig4a.pdf'), width=6, height=6)

## PFS ~ ChrY Del in UVM
cli.karyo.uvm.chrY.del <- cli.karyo %>%
  filter(project=='UVM') %>%
  filter(chr=='Y') %>%
  filter(!is.na(if.del)) %>%
  filter(!is.na(PFS)) %>%
  rename(if.alt=if.del)

surv.pfs.uvm.chrY.del <- ggsurvfit::survfit2(Surv(PFS.time, PFS) ~ if.alt, data=cli.karyo.uvm.chrY.del)
survdiff.pfs.uvm.chrY.del <- survdiff(Surv(PFS.time, PFS) ~ if.alt, data=cli.karyo.uvm.chrY.del)
g <- ggsurvfit(surv.pfs.uvm.chrY.del) +
  add_censor_mark() +
  add_pvalue(location='annotation', caption=paste('p=', as.character(format(survdiff.pfs.uvm.chrY.del$pvalue, digit=3))), size=5, hjust=1, y=1) +
  scale_color_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
  scale_fill_manual(values=c('No.Alt'='black', 'Whole.Amp'='red', 'Whole.Del'='blue')) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
  labs(title='UVM ~ ChrY deletion', x='Progression-free survival (Days)') +
  theme_classic(base_size=20) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5)) +
  theme(legend.position = c(0.8, 0.85)) +
  theme(legend.background=element_blank())
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis', 'Fig4b.png'), dpi=100, width=6, height=6)
ggsave(g, file=here('09_TCGA_SCNA_vs_prognosis/output/01_SCNA_vs_prognosis', 'Fig4b.pdf'), width=6, height=6)

