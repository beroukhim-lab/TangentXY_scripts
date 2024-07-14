library(tidyverse)
library(here)

sif <- read.csv(here('05_CCLE_data_preparation/output/02_sample_annotation', 'SampleInfo.csv'), na.strings='') %>%
  mutate(ModelID=gsub('-', '.', ModelID))

summary.data.file <- here('06_CCLE_TangentXY/data', 'TangentXY_accuracy.xlsx')
summary.data <- readxl::read_xlsx(summary.data.file, skip=2) %>%
  rename(correct_Tangent_ChrX=Correct...5) %>%
  rename(wrong_Tangent_ChrX=Wrong...6) %>%
  rename(unclear_Tangent_ChrX=Unclear...7) %>%
  rename(correct_Tangent_ChrY=Correct...8) %>%
  rename(wrong_Tangent_ChrY=Wrong...9) %>%
  rename(unclear_Tangent_ChrY=Unclear...10) %>%
  rename(correct_TangentXY_ChrX=Correct...11) %>%
  rename(wrong_TangentXY_ChrX=Wrong...12) %>%
  rename(unclear_TangentXY_ChrX=Unclear...13) %>%
  rename(correct_TangentXY_ChrY=Correct...14) %>%
  rename(wrong_TangentXY_ChrY=Wrong...15) %>%
  rename(unclear_TangentXY_ChrY=Unclear...16)

summary.longer <- summary.data %>%
  filter(karyo.info=='Y') %>%
  select(-karyo.info) %>%
  pivot_longer(names_to='class', values_to='value', cols=-c('ModelID', 'CellLineName', 'Sex', 'OncotreePrimaryDisease', 'karyo.info_ChrX', 'karyo.info_ChrY',)) %>%
  filter(!is.na(value)) %>%
  separate(col=class, into=c('result', 'method', 'chr'), sep='_') %>%
  mutate(karyo.info=case_when(chr=='ChrX' ~ karyo.info_ChrX, chr=='ChrY' ~ karyo.info_ChrY))

correct.rate <- summary.longer %>%
  filter(result!='unclear') %>%
  group_by(method) %>%
  summarize(total=n(), correct=sum(result=='correct'), wrong=sum(result=='wrong')) %>%
  mutate(correct.rate=correct/total)

correct.rate.XY.FM <- summary.longer %>%
  filter(result!='unclear') %>%
  group_by(Sex, method, chr) %>%
  summarize(total=n(), correct=sum(result=='correct'), wrong=sum(result=='wrong')) %>%
  mutate(correct.rate=correct/total)

g <- ggplot(correct.rate.XY.FM, aes(x=Sex, y=correct.rate)) +
  geom_bar(stat='identity', aes(fill=method), position='dodge') +
  scale_fill_manual(values=c('Tangent'='#AF38EB', 'TangentXY'='#7ECC49')) +
  lemon::facet_rep_grid(~chr, scales='free_x', space='free') +
  labs(y='Correct rate', fill='Method') +
  ylim(0, 1) +
  theme_classic(base_size=20) +
  theme(panel.grid=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/09_TangentXY_performance', 'correct_rate_XY_FM.png'), bg='white', dpi=100, width=8, height=4)

g <- ggplot(correct.rate.XY.FM, aes(x=chr, y=correct.rate)) +
  geom_bar(stat='identity', aes(fill=method), position='dodge') +
  scale_fill_manual(values=c('Tangent'='#AF38EB', 'TangentXY'='#7ECC49')) +
  lemon::facet_rep_grid(~chr, scales='free_x', space='free') +
  labs(y='Correct rate', fill='Method') +
  ylim(0, 1) +
  theme_classic(base_size=20) +
  theme(panel.grid=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/09_TangentXY_performance', 'correct_rate_XY.png'), bg='white', dpi=100, width=8, height=4)

correct.rate.detail <- summary.longer %>%
  filter(result!='unclear') %>%
  group_by(Sex, method, chr, karyo.info) %>%
  summarize(correct=sum(result=='correct'), wrong=sum(result=='wrong')) %>%
  pivot_longer(names_to='result', values_to='count', cols=c('correct', 'wrong')) %>%
  mutate(total=sum(count)) %>%
  mutate(rate=count/total) %>%
  ungroup() %>%
  mutate(result=factor(.$result, levels=c('wrong', 'correct')))

g <- ggplot(correct.rate.detail, aes(x=method, y=rate)) +
  geom_bar(stat='identity', aes(fill=result, color=karyo.info), linewidth=1) +
  ggh4x::facet_nested(. ~ Sex + chr) +
  ylim(0, 1) +
  theme_classic(base_size=20) +
  theme(panel.grid=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/09_TangentXY_performance', 'correct_rate_XY_FM_detail.png'), bg='white', dpi=100, width=8, height=6)

g <- ggplot(correct.rate.detail, aes(x=method, y=count)) +
  geom_bar(stat='identity', aes(fill=result, color=karyo.info), position='stack', linewidth=1) +
  ggh4x::facet_nested(. ~ Sex + chr) +
  theme_classic(base_size=20) +
  theme(panel.grid=element_blank()) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(axis.ticks.x=element_blank()) +
  theme(strip.background=element_blank()) +
  theme(axis.line.x=element_line(linewidth=0.5)) +
  theme(axis.line.y=element_line(linewidth=0.5))
ggsave(g, file=here('06_CCLE_TangentXY/output/09_TangentXY_performance', 'correct_number_XY_FM_detail.png'), bg='white', dpi=100, width=8, height=6)
