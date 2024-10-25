library(data.table)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(janitor)
library(ggbreak)
library(scales)
library(ggpubr)
library(bayestestR)

options(bitmapType='cairo')
options(scipen=10000)

##########################################################################################################################################
#####	SBS Count CI Plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i$MutationType = 'Mutation'
i = aggregate(.~MutationType, i, sum)

i = as.data.frame(i)

df = melt(i, id.vars = 'MutationType')
colnames(df) = c('MutationType', 'Sample', 'Count')

meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')
meta2 = meta
meta2$Sample = gsub('\\.\\S+', '', meta2$Sample)

df = merge(df, meta2, by='Sample')

df$Gene = gsub('\\-\\S+', '', df$Group)
df$Drug = gsub('\\S+\\-', '', df$Group)
df$Dose = sapply(strsplit(as.character(unlist(df$Group)), '-'), '[[', 3)
df$Dose = gsub('nodrug', 0, df$Dose)
df$Tissue = sapply(strsplit(as.character(unlist(df$Group)), '-'), '[[', 2)

df2 = df
df2$Group = gsub('-cisplatin', '-Cisplatin', df2$Group)
df2$Group = gsub('\\-c\\S+', '', df2$Group)

df2 = df2[df2$Tissue=='TK6',]
df2 = df2[df2$Drug %in% c('nodrug', 'TMZ'),]
#df3 = df2[grepl('^MSH2-', df2$Group) | grepl('WT', df2$Group),]
df3 = df2[df2$Sample!='D6_270k',]

############################################################################################################################################
#####	Group level
g.df = df3[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
g.df2 = g.df[!grepl('-p', g.df$Group),]

gene.level = c('WT', c('MGMT', 'ALKBH2', 'ALKBH3'), c('ERCC1', 'XPA'), 'XPA_XRCC1', c('XRCC1', 'POLB', 'MPG'), c('RAD18', 'REV1', 'REV3', 'REV7', 'POLH', 'POLK', 'POLI'), c('FANCD2', 'FANCC', 'FANCM'), c('EXO1'), 'p53', 'ATAD5', c('ESCO1', 'WAPL'))
gene.level = c(gene.level, 'LIG4', 'MUS81', 'POLQ', 'POLHKI', 'RAD54L_54B', 'MSH2', 'ATAD5_MSH2', 'FANCD2_MSH2', 'MSH2_p53', 'MSH2_XRCC1', 'MSH2_ALKBH3', 'MSH2_RAD18', 'MSH2_ALKBH2', 'MSH2_MPG', 'MSH2_REV1', 'MSH2_REV7')

tissue.level = c('TK6', 'HAP1')

gen.g.df = gen.g.df[order(factor(gen.g.df$Tissue, levels=tissue.level), factor(gen.g.df$Gene, levels=gene.level)),]
gen.group.level = gen.g.df$Group

drug.level = c('nodrug', 'TMZ', 'cisplatin')

g.df2$Prefix = ifelse(g.df2$Dose == 0, 0, 
				ifelse(grepl('nanoM', g.df2$Dose), 'nano',
				ifelse(grepl('microM', g.df2$Dose), 'micro', 'error')))

g.df2$Dose2 = ifelse(g.df2$Dose == 0, 0, 
				ifelse(grepl('nanoM', g.df2$Dose), gsub('nanoM','',g.df2$Dose),
				ifelse(grepl('microM', g.df2$Dose), gsub('microM','',g.df2$Dose), 'error')))

g.df2$Dose2 = as.numeric(g.df2$Dose2)

g.df2$Dose2 = ifelse(g.df2$Prefix == 0, 0, 
				ifelse(g.df2$Prefix == 'nano', g.df2$Dose2/1000,
				ifelse(g.df2$Prefix == 'micro', g.df2$Dose2, 'error')))

g.df2$Dose2 = as.numeric(g.df2$Dose2)

g.df2 = g.df2[order(factor(g.df2$Gene, levels=gene.level), factor(g.df2$Tissue, levels=tissue.level), factor(g.df2$Drug, levels=drug.level), g.df2$Dose2),]
group.level = g.df2$Group

group.level = c(group.level, gen.group.level)

group.level = unique(group.level)

############################################################################################################################################
#####	Draw plot
df4 = df3
df4$Dose = gsub('microM', '', df4$Dose)
df4$Dose = as.numeric(df4$Dose)

##### Pearson
r = round(cor(df4$Dose, df4$Count), 2)
p = cor.test(df4$Dose, df4$Count)$p.value

pdf(file=file.path(output.dir, '03.review2.10_sbs.tmz.abs.dot.corr.dose.plot.pdf'), width=5, height=5)
ggplot(df4, aes(x=Dose, y=Count)) + 
	geom_point(size=1) +
	geom_smooth(method='lm', col='black') +
	annotate('text', x=50, y=111000, label=paste0('R = ', r), hjust=0) +
	annotate('text', x=50, y=105000, label='P-value < 0.00000000000000022', hjust=0) +
#	annotate('text', x=50, y=105000, label=paste0('P-value = ', round(p, 3)), hjust=0) +
	theme_light(base_size=10, base_family='sans') +
	ylab('Mutation Count') +
	xlab('Temozolomide Dose')
dev.off()

#2.2e-16

