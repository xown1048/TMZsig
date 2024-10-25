library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)
library(ggh4x)

options(bitmapType='cairo')

##########################################################################################################################################
#####	INDEL 83 Matrix Bar Count Plot for MSH2
rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/ID/sigProfilerExtractor.ID83.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure2'
if(!dir.exists(output.dir)){dir.create(output.dir)}

df = melt(i, id.vars = 'MutationType')
colnames(df) = c('MutationType', 'Sample', 'Count')

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
df2 = df2[!grepl('-p0', df2$Group),]

df2 = df2[df2$Tissue=='TK6',]
df2 = df2[df2$Drug=='nodrug',]
df3 = df2[grepl('MSH2', df2$Group),]

##########################################################################################################################################
#####	Group level
g.df = df3[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
g.df2 = g.df[!grepl('-p', g.df$Group),]

gene.level = c('WT',  c('p53'), c('ALKBH2', 'ALKBH3', 'MGMT'), c('MPG', 'POLB', 'XRCC1'), c('ERCC1', 'XPA'), 
	c('XPA_XRCC1'), c('POLQ'), c('LIG4'), c('RAD54L_54B'), c('ATAD5'), c('FANCC', 'FANCD2', 'FANCM', 'MUS81'),
	c('POLH', 'POLK', 'POLI', 'POLHKI', 'RAD18', 'REV1', 'REV7'), c('EXO1', 'MSH2'),
	c('MSH2_p53', 'MSH2_MPG', 'MSH2_XRCC1', 'ATAD5_MSH2', 'MSH2_ALKBH2', 'MSH2_ALKBH3', 'FANCD2_MSH2', 'MSH2_RAD18', 'MSH2_REV1', 'MSH2_REV7'))

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

##########################################################################################################################################
#####	Draw plot
temp.df = df3[,c('Sample', 'MutationType', 'Gene', 'Count')]
temp.df = temp.df[temp.df$MutationType %in% c('1:Del:T:5', '1:Ins:T:5', '2:Del:R:5'),]

test.df = compare_means(Count~Gene, method='t.test', group='MutationType', data=temp.df)
test.df2 = as.data.frame(subset(test.df, group1=='MSH2'))
test.df3 = as.data.frame(subset(test.df, group2=='MSH2'))
test.df3 = test.df3[,c(1,2,4,3,5:9)]
colnames(test.df3) = colnames(test.df2)
test.df4 = rbind(test.df2, test.df3)

test.df5 = test.df4[test.df4$p<=0.15,]
print(test.df5)

mean.df = aggregate(Count ~ MutationType+Group, df3, mean)
ci.high.df = aggregate(Count ~ MutationType+Group, df3, function(x){ci(x, ci=0.95, method='ETI')$CI_high})
ci.low.df = aggregate(Count ~ MutationType+Group, df3, function(x){ci(x, ci=0.95, method='ETI')$CI_low})

colnames(mean.df)[3] = 'Mean'
colnames(ci.high.df)[3] = 'CI_high'
colnames(ci.low.df)[3] = 'CI_low'

head(df3)
head(mean.df)
head(ci.high.df)
head(ci.low.df)

merge.df = merge(mean.df, ci.high.df, by=c('MutationType', 'Group'))
merge.df = merge(merge.df, ci.low.df, by=c('MutationType', 'Group'))

merge.df2 = merge.df[merge.df$MutationType %in% c('1:Del:T:5', '1:Ins:T:5', '2:Del:R:5'),]
merge.df2$Gene = gsub('\\-\\S+', '', merge.df2$Group)
merge.df2$Group = factor(merge.df2$Group, levels=group.level)
merge.df2$Gene = factor(merge.df2$Gene, levels=gene.level)

merge.df2$MutationType = factor(merge.df2$MutationType, levels=c('1:Del:T:5', '1:Ins:T:5', '2:Del:R:5'))

facet.name = c(gene.level, c('1:Del:T:5', '1:Ins:T:5', '2:Del:R:5'))
facet.label = c(gsub('_', '\n', gene.level), c('1:Del:T:5', '1:Ins:T:5', '2:Del:R:5'))
names(facet.label) = facet.name

colors = c('#FD8002', '#35A02E', '#FCC9B4')

pdf(file=file.path(output.dir, 'figure2.e.nodrug.msh2.indel83.pdf'), width=50, height=50)
ggplot(merge.df2, aes(x=Group, y=Mean)) +
	geom_bar(aes(fill=MutationType), stat='identity') +
	geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.2, position=position_dodge(.9)) +
	facet_grid(MutationType~Gene, scales='free', space='free_x', labeller=as_labeller(facet.label)) + 
	scale_fill_manual(values=colors) +
	ylab(NULL) + xlab(NULL) +
	theme_light(base_size=80, base_family='sans') +
	theme(axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.text.x=element_text(colour='black', face='bold.italic', size=60),
		strip.text.y=element_text(colour='black', size=100),
		strip.background=element_rect(colour='white', fill='white')) +
	theme(legend.position='none')
dev.off()
