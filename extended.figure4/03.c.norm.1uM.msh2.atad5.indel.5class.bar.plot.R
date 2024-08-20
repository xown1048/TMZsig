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

rm(list=ls())

##########################################################################################################################################
#####	ID bar absolute number for all samples
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/ID/sigProfilerExtractor.ID83.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/extended.figure4'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i = as.data.frame(i)
i$Len = substr(i$MutationType, 9, 9)
i$Len = with(i, ifelse(substr(MutationType, 3, 5) == 'Del', as.numeric(Len)+1, Len))
i$MutationType = substr(i$MutationType, 1, 8)
i$MutationType = paste0(i$MutationType, i$Len)
i$Len = NULL

i$MutationType = c(rep('1bp_Deletion',12), rep('1bp_Insertion',12), rep('Repeat_Deletion',24), rep('Repeat_Insertion',24), rep('Microhomology',11))
i = aggregate(.~MutationType, i, sum)

meta$Sample = gsub('\\.\\S+', '', meta$Sample)

df = melt(i, id.vars = 'MutationType')
colnames(df) = c('MutationType', 'Sample', 'Count')

df = merge(df, meta, by='Sample')

df$Gene = gsub('\\-\\S+', '', df$Group)
df$Drug = gsub('\\S+\\-', '', df$Group)
df$Dose = sapply(strsplit(as.character(unlist(df$Group)), '-'), '[[', 3)
df$Dose = gsub('nodrug', 0, df$Dose)
df$Tissue = sapply(strsplit(as.character(unlist(df$Group)), '-'), '[[', 2)

df2 = df
df2$Group = gsub('-cisplatin', '-Cisplatin', df2$Group)
df2$Group = gsub('\\-c\\S+', '', df2$Group)

df2 = df2[df2$Tissue == 'TK6',]
df2 = df2[df2$Drug %in% c('nodrug', 'TMZ'),]
df3 = df2[grepl('MSH2', df2$Group),]
df3 = df3[df3$Gene %in% c('MSH2', 'ATAD5_MSH2'),]

############################################################################################################################################
#####	Group level
g.df = df3[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
g.df2 = g.df[!grepl('-p', g.df$Group),]

gene.level = c('WT', c('MGMT', 'ALKBH2', 'ALKBH3'), c('ERCC1', 'XPA'), 'XPA_XRCC1', c('XRCC1', 'POLB', 'MPG'),
	c('RAD18', 'REV1', 'REV3', 'REV7', 'POLH', 'POLK', 'POLI'), c('FANCD2', 'FANCC', 'FANCM'), c('EXO1'), 'p53', 'ATAD5', c('ESCO1', 'WAPL'))
gene.level = c(gene.level, 'LIG4', 'MUS81', 'POLQ', 'POLHKI', 'RAD54L_54B', 'MSH2', 'MSH2_ALKBH3', 'MSH2_ALKBH2',
	'MSH2_XRCC1', 'MSH2_MPG', 'MSH2_RAD18', 'MSH2_REV1', 'MSH2_REV7', 'FANCD2_MSH2', 'MSH2_p53', 'ATAD5_MSH2')

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
df3$Prefix = ifelse(df3$Dose == 0, 0, 
				ifelse(grepl('nanoM', df3$Dose), 'nano',
				ifelse(grepl('microM', df3$Dose), 'micro', 'error')))

df3$Dose2 = ifelse(df3$Dose == 0, 0, 
				ifelse(grepl('nanoM', df3$Dose), gsub('nanoM','',df3$Dose),
				ifelse(grepl('microM', df3$Dose), gsub('microM','',df3$Dose), 'error')))

df3$Dose2 = as.numeric(df3$Dose2)
df3$Dose2 = ifelse(df3$Prefix == 0, 0, 
				ifelse(df3$Prefix == 'nano', df3$Dose2/1000,
				ifelse(df3$Prefix == 'micro', df3$Dose2, 'error')))

df4 = subset(df3, Drug=='nodrug')
df5 = subset(df3, Drug!='nodrug')
df5$Category = paste0(df5$Gene, '-', df5$Drug)
df5$Dose2 = as.numeric(df5$Dose2)
df4$Category = paste0(df4$Gene, '-', df4$Drug)
df4$Count2 = df4$Count

cates = sort(unique(df5$Category))

merge.df = data.frame()
for (cate in cates){
	temp.df = subset(df5, Category==cate)

	nodrug.temp.df = subset(df4, df4$Gene == unique(temp.df$Gene))
	nodrug.temp.df2 = nodrug.temp.df[,c(2,3)]
	nodrug.temp.df3 = aggregate(Count ~ MutationType, nodrug.temp.df2, mean)

	temp.df$Count2 = temp.df$Count - nodrug.temp.df3$Count[match(temp.df$MutationType, nodrug.temp.df3$MutationType)]
	temp.df$Count2 = temp.df$Count2/temp.df$Dose2
	#temp.df$Count2 = temp.df$Count2*500/temp.df$Dose2

	temp.df$Count2 = ifelse(temp.df$Count2 < 0, 0, temp.df$Count2)

	merge.df = rbind(merge.df, temp.df)
}

merge.df$Gene = factor(merge.df$Gene, levels=gene.level)
merge.df$Group = factor(merge.df$Group, levels=group.level)

group.level2 = group.level[!grepl('nodrug', group.level)]
facet.name = group.level2

temp.label = gsub('-TK6-', '-', group.level2)
temp.label = gsub('-TMZ', '', temp.label)
temp.label = gsub('-', '\n', temp.label)
temp.label = gsub('_', '\n', temp.label)
temp.label = gsub('micro', 'u', temp.label)

facet.label = temp.label
names(facet.label) = facet.name

colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#264653')
merge.df$MutationType = factor(merge.df$MutationType, levels=c('1bp_Deletion', '1bp_Insertion', 'Repeat_Deletion', 'Repeat_Insertion', 'Microhomology'))

pdf(file=file.path(output.dir, 'extended.figure4.a.norm.1uM.msh2.indel.5class.barplot.pdf'), width=100, height=60)
ggplot(merge.df, aes(x=Sample, y=Count2)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=200, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(legend.key.size=unit(3, 'cm'), strip.text=element_text(colour='black', size=200, face='bold')) +
	theme(axis.text.x=element_blank(), strip.background=element_rect(colour='white', fill='white'),
		axis.title=element_text(size=250, face='bold'))
dev.off()

