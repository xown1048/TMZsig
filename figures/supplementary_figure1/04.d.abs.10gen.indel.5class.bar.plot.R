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

rm(list=ls())

######################################################################
##### ID bar absolute number for all samples
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/ID/sigProfilerExtractor.ID83.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/extended.figure1'
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

df2 = df2[df2$Tissue=='TK6',]
gene.level = c('WT', 'XRCC1', 'POLHKI', 'REV1', 'REV7')
df2 = df2[df2$Gene %in% gene.level,]

df3 = df2[grepl('-p10', df2$Group),]

######################################################################
##### Group level
g.df = df3[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
tissue.level = c('TK6', 'HAP1')

gen.g.df = gen.g.df[order(factor(gen.g.df$Tissue, levels=tissue.level), factor(gen.g.df$Gene, levels=gene.level)),]
gen.group.level = gen.g.df$Group

group.level = unique(gen.group.level)

######################################################################
#####	Draw plots
colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#264653')

df3$Group = factor(df3$Group, levels=group.level)
df3$Gene = factor(df3$Gene, levels=gene.level)
df3$MutationType = factor(df3$MutationType, levels=c('1bp_Deletion', '1bp_Insertion', 'Repeat_Deletion', 'Repeat_Insertion', 'Microhomology'))

pdf(file=file.path(output.dir, 'extended.figure1.d.abs.tk6.p10.10gen.indel.5class.barplot.pdf'), width=100, height=65)
ggplot(df3, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Gene, scales='free_x', space='free_x') +
	theme_light(base_size=200, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_text(colour='black', size=200, face='bold')) +
	theme(legend.position='none') +
	theme(axis.text.x=element_blank(), strip.background=element_rect(colour='white', fill='white'),
		axis.title=element_text(size=250, face='bold'), axis.ticks.x=element_blank())
dev.off()

######################################################################
##### Count write table
output.dir2 = paste0(output.dir, '/table')
if(!dir.exists(output.dir2)){dir.create(output.dir2)}

head(df3)

df4 = df3[,c('Count', 'MutationType', 'Gene')]
mean.df = aggregate(Count~MutationType+Gene, df4, mean)
sd.df = aggregate(Count~MutationType+Gene, df4, sd)
colnames(mean.df)[3] = 'Mean'
colnames(sd.df)[3] = 'Sd'

df5 = merge(mean.df, sd.df, by=c('Gene', 'MutationType'))
df5$Sd = round(df5$Sd, 1)
df5$Mean = round(df5$Mean, 1)

df5 = df5[order(factor(df5$Gene, levels=group.level), df5$MutationType),]

sum.df = aggregate(Count~Sample+Gene, df3, sum)
mean.df = aggregate(Count~Gene, sum.df, mean)
sd.df = aggregate(Count~Gene, sum.df, sd)
colnames(mean.df)[2] = 'Mean'
colnames(sd.df)[2] = 'Sd'

df6 = merge(mean.df, sd.df, by='Gene')
df6$Sd = round(df6$Sd, 1)
df6$Mean = round(df6$Mean, 1)

df6 = df6[order(factor(df6$Gene, levels=group.level)),]

write.table(df5, paste0(output.dir2, '/extended.figure1.d.5class.indel.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)
write.table(df6, paste0(output.dir2, '/extended.figure1.d.total.indel.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)

head(df5)
head(df6)

df6[df6$Gene == 'WT',]
103/180
0.5722222
