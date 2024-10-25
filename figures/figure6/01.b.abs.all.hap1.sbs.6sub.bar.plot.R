library(data.table)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(janitor)
library(ggbreak)
library(scales)
library(ggpubr)
library(bayestestR)
library(patchwork)
library(ggh4x)

options(bitmapType='cairo')

rm(list=ls())

######################################################################
##### SBS bar absolute number for all samples
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/figure4'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i$MutationType = substr(i$MutationType, 3, 5)
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
df2 = df2[!grepl('-p0', df2$Group),]
df2 = df2[df2$Tissue=='HAP1',]

######################################################################
##### Group level
g.df = df2[,4:8]
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

######################################################################
#####	Cosmic signature colors
colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')

df2$Group = factor(df2$Group, levels=group.level)
df3 = df2[df2$Gene %in% c('WT', 'MGMT'),]

p1 = ggplot(df3, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Group, scales='free_x', space='free_x') +
	theme_light(base_size=200, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_blank(), axis.ticks.x=element_blank()) +
	theme(legend.position='none') +
	theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), strip.background=element_rect(colour='white', fill='white'),
		axis.title=element_text(size=250, face='bold'))
dev.off()

df3 = df2[df2$Gene %in% c('WT', 'MSH2'),]

p2 = ggplot(df3, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Group, scales='free_x', space='free_x') +
	theme_light(base_size=200, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_blank(), axis.ticks.x=element_blank()) +
	theme(legend.position='none') +
	theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), strip.background=element_rect(colour='white', fill='white'),
		axis.title=element_text(size=250, face='bold'))


pdf(file=file.path(output.dir, 'figure4.d.abs.all.hap1.sbs.6sub.barplot.pdf'), width=120, height=60)
p1 + p2 + plot_layout(widths=unit(c(1, 1), rep('null', 2)))
dev.off()

######################################################################
##### Count write table
output.dir2 = paste0(output.dir, '/table')
if(!dir.exists(output.dir2)){dir.create(output.dir2)}

table.df0 = df2
head(table.df0)

table.df = table.df0[,c('Count', 'MutationType', 'Gene')]
mean.df = aggregate(Count~MutationType+Gene, table.df, mean)
sd.df = aggregate(Count~MutationType+Gene, table.df, sd)
colnames(mean.df)[3] = 'Mean'
colnames(sd.df)[3] = 'Sd'

table.df2 = merge(mean.df, sd.df, by=c('Gene', 'MutationType'))
table.df2$Sd = round(table.df2$Sd, 1)
table.df2$Mean = round(table.df2$Mean, 1)

table.df2 = table.df2[order(factor(table.df2$Gene, levels=group.level), table.df2$MutationType),]

sum.df = aggregate(Count~Sample+Gene, table.df0, sum)
mean.df = aggregate(Count~Gene, sum.df, mean)
sd.df = aggregate(Count~Gene, sum.df, sd)
colnames(mean.df)[2] = 'Mean'
colnames(sd.df)[2] = 'Sd'

table.df3 = merge(mean.df, sd.df, by='Gene')
table.df3$Sd = round(table.df3$Sd, 1)
table.df3$Mean = round(table.df3$Mean, 1)

table.df3 = table.df3[order(factor(table.df3$Gene, levels=group.level)),]

head(table.df2)
head(table.df3)

write.table(table.df2, paste0(output.dir2, '/figure4.d.6sub.sbs.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)
write.table(table.df3, paste0(output.dir2, '/figure4.d.total.sbs.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)
write.table(sum.df, paste0(output.dir2, '/figure4.d.total.sample.sbs.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)

(766+697)/2