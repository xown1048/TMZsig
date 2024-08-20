library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
library(janitor)
library(ggh4x)
library(plotly)
library(patchwork)
library(ggplotify)
library(grid)
library(gridExtra)

options(bitmapType='cairo')
options(scipen=10000)

rm(list=ls())

##########################################################################################################################################
#####	ID Count and Proportion Box Plot
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/ID/sigProfilerExtractor.ID83.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.3/output/figure3'
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
df2 = df2[!grepl('-p', df2$Group),]
df2 = df2[df2$Drug %in% c('nodrug', 'TMZ'),]
df2 = df2[df2$Sample!='D6_270k',]

'%notin%' = Negate('%in%')
df3 = df2[df2$Gene %notin% c('WAPL', 'ESCO1', 'REV3'),]

############################################################################################################################################
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

############################################################################################################################################
#####	Draw plot 1 (WT)
df4 = df3[!grepl('MSH2', df3$Group),]
gene.level2 = gene.level[!grepl('MSH2', gene.level)]
group.level2 = group.level[!grepl('MSH2', group.level)]

pathways = c(rep('WT', 1), rep('Apoptosis', 1), rep('DR', 3), rep('BER', 3), rep('NER', 2),
	rep('NER/BER', 1), rep('TMEJ', 1), rep('NHEJ', 1), rep('HR', 1), rep('PCNA', 1), rep('ICLR', 4),
	rep('TLS', 7), rep('MMR', 1))

path = data.frame(Gene=gene.level2, Path=pathways)
df4 = merge(df4, path, 'Gene')

pathway.level = unique(pathways)

df4$Gene = factor(df4$Gene, levels=gene.level2)
df4$Group = factor(df4$Group, levels=group.level2)
df4$Path = factor(df4$Path, levels=pathway.level)

df4 = df4[order(factor(df4$Group, levels=group.level2)),]
sample.level = unique(df4$Sample)
df4$Sample = factor(df4$Sample, levels=sample.level)

facet.name = c(pathway.level, gene.level2, group.level2)
temp.label = gsub('\\S+\\-nodrug', '0', facet.name)
temp.label = gsub('-TMZ', '', temp.label)
temp.label = gsub('\\S+\\-', '', temp.label)
temp.label = gsub('microM', '', temp.label)
temp.label = gsub('_', '\n', temp.label)
temp.label = ifelse(grepl('\\.', temp.label)==TRUE, as.character(round(as.numeric(temp.label),1)), temp.label)

temp.label[2] = 'Apop.'
facet.label = temp.label
names(facet.label) = facet.name

df4$MutationType = factor(df4$MutationType, levels=c('1bp_Deletion', '1bp_Insertion', 'Repeat_Deletion', 'Repeat_Insertion', 'Microhomology'))

colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#264653')

strip.colors = c(rep(c('#9EDDFF', '#B0D9B1'), length(pathway.level)/2), '#9EDDFF', rep('white', length(gene.level2)), rep(c('white', '#F4CE14', '#C70039'), length(group.level2)/3))

strip = strip_nested(background_x=elem_list_rect(fill=strip.colors))

#####	plot.margin=margin(Top, Right, Bottom, Left, 'cm')
left.sample = length(unique(df4$Sample))
#####	Left
p1 = ggplot(df4, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_nested(~Path+Gene+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip) +
	theme_light(base_size=120, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_text(colour='black', face='bold', size=80)) +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
	theme(legend.position='none') +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'), panel.border=element_rect(color='grey60'), plot.margin=margin(1, 0.5, 1, 1, 'cm'))

############################################################################################################################################
#####	Draw plot 1 (MSH2)
df4 = df3[grepl('MSH2', df3$Group),]
gene.level2 = gene.level[grepl('MSH2', gene.level)]
group.level2 = group.level[grepl('MSH2', group.level)]

pathways = c(rep('MMR', 1), rep('MSH2 Double Knockouts', 10))

path = data.frame(Gene=gene.level2, Path=pathways)
df4 = merge(df4, path, 'Gene')

pathway.level = unique(pathways)

df4$Gene = factor(df4$Gene, levels=gene.level2)
df4$Group = factor(df4$Group, levels=group.level2)
df4$Path = factor(df4$Path, levels=pathway.level)

df4 = df4[order(factor(df4$Group, levels=group.level2)),]
sample.level = unique(df4$Sample)
df4$Sample = factor(df4$Sample, levels=sample.level)

facet.name = c(pathway.level, gene.level2, group.level2)
temp.label = gsub('\\S+\\-nodrug', '0', facet.name)
temp.label = gsub('-TMZ', '', temp.label)
temp.label = gsub('\\S+\\-', '', temp.label)
temp.label = gsub('microM', '', temp.label)
temp.label = gsub('_', '\n', temp.label)
temp.label = ifelse(grepl('\\.', temp.label)==TRUE, as.character(round(as.numeric(temp.label),1)), temp.label)

facet.label = temp.label
names(facet.label) = facet.name

df4$MutationType = factor(df4$MutationType, levels=c('1bp_Deletion', '1bp_Insertion', 'Repeat_Deletion', 'Repeat_Insertion', 'Microhomology'))

colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#264653')

strip.colors = c(rep(c('#9EDDFF', '#B0D9B1'), length(pathway.level)/2), rep('white', length(gene.level2)), rep(c('white', '#F4CE14', '#C70039'), length(group.level2)/3))

strip = strip_nested(background_x=elem_list_rect(fill=strip.colors))

right.sample = length(unique(df4$Sample))

#####	Right
p2 = ggplot(df4, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_nested(~Path+Gene+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip) +
	theme_light(base_size=120, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_text(colour='black', face='bold', size=80)) +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
	theme(legend.position='none') +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'), panel.border=element_rect(color='grey60'), plot.margin=margin(1, 1, 1, 0.01, 'cm'))

pdf(file=file.path(output.dir, 'figure3.b.nodrug.indel.5class.bar.plot.pdf'), width=200, height=40)
p1 + p2 + plot_layout(widths=unit(c(left.sample, right.sample), rep('null', 2)))
dev.off()

p2 = ggplot(df4, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_nested(~Path+Gene+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip) +
	theme_light(base_size=10, base_family='sans') +
	ylab(NULL) + xlab(NULL)

legend = cowplot::get_legend(p2)

pdf(file=file.path(output.dir, 'figure3.b.nodrug.indel.5class.legend.bar.plot.pdf'), width=5, height=5)
grid.newpage()
grid.draw(legend)
dev.off()
