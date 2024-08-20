library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
library(janitor)
library(patchwork)

options(bitmapType='cairo')
options(scipen=10000)

rm(list=ls())

##########################################################################################################################################
#####	SBS Count and Proportion Box Plot
rm(list=ls())
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure5'
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

df2 = df2[df2$Tissue=='TK6',]
df2 = df2[df2$Drug %in% c('nodrug', 'TMZ'),]
df3 = df2[grepl('MSH2', df2$Group),]
df3 = df3[df3$Gene %in% c('MSH2', 'MSH2_MPG', 'MSH2_RAD18', 'MSH2_REV1', 'MSH2_REV7'),]

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
	#temp.df$Count2 = temp.df$Count2*500/temp.df$Dose2
	temp.df$Count2 = temp.df$Count2/temp.df$Dose2

	temp.df$Count2 = ifelse(temp.df$Count2 < 0, 0, temp.df$Count2)

	merge.df = rbind(merge.df, temp.df)
}

merge.df$Gene = factor(merge.df$Gene, levels=gene.level)
merge.df$Group = factor(merge.df$Group, levels=group.level)

gene.level2 = gene.level[gene.level %in% unique(df3$Gene)]
facet.name = gene.level2

temp.label = gsub('_', '\n', facet.name)

facet.label = temp.label
names(facet.label) = facet.name

sum.df = merge.df[,c('Sample', 'Count2')]
sum.df = aggregate(.~Sample, sum.df, sum)
sum.df = merge(sum.df, meta, by='Sample')
sum.df$Gene = gsub('\\-\\S+', '', sum.df$Group)

sum.df2 = sum.df[,c(4,2)]
mean.df = aggregate(.~Gene, sum.df2, mean)

control.mean.df = mean.df[mean.df$Gene=='MSH2',]
control.mean = control.mean.df[1,2]

#####	Cosmic signature colors
colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')
#colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#2A9D8F', '#264653')

#####	plot.margin=margin(Top, Right, Bottom, Left, 'cm')
p1 = ggplot(merge.df, aes(x=Sample, y=Count2)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	geom_hline(yintercept=control.mean, linetype='dashed', colour='black', linewidth=5) +
	scale_fill_manual(values=colors) +
	facet_grid(~Gene, scales='free_x', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=250, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(legend.position='none',
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.text=element_text(colour='black', size=200, face='bold.italic'),
		strip.background=element_rect(colour='white', fill='white'),
		panel.spacing=unit(0,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.3, 1, 0.01, 1, 'cm'))

merge.df2 = merge.df[merge.df$MutationType!='C>T',]
sum.df = merge.df2[,c('Sample', 'Count2')]
sum.df = aggregate(.~Sample, sum.df, sum)
sum.df = merge(sum.df, meta, by='Sample')
sum.df$Gene = gsub('\\-\\S+', '', sum.df$Group)

sum.df2 = sum.df[,c(4,2)]
mean.df = aggregate(.~Gene, sum.df2, mean)

control.mean.df = mean.df[mean.df$Gene=='MSH2',]
control.mean = control.mean.df[1,2]

colors = c('#01bbed', '#121212', '#ccc8c9', '#a2cd61', '#ecc6c4')

p2 = ggplot(merge.df2, aes(x=Sample, y=Count2)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	geom_hline(yintercept=control.mean, linetype='dashed', colour='red', linewidth=5) +
	scale_fill_manual(values=colors) +
	facet_grid(~Gene, scales='free_x', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=250, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(legend.position='none',
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.text=element_blank(),
		strip.background=element_blank(),
		panel.spacing=unit(0,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.01, 1, 0.3, 1, 'cm'))

pdf(file=file.path(output.dir, 'figure5.c.norm.1uM.msh2.sbs.6sub.barplot.pdf'), width=120, height=120)
p1/p2  + plot_layout(heights=unit(rep(c(1, 1), 1), rep('null', 2)))
dev.off()
