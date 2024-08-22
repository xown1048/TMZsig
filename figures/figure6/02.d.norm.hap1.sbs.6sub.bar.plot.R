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

options(bitmapType='cairo')
options(scipen=10000)

rm(list=ls())

##########################################################################################################################################
#####	SBS Count and Proportion Box Plot
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

df2 = df2[df2$Tissue == 'HAP1',]
df2 = df2[!grepl('-p', df2$Group),]
df2 = df2[df2$Drug %in% c('nodrug', 'TMZ'),]

df3 = df2[df2$Gene %in% c('WT', 'MGMT', 'MSH2'),]

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

pathways = c('WT', 'DR', 'MMR')
gene.level = c('WT', 'MGMT', 'MSH2')

path = data.frame(Gene=gene.level, Path=pathways)
merge.df = merge(merge.df, path, 'Gene')

pathway.level = unique(pathways)

merge.df$Gene = factor(merge.df$Gene, levels=gene.level)
merge.df$Group = factor(merge.df$Group, levels=group.level)
merge.df$Path = factor(merge.df$Path, levels=pathway.level)

merge.df = merge.df[order(factor(merge.df$Group, levels=group.level)),]
sample.level = unique(merge.df$Sample)
merge.df$Sample = factor(merge.df$Sample, levels=sample.level)

sum.df = merge.df[,c('Sample', 'Count2')]
sum.df = aggregate(.~Sample, sum.df, sum)
sum.df = merge(sum.df, meta, by='Sample')
sum.df$Gene = gsub('\\-\\S+', '', sum.df$Group)

test.df = compare_means(Count2~Gene, method='t.test', data=sum.df)
test.df2 = as.data.frame(subset(test.df, group1=='WT'))
test.df3 = as.data.frame(subset(test.df, group2=='WT'))
test.df3 = test.df3[,c(1,3,2,4:8)]
colnames(test.df3) = colnames(test.df2)
test.df4 = rbind(test.df2, test.df3)

row_genes = gene.level[gene.level!='WT']
test.df4 = test.df4[match(row_genes, test.df4$group2),]

sum.df2 = sum.df[,c(4,2)]
mean.df = aggregate(.~Gene, sum.df2, mean)

wt.mean.df = mean.df[mean.df$Gene=='WT',]
mean.df = mean.df[mean.df$Gene!='WT',]

mean.df = mean.df[match(row_genes, mean.df$Gene),]

mean.df$Count2 = mean.df$Count2/wt.mean.df[1,2]

test.df5 = test.df4
test.df5$Ratio = mean.df$Count2
test.df5$p.format = as.numeric(test.df5$p.format)
test.df5$colors = ifelse(test.df5$p.format>0.1, 'white',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.05, '#e899af',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.01, '#dd6688',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.001, '#d23260',
					ifelse(test.df5$Ratio>1 & test.df5$p.format<=0.001, '#C70039',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.05, '#DDF2FD',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.01, '#9BBEC8',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.001, '#427D9D', '#164863'))))))))

#####	Down
#####	0.1, 0.05, 0.01, 0.001
#####	#DDF2FD, #9BBEC8, #427D9D, #164863

#####	UP
#####	0.1, 0.05, 0.01, 0.001
#####	#e899af, #dd6688, #d23260, #C70039

pval.strips.colors = c('white', test.df5$colors)

pval.facet.label = ifelse(test.df5$p.format<0.05, paste0(test.df5$p.signif, '\n', test.df5$group2), 
				ifelse(test.df5$p.format<0.1, paste0('#', '\n', test.df5$group2), test.df5$group2))

pval.facet.label = c('WT', pval.facet.label)

facet.name = c(pathway.level, gene.level)
temp.label = c(pathway.level, pval.facet.label)
temp.label = gsub('_', '\n', temp.label)

facet.label = temp.label
names(facet.label) = facet.name

#####	Cosmic signature colors
colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')

strip.colors = c(c('#9EDDFF', '#B0D9B1', '#9EDDFF'), pval.strips.colors)
strip = strip_nested(background_x=elem_list_rect(fill=strip.colors), size='variable')

pdf(file=file.path(output.dir, 'figure4.f.norm.sbs.hap1.6sub.bar.plot.pdf'), width=30, height=28)
ggplot(merge.df, aes(x=Sample, y=Count2)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	geom_hline(yintercept=wt.mean.df[1,2], linetype='dashed', colour='black', linewidth=3) +
	scale_fill_manual(values=colors) +
	facet_nested(~Path+Gene, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip) +
	theme_light(base_size=100, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(legend.position='none') +
	theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1), axis.title=element_text(size=150, face='bold')) +
	theme(strip.text=element_text(colour='black', size=50, face='bold')) +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'),
		panel.border=element_rect(color='grey60'))
dev.off()
