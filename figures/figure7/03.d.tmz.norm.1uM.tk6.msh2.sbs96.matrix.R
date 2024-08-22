library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)
library(bayestestR)
library(patchwork)
library(patchwork)
library(ggh4x)

options(bitmapType='cairo')

############################################################################################################################################
### SBS 96 Matrix Bar Count Plot
rm(list=ls())
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure5'
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
df2 = df2[!grepl('-p', df2$Group),]

df2 = df2[df2$Tissue == 'TK6',]
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
	#temp.df$Count2 = temp.df$Count2/temp.df$Dose2*500
	temp.df$Count2 = temp.df$Count2/temp.df$Dose2

	temp.df$Count2 = ifelse(temp.df$Count2 < 0, 0, temp.df$Count2)

	merge.df = rbind(merge.df, temp.df)
}

mean.df = aggregate(Count2 ~ MutationType+Gene, merge.df, mean)
ci.high.df = aggregate(Count2 ~ MutationType+Gene, merge.df, function(x){ci(x, ci=0.95, method='ETI')$CI_high})
ci.low.df = aggregate(Count2 ~ MutationType+Gene, merge.df, function(x){ci(x, ci=0.95, method='ETI')$CI_low})

colnames(mean.df)[3] = 'Mean'
colnames(ci.high.df)[3] = 'CI_high'
colnames(ci.low.df)[3] = 'CI_low'

colnames(mean.df)[2] = 'Group'
colnames(ci.high.df)[2] = 'Group'
colnames(ci.low.df)[2] = 'Group'

head(merge.df)
head(mean.df)
head(ci.high.df)
head(ci.low.df)

merge.df = merge(mean.df, ci.high.df, by=c('MutationType', 'Group'))
merge.df = merge(merge.df, ci.low.df, by=c('MutationType', 'Group'))

merge.df$Substitution = substr(merge.df$MutationType, 3, 5)
merge.df$Side = paste0(substr(merge.df$MutationType, 1, 1), substr(merge.df$MutationType, 7, 7))
merge.df$Substitution = factor(merge.df$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))

merge.df$Gene = gsub('\\-\\S+', '', merge.df$Group)

head(merge.df)

merge.df$Gene = factor(merge.df$Gene, levels=gene.level)

gene.level2 = gene.level[gene.level %in% unique(df3$Gene)]

merge.df2 = merge.df[merge.df$Substitution!='C>T',]

colors = c('#01bbed', '#121212', '#ccc8c9', '#a2cd61', '#ecc6c4')
strip = strip_themed(background_x=elem_list_rect(fill=colors), text_x=elem_list_text(color=colors))
facet.name = c(c('C>A', 'C>G', 'T>A', 'T>C', 'T>G'), gene.level2)
facet.label = c(c('C>A', 'C>G', 'T>A', 'T>C', 'T>G'), gsub('_', '\n', gene.level2))
names(facet.label) = facet.name

p1 = ggplot(merge.df2, aes(x=MutationType, y=Mean)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	#geom_errorbar(aes(ymin=CI_low, ymax=CI_high)) +
	facet_grid2(Gene~Substitution, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip) + 
	theme_light(base_size=70, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_text(colour='white', size=80),
		strip.text.y=element_blank(),
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
		axis.text.y=element_text(size=85)) +
	theme(panel.spacing=unit(0,'lines'), axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60')) +
	xlab(NULL) + ylab(NULL)

merge.df2 = merge.df[merge.df$Substitution=='C>T',]

colors = '#db2d2a'
strip = strip_themed(background_x=elem_list_rect(fill=colors), text_x=elem_list_text(color=colors),
					background_y=elem_list_rect(fill=rep('white', length(gene.level2))), text_y=elem_list_text(color=rep('black', length(gene.level2))))

facet.name = c('C>T', gene.level2)
facet.label = c('C>T', gsub('_', '\n', gene.level2))
names(facet.label) = facet.name

p2 = ggplot(merge.df2, aes(x=MutationType, y=Mean)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	#geom_errorbar(aes(ymin=CI_low, ymax=CI_high)) +
	facet_grid2(Gene~Substitution, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip) + 
	theme_light(base_size=70, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_text(colour='white', size=80),
		strip.text.y=element_text(colour='black', size=100, face='bold.italic'),
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
		axis.text.y=element_text(size=85)) +
	theme(panel.spacing=unit(0,'lines'), axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60')) +
	xlab(NULL) + ylab(NULL)

pdf(file=file.path(output.dir, 'figure5.d.sbs96.tmz.norm.1uM.tk6.msh2.pdf'), width=90, height=75)
p1 + p2 + plot_layout(widths=unit(c(5, 1), rep('null', 2)))
dev.off()
