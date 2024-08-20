library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)
library(patchwork)
library(ggh4x)

options(bitmapType='cairo')

######################################################################
### SBS 96 Matrix Bar Count Plot
rm(list=ls())
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure4'
if(!dir.exists(output.dir)){dir.create(output.dir)}

df = melt(i, id.vars='MutationType')
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
df2 = df2[df2$Drug %in% c('nodrug', 'TMZ'),]
df3 = df2[df2$Gene %in% c('WT', 'MGMT', 'EXO1', 'MSH2', 'MSH2_REV7'),]

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

##################################################################################
##### COSMIC
cosmic.df = as.data.frame(fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/COSMIC_v3.3.1_SBS_GRCh38.txt'))
colnames(cosmic.df)[1] = 'MutationType'

cosmic.df2 = melt(cosmic.df, id.vars='MutationType')
colnames(cosmic.df2) = c('MutationType', 'Group', 'Proportion')

cosmic.df3 = cosmic.df2[cosmic.df2$Group %in% c('SBS5', 'SBS11'),]
head(cosmic.df3)

##################################################################################
##### Serena hiPSC_TMZ
serena.df = as.data.frame(fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/serena.data/MutagenMutationalSignature/Mutagen53_sub_signature.txt'))
serena.df2 = melt(serena.df, id.vars='MutationType')
colnames(serena.df2) = c('MutationType', 'Group', 'Proportion')

serena.df3 = serena.df2[serena.df2$Group=='Temozolomide (200 uM)',]
serena.df3$Group = 'hiPSC_TMZ'
head(serena.df3)

##################################################################################
##### Count Plot
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

	temp.df$Count2 = ifelse(temp.df$Count2 < 0, 0, temp.df$Count2)

	merge.df = rbind(merge.df, temp.df)
}

mean.df = aggregate(Count2 ~ MutationType+Gene, merge.df, mean)

mean.df$Substitution = substr(mean.df$MutationType, 3, 5)
mean.df$Side = paste0(substr(mean.df$MutationType, 1, 1), substr(mean.df$MutationType, 7, 7))

result.df = rbind(cosmic.df3, serena.df3)

result.df$Substitution = substr(result.df$MutationType, 3, 5)
result.df$Side = paste0(substr(result.df$MutationType, 1, 1), substr(result.df$MutationType, 7, 7))

mean.df = mean.df[mean.df$Substitution!='C>T',]
mean.df$Substitution = factor(mean.df$Substitution, levels=c('C>A', 'C>G', 'T>A', 'T>C', 'T>G'))

result.df = result.df[result.df$Substitution!='C>T',]
result.df$Substitution = factor(result.df$Substitution, levels=c('C>A', 'C>G', 'T>A', 'T>C', 'T>G'))

mean.df$Gene = factor(mean.df$Gene, levels=c('WT', 'MGMT', 'EXO1', 'MSH2', 'MSH2_REV7'))
result.df$Group = factor(result.df$Group, levels=c('hiPSC_TMZ', 'SBS5', 'SBS11'))

colors = c('#01bbed', '#121212', '#ccc8c9', '#a2cd61', '#ecc6c4')
strip = strip_themed(background_x=elem_list_rect(fill=colors), text_x=elem_list_text(color=colors),
					background_y=elem_list_rect(fill=rep('white', length(c('WT', 'MGMT', 'EXO1', 'MSH2', 'MSH2_REV7')))), text_y=elem_list_text(color=rep('black', length(c('WT', 'MGMT', 'EXO1', 'MSH2', 'MSH2_REV7')))))

facet.name = c(c('C>A', 'C>G', 'T>A', 'T>C', 'T>G'), c('WT', 'MGMT', 'EXO1', 'MSH2', 'MSH2_REV7'))
facet.label = gsub('_', '\n', facet.name)
names(facet.label) = facet.name

#####	plot.margin=margin(Top, Right, Bottom, Left, 'cm')
#####	Top
p1 = ggplot(mean.df, aes(x=MutationType, y=Count2)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Gene~Substitution, scales='free_x', space='free_x', strip=strip, labeller=as_labeller(facet.label)) +
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_text(colour='white', size=60),
		strip.text.y=element_text(colour='black', size=100),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		panel.spacing=unit(2,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.3, 1, 0.01, 1, 'cm')) +
	xlab(NULL) + ylab(NULL)

facet.name = c(c('C>A', 'C>G', 'T>A', 'T>C', 'T>G'), c('hiPSC_TMZ', 'SBS5', 'SBS11'))
facet.label = gsub('_', '\n', facet.name)
names(facet.label) = facet.name

#####	Bottom
p2 = ggplot(result.df, aes(x=MutationType, y=Proportion)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Group~Substitution, scales='free', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_blank(),
		strip.text.y=element_text(colour='black', size=100),
		axis.text.x=element_text(colour='black', angle=90, vjust=0.5, hjust=1, size=60),
		strip.background=element_blank(),
		panel.spacing=unit(2,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.01, 1, 0.3, 1, 'cm')) +
	xlab(NULL) + ylab(NULL)

pdf(file=file.path(output.dir, 'figure4.b.tmz.nrom.sbs96.matrix.plot.pdf'), width=125, height=100)
p1/p2 + plot_layout(heights=unit(c(5, 3), rep('null', 2)))
dev.off()

