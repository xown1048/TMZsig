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

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output'
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

df3 = df2[grepl('MSH2', df2$Group) | grepl('WT', df2$Group) | grepl('EXO1', df2$Group),]
df3 = df3[grepl('nodrug', df3$Group),]
df3 = df3[grepl('TK6', df3$Group),]
df3 = df3[!grepl('-p', df3$Group),]

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
colnames(cosmic.df2) = c('MutationType', 'Group', 'Count')

cosmic.df3 = cosmic.df2[cosmic.df2$Group %in% c('SBS26', 'SBS44'),]

##################################################################################
##### RefSig_MMR1
refsig.df = as.data.frame(fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/refsig.mmr.txt'))

refsig.df2 = melt(refsig.df, id.vars='Group')
refsig.df2 = refsig.df2[,c('variable', 'Group', 'value')]
colnames(refsig.df2) = c('MutationType', 'Group', 'Count')

refsig.df3 = refsig.df2[refsig.df2$Group=='RefSig_MMR1',]
refsig.df3$Group = 'SigMMR1'

##################################################################################
##### sigOX
oxsig.df = as.data.frame(fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/ox.sig.txt'))

oxsig.df2 = melt(oxsig.df, id.vars='Group')
oxsig.df2 = oxsig.df2[,c('variable', 'Group', 'value')]

colnames(oxsig.df2) = c('MutationType', 'Group', 'Count')

oxsig.df3 = oxsig.df2[oxsig.df2$Group=='sigOX',]
oxsig.df3$Group = 'SigOX'

##################################################################################
##### Count Plot
mean.df = aggregate(Count ~ MutationType+Group, df3, mean)

mean.df$Substitution = substr(mean.df$MutationType, 3, 5)
mean.df$Side = paste0(substr(mean.df$MutationType, 1, 1), substr(mean.df$MutationType, 7, 7))
mean.df$Substitution = factor(mean.df$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))
mean.df$Gene = gsub('\\-\\S+', '', mean.df$Group)

merge.df = rbind(cosmic.df3, refsig.df3)
merge.df = rbind(merge.df, oxsig.df3)

merge.df$Substitution = substr(merge.df$MutationType, 3, 5)
merge.df$Side = paste0(substr(merge.df$MutationType, 1, 1), substr(merge.df$MutationType, 7, 7))
merge.df$Substitution = factor(merge.df$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))

mean.df$Group = gsub('-TK6-p10', '', mean.df$Group)

mean.df$Group = factor(mean.df$Group, levels=group.level)
mean.df$Gene = factor(mean.df$Gene, levels=gene.level)
merge.df$Group = factor(merge.df$Group, levels=rev(c('SBS44', 'SBS26', 'SigOX', 'SigMMR1')))


facet.name = c(c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'), gene.level)
facet.label = c(c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'), gsub('_', '\n', gene.level))
names(facet.label) = facet.name

colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')
strip = strip_themed(background_x=elem_list_rect(fill=colors), text_x=elem_list_text(color=colors),
					background_y=elem_list_rect(fill=rep('white', length(group.level))), text_y=elem_list_text(color=rep('black', length(group.level))))

#####	plot.margin=margin(Top, Right, Bottom, Left, 'cm')
#####	Top
p1 = ggplot(mean.df, aes(x=MutationType, y=Count)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Gene~Substitution, scales='free_x', space='free_x', strip=strip, labeller=as_labeller(facet.label)) +
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_text(colour='white', size=60),
		strip.text.y=element_text(colour='black', angle=360, hjust=0.1, size=150),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		panel.spacing=unit(2,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.3, 1, 0.01, 1, 'cm')) +
	xlab(NULL) + ylab(NULL)

#####	Bottom
p2 = ggplot(merge.df, aes(x=MutationType, y=Count)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Group~Substitution, scales='free', space='free_x') +
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_blank(),
		strip.text.y=element_text(colour='black', angle=360, hjust=0.1, size=150),
		axis.text.x=element_text(colour='black', angle=90, vjust=0.5, hjust=1, size=60),
		strip.background=element_blank(),
		panel.spacing=unit(2,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.01, 1, 0.3, 1, 'cm')) +
	xlab(NULL) + ylab(NULL)

pdf(file=file.path(output.dir, 'figure2.f.abs.nodrug.msh2.sbs96.matrix.plot.pdf'), width=125, height=100)
p1/p2 + plot_layout(heights=unit(c(length(group.level), 4), rep('null', 2)))
dev.off()
