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

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure1'
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

df3 = df2[grepl('-p10', df2$Group),]

##################################################################################
##### COSMIC
cosmic.df = as.data.frame(fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/COSMIC_v3.3.1_SBS_GRCh38.txt'))
colnames(cosmic.df)[1] = 'MutationType'

cosmic.df2 = melt(cosmic.df, id.vars='MutationType')
colnames(cosmic.df2) = c('MutationType', 'Group', 'Count')

cosmic.df3 = cosmic.df2[cosmic.df2$Group %in% c('SBS5', 'SBS40'),]

##################################################################################
##### P10 Count Plot
mean.df = aggregate(Count ~ MutationType+Group, df3, mean)

mean.df$Substitution = substr(mean.df$MutationType, 3, 5)
mean.df$Side = paste0(substr(mean.df$MutationType, 1, 1), substr(mean.df$MutationType, 7, 7))
mean.df$Substitution = factor(mean.df$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))

cosmic.df3$Substitution = substr(cosmic.df3$MutationType, 3, 5)
cosmic.df3$Side = paste0(substr(cosmic.df3$MutationType, 1, 1), substr(cosmic.df3$MutationType, 7, 7))
cosmic.df3$Substitution = factor(cosmic.df3$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))

mean.df$Group = gsub('-TK6-p10', '', mean.df$Group)

mean.df$Group = factor(mean.df$Group, levels=c('WT', 'XRCC1', 'POLHKI', 'REV1', 'REV7'))
cosmic.df3$Group = factor(cosmic.df3$Group, levels=c('SBS5', 'SBS40'))

colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')
strip = strip_themed(background_x=elem_list_rect(fill=colors), text_x=elem_list_text(color=colors),
					background_y=elem_list_rect(fill=rep('white', 5)), text_y=elem_list_text(color=rep('black', 5)))

#####	plot.margin=margin(Top, Right, Bottom, Left, 'cm')
#####	Top
p1 = ggplot(mean.df, aes(x=MutationType, y=Count)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Group~Substitution, scales='free_x', space='free_x', strip=strip) + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_text(colour='white', size=60),
		strip.text.y=element_text(colour='black', size=100),
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		panel.spacing=unit(0,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.3, 1, 0.01, 1, 'cm')) +
	xlab(NULL) + ylab(NULL)

#####	Bottom
p2 = ggplot(cosmic.df3, aes(x=MutationType, y=Count)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Group~Substitution, scales='free', space='free_x') + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_blank(),
		strip.text.y=element_text(colour='black', size=100),
		axis.text.x=element_text(colour='black', angle=90, vjust=0.5, hjust=1, size=60),
		strip.background=element_blank(),
		panel.spacing=unit(0,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.01, 1, 0.3, 1, 'cm')) +
	xlab(NULL) + ylab(NULL)

pdf(file=file.path(output.dir, 'figure1.f.abs.10gen.sbs96.matrix.plot.pdf'), width=100, height=65)
p1/p2 + plot_layout(heights=unit(c(2, 0.8), rep('null', 2)))
dev.off()
