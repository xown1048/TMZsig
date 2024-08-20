library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)
library(bayestestR)

options(bitmapType='cairo')

############################################################################################################################################
### SBS 96 Matrix Bar Count Plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.3/output/figure1'
if(!dir.exists(output.dir)){dir.create(output.dir)}

df <- melt(i, id.vars = 'MutationType')
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
df2 = df2[df2$Sample!='D6_270k',]

df3 = df2[df2$Group %in% c('WT-TK6-nodrug', 'WT-HAP1-nodrug'),]

#############################################################################################
#####	Draw plot
df3$Count = as.numeric(df3$Count)
df3$Count = df3$Count/18

mean.df = aggregate(Count ~ MutationType+Group, df3, mean)
ci.high.df = aggregate(Count ~ MutationType+Group, df3, function(x){ci(x, ci=0.95, method='ETI')$CI_high})
ci.low.df = aggregate(Count ~ MutationType+Group, df3, function(x){ci(x, ci=0.95, method='ETI')$CI_low})

colnames(mean.df)[3] = 'Mean'
colnames(ci.high.df)[3] = 'CI_high'
colnames(ci.low.df)[3] = 'CI_low'

head(df3)
head(mean.df)
head(ci.high.df)
head(ci.low.df)

merge.df = merge(mean.df, ci.high.df, by=c('MutationType', 'Group'))
merge.df = merge(merge.df, ci.low.df, by=c('MutationType', 'Group'))

merge.df$Substitution = substr(merge.df$MutationType, 3, 5)
merge.df$Side = paste0(substr(merge.df$MutationType, 1, 1), substr(merge.df$MutationType, 7, 7))
merge.df$Substitution = factor(merge.df$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))

head(merge.df)

group.level = c('WT-TK6-nodrug', 'WT-HAP1-nodrug')
merge.df$Group = factor(merge.df$Group, levels=group.level)

colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')

p = ggplot(merge.df, aes(x=MutationType, y=Mean)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	geom_errorbar(aes(ymin=CI_low, ymax=CI_high)) +
	facet_grid(Group~Substitution, scales='free_x', space='free_x') + 
	theme_light(base_size=70, base_family='sans') +
	ylab('Mutation Count') +
	theme(legend.position='none',
		axis.text.y=element_text(colour='black', size=100),
		strip.text.x=element_text(colour='white', size=50),
		strip.text.y=element_text(colour='black', size=60),
		axis.text.x=element_text(colour='black', angle=90, vjust=0.5, hjust=1)) +
	xlab(NULL) + ylab(NULL)

#https://stackoverflow.com/questions/53455092/r-ggplot2-change-colour-of-font-and-background-in-facet-strip
pdf(file=file.path(output.dir, 'figure1.c.wt.tk6.hap1.sbs.96bar.pdf'), width=90, height=30)
g <- ggplot_gtable(ggplot_build(p))
dev.off()

strips <- which(grepl('strip-', g$layout$name))

for (i in seq_along(strips)) {
	k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))

	g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- colors[i]
	g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- colors[i]
}	

pdf(file=file.path(output.dir, 'figure1.c.wt.tk6.hap1.sbs.96bar.pdf'), width=90, height=30)
plot(g)
dev.off()

