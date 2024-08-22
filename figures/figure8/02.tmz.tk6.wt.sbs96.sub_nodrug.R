library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)
library(bayestestR)

options(bitmapType='cairo')

###########################################################################################################################
rm(list=ls())
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure6'
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

df3 = df2[df2$Group %in% c('WT-TK6-nodrug', 'WT-TK6-5microM-TMZ'),]

###########################################################################################################################
groups = unique(df3$Group)
wt.group = 'WT-TK6-nodrug'
case.groups = groups[groups!=wt.group]

merge.df = data.frame()
for (a.group in case.groups){
	temp.df = subset(df3, Group==a.group)

	wt.temp.df = subset(df3, Group==wt.group)
	wt.temp.df2 = wt.temp.df[,c(2,3)]
	wt.temp.df3 = aggregate(Count ~ MutationType, wt.temp.df2, mean)

	temp.df$Count2 = temp.df$Count - wt.temp.df3$Count[match(temp.df$MutationType, wt.temp.df3$MutationType)]
	temp.df$Count2 = ifelse(temp.df$Count2 < 0, 0, temp.df$Count2)

	merge.df = rbind(merge.df, temp.df)
}

mean.df = aggregate(Count2 ~ MutationType+Group, merge.df, mean)
ci.high.df = aggregate(Count2 ~ MutationType+Group, merge.df, function(x){ci(x, ci=0.95, method='ETI')$CI_high})
ci.low.df = aggregate(Count2 ~ MutationType+Group, merge.df, function(x){ci(x, ci=0.95, method='ETI')$CI_low})

colnames(mean.df)[3] = 'Mean'
colnames(ci.high.df)[3] = 'CI_high'
colnames(ci.low.df)[3] = 'CI_low'

head(mean.df)
head(ci.high.df)
head(ci.low.df)

merge.df = merge(mean.df, ci.high.df, by=c('MutationType', 'Group'))
merge.df = merge(merge.df, ci.low.df, by=c('MutationType', 'Group'))

head(merge.df)

merge.df$Substitution = substr(merge.df$MutationType, 3, 5)
merge.df$Side = paste0(substr(merge.df$MutationType, 1, 1), substr(merge.df$MutationType, 7, 7))
merge.df$Substitution = factor(merge.df$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))
merge.df$Gene = gsub('\\-\\S+', '', merge.df$Group)

head(merge.df)

colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')

p = ggplot(merge.df, aes(x=MutationType, y=Mean)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(Gene~Substitution, scales='free', space='free_x') + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank()) +
	xlab(NULL) + ylab(NULL)

pdf(file=file.path(output.dir, 'figure6.2.tmz.tk6.wt.sbs96.sub_nodrug.pdf'), width=90, height=30)
g <- ggplot_gtable(ggplot_build(p))
dev.off()

strips <- which(grepl('strip-', g$layout$name))

for (i in seq_along(strips)) {
	k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))

	g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- colors[i]
	g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- colors[i]
}	

pdf(file=file.path(output.dir, 'figure6.2.tmz.tk6.wt.sbs96.sub_nodrug.pdf'), width=90, height=30)
plot(g)
dev.off()
