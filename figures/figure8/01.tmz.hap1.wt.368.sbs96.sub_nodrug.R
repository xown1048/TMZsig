library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)

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

df3 = df2[df2$Group %in% c('WT-HAP1-nodrug', 'WT-HAP1-150microM-TMZ', 'WT-HAP1-300microM-TMZ'),]

###########################################################################################################################
groups = unique(df3$Group)
wt.group = 'WT-HAP1-nodrug'
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

merge.df$Substitution = substr(merge.df$MutationType, 3, 5)
merge.df$Side = paste0(substr(merge.df$MutationType, 1, 1), substr(merge.df$MutationType, 7, 7))
merge.df$Substitution = factor(merge.df$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))

merge.df2 = merge.df[merge.df$Sample==368,]

colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')

p = ggplot(merge.df2, aes(x=MutationType, y=Count2)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(Sample~Substitution, scales='free', space='free_x') + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank()) +
	xlab(NULL) + ylab(NULL)

pdf(file=file.path(output.dir, 'figure6.1.tmz.hap1.wt.368.sbs96.sub_nodrug.pdf'), width=90, height=30)
g <- ggplot_gtable(ggplot_build(p))
dev.off()

strips <- which(grepl('strip-', g$layout$name))

for (i in seq_along(strips)) {
	k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))

	g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- colors[i]
	g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- colors[i]
}	

pdf(file=file.path(output.dir, 'figure6.1.tmz.hap1.wt.368.sbs96.sub_nodrug.pdf'), width=90, height=30)
plot(g)
dev.off()

