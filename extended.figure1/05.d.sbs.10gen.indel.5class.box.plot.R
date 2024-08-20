library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
library(janitor)
library(ggh4x)

options(bitmapType='cairo')
options(scipen=10000)

rm(list=ls())

##########################################################################################################################################
##### SBS Count and Proportion Box Plot
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/ID/sigProfilerExtractor.ID83.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/extended.figure1'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i = as.data.frame(i)
i$Len = substr(i$MutationType, 9, 9)
i$Len = with(i, ifelse(substr(MutationType, 3, 5) == 'Del', as.numeric(Len)+1, Len))
i$MutationType = substr(i$MutationType, 1, 8)
i$MutationType = paste0(i$MutationType, i$Len)
i$Len = NULL

i$MutationType = c(rep('1bp_Deletion',12), rep('1bp_Insertion',12), rep('Repeat_Deletion',24), rep('Repeat_Insertion',24), rep('Microhomology',11))
i = aggregate(.~MutationType, i, sum)
i = adorn_totals(i)

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

df2 = df2[df2$Tissue=='TK6',]
gene.level = c('WT', 'XRCC1', 'POLHKI', 'REV1', 'REV7')
df2 = df2[df2$Gene %in% gene.level,]

df3 = df2[grepl('-p10', df2$Group),]

######################################################################
##### Group level
g.df = df3[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
tissue.level = c('TK6', 'HAP1')

gen.g.df = gen.g.df[order(factor(gen.g.df$Tissue, levels=tissue.level), factor(gen.g.df$Gene, levels=gene.level)),]
gen.group.level = gen.g.df$Group

group.level = unique(gen.group.level)

######################################################################
df3$Gene = factor(df3$Gene, levels=gene.level)
df3$Group = factor(df3$Group, levels=group.level)

my_comparisons = combn(levels(df3$Gene), 2, simplify=F)
t = compare_means(Count~Gene, comparisons=my_comparisons,
		group.by='MutationType', method='t.test', data=df3)
t = subset(t, p<= 0.1)
t.df = as.data.frame(t)

merge.t = data.frame()
for (c in unique(t$MutationType)) {
	tmp = subset(df3, MutationType == c)
	max = max(as.numeric(tmp$Count))
	tmp.t = subset(t.df, MutationType == c)

	if(nrow(tmp.t) == 1) {
		tmp.t$'y.position' = max + max/10
	} else {
	interval = max/nrow(tmp.t)
	max2 = max + interval
	to = max2 + ((interval*(nrow(tmp.t)-1))*0.9)
	tmp.t$'y.position' = seq(from=max2, to=to, by=interval*0.9)
	}

	merge.t = rbind(merge.t, tmp.t)
}

merge.t$p = as.numeric(merge.t$p)
merge.t$p.format = as.numeric(merge.t$p.format)
merge.t$symbol = ifelse(merge.t$p <= 0.05, merge.t$p.signif, '#')

merge.t = tibble::as_tibble(merge.t)

merge.t$MutationType = factor(merge.t$MutationType, levels=c('1bp_Deletion', '1bp_Insertion', 'Repeat_Deletion', 'Repeat_Insertion', 'Microhomology', 'Total'))
df3$MutationType = factor(df3$MutationType, levels=c('1bp_Deletion', '1bp_Insertion', 'Repeat_Deletion', 'Repeat_Insertion', 'Microhomology', 'Total'))
colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#264653', 'lightgrey')

strip = strip_themed(background_x=elem_list_rect(fill=colors))

pdf(file=file.path(output.dir, 'extended.figure1.d.abs.tk6.p10.10gen.indel.5class.boxplot.pdf'), width=60, height=40)
ggplot(df3, aes(x=Gene, y=Count)) +
	geom_boxplot(aes(fill=Gene), colour='#264653', fill='white', size=2) +
	xlab(NULL) + ylab(NULL) +
	facet_wrap2(~MutationType, nrow=2, scales='free', strip=strip) + 
	stat_pvalue_manual(merge.t, size=20, label='symbol') +
	theme_light(base_size=100, base_family='sans') +
	theme(strip.text=element_text(colour='white', size=80, face='bold'), 
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face='bold.italic', colour='black'),
		strip.background=element_rect(colour='grey', fill='white'))
dev.off()

