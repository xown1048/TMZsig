library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
library(ggh4x)

options(bitmapType='cairo')
options(scipen=10000)

##########################################################################################################################################
##### SBS Signature Count and Proportion Box Plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part10/output/sigProfilerExtractor.p0wt.with.all/output/signature/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/extended.figure1'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i = as.data.frame(i)
df = melt(i, id.vars='Samples')
colnames(df) = c('Sample', 'Signature', 'Count')

meta2 = meta
meta2$Sample = gsub('\\.\\S+', '', meta2$Sample)

df = merge(df, meta2, by='Sample')

df$Gene = gsub('\\-\\S+', '', df$Group)
df$Drug = gsub('\\S+\\-', '', df$Group)
df$Dose = sapply(strsplit(as.character(unlist(df$Group)), '-'), '[[', 3)
df$Dose = gsub('nodrug', 0, df$Dose)
df$Tissue = sapply(strsplit(as.character(unlist(df$Group)), '-'), '[[', 2)

prop.df = cbind(Samples=i[,1], i[,-1]/rowSums(i[,-1]))
prop.df = melt(prop.df, id.vars='Samples')
colnames(prop.df) = c('Sample', 'Signature', 'Proportion')

df2 = merge(df, prop.df, by=c('Sample', 'Signature'))
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
#df3 = df3[df3$Signature %in% c('SBS40', 'SBS18'),]

df3$Gene = factor(df3$Gene, levels=gene.level)
df3$Group = factor(df3$Group, levels=group.level)
#df3$Signature = factor(df3$Signature, levels=c('SBS40', 'SBS18'))

sum.sig.df = df3[,c('Signature', 'Count')]
sum.sig.df = aggregate(Count~Signature, sum.sig.df, sum)
sum.sig.df = sum.sig.df[sum.sig.df$Count!=0,]

df3 = df3[df3$Signature %in% sum.sig.df$Signature,]

my_comparisons = combn(levels(df3$Gene), 2, simplify=F)
t = compare_means(Count~Gene, comparisons=my_comparisons,
		group.by='Signature', method='t.test', data=df3)
t = subset(t, p<= 0.1)
t.df = as.data.frame(t)

merge.t = data.frame()
for (c in unique(t$Signature)) {
	tmp = subset(df3, Signature == c)
	max = max(as.numeric(tmp$Count))
	tmp.t = subset(t.df, Signature == c)

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

colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#2A9D8F', '#118AB2', '#264653', '#613F75')
#colors = c('#613F75', '#2A9D8F')

strip = strip_themed(background_x=elem_list_rect(fill=colors))

pdf(file=file.path(output.dir, 'extended.figure1.b.sig.sbs.10gen.p10.box.plot.pdf'), width=60, height=40)
ggplot(df3, aes(x=Gene, y=Count)) +
	geom_boxplot(aes(fill=Gene), colour='#264653', fill='white', size=2) +
	xlab(NULL) + ylab(NULL) +
	facet_wrap2(~Signature, nrow=2, scales='free', strip=strip) + 
	stat_pvalue_manual(merge.t, size=25, label='symbol') +
	theme_light(base_size=100, base_family='sans') +
	theme(legend.key.size=unit(3,'cm'), legend.text=element_text(size=50), axis.title=element_text(size=100),
		legend.title=element_text(size=70), strip.text=element_text(colour='black', size=80, face='bold'), 
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face='bold.italic', colour='black'), strip.background=element_rect(colour='grey', fill='white'))
dev.off()
