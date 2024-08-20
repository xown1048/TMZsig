library(data.table)
library(ggplot2)
library(reshape2)
library(ggsignif)
library(Rsamtools)
library(RecordLinkage)

options(bitmapGroup='cairo')
options(scipen=10000)

rm(list=ls())

##############################################################################################################
#####	Directory
output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/extended.figure2'
if(!dir.exists(output.dir)) {dir.create(output.dir)}

fasta.file = FaFile(file='/BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230313.sample.list.txt')

##############################################################################################################
#####	Merge INDEL + SV and 1-100, over 100 + MH 20 bp (48% similarity 12 match/ 8 not-match)
merge.df = as.data.frame(fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part9/output/comb.mh.not.perfect.match/microhomology.txt'))

merge.df2 = merge.df
merge.df2$cate = ifelse(merge.df2$MH_LEN==0, 'no_mh', 'mh')

table(merge.df2$cate)
#   0     1     2     3     4     5     6     7     8     9    10 11-20   20+
#1639   231   138    80    37    20     8     8     3     3     2     8     7

merge.df2$cate2 = ifelse(merge.df2$DEL_LEN<=100, '1-100', 'over100')
merge.df2$cate3 = paste0('Del:M:', merge.df2$cate2, ':', merge.df2$cate)

head(merge.df2)

merge.df3 = merge.df2[,c('Sample', 'cate3')]

head(merge.df3)

df = as.data.frame.matrix(table(merge.df3$cate3, merge.df3$Sample))
df2 = tibble::rownames_to_column(df, 'MutationType')

#######################################################################################################################################
### Dot plot
meta$Sample = gsub('\\.\\S+', '', meta$Sample)

df3 = melt(df2, id.vars='MutationType')
colnames(df3) = c('MutationType', 'Sample', 'Count')
head(df3)

df4 = merge(df3, meta, by='Sample')

df4$Gene = gsub('\\-\\S+', '', df4$Group)
df4$Drug = gsub('\\S+\\-', '', df4$Group)
df4$Dose = sapply(strsplit(as.character(unlist(df4$Group)), '-'), '[[', 3)
df4$Dose = gsub('nodrug', 0, df4$Dose)
df4$Tissue = sapply(strsplit(as.character(unlist(df4$Group)), '-'), '[[', 2)

df5 = df4
df5$Group = gsub('-cisplatin', '-Cisplatin', df5$Group)
df5$Group = gsub('\\-c\\S+', '', df5$Group)
df5 = df5[!grepl('-p0', df5$Group),]

df5 = df5[df5$Drug=='nodrug',]
df5 = df5[df5$Tissue=='TK6',]
df5 = df5[df5$Gene %in% c('WT', 'ERCC1', 'FANCD2', 'FANCC'),]
df5 = df5[df5$MutationType=='Del:M:1-100:mh',]

######################################################################
##### Group level
g.df = df5[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
g.df2 = g.df[!grepl('-p', g.df$Group),]

gene.level = c('WT', c('MGMT', 'ALKBH2', 'ALKBH3'), c('ERCC1', 'XPA'), 'XPA_XRCC1', c('XRCC1', 'POLB', 'MPG'), c('RAD18', 'REV1', 'REV3', 'REV7', 'POLH', 'POLK', 'POLI', 'POLHKI'),
	c('FANCD2', 'FANCC', 'FANCM'), c('EXO1', 'MSH2'), 'p53', 'ATAD5', c('ESCO1', 'WAPL'), 'LIG4', 'MUS81', 'POLQ', 'RAD54L_54B',
	'ATAD5_MSH2', 'FANCD2_MSH2', 'MSH2_p53', 'MSH2_XRCC1', 'MSH2_ALKBH3', 'MSH2_RAD18', 'MSH2_ALKBH2', 'MSH2_MPG')

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

group.level = unique(group.level)

######################################################################
df5$Gene = gsub('\\-\\S+', '', df5$Group)
df5$Group = factor(df5$Group, levels=group.level)
df5$Gene = factor(df5$Gene, levels=gene.level)

colors = '#E76F51'

pdf(file=file.path(output.dir, 'extende.figure2.b.nodrug.tk6.t.test.mh.bar.plot.pdf'), width=60, height=30)
ggplot(df5, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Gene, scales='free') +
	xlab(NULL) + ylab(NULL) + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.key.size=unit(5,'cm'), legend.text=element_text(size=100),
		legend.title=element_text(size=100), strip.text=element_text(colour='black'),
		strip.background=element_rect(colour='grey', fill='white'),
		axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

my_comparisons = combn(levels(df5$Gene), 2, simplify=F)
t = compare_means(Count~Gene, comparisons=my_comparisons,
		group.by='MutationType', method='t.test', data=df5)
t = subset(t, p<= 0.1)
t.df = as.data.frame(t)

merge.t = data.frame()
for (c in unique(t$MutationType)) {
	tmp = subset(df5, MutationType == c)
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

colors = '#E76F51'

strip = strip_themed(background_x=elem_list_rect(fill=colors))

pdf(file=file.path(output.dir, 'extende.figure2.b.nodrug.tk6.t.test.mh.box.plot.pdf'), width=20, height=20)
ggplot(df5, aes(x=Gene, y=Count)) +
	geom_boxplot(aes(fill=Gene), colour='#264653', fill='white', size=2) +
	xlab(NULL) + ylab(NULL) +
	facet_wrap2(~MutationType, nrow=1, scales='free', strip=strip) + 
	stat_pvalue_manual(merge.t, size=20, label='symbol') +
	theme_light(base_size=80, base_family='sans') +
	theme(strip.text=element_text(colour='white', size=70, face='bold'), 
		axis.text.x=element_text(colour='black'),
		strip.background=element_rect(colour='grey', fill='white'))
dev.off()






