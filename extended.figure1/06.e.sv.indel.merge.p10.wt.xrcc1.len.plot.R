library(data.table)
library(ggplot2)
library(reshape2)
library(ggsignif)
library(Rsamtools)
library(bayestestR)
library(ggpubr)

options(bitmapGroup='cairo')
options(scipen=10000)

rm(list=ls())

##############################################################################################################
#####	Directory
output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/extended.figure1'
if(!dir.exists(output.dir)) {dir.create(output.dir)}

meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

##############################################################################################################
#####	INDEL
input.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/vcf_files/ID'
chrs = c(1:22, 'X', 'Y')
df = data.frame()
for (chr in chrs) {
	vcf = paste0(input.dir, '/', chr, '_seqinfo.txt')
	i = as.data.frame(fread(vcf))
	i2 = i[grepl('Del', i$V4),]
	df = rbind(df, i2)
}

head(df)

df = df[,-ncol(df)]
colnames(df) = c('SAMPLE', 'CHROM', 'POS', 'TYPE', 'REF', 'ALT')
df$MH_LEN = substr(df$TYPE, nchar(df$TYPE), nchar(df$TYPE))
df$MH_LEN = as.numeric(df$MH_LEN)
df$DEL_LEN = nchar(df$REF) - nchar(df$ALT)
df$DEL_LEN = as.numeric(df$DEL_LEN)
df$TYPE = substr(df$TYPE, 3, nchar(df$TYPE))

df2 = df

df2$CHROM = paste0('chr', df2$CHROM)
df2$START = df2$POS - df2$DEL_LEN + 1
df2$END = df2$POS + df2$DEL_LEN + df2$DEL_LEN

indel.df = df2
indel.df = indel.df[,c('SAMPLE', 'CHROM', 'START', 'END', 'DEL_LEN')]

head(indel.df)

##############################################################################################################
#####	SV
vcfs = Sys.glob('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/delly.merge/no.read.filter/all/*.all.2pass.filtered.somatic.vcf')

df = data.frame()
for (vcf in vcfs) {
	data.dir = paste0('zgrep -v "^##" ', vcf)
	j = fread(data.dir)

	colnames(j)[1] = 'CHROM'
	i = subset(j, FILTER == 'PASS')
	i = subset(i, ALT == '<DEL>')

	i2 = i[,c(1,2,4,5,8)]

	i2$END = gsub('\\S+\\;END=', '', i2$INFO)
	i2$END = gsub('\\;\\S+', '', i2$END)
	i2$END = as.integer(i2$END)
	i2$END = i2$END - 1
	i2$START = i2$POS + 1
	i2$SAMPLE = colnames(i)[10]

	i3 = i2[,c('SAMPLE', 'CHROM', 'START', 'END', 'ALT')]
	df = rbind(df, i3)
}

df = as.data.frame(df)

df2 = df
df2$DEL_LEN = df2$END - df2$START + 1 

sv.df = df2
sv.df = sv.df[,c('SAMPLE', 'CHROM', 'START', 'END', 'DEL_LEN')]

head(sv.df)

##############################################################################################################
#####	Merge INDEL + SV
merge.df = rbind(indel.df, sv.df)

table(duplicated(merge.df[,c(2,3,4)]))
#FALSE
#210211

merge.df = merge.df[!duplicated(merge.df[,c(2,3,4)]),]

merge.df2 = merge.df
colnames(merge.df2)[1] = 'Sample'

meta2 = meta
meta2$Sample = gsub('\\.\\S+', '', meta2$Sample)
merge.df2$Sample = gsub('\\.\\S+', '', merge.df2$Sample)

merge.df2 = merge(merge.df2, meta2, by='Sample')

dim(merge.df)
#[1] 210211     5
dim(merge.df2)
#[1] 210211     6

head(merge.df2)

merge.df3 = merge.df2[grepl('WT-TK6-p10', merge.df2$Group) | grepl('XRCC1-TK6-p10', merge.df2$Group),]

head(merge.df3)

merge.df3$Range = ifelse(merge.df3$DEL_LEN<101, '1-100bp',
					ifelse(merge.df3$DEL_LEN<501, '101-500bp',
					ifelse(merge.df3$DEL_LEN<1001, '501-1000bp', 'over1000bp')))

merge.df3$MutationType = paste0(merge.df3$Range, '_Deletions')

merge.df4 = as.data.frame(table(merge.df3$Sample, merge.df3$MutationType))
colnames(merge.df4) = c('Sample', 'MutationType', 'Count')

merge.df4 = merge(merge.df4, meta, by='Sample')

merge.df4$Group = gsub('-c\\S+', '', merge.df4$Group)
merge.df4$Group = gsub('TK6-', '', merge.df4$Group)

merge.df4$Group = factor(merge.df4$Group, levels=c('WT-p10', 'XRCC1-p10'))
merge.df4$MutationType = factor(merge.df4$MutationType, levels=c('1-100bp_Deletions', '101-500bp_Deletions', '501-1000bp_Deletions', 'over1000bp_Deletions'))

colors = c('#E76F51', '#E9C46A', '#6A994E', '#264653')

pdf(file=file.path(output.dir, 'extended.figure1.e.len.group.sv.bar.plot.pdf'), width=60, height=30)
ggplot(merge.df4, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Group, scales='free') +
	xlab(NULL) + ylab(NULL) + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.key.size=unit(5,'cm'), legend.text=element_text(size=100),
		legend.title=element_text(size=100), strip.text=element_text(colour='black'),
		strip.background=element_rect(colour='grey', fill='white'))
dev.off()

my_comparisons = combn(levels(merge.df4$Group), 2, simplify=F)
t = compare_means(Count~Group, comparisons=my_comparisons,
		group.by='MutationType', method='t.test', data=merge.df4)
t = subset(t, p<= 0.2)
t.df = as.data.frame(t.2)

merge.t = data.frame()
for (c in unique(t$MutationType)) {
	tmp = subset(merge.df4, MutationType == c)
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

merge.t$MutationType = factor(merge.t$MutationType, levels=c('1-100bp_Deletions', '101-500bp_Deletions', '501-1000bp_Deletions', 'over1000bp_Deletions'))
merge.df4$MutationType = factor(merge.df4$MutationType, levels=c('1-100bp_Deletions', '101-500bp_Deletions', '501-1000bp_Deletions', 'over1000bp_Deletions'))

colors = c('#E76F51', '#E9C46A', '#6A994E', '#264653')

strip = strip_themed(background_x=elem_list_rect(fill=colors))

pdf(file=file.path(output.dir, 'extended.figure1.e.len.group.sv.box.plot.pdf'), width=60, height=20)
ggplot(merge.df4, aes(x=Group, y=Count)) +
	geom_boxplot(aes(fill=Group), colour='#264653', fill='white', size=2) +
	xlab(NULL) + ylab(NULL) +
	facet_wrap2(~MutationType, nrow=1, scales='free', strip=strip) + 
	stat_pvalue_manual(merge.t, size=20, label='symbol') +
	theme_light(base_size=80, base_family='sans') +
	theme(strip.text=element_text(colour='white', size=70, face='bold'), 
		axis.text.x=element_text(colour='black'),
		strip.background=element_rect(colour='grey', fill='white'))
dev.off()



