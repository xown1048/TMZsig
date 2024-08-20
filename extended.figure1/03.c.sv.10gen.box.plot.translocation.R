library(data.table)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(janitor)
library(ggbreak)
library(scales)
library(ggpubr)
library(bayestestR)
library(ggh4x)

options(bitmapType='cairo')

rm(list=ls())

######################################################################
##### SV bar count for all samples
#vcfs = Sys.glob('/BiO/Research/UNIST-Toni-Sig-2020-1026/part10/output/delly.merge/10read.filter/p0wt.all/*.all.2pass.filtered.somatic.vcf')
vcfs = Sys.glob('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/delly.merge/no.read.filter/all/*.all.2pass.filtered.somatic.vcf')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/extended.figure1'
if(!dir.exists(output.dir)){dir.create(output.dir)}

no.samples = c()
samples = c()

df = data.frame()
for (vcf in vcfs) {
	data.dir = paste0('zgrep -v "^##" ', vcf)
	i = fread(data.dir, header=TRUE)
	colnames(i)[1] = 'CHROM'
	samples = append(colnames(i)[10], samples)

	if (nrow(i) == 0) {
		no.samples = append(colnames(i)[10], no.samples)
		next;
	}

	i2 = i[,c(1,2,3,4,5,8)]
	i2$Sample = colnames(i)[10]
	df = rbind(df, i2)
}

print(no.samples)
length(no.samples)

df2 = df

df2$END = gsub('\\S+\\;END=', '', df2$INFO)
df2$END = gsub('\\;\\S+', '', df2$END)

df2$LEN = as.numeric(df2$END) - as.numeric(df2$POS)

df2$INSLEN = gsub('\\S+\\;INSLEN=', '', df2$INFO)
df2$INSLEN = gsub('\\;\\S+', '', df2$INSLEN)

df2$INFO = NULL

df2$ALT = ifelse(grepl('[', df2$ALT, fixed=TRUE), '<TRA>', df2$ALT)
df2$ALT = ifelse(grepl(']', df2$ALT, fixed=TRUE), '<TRA>', df2$ALT)
df2$ALT = gsub('<', '', df2$ALT)
df2$ALT = gsub('>', '', df2$ALT)
df2$ALT = ifelse(df2$ALT=='INV', 'Inversions', 
			ifelse(df2$ALT=='DUP', 'Duplications',
			ifelse(df2$ALT=='TRA', 'Translocations',
			ifelse(df2$ALT=='DEL', 'Deletions',
			ifelse(df2$ALT=='INS', 'Insertions', 0)))))

df2$LEN2 = ifelse(df2$ALT=='Insertions', df2$INSLEN, df2$LEN)
df2$LEN2 = as.numeric(df2$LEN2)

df2[df2$ALT=='Insertions',]

temp.df = df2[df2$LEN2>200 & df2$ALT=='Insertions',]
temp.df2 = df2[df2$LEN2>200 & df2$ALT=='Deletions',]

'%notin%' = Negate('%in%')
temp.df3 = df2[df2$ALT %notin% c('Insertions', 'Deletions'),]

df3 = rbind(temp.df, temp.df2)
df3 = rbind(df3, temp.df3)

table(df3$ALT)
#Deletions Inversions Duplications Translocations
#964       136        103                 68
#3708      1883       1553                164

df4 = as.data.frame(table(df3$Sample, df3$ALT))
colnames(df4) = c('Sample', 'MutationType', 'Count')

df4$Sample = as.character(df4$Sample)
df4$MutationType = as.character(df4$MutationType)

##### no over 200 insertions
#all.types.order = c('Insertions', 'Deletions', 'Inversions', 'Duplications', 'Translocations')
all.types.order = c('Deletions', 'Inversions', 'Duplications', 'Translocations')

df5 = df4
df5$Sample = as.character(df5$Sample)
no.samples = setdiff(samples, unique(df5$Sample))

sv.types = all.types.order

if(length(no.samples)!=0) {
	for (i in 1:length(no.samples)) {
		for (j in 1:length(sv.types)) {
			df5[nrow(df5)+1,] = c(no.samples[i], sv.types[j], 0)
		}
	}
}

df6 = merge(df5, meta, by='Sample')
df6$Count = as.numeric(df6$Count)

df6$Gene = gsub('\\-\\S+', '', df6$Group)
df6$Drug = gsub('\\S+\\-', '', df6$Group)
df6$Dose = sapply(strsplit(as.character(unlist(df6$Group)), '-'), '[[', 3)
df6$Dose = gsub('nodrug', 0, df6$Dose)
df6$Tissue = sapply(strsplit(as.character(unlist(df6$Group)), '-'), '[[', 2)

df7 = df6
df7$Group = gsub('-cisplatin', '-Cisplatin', df7$Group)
df7$Group = gsub('\\-c\\S+', '', df7$Group)

df7 = df7[df7$Tissue=='TK6',]

gene.level = c('WT', 'XRCC1', 'POLHKI', 'REV1', 'REV7')
df7 = df7[df7$Gene %in% gene.level,]

df7 = df7[grepl('-p10', df7$Group),]

######################################################################
##### Group level
g.df = df7[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
tissue.level = c('TK6', 'HAP1')

gen.g.df = gen.g.df[order(factor(gen.g.df$Tissue, levels=tissue.level), factor(gen.g.df$Gene, levels=gene.level)),]
gen.group.level = gen.g.df$Group

group.level = unique(gen.group.level)

######################################################################
df7$Gene = factor(df7$Gene, levels=gene.level)
df7$Group = factor(df7$Group, levels=group.level)
df7$MutationType = factor(df7$MutationType, levels=c('Deletions', 'Duplications', 'Inversions', 'Translocations'))

sum.df = df7[,c('Sample', 'Count')]
sum.df = aggregate(.~Sample, sum.df, sum)
sum.df = merge(sum.df, meta, by='Sample')
sum.df$Gene = gsub('\\-\\S+', '', sum.df$Group)

test.df = compare_means(Count~Gene, method='t.test', data=sum.df)
test.df2 = as.data.frame(subset(test.df, group1=='WT'))

my_comparisons = combn(levels(df7$Gene), 2, simplify=F)
t = compare_means(Count~Gene, comparisons=my_comparisons,
		group.by='MutationType', method='t.test', data=df7)
t = subset(t, p<= 0.1)
t.df = as.data.frame(t)

merge.t = data.frame()
for (c in unique(t$MutationType)) {
	tmp = subset(df7, MutationType == c)
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

colors = c('#E76F51', '#E9C46A', '#6A994E', '#264653')

strip = strip_themed(background_x=elem_list_rect(fill=colors))

pdf(file=file.path(output.dir, 'extended.figure1.c.no.read.sv.p10.10gen.boxplot.pdf'), width=30, height=40)
ggplot(df7, aes(x=Gene, y=Count)) +
	geom_boxplot(aes(fill=Gene), colour='#264653', fill='white', size=2) +
	xlab(NULL) + ylab(NULL) +
	facet_wrap2(~MutationType, nrow=2, scales='free', strip=strip) + 
	stat_pvalue_manual(merge.t, size=25, label='symbol') +
	theme_light(base_size=100, base_family='sans') +
	theme(legend.key.size=unit(3,'cm'), legend.text=element_text(size=50), axis.title=element_text(size=100),
		legend.title=element_text(size=70), strip.text=element_text(colour='black', size=70, face='bold'), 
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face='bold.italic', colour='black', size=60),
		strip.background=element_rect(colour='grey', fill='white'))
dev.off()


