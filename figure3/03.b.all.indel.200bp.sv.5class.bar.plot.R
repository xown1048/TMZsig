library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
library(janitor)
library(ggh4x)
library(plotly)
library(patchwork)
library(ggplotify)

options(bitmapType='cairo')
options(scipen=10000)

rm(list=ls())

##########################################################################################################################################
#####	SV Count and Proportion Box Plot
vcfs = Sys.glob('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/delly.merge/no.read.filter/all/*.all.2pass.filtered.somatic.vcf')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/figure3'
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
#[1] 163

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
df2[df2$ALT=='Duplications',]
df2[df2$ALT=='Inversions',]

summary(df2[df2$ALT=='Insertions', LEN2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#24.00   27.00   28.00   28.31   29.00   36.00
summary(df2[df2$ALT=='Deletions', LEN2])
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#18        42       361    192695      1907 135710121
summary(df2[df2$ALT=='Translocations', LEN2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1       1       1       1       1       1
summary(df2[df2$ALT=='Duplications', LEN2])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#136      604     9113  1423638    87911 50865216
summary(df2[df2$ALT=='Inversions', LEN2])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#123     1371     7701  2959983    26728 98587167


#####	Only indel 200bp, 다른 것들은 short indel에 없는 event이기 때문에 
temp.df = df2[df2$LEN2>200 & df2$ALT=='Insertions',]
temp.df2 = df2[df2$LEN2>200 & df2$ALT=='Deletions',]

'%notin%' = Negate('%in%')
temp.df3 = df2[df2$ALT %notin% c('Insertions', 'Deletions'),]

df3 = rbind(temp.df, temp.df2)
df3 = rbind(df3, temp.df3)

table(df3$ALT)
#Deletions Inversions Duplications Translocations
#607       100        60                  47

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

df7 = df7[df7$Tissue == 'TK6',]
df7 = df7[!grepl('-p', df7$Group),]
df7 = df7[df7$Drug %in% c('nodrug', 'TMZ'),]

'%notin%' = Negate('%in%')
df8 = df7[df7$Gene %notin% c('WAPL', 'ESCO1', 'REV3'),]
#df8 = df8[df8$Sample %notin% c('D6_270k', 'E6_2.5Tmz_204k'),]
df8 = df8[df8$Sample!='D6_270k',]

############################################################################################################################################
#####	Group level
g.df = df8[,4:8]
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

############################################################################################################################################
#####	Draw plot 1 (WT)
pathways = c(rep('WT', 1), rep('Apoptosis', 1), rep('DR', 3), rep('BER', 3), rep('NER', 2),
	rep('NER/BER', 1), rep('TMEJ', 1), rep('NHEJ', 1), rep('HR', 1), rep('PCNA', 1), rep('ICLR', 4),
	rep('TLS', 7), rep('MMR', 2), rep('MSH2 Double Knockouts', 10))

path = data.frame(Gene=gene.level, Path=pathways)
df9 = merge(df8, path, 'Gene')

pathway.level = unique(pathways)

df9$Gene = factor(df9$Gene, levels=gene.level)
df9$Group = factor(df9$Group, levels=group.level)
df9$Path = factor(df9$Path, levels=pathway.level)
df9$MutationType = factor(df9$MutationType, levels=c('Deletions', 'Duplications', 'Inversions', 'Translocations'))

df9 = df9[order(factor(df9$Group, levels=group.level)),]
sample.level = unique(df9$Sample)
df9$Sample = factor(df9$Sample, levels=sample.level)

facet.name = c(pathway.level, gene.level, group.level)
temp.label = gsub('\\S+\\-nodrug', '0', facet.name)
temp.label = gsub('-TMZ', '', temp.label)
temp.label = gsub('\\S+\\-', '', temp.label)
temp.label = gsub('microM', '', temp.label)
temp.label = gsub('_', '\n', temp.label)
temp.label = ifelse(grepl('\\.', temp.label)==TRUE, as.character(round(as.numeric(temp.label),1)), temp.label)

temp.label[2] = 'Apop.'
facet.label = temp.label
names(facet.label) = facet.name

colors = c('#E76F51', '#E9C46A', '#6A994E', '#264653')
df9$MutationType = factor(df9$MutationType, levels=c('Deletions', 'Duplications', 'Inversions', 'Translocations'))

strip.colors = c(rep(c('#9EDDFF', '#B0D9B1'), length(pathway.level)/2), rep('white', length(gene.level)), rep(c('white', '#F4CE14', '#C70039'), length(group.level)/3))

strip = strip_nested(background_x=elem_list_rect(fill=strip.colors))

#####	plot.margin=margin(Top, Right, Bottom, Left, 'cm')
left.sample = length(unique(df9$Sample))
#####	Left
pdf(file=file.path(output.dir, 'figure3.b.nodrug.sv.5class.bar.plot.pdf'), width=200, height=40)
ggplot(df9, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_nested(~Path+Gene+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip) +
	theme_light(base_size=120, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_text(colour='black', face='bold', size=80)) +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
	theme(legend.position='none') +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'), panel.border=element_rect(color='grey60'), plot.margin=margin(1, 0.5, 1, 1, 'cm'))
dev.off()


