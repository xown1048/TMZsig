library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
library(ggh4x)

options(bitmapType='cairo')
options(scipen=10000)

##########################################################################################################################################
#### SBS signature point prop and count plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/signature/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/extended.figure3'
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
rowSums(prop.df[,-1])

prop.df = melt(prop.df, id.vars='Samples')
colnames(prop.df) = c('Sample', 'Signature', 'Proportion')

df2 = merge(df, prop.df, by=c('Sample', 'Signature'))
df3 = aggregate(Proportion~Signature+Group, df2, mean)
df4 = aggregate(Count~Signature+Group, df2, mean)

df4 = merge(df4, df3, by=c('Group', 'Signature'))
df4$Count = log2(df4$Count+1)

df4$Gene = gsub('\\-\\S+', '', df4$Group)
df4$Drug = gsub('\\S+\\-', '', df4$Group)
df4$Dose = sapply(strsplit(as.character(unlist(df4$Group)), '-'), '[[', 3)
df4$Dose = gsub('nodrug', 0, df4$Dose)
df4$Tissue = sapply(strsplit(as.character(unlist(df4$Group)), '-'), '[[', 2)

df5 = df4
df5$Group = gsub('-cisplatin', '-Cisplatin', df5$Group)
df5$Group = gsub('\\-c\\S+', '', df5$Group)

df5 = df5[df5$Tissue == 'TK6',]
df5 = df5[!grepl('-p', df5$Group),]
df5 = df5[df5$Drug %in% c('nodrug', 'TMZ'),]

'%notin%' = Negate('%in%')
df6 = df5[df5$Gene %notin% c('WAPL', 'ESCO1', 'REV3'),]

############################################################################################################################################
#####	Group level
g.df = df6[,c('Group', 'Gene', 'Drug', 'Dose', 'Tissue')]
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
#####	Draw plot
pathways = c(rep('WT', 1), rep('Apoptosis', 1), rep('DR', 3), rep('BER', 3), rep('NER', 2),
	rep('NER/BER', 1), rep('TMEJ', 1), rep('NHEJ', 1), rep('HR', 1), rep('PCNA', 1), rep('ICLR', 4),
	rep('TLS', 7), rep('MMR', 2), rep('MSH2 Double Knockouts', 10))

path = data.frame(Gene=gene.level, Path=pathways)
df7 = merge(df6, path, 'Gene')

sum.sig.df = df7[,c('Signature', 'Count')]
sum.sig.df = aggregate(Count~Signature, sum.sig.df, sum)
sum.sig.df = sum.sig.df[sum.sig.df$Count!=0,]

df7 = df7[df7$Signature %in% sum.sig.df$Signature,]
df7 = df7[df7$Signature %in% c('SBS5', 'SBS11'),]

pathway.level = unique(pathways)

df7$Gene = factor(df7$Gene, levels=gene.level)
df7$Group = factor(df7$Group, levels=group.level)
df7$Path = factor(df7$Path, levels=pathway.level)

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

strip.colors = c(rep(c('#9EDDFF', '#B0D9B1'), length(pathway.level)/2), rep('white', length(gene.level)), rep(c('white', '#F4CE14', '#C70039'), length(group.level)/3))
strip = strip_nested(background_x=elem_list_rect(fill=strip.colors), size='variable')

mid = mean(df7$Count)
print(mid)

mid = 8

pdf(file=paste0(output.dir, '/extended.figure3.b.sig.sbs.point.plot.pdf'), width=20, height=4)
ggplot(data=df7, aes(Group, Signature, color=Count)) + 
	geom_point(aes(size=Proportion)) +
	scale_color_gradient2(midpoint=mid, low='blue', mid='yellow', high='red', space ='Lab') +
	facet_nested(.~Path+Gene+Group, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip) +
	theme_light(base_size=10, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(axis.text.y=element_text(size=15), axis.text.x=element_blank(), axis.title=element_blank()) +
	theme(strip.text=element_text(colour='black', face='bold')) +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'),
		axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60'))
dev.off()

