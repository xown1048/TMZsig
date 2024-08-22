library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)
library(lsa)

options(bitmapType='cairo')

##########################################################################################################################################
#####	SBS 96 Matrix Bar Count Plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure2'
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

df2 = df2[df2$Tissue=='TK6',]
df2 = df2[df2$Drug=='nodrug',]
df3 = df2[grepl('MSH2', df2$Group),]

############################################################################################################################################
#####	Group level
g.df = df3[,4:8]
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

##################################################################################
##### MSH2 Subtraction
groups = unique(df3$Group)
wt.group = 'MSH2-TK6-nodrug'
case.groups = groups[groups!='MSH2-TK6-nodrug']

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
colnames(mean.df)[3] = 'Mean'

##################################################################################
##### To proportion
mean.df2 = dcast(mean.df, MutationType~Group)
mean.df3 = as.data.frame(t(t(mean.df2[,-1])/colSums(mean.df2[,-1])))
colSums(mean.df3)

mean.df4 = cbind(mean.df2[,1], mean.df3)
colnames(mean.df4)[1] = 'MutationType'

mean.df5 = melt(mean.df4, id.vars='MutationType')
colnames(mean.df5) = c('MutationType', 'Group', 'Mean')

##################################################################################
##### Read COSMIC data
cosmic.df = as.data.frame(fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/COSMIC_v3.3.1_SBS_GRCh38.txt'))
colnames(cosmic.df)[1] = 'MutationType'

cosmic.df2 = melt(cosmic.df, id.vars='MutationType')
colnames(cosmic.df2) = c('MutationType', 'Group', 'Mean')

cosmic.df3 = cosmic.df2[cosmic.df2$Group %in% c('SBS26', 'SBS44'),]

##################################################################################
##### Merge and draw cosine correlation plot
merge.df = rbind(mean.df5, cosmic.df3)
merge.df2 = dcast(merge.df, MutationType~Group)
merge.df2 = merge.df2[,c('MutationType', group.level[-1], 'SBS26', 'SBS44')]
colnames(merge.df2) = gsub('\\-TK6\\S+', '', colnames(merge.df2)) 

merge.df3 = merge.df2
rownames(merge.df3) = merge.df3[,1]
merge.df3[,1] = NULL

merge.mat = as.matrix(merge.df3)
cormat = round(cosine(merge.mat), 2)

get_upper_tri <- function(cormat){
	cormat[lower.tri(cormat)]<- NA
	return(cormat)
}

cormat2 = cormat
upper_tri = get_upper_tri(cormat2)
melted_cormat = melt(upper_tri, na.rm=TRUE)

pdf(file=file.path(output.dir, 'figure2.c.cosine.corr.heatmap.sbs96.nodrug.mshs2.sub_msh2.cosmic.pdf'), width=10, height=10)
ggplot(data=melted_cormat, aes(Var2, Var1, fill=value)) +
	geom_tile(color='black') +
	geom_text(data=melted_cormat, aes(Var2, Var1, label=value)) +
	scale_fill_gradient2(low='blue', high='red', mid='white', midpoint=0.5,
		limit=c(0,1), space='Lab', name='Cosine\nCorrelation') +
	theme_light(base_size=23, base_family='sans') +
	theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=15, face='bold', color='black')) +
	theme(axis.text.y=element_text(size=15, face='bold', color='black')) +
	coord_fixed() +
	xlab(NULL) + ylab(NULL)
dev.off()













