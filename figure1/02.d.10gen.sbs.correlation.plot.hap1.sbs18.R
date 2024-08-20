library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)

options(bitmapType='cairo')

######################################################################
### SBS 96 Matrix Bar Count Plot
rm(list=ls())
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure1'

df = melt(i, id.vars='MutationType')
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

df3 = df2[grepl('-p10', df2$Group) | df2$Group=='WT-HAP1-nodrug',]

##################################################################################
mean.df = aggregate(Count ~ MutationType+Group, df3, mean)

cosmic.df = as.data.frame(fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/COSMIC_v3.3.1_SBS_GRCh38.txt'))
colnames(cosmic.df)[1] = 'MutationType'

cosmic.df2 = melt(cosmic.df, id.vars='MutationType')
colnames(cosmic.df2) = c('MutationType', 'Group', 'Mean')

cosmic.df3 = cosmic.df2[cosmic.df2$Group %in% c('SBS5', 'SBS18', 'SBS40'),]

##################################################################################
##### To proportion
mean.df2 = dcast(mean.df, MutationType~Group)
mean.df3 = as.data.frame(t(t(mean.df2[,-1])/colSums(mean.df2[,-1])))

colSums(mean.df3)

mean.df4 = cbind(mean.df2[,1], mean.df3)
colnames(mean.df4)[1] = 'MutationType'

mean.df5 = melt(mean.df4, id.vars='MutationType')
colnames(mean.df5) = c('MutationType', 'Group', 'Mean')

merge.df = rbind(mean.df5, cosmic.df3)

########################################################################################################################################################################
##### Signature Cosine Correlation
library(lsa)

merge.df2 = dcast(merge.df, MutationType~Group)

merge.df3 = merge.df2
rownames(merge.df3) = merge.df3[,1]
merge.df3[,1] = NULL

merge.df3 = merge.df3[,c('WT-TK6-p10', 'WT-HAP1-nodrug', 'XRCC1-TK6-p10', 'POLHKI-TK6-p10', 'REV1-TK6-p10', 'REV7-TK6-p10', 'SBS5', 'SBS18', 'SBS40')]
colnames(merge.df3) = c('WT-TK6', 'WT-HAP1', 'XRCC1', 'POLHKI', 'REV1', 'REV7', 'SBS5', 'SBS18', 'SBS40')

merge.mat = as.matrix(merge.df3)
cormat = round(cosine(merge.mat), 2)

get_upper_tri <- function(cormat){
	cormat[lower.tri(cormat)]<- NA
	return(cormat)
}

cormat2 = cormat
upper_tri = get_upper_tri(cormat2)
melted_cormat = melt(upper_tri, na.rm=TRUE)

pdf(file=file.path(output.dir, 'figure1.d.cosine.corr.heatmap.p10.cosmic.sbs.spectrum.pdf'), width=10, height=10)
ggplot(data=melted_cormat, aes(Var2, Var1, fill=value)) +
	geom_tile(color='black') +
	geom_text(data=melted_cormat, aes(Var2, Var1, label=value), size=7) +
	scale_fill_gradient2(low='blue', high='red', mid='white', midpoint=0.5,
		limit=c(0,1), space='Lab', name='Cosine\nCorrelation') +
	theme_light(base_size=23, base_family='sans') +
	theme(axis.text.x=element_text(angle=90, vjust=1, hjust=1, color='black', face='bold.italic')) +
	theme(axis.text.y=element_text(color='black', face='bold.italic')) +
	coord_fixed() +
	xlab(NULL) + ylab(NULL)
dev.off()

