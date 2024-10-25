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

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output'
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
df2 = df2[df2$Sample!='D6_270k',]

df3 = df3[df3$Tissue=='TK6',]
df3 = df3[df3$Drug %in% c('nodrug','TMZ'),]
df3 = df3[df3$Gene %in% c('MSH2', 'WT'),]
df3 = df3[df3$Dose %in% c('0', '5microM', '500microM'),]

df3 = df3[df3$Group!='WT-TK6-nodrug',]

df3 = rbind(df3, df2[df2$Group=='WT-TK6-p10',])

##################################################################################
##### Draw cosine correlation plot
meta$Sample = gsub('\\.\\S+', '', meta$Sample)
meta2 = meta[meta$Sample %in% unique(df3$Sample),]
meta2 = meta2[order(-meta2$Group),]

meta2$cate = gsub('-TK6-', '-', meta2$Group)
meta2$cate = gsub('5microM-', '', meta2$cate)
meta2$cate = gsub('500microM-', '', meta2$cate)
meta2$cate = gsub('\\-p10\\S+', '-nodrug', meta2$cate)

meta2$cate = paste0(meta2$cate, '_', rep(c(1,2,3),4))

meta2 = meta2[,c(2,3)]

df3 = df3[,c(1,2,3)]
df4 = merge(df3, meta2, 'Sample')
df4 = df4[,c(4,2,3)]

df5 = dcast(df4, MutationType~cate)
df5 = df5[,meta2$cate]

mat = as.matrix(df5)
cormat = round(cosine(mat), 2)

get_upper_tri <- function(cormat){
	cormat[lower.tri(cormat)]<- NA
	return(cormat)
}

cormat2 = cormat
upper_tri = get_upper_tri(cormat2)
melted_cormat = melt(upper_tri, na.rm=TRUE)

pdf(file=file.path(output.dir, '01.triplicate.cosine.corr.v2.pdf'), width=10, height=10)
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




