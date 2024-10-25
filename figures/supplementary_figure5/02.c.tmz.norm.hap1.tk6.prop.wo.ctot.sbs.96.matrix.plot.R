library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggsignif)
library(reshape2)
library(janitor)
library(patchwork)
library(ggh4x)

options(bitmapType='cairo')

######################################################################
### SBS 96 Matrix Bar Count Plot
rm(list=ls())
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output'
if(!dir.exists(output.dir)){dir.create(output.dir)}

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
df2 = df2[!grepl('-p0', df2$Group),]

df2 = df2[df2$Tissue %in% c('HAP1', 'TK6'),]
df2 = df2[df2$Drug %in% c('nodrug', 'TMZ'),]
df3 = df2[df2$Gene %in% c('WT', 'MGMT', 'MSH2'),]
df3 = df3[df3$Group %in% c('WT-TK6-nodrug', 'MGMT-HAP1-nodrug', 'WT-TK6-5microM-TMZ', 'MGMT-HAP1-5microM-TMZ'),]

group.level = c('WT-TK6-5microM-TMZ', 'MGMT-HAP1-5microM-TMZ')

##################################################################################
##### Count Plot
df3$Prefix = ifelse(df3$Dose == 0, 0, 
				ifelse(grepl('nanoM', df3$Dose), 'nano',
				ifelse(grepl('microM', df3$Dose), 'micro', 'error')))

df3$Dose2 = ifelse(df3$Dose == 0, 0, 
				ifelse(grepl('nanoM', df3$Dose), gsub('nanoM','',df3$Dose),
				ifelse(grepl('microM', df3$Dose), gsub('microM','',df3$Dose), 'error')))

df3$Dose2 = as.numeric(df3$Dose2)
df3$Dose2 = ifelse(df3$Prefix == 0, 0, 
				ifelse(df3$Prefix == 'nano', df3$Dose2/1000,
				ifelse(df3$Prefix == 'micro', df3$Dose2, 'error')))

df4 = subset(df3, Drug=='nodrug')
df5 = subset(df3, Drug!='nodrug')
df5$Category = paste0(df5$Gene, '-', df5$Drug)
df5$Dose2 = as.numeric(df5$Dose2)
df4$Category = paste0(df4$Gene, '-', df4$Drug)
df4$Count2 = df4$Count

cates = sort(unique(df5$Category))

merge.df = data.frame()
for (cate in cates){
	temp.df = subset(df5, Category==cate)

	nodrug.temp.df = subset(df4, df4$Gene == unique(temp.df$Gene))
	nodrug.temp.df2 = nodrug.temp.df[,c(2,3)]
	nodrug.temp.df3 = aggregate(Count ~ MutationType, nodrug.temp.df2, mean)

	temp.df$Count2 = temp.df$Count - nodrug.temp.df3$Count[match(temp.df$MutationType, nodrug.temp.df3$MutationType)]
	temp.df$Count2 = temp.df$Count2/temp.df$Dose2

	temp.df$Count2 = ifelse(temp.df$Count2 < 0, 0, temp.df$Count2)

	merge.df = rbind(merge.df, temp.df)
}

merge.df2 = merge.df
merge.df2 = merge.df2[,c('Group', 'MutationType', 'Count2')]

mean.df = dcast(merge.df2, MutationType~Group, mean)
head(mean.df)

mean.df2 = as.data.frame(t(t(mean.df[,-1])/colSums(mean.df[,-1])))
head(mean.df2)
colSums(mean.df2)

mean.df3 = cbind(mean.df[,1], mean.df2)
colnames(mean.df3)[1] = 'MutationType'
head(mean.df3)

mean.df4 = melt(mean.df3, id.vars='MutationType')
colnames(mean.df4) = c('MutationType', 'Group', 'Proportion')
head(mean.df4)

result.df = mean.df4

result.df$Substitution = substr(result.df$MutationType, 3, 5)
result.df$Side = paste0(substr(result.df$MutationType, 1, 1), substr(result.df$MutationType, 7, 7))

head(result.df)

result.df$Group = gsub('-5microM-TMZ', '', result.df$Group)

result.df$Substitution = factor(result.df$Substitution, levels=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))
group.level = c('WT-TK6', 'MGMT-HAP1')
result.df$Group = factor(result.df$Group, levels=group.level)

head(result.df)

result.df2 = result.df[result.df$Substitution!='C>T',]
facet.name = c(c('C>A', 'C>G', 'T>A', 'T>C', 'T>G'), group.level)
temp.label = gsub('_', '\n', group.level)
facet.label = c(c('C>A', 'C>G', 'T>A', 'T>C', 'T>G'), temp.label)
names(facet.label) = facet.name

colors = c('#01bbed', '#121212', '#ccc8c9', '#a2cd61', '#ecc6c4')
strip = strip_themed(background_x=elem_list_rect(fill=colors), text_x=elem_list_text(color=colors))

p1 = ggplot(result.df2, aes(x=MutationType, y=Proportion)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Group~Substitution, scales='free', space='free_x', strip=strip, labeller=as_labeller(facet.label)) + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_text(colour='white', size=60),
		strip.text.y=element_blank(),
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
		panel.spacing=unit(1.5,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.3, 1, 0.01, 1, 'cm')) +
	xlab(NULL) + ylab(NULL)

result.df2 = result.df[result.df$Substitution=='C>T',]
facet.name = c('C>T', group.level)
temp.label = gsub('_', '\n', group.level)
facet.label = c('C>T', temp.label)
names(facet.label) = facet.name

colors = '#db2d2a'
strip = strip_themed(background_x=elem_list_rect(fill=colors), text_x=elem_list_text(color=colors),
					background_y=elem_list_rect(fill=rep('white', length(group.level))), text_y=elem_list_text(color=rep('black', length(group.level))))

p2 = ggplot(result.df2, aes(x=MutationType, y=Proportion)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Group~Substitution, scales='free_x', space='free_x', strip=strip, labeller=as_labeller(facet.label)) + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_text(colour='white', size=60),
		strip.text.y=element_text(colour='black', face='bold', size=150),
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
		panel.spacing=unit(1.5,'lines'),
		panel.border=element_rect(color='grey60'),
		plot.margin=margin(0.3, 1, 0.01, 1, 'cm')) +
	xlab(NULL) + ylab(NULL)

pdf(file=file.path(output.dir, '02.review2.16_tmz.norm.hap1.tk6.prop.wo.ctot.sbs.96.matrix.plot.pdf'), width=150, height=60)
p1 + p2 + plot_layout(widths=unit(c(5, 1), rep('null', 2)))
dev.off()

facet.name = c(c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'), group.level)
temp.label = gsub('_', '\n', group.level)
facet.label = c(c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'), temp.label)
names(facet.label) = facet.name

colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')
strip = strip_themed(background_x=elem_list_rect(fill=colors), text_x=elem_list_text(color=colors),
					background_y=elem_list_rect(fill=rep('white', length(group.level))), text_y=elem_list_text(color=rep('black', length(group.level))))

pdf(file=file.path(output.dir, '02.review2.16_tmz.norm.hap1.tk6.prop.wo.ctot.sbs.96.matrix.plot.pdf'), width=150, height=60)
ggplot(result.df, aes(x=MutationType, y=Proportion)) +
	geom_bar(aes(fill=Substitution), stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid2(Group~Substitution, scales='free', space='free_x', strip=strip, labeller=as_labeller(facet.label)) + 
	theme_light(base_size=100, base_family='sans') +
	theme(legend.position='none',
		strip.text.x=element_text(colour='white', size=60),
		strip.text.y=element_text(colour='black', face='bold', size=150),
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
		panel.spacing=unit(1.5,'lines'),
		panel.border=element_rect(color='grey60')) +
	xlab(NULL) + ylab(NULL)
dev.off()


library(lsa)
cosine(mean.df3[,2], mean.df3[,3])
#0.95