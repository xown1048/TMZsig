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
#####	SBS Count and Proportion Box Plot
i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/figure4'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i$MutationType = substr(i$MutationType, 3, 5)
i = aggregate(.~MutationType, i, sum)

meta$Sample = gsub('\\.\\S+', '', meta$Sample)

df = melt(i, id.vars = 'MutationType')
colnames(df) = c('MutationType', 'Sample', 'Count')

df = merge(df, meta, by='Sample')

df$Gene = gsub('\\-\\S+', '', df$Group)
df$Drug = gsub('\\S+\\-', '', df$Group)
df$Dose = sapply(strsplit(as.character(unlist(df$Group)), '-'), '[[', 3)
df$Dose = gsub('nodrug', 0, df$Dose)
df$Tissue = sapply(strsplit(as.character(unlist(df$Group)), '-'), '[[', 2)

df2 = df
df2$Group = gsub('-cisplatin', '-Cisplatin', df2$Group)
df2$Group = gsub('\\-c\\S+', '', df2$Group)

df2 = df2[df2$Tissue == 'TK6',]
df2 = df2[!grepl('-p', df2$Group),]
df2 = df2[df2$Drug %in% c('nodrug', 'TMZ'),]

'%notin%' = Negate('%in%')
df3 = df2[df2$Gene %notin% c('WAPL', 'ESCO1', 'REV3'),]

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

############################################################################################################################################
#####	Draw plot
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
	#temp.df$Count2 = temp.df$Count2*500/temp.df$Dose2
	temp.df$Count2 = temp.df$Count2/temp.df$Dose2

	temp.df$Count2 = ifelse(temp.df$Count2 < 0, 0, temp.df$Count2)

	merge.df = rbind(merge.df, temp.df)
}

pathways = c(rep('WT', 1), rep('Apoptosis', 1), rep('DR', 3), rep('BER', 3), rep('NER', 2),
	rep('NER/BER', 1), rep('TMEJ', 1), rep('NHEJ', 1), rep('HR', 1), rep('PCNA', 1), rep('ICLR', 4),
	rep('TLS', 7), rep('MMR', 2), rep('MSH2 Double Knockouts', 10))

path = data.frame(Gene=gene.level, Path=pathways)
merge.df = merge(merge.df, path, 'Gene')

pathway.level = unique(pathways)

merge.df$Gene = factor(merge.df$Gene, levels=gene.level)
merge.df$Group = factor(merge.df$Group, levels=group.level)
merge.df$Path = factor(merge.df$Path, levels=pathway.level)

merge.df = merge.df[order(factor(merge.df$Group, levels=group.level)),]
sample.level = unique(merge.df$Sample)
merge.df$Sample = factor(merge.df$Sample, levels=sample.level)

sum.df = merge.df[,c('Sample', 'Count2')]
sum.df = aggregate(.~Sample, sum.df, sum)
sum.df = merge(sum.df, meta, by='Sample')
sum.df$Gene = gsub('\\-\\S+', '', sum.df$Group)

test.df = compare_means(Count2~Gene, method='t.test', data=sum.df)
test.df2 = as.data.frame(subset(test.df, group1=='WT'))
test.df3 = as.data.frame(subset(test.df, group2=='WT'))
test.df3 = test.df3[,c(1,3,2,4:8)]
colnames(test.df3) = colnames(test.df2)
test.df4 = rbind(test.df2, test.df3)

row_genes = gene.level[gene.level!='WT']
test.df4 = test.df4[match(row_genes, test.df4$group2),]

sum.df2 = sum.df[,c(4,2)]
mean.df = aggregate(.~Gene, sum.df2, mean)

wt.mean.df = mean.df[mean.df$Gene=='WT',]
mean.df = mean.df[mean.df$Gene!='WT',]

mean.df = mean.df[match(row_genes, mean.df$Gene),]

mean.df$Count2 = mean.df$Count2/wt.mean.df[1,2]

test.df5 = test.df4
test.df5$Ratio = mean.df$Count2
test.df5$p.format = as.numeric(test.df5$p.format)
test.df5$colors = ifelse(test.df5$p.format>0.1, 'white',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.05, '#e899af',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.01, '#dd6688',
					ifelse(test.df5$Ratio>1 & test.df5$p.format>0.001, '#d23260',
					ifelse(test.df5$Ratio>1 & test.df5$p.format<=0.001, '#C70039',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.05, '#DDF2FD',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.01, '#9BBEC8',
					ifelse(test.df5$Ratio<1 & test.df5$p.format>0.001, '#427D9D', '#164863'))))))))

#####	Down
#####	0.1, 0.05, 0.01, 0.001
#####	#DDF2FD, #9BBEC8, #427D9D, #164863

#####	UP
#####	0.1, 0.05, 0.01, 0.001
#####	#e899af, #dd6688, #d23260, #C70039

pval.strips.colors = c('white', test.df5$colors)

pval.facet.label = ifelse(test.df5$p.format<0.05, paste0(test.df5$p.signif, '\n', test.df5$group2), 
				ifelse(test.df5$p.format<0.1, paste0('#', '\n', test.df5$group2), test.df5$group2))

pval.facet.label = c('WT', pval.facet.label)

facet.name = c(pathway.level, gene.level)
temp.label = c(pathway.level, pval.facet.label)
temp.label = gsub('_', '\n', temp.label)

temp.label[2] = 'Apop.'
facet.label = temp.label
names(facet.label) = facet.name

#####	Cosmic signature colors
colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')
#colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#2A9D8F', '#264653')

pdf(file=file.path(output.dir, 'figure4.a.norm.sbs.6sub.bar.plot.pdf'), width=200, height=30)
ggplot(merge.df, aes(x=Sample, y=Count2)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	geom_hline(yintercept=wt.mean.df[1,2], linetype='solid', colour='#192655', linewidth=7) +
	scale_fill_manual(values=colors) +
	facet_nested(~Path+Gene, scales='free_x', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=100, base_family='sans') +
	ylab('Mutation Count per 1uM') + xlab(NULL) +
	theme(axis.text.x=element_blank(), axis.title=element_text(size=150, face='bold')) +
	theme(strip.text=element_text(colour='black', size=50)) +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'),
		axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60'), legend.key.size=unit(3, 'cm'))
dev.off()

p = ggplot(merge.df, aes(x=Sample, y=Count2)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	geom_hline(yintercept=wt.mean.df[1,2], linetype='dashed', colour='black', linewidth=5) +
	scale_fill_manual(values=colors) +
	facet_nested(~Path+Gene, scales='free_x', space='free_x', labeller=as_labeller(facet.label), strip=strip_nested(size='variable')) +
	theme_light(base_size=120, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_text(colour='black', face='bold', size=80)) +
	theme(axis.text.x=element_blank(), axis.title=element_blank()) +
	theme(legend.position='none') +
	theme(panel.spacing=unit(0,'lines'), strip.background=element_rect(colour='grey60', fill='white'),
		axis.ticks.x=element_blank(), panel.border=element_rect(color='grey60')) +
	guides(x='none')

strip.colors = c(rep(c('#9EDDFF', '#B0D9B1'), length(pathway.level)/2), pval.strips.colors)

pdf(file=file.path(output.dir, 'figure4.a.norm.sbs.6sub.bar.plot.pdf'), width=200, height=40)
g = ggplot_gtable(ggplot_build(p))
dev.off()

strips <- which(grepl('strip-', g$layout$name))

for (i in seq_along(strips)) {
	k = which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- strip.colors[i]
}

pdf(file=file.path(output.dir, 'figure4.a.norm.sbs.6sub.bar.plot.pdf'), width=200, height=40)
plot(g)
dev.off()

######################################################################
##### Count write table
output.dir2 = paste0(output.dir, '/table')
if(!dir.exists(output.dir2)){dir.create(output.dir2)}

merge.df2 = merge.df[,c('Count2','MutationType','Gene')]
mean.df = aggregate(Count2~MutationType+Gene, merge.df2, mean)
sd.df = aggregate(Count2~MutationType+Gene, merge.df2, sd)
colnames(mean.df)[3] = 'Mean'
colnames(sd.df)[3] = 'Sd'

merge.df3 = merge(mean.df, sd.df, by=c('Gene', 'MutationType'))
merge.df3$Sd = round(merge.df3$Sd, 1)
merge.df3$Mean = round(merge.df3$Mean, 1)

merge.df3 = merge.df3[order(factor(merge.df3$Gene, levels=group.level), merge.df3$MutationType),]

sum.df = aggregate(Count2~Sample+Gene, merge.df, sum)
mean.df = aggregate(Count2~Gene, sum.df, mean)
sd.df = aggregate(Count2~Gene, sum.df, sd)
colnames(mean.df)[2] = 'Mean'
colnames(sd.df)[2] = 'Sd'

merge.df4 = merge(mean.df, sd.df, by='Gene')
merge.df4$Sd = round(merge.df4$Sd, 1)
merge.df4$Mean = round(merge.df4$Mean, 1)

merge.df4 = merge.df4[order(factor(merge.df4$Gene, levels=group.level)),]

write.table(merge.df3, paste0(output.dir2, '/figure4.a.6sub.sbs.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)
write.table(merge.df4, paste0(output.dir2, '/figure4.a.total.sbs.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)

head(merge.df3)
head(merge.df4)

merge.df3[merge.df3$MutationType=='C>T' & merge.df3$Gene=='MSH2',]
merge.df3[merge.df3$MutationType=='C>T' & merge.df3$Gene=='EXO1',]

merge.df5 = merge.df3[merge.df3$Gene %in% c('WT', 'MSH2'),]
merge.df6 = dcast(merge.df5, Gene~MutationType, value.var='Mean')
merge.df7 = merge.df6[,-1]/rowSums(merge.df6[,-1])
rowSums(merge.df7)
merge.df8 = cbind(merge.df6[,1], merge.df7)
colnames(merge.df8)[1] = 'Gene'
#  Gene         C>A          C>G       C>T         T>A         T>C         T>G
#1   WT 0.074235808 0.0203784571 0.7714702 0.039301310 0.066957787 0.027656477
#2 MSH2 0.002257336 0.0004514673 0.9832957 0.003160271 0.009480813 0.001354402

mean(merge.df3[merge.df3$MutationType=='C>T' & merge.df3$Gene!='EXO1' & !grepl('MSH2', merge.df3$Gene), 'Mean'])
55.26154

merge.df3[merge.df3$MutationType=='C>T' & merge.df3$Gene %in% c('RAD18', 'ATAD5', 'MSH2', 'EXO1'),]
#     Gene MutationType  Mean   Sd
#15  ATAD5          C>T  75.5  8.5
#33   EXO1          C>T 176.4 16.0
#81   MSH2          C>T 217.8 20.6
#183 RAD18          C>T  82.3  5.6
#ATAD5
75.5/55.3
1.36528
#ATAD5
82.3/55.3
1.488246

(1.488246+1.36528)/2
1.426763