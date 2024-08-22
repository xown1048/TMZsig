library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
options(bitmapType='cairo')
options(scipen=10000)

##########################################################################################################################################
##### SBS Signature Count and Proportion Bar Plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part10/output/sigProfilerExtractor.p0wt.with.all/output/signature/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities_refit.txt')

meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2.2/output/figure1'
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
prop.df = melt(prop.df, id.vars='Samples')
colnames(prop.df) = c('Sample', 'Signature', 'Proportion')

df2 = merge(df, prop.df, by=c('Sample', 'Signature'))
df2$Group = gsub('-cisplatin', '-Cisplatin', df2$Group)
df2$Group = gsub('\\-c\\S+', '', df2$Group)

df2 = df2[df2$Tissue=='TK6',]
gene.level = c('WT', 'XRCC1', 'POLHKI', 'REV1', 'REV7')
df2 = df2[df2$Gene %in% gene.level,]
df3 = df2[grepl('-p10', df2$Group),]

######################################################################
##### Group level
g.df = df3[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
tissue.level = c('TK6', 'HAP1')

gen.g.df = gen.g.df[order(factor(gen.g.df$Tissue, levels=tissue.level), factor(gen.g.df$Gene, levels=gene.level)),]
gen.group.level = gen.g.df$Group

group.level = unique(gen.group.level)

######################################################################
df3$Gene = factor(df3$Gene, levels=gene.level)
df3$Group = factor(df3$Group, levels=group.level)

sum.sig.df = df3[,c('Signature', 'Count')]
sum.sig.df = aggregate(Count~Signature, sum.sig.df, sum)
sum.sig.df = sum.sig.df[sum.sig.df$Count!=0,]

df3 = df3[df3$Signature %in% sum.sig.df$Signature,]

length(unique(df3$Signature))
#[1] 8

colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#2A9D8F', '#118AB2', '#264653', '#613F75')

pdf(file=file.path(output.dir, 'figure1.g.sig.sbs.10gen.p10.bar.count.pdf'), width=100, height=65)
ggplot(df3, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=Signature), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Gene, scales='free_x', space='free_x') +
	theme_light(base_size=200, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(legend.position='none') +
	theme(strip.text=element_text(colour='black', size=200, face='bold')) +
	theme(legend.key.size=unit(3, 'cm'), axis.text.x=element_blank(),
		strip.background=element_rect(colour='white', fill='white'),
		axis.title=element_text(size=250, face='bold'), axis.ticks.x=element_blank())
dev.off()

######################################################################
##### Count write table
output.dir2 = paste0(output.dir, '/table')
if(!dir.exists(output.dir2)){dir.create(output.dir2)}

table.df0 = df3
head(table.df0)

table.df = table.df0[,c('Count', 'Signature', 'Gene')]
mean.df = aggregate(Count~Signature+Gene, table.df, mean)
sd.df = aggregate(Count~Signature+Gene, table.df, sd)
colnames(mean.df)[3] = 'Mean'
colnames(sd.df)[3] = 'Sd'

table.df2 = merge(mean.df, sd.df, by=c('Gene', 'Signature'))
table.df2$Sd = round(table.df2$Sd, 1)
table.df2$Mean = round(table.df2$Mean, 1)

table.df2 = table.df2[order(factor(table.df2$Gene, levels=group.level), table.df2$Signature),]

sum.df = aggregate(Count~Sample+Gene, table.df0, sum)
mean.df = aggregate(Count~Gene, sum.df, mean)
sd.df = aggregate(Count~Gene, sum.df, sd)
colnames(mean.df)[2] = 'Mean'
colnames(sd.df)[2] = 'Sd'

table.df3 = merge(mean.df, sd.df, by='Gene')
table.df3$Sd = round(table.df3$Sd, 1)
table.df3$Mean = round(table.df3$Mean, 1)

table.df3 = table.df3[order(factor(table.df3$Gene, levels=group.level)),]

head(table.df2)
head(table.df3)

write.table(table.df2, paste0(output.dir2, '/figure1.h.5class.sv.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)
write.table(table.df3, paste0(output.dir2, '/figure1.h.total.sv.count.txt'), sep='\t', quote=F, col.names=T, row.names=F)

table.df3 = table.df2[table.df2$Gene %in% c('WT', 'REV7'),]
table.df4 = dcast(table.df3, Gene~Signature, value.var='Mean')
table.df5 = table.df4[,-1]/rowSums(table.df4[,-1])
rowSums(table.df5)

table.df4
table.df5
#  Gene SBS1  SBS5 SBS7a SBS8 SBS18 SBS26 SBS31 SBS40
#1   WT 10.7 297.3     0    0   0.0  39.7   0.0 637.7
#2 REV7 24.0 181.7    28    0 139.3   0.0  25.7  37.7
#        SBS1      SBS5      SBS7a SBS8     SBS18      SBS26      SBS31
#1 0.01085853 0.3017049 0.00000000    0 0.0000000 0.04028821 0.00000000
#2 0.05499542 0.4163611 0.06416132    0 0.3192026 0.00000000 0.05889093
#       SBS40
#1 0.64714837
#2 0.08638863

