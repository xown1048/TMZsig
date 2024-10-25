library(data.table)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(janitor)
library(ggbreak)
library(scales)
library(ggpubr)
library(bayestestR)

options(bitmapType='cairo')

rm(list=ls())

############################################################################################################################################
#####	SV bar count for all samples
vcfs = Sys.glob('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/delly.merge/no.read.filter/all/*.all.2pass.filtered.somatic.vcf')

meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure5'
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

temp.df = df2[df2$LEN2>200 & df2$ALT=='Insertions',]
temp.df2 = df2[df2$LEN2>200 & df2$ALT=='Deletions',]

'%notin%' = Negate('%in%')
temp.df3 = df2[df2$ALT %notin% c('Insertions', 'Deletions'),]

df3 = rbind(temp.df, temp.df2)
df3 = rbind(df3, temp.df3)

table(df3$ALT)
#Deletions Inversions Duplications Translocations
#607       101        66                  47

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
df7 = df7[!grepl('-p0', df7$Group),]

df7 = df7[df7$Tissue=='TK6',]
df7 = df7[df7$Drug %in% c('nodrug', 'TMZ'),]
df8 = df7[grepl('MSH2', df7$Group),]
df8 = df8[df8$Gene %in% c('MSH2', 'MSH2_MPG', 'MSH2_RAD18', 'MSH2_FANCD2', 'ATAD5_MSH2', 'MSH2_REV1', 'MSH2_REV7'),]

############################################################################################################################################
#####	Group level
g.df = df8[,c('Group', 'Gene', 'Drug', 'Dose', 'Tissue')]
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
df8$Prefix = ifelse(df8$Dose == 0, 0, 
				ifelse(grepl('nanoM', df8$Dose), 'nano',
				ifelse(grepl('microM', df8$Dose), 'micro', 'error')))

df8$Dose2 = ifelse(df8$Dose == 0, 0, 
				ifelse(grepl('nanoM', df8$Dose), gsub('nanoM','',df8$Dose),
				ifelse(grepl('microM', df8$Dose), gsub('microM','',df8$Dose), 'error')))

df8$Dose2 = as.numeric(df8$Dose2)
df8$Dose2 = ifelse(df8$Prefix == 0, 0, 
				ifelse(df8$Prefix == 'nano', df8$Dose2/1000,
				ifelse(df8$Prefix == 'micro', df8$Dose2, 'error')))

df9 = subset(df8, Drug=='nodrug')
df10 = subset(df8, Drug!='nodrug')
df10$Category = paste0(df10$Gene, '-', df10$Drug)
df10$Dose2 = as.numeric(df10$Dose2)
df9$Category = paste0(df9$Gene, '-', df9$Drug)
df9$Count2 = df9$Count

cates = sort(unique(df10$Category))

merge.df = data.frame()
for (cate in cates){
	temp.df = subset(df10, Category==cate)

	nodrug.temp.df = subset(df9, df9$Gene == unique(temp.df$Gene))
	nodrug.temp.df2 = nodrug.temp.df[,c(2,3)]
	nodrug.temp.df8 = aggregate(Count ~ MutationType, nodrug.temp.df2, mean)

	temp.df$Count2 = temp.df$Count - nodrug.temp.df8$Count[match(temp.df$MutationType, nodrug.temp.df8$MutationType)]
	#temp.df$Count2 = temp.df$Count2*500/temp.df$Dose2
	temp.df$Count2 = temp.df$Count2/temp.df$Dose2

	temp.df$Count2 = ifelse(temp.df$Count2 < 0, 0, temp.df$Count2)

	merge.df = rbind(merge.df, temp.df)
}

merge.df$Gene = factor(merge.df$Gene, levels=gene.level)
merge.df$Group = factor(merge.df$Group, levels=group.level)

merge.df = merge.df[order(factor(merge.df$Group, levels=group.level)),]
sample.level = unique(merge.df$Sample)
merge.df$Sample = factor(merge.df$Sample, levels=sample.level)

sum.df = merge.df[,c('Sample', 'Count2')]
sum.df = aggregate(.~Sample, sum.df, sum)
sum.df = merge(sum.df, meta, by='Sample')
sum.df$Gene = gsub('\\-\\S+', '', sum.df$Group)

test.df = compare_means(Count2~Gene, method='t.test', data=sum.df)
test.df2 = as.data.frame(subset(test.df, group1=='MSH2'))
test.df3 = as.data.frame(subset(test.df, group2=='MSH2'))
test.df3 = test.df3[,c(1,3,2,4:8)]
colnames(test.df3) = colnames(test.df2)
test.df4 = rbind(test.df2, test.df3)

gene.level = gene.level[gene.level %in% unique(df8$Gene)]

row_genes = gene.level[gene.level!='MSH2']
test.df4 = test.df4[match(row_genes, test.df4$group2),]

sum.df2 = sum.df[,c(4,2)]
mean.df = aggregate(.~Gene, sum.df2, mean)

wt.mean.df = mean.df[mean.df$Gene=='MSH2',]
mean.df = mean.df[mean.df$Gene!='MSH2',]

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

pval.facet.label = c('MSH2', pval.facet.label)


facet.name =  gene.level
temp.label = pval.facet.label
facet.label = gsub('_', '\n', temp.label)

names(facet.label) = facet.name

colors = c('#E76F51', '#E9C46A', '#6A994E', '#264653')
merge.df$MutationType = factor(merge.df$MutationType, levels=c('Deletions', 'Duplications', 'Inversions', 'Translocations'))

pdf(file=file.path(output.dir, 'figure5.b.no.read.filter.norm.500uM.msh2.sv.5class.barplot.pdf'), width=100, height=60)
ggplot(merge.df, aes(x=Sample, y=Count2)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Gene, scales='free_x', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=200, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(legend.position='none',
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.text=element_text(colour='black', size=120, face='bold.italic'),
		panel.spacing=unit(0,'lines'),
		strip.background=element_rect(colour='grey60', fill='white'),
		panel.border=element_rect(color='grey60'))
dev.off()

p = ggplot(merge.df, aes(x=Sample, y=Count2)) +
	geom_bar(aes(fill=MutationType), position='stack', stat='identity') +
	geom_hline(yintercept=wt.mean.df[1,2], linetype='dashed', colour='red', linewidth=5) +
	scale_fill_manual(values=colors) +
	facet_grid(~Gene, scales='free_x', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=200, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(legend.position='none',
		axis.text.x=element_blank(),
		axis.ticks.x=element_blank(),
		strip.text=element_text(colour='black', size=120, face='bold.italic'),
		panel.spacing=unit(0,'lines'),
		strip.background=element_rect(colour='grey60', fill='white'),
		panel.border=element_rect(color='grey60'))

strip.colors =  pval.strips.colors

pdf(file=file.path(output.dir, 'figure5.b.no.read.filter.norm.500uM.msh2.sv.5class.barplot.pdf'), width=100, height=60)
g = ggplot_gtable(ggplot_build(p))
dev.off()

strips <- which(grepl('strip-', g$layout$name))

for (i in seq_along(strips)) {
	k = which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
	g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- strip.colors[i]
}

pdf(file=file.path(output.dir, 'figure5.b.no.read.filter.norm.500uM.msh2.sv.5class.barplot.pdf'), width=100, height=60)
plot(g)
dev.off()
