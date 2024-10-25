library(data.table)
library(ggplot2)
library(ggsignif)
library(bayestestR)

rm(list=ls())

is = Sys.glob('/BiO/Research/UNIST-Toni-Sig-2020-1026/part7/output/mantis/*/*.msi.mantis.txt.status')
is = c(is, Sys.glob('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/mantis/*/*.msi.mantis.txt.status'))

meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/figures/part2/output/figure2'
if(!dir.exists(output.dir)){dir.create(output.dir)}

meta$Sample = gsub('\\.\\S+', '', meta$Sample)
colnames(meta) = c('Group', 'Sample')

df = data.frame()
for (i in is) {
	sample = basename(dirname(i))
	temp.df = as.data.frame(fread(i))

	temp.df = temp.df[,c(2,4)]
	temp.df = cbind(matrix=c('dif', 'euc', 'cos'), Sample=c(sample, sample, sample), temp.df)
	colnames(temp.df) = c('Matrix', 'Sample', 'Value', 'Status')

	df = rbind(df, temp.df)
}

df2 = merge(df, meta, by='Sample')

df2$Gene = gsub('\\-\\S+', '', df2$Group)
df2$Drug = gsub('\\S+\\-', '', df2$Group)
df2$Dose = sapply(strsplit(as.character(unlist(df2$Group)), '-'), '[[', 3)
df2$Dose = gsub('nodrug', 0, df2$Dose)
df2$Tissue = sapply(strsplit(as.character(unlist(df2$Group)), '-'), '[[', 2)

df2$Group = gsub('-cisplatin', '-Cisplatin', df2$Group)
df2$Group = gsub('\\-c\\S+', '', df2$Group)

df3 = df2[df2$Tissue == 'TK6',]

df3 = df3[!grepl('-p', df3$Group),]
df3 = df3[df3$Drug == 'nodrug',]

df3 = df3[df3$Gene %in% c('WT', 'EXO1') | grepl('MSH2', df3$Group),]

############################################################################################################################################
#####	Group level
g.df = df3[,c('Group', 'Gene', 'Drug', 'Dose', 'Tissue')]
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
df3$Gene = factor(df3$Gene, levels=gene.level)
df3$Group = factor(df3$Group, levels=group.level)

wt.df = subset(df3, Group=='MSH2-TK6-nodrug')

sub.groups = group.level[-which(group.level=='MSH2-TK6-nodrug')]
result.df = data.frame()
for (gro in sub.groups){
	temp.df = subset(df3, Group==gro)

	mats = unique(df3$Matrix)
	for(mat in mats){
		temp.df3 = subset(temp.df, Matrix==mat)
		temp.wt.df = subset(wt.df, Matrix==mat)

		counts = temp.df3$Value
		wt.counts = temp.wt.df$Value

		if(all(counts == '1') & all(wt.counts == '1')){
			p.val = 1
		} else if(length(unique(counts)) == 1 & length(unique(wt.counts)) == 1){
			if(mean(counts) == mean(wt.counts)) {
				p.val = 1
			} else{
				p.val = 0
			}
		} else{ 
			p.val = t.test(wt.counts, counts)$'p.value'
		}

		min.ci = ci(counts, ci=0.95, method='ETI')$CI_low
		max.ci = ci(counts, ci=0.95, method='ETI')$CI_high
		amean = mean(counts)

		afc = ifelse(mean(wt.counts)==0 & mean(counts)==0, 1, 
			ifelse(mean(wt.counts)==0, mean(counts)/mean(wt.counts+1), mean(counts)/mean(wt.counts)))

		temp.df3 = data.frame(Group=gro, Matrix=mat, Mean=amean, Minci=min.ci, Maxci=max.ci, Pval=p.val, Fc=afc)

		result.df = rbind(result.df, temp.df3)
	}
}

wt.mean.df = aggregate(Value ~ Matrix+Group, wt.df, function(x){mean(x)})
wt.ci.high.df = aggregate(Value ~ Matrix+Group, wt.df, function(x){ci(x, ci=0.95, method='ETI')$CI_high})
wt.ci.low.df = aggregate(Value ~ Matrix+Group, wt.df, function(x){ci(x, ci=0.95, method='ETI')$CI_low})

colnames(wt.mean.df)[3] = 'Mean'
colnames(wt.ci.high.df)[3] = 'Maxci'
colnames(wt.ci.low.df)[3] = 'Minci'

wt.result.df = merge(wt.mean.df, wt.ci.high.df, by=c('Matrix', 'Group'))
wt.result.df = merge(wt.result.df, wt.ci.low.df, by=c('Matrix', 'Group'))

wt.result.df = wt.result.df[,c('Group', 'Matrix', 'Mean', 'Minci', 'Maxci')]
wt.result.df = cbind(wt.result.df, Pval=1)
wt.result.df = cbind(wt.result.df, Fc=1)

result.df3 = rbind(wt.result.df, result.df)

s.group.level = c('MSH2-TK6-nodrug', sub.groups)

result.df3$Group = factor(result.df3$Group, levels=s.group.level)
result.df3$Gene = gsub('\\-\\S+', '', result.df3$Group)

s.gene.level = gsub('\\-\\S+', '', s.group.level)

result.df3$Gene = factor(result.df3$Gene, levels=s.gene.level)

facet.name = c(s.gene.level, 'cos', 'dif', 'euc')
facet.label = gsub('_', '\n', facet.name)

names(facet.label) = facet.name


threshold.df = data.frame(Matrix=c('cos', 'dif', 'euc'), Mean=c(0.07, 0.4, 0.187))


pdf(file=file.path(output.dir, 'figure2.b.wt.msh2.t.test.ci.dot.msi.pdf'), width=40, height=40)
ggplot(result.df3, aes(x=Group, y=Mean)) + 
	geom_point(size=20, shape=19, aes(alpha=ifelse(Pval<=0.1, 1, 0), colour=ifelse(Fc==1, 'black', ifelse(Fc>1, 'darkred', 'darkblue')))) +
	geom_errorbar(width=0.6, size=2, aes(ymin=Minci, ymax=Maxci, alpha=ifelse(Pval<=0.1, 1, 0), colour=ifelse(Fc==1, 'black', ifelse(Fc>1, 'darkred', 'darkblue')))) +
	scale_alpha(range=c(0.3, 1)) +
	scale_color_identity() +
	geom_hline(data=threshold.df, aes(yintercept=Mean), linetype='dashed', colour='red', size=2) +
	facet_grid(Matrix~Gene, scales='free', space='free_x', labeller=as_labeller(facet.label)) +
	theme_light(base_size=60, base_family='sans') +
	ylab(NULL) +
	xlab(NULL) +
	theme(axis.title=element_text(size=60), strip.background=element_rect(colour='white', fill='white'), 
		legend.position='none', strip.text=element_text(colour='black'), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
	theme(strip.text.x=element_text(size=40, face='bold.italic'), strip.text.y=element_text(size=80))
dev.off()

