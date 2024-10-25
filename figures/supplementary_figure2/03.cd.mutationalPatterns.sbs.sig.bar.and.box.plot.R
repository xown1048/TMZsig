#####	https://bioconductor.org/packages/release/bioc/vignettes/MutationalPatterns/inst/doc/Introduction_to_MutationalPatterns.html

library(MutationalPatterns)
library(data.table)

rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/part11/output/sigProfilerExtractor/output/SBS/sigProfilerExtractor.SBS96.all')
meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output/mutationalPatterns'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i = as.data.frame(i)

vcf_files <- list.files(system.file("extdata", package = "MutationalPatterns"),
	pattern = "sample.vcf", full.names = TRUE
)

sample_names <- c(
	"colon1", "colon2", "colon3",
	"intestine1", "intestine2", "intestine3",
	"liver1", "liver2", "liver3"
)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

row.order = rownames(mut_mat)

meta2 = meta
meta2 = meta2[grepl('-p10', meta2$Group),]
samples = meta2$Sample

i2 = i[,which(colnames(i) %in% c('MutationType', samples))]
i3 = i2[order(factor(i2$MutationType, levels=row.order)),]
rownames(i3) = i3[,1]
i3$MutationType = NULL

mut.mat = as.matrix(i3)

#####	COSMIC (v3.2) (Alexandrov et al. 2020)
signatures = get_known_signatures()

selected.signatures = signatures[,c('SBS1', 'SBS5', 'SBS7a', 'SBS8', 'SBS18', 'SBS26', 'SBS31', 'SBS40')]
best_subset_refit <- fit_to_signatures_strict(mut.mat, selected.signatures, max_delta = 0.002)

fit_best_res_strict <- best_subset_refit$fit_res

pdf(file=file.path(output.dir, '01.review1.1_mutationalPatterns.best.fit_res_strict.sbs.sig.bar.plot.pdf'), width=10, height=5)
plot_contribution(fit_best_res_strict$contribution,
	coord_flip = FALSE,
	mode = "absolute"
)
dev.off()

res.best.fit_res = plot_contribution(fit_best_res_strict$contribution,
	coord_flip = FALSE,
	mode = "absolute"
)

res.best.fit_res.df = as.data.frame(res.best.fit_res$data)

output.dir2 = paste0(output.dir, '/01.review1.1_mutationalPatterns.best.sbs.sig.fit_res.txt')
write.table(x=res.best.fit_res.df, file=output.dir2, quote=F, sep='\t', row.names=F, col.names=T)

library(data.table)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(bayestestR)
library(ggh4x)

options(bitmapType='cairo')
options(scipen=10000)

##########################################################################################################################################
#####	SBS Signature Count and Proportion Box Plot
rm(list=ls())

i = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output/mutationalPatterns/01.review1.1_mutationalPatterns.best.sbs.sig.fit_res.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output/mutationalPatterns'
if(!dir.exists(output.dir)){dir.create(output.dir)}

i = as.data.frame(i)
df = i
colnames(df) = c('Signature', 'Sample', 'Count')

meta = fread('/BiO/Research/UNIST-Toni-Sig-2020-1026/20230914.sample.list.txt')
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

df2 = df2[df2$Tissue=='TK6',]
gene.level = c('WT', 'XRCC1', 'POLHKI', 'REV1', 'REV7')
df2 = df2[df2$Gene %in% gene.level,]
df3 = df2[grepl('-p10', df2$Group),]

############################################################################################################################################
#####	Group level
g.df = df3[,4:8]
g.df = unique(g.df)

gen.g.df = g.df[grepl('-p', g.df$Group),]
tissue.level = c('TK6', 'HAP1')

gen.g.df = gen.g.df[order(factor(gen.g.df$Tissue, levels=tissue.level), factor(gen.g.df$Gene, levels=gene.level)),]
gen.group.level = gen.g.df$Group

group.level = unique(gen.group.level)

############################################################################################################################################
#####	Draw bar plot
df3$Gene = factor(df3$Gene, levels=gene.level)
df3$Group = factor(df3$Group, levels=group.level)
df3$Signature = factor(df3$Signature, levels=c('SBS1', 'SBS5', 'SBS7a', 'SBS8', 'SBS18', 'SBS26', 'SBS31', 'SBS40'))

colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#2A9D8F', '#118AB2', '#264653', '#613F75')

pdf(file=file.path(output.dir, 'supp.fig2c.mutationalPatterns.p10.sbs.sig.bar.plot.pdf'), width=80, height=65)
ggplot(df3, aes(x=Sample, y=Count)) +
	geom_bar(aes(fill=Signature), position='stack', stat='identity') +
	scale_fill_manual(values=colors) +
	facet_grid(~Gene, scales='free_x', space='free_x') +
	theme_light(base_size=200, base_family='sans') +
	ylab(NULL) + xlab(NULL) +
	theme(strip.text=element_text(colour='black', size=150, face='bold.italic')) +
	theme(legend.key.size=unit(3, 'cm'), axis.text.x=element_blank(),
		strip.background=element_rect(colour='white', fill='white'),
		axis.title=element_text(size=250, face='bold'), axis.ticks.x=element_blank())
dev.off()

############################################################################################################################################
#####	Draw box plot
my_comparisons = combn(levels(df3$Gene), 2, simplify=F)
t = compare_means(Count~Gene, comparisons=my_comparisons,
		group.by='Signature', method='t.test', data=df3)
t2 = subset(t, p<= 0.1)
t.df = as.data.frame(t2)

merge.t = data.frame()
for (c in unique(t2$Signature)) {
	tmp = subset(df3, Signature == c)
	max = max(as.numeric(tmp$Count))
	tmp.t = subset(t.df, Signature == c)

	if(nrow(tmp.t) == 1) {
		tmp.t$'y.position' = max + max/10
	} else {
	interval = max/nrow(tmp.t)
	max2 = max + interval
	to = max2 + ((interval*(nrow(tmp.t)-1))*0.9)
	tmp.t$'y.position' = seq(from=max2, to=to, by=interval*0.9)
	}

	merge.t = rbind(merge.t, tmp.t)
}

merge.t$p = as.numeric(merge.t$p)
merge.t$p.format = as.numeric(merge.t$p.format)
merge.t$symbol = ifelse(merge.t$p <= 0.05, merge.t$p.signif, '#')

merge.t = tibble::as_tibble(merge.t)

colors = c('#E76F51', '#F4A261', '#E9C46A', '#6A994E', '#2A9D8F', '#118AB2', '#264653', '#613F75')

strip = strip_themed(background_x=elem_list_rect(fill=colors))

pdf(file=file.path(output.dir, 'supp.fig2d.mutationalPatterns.p10.sbs.sig.box.plot.pdf'), width=60, height=40)
ggplot(df3, aes(x=Gene, y=Count)) +
	geom_boxplot(aes(fill=Gene), colour='#264653', fill='white', size=2) +
	xlab(NULL) + ylab(NULL) +
	facet_wrap2(~Signature, nrow=2, scales='free', strip=strip) + 
	stat_pvalue_manual(merge.t, size=25, label='symbol') +
	theme_light(base_size=100, base_family='sans') +
	theme(legend.key.size=unit(3,'cm'), legend.text=element_text(size=50), axis.title=element_text(size=100),
		legend.title=element_text(size=70), strip.text=element_text(colour='black', size=80, face='bold'), 
		axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face='bold.italic', colour='black'),
		strip.background=element_rect(colour='grey', fill='white'))
dev.off()
