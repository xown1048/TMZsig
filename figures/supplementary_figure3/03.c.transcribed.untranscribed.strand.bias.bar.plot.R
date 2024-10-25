library(data.table)
library(reshape2)
library(ggplot2)
library(ggpattern)

rm(list=ls())

#####	ggpattern
#####	https://stackoverflow.com/questions/62393159/how-can-i-add-hatches-stripes-or-another-pattern-or-texture-to-a-barplot-in-ggp

######################################################################
#####	Transcribed strand bias (transcribed vs untranscribed)bar plot
inputs = Sys.glob('/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output/sigProfilerTopography/part5/*/output/result/data/transcription_strand_bias/text_files/Type_Transcribed_Versus_Untranscribed.txt')

output.dir = '/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output/sigProfilerTopography/part5/merge'
if(!dir.exists(output.dir)){dir.create(output.dir)}

df = data.frame()
for (temp in inputs){
	i = fread(temp)
	i2 = i[i$type %in% c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'),]
	i2 = i2[,c('type', 'Transcribed_real_count', 'UnTranscribed_real_count', 'Transcribed_mean_sims_count', 'UnTranscribed_mean_sims_count', 'transcribed_versus_untranscribed_p_value')]

	i2$real_ratio = i2$Transcribed_real_count/i2$UnTranscribed_real_count
	i2$sims_ratio = i2$Transcribed_mean_sims_count/i2$UnTranscribed_mean_sims_count
	i2$odds_ratio = i2$real_ratio/i2$sims_ratio

	group = temp
	group = gsub('/BiO/Research/UNIST-Toni-Sig-2020-1026/peer.review/nar/part1/output/sigProfilerTopography/part5/', '', group)
	group = gsub('/output/result/data/transcription_strand_bias/text_files/Type_Transcribed_Versus_Untranscribed.txt', '', group)
	i2$group = group
	df = rbind(df, i2)
}

write.table(df, paste0(output.dir, '/transcribed.untranscribed.txt'), sep='\t', quote=F, col.names=T, row.names=F)

df2 = df
df2 = melt(df2, id.vars=c('group', 'type'))
df3 = df2[df2$variable == 'odds_ratio',]

df3$group = gsub('-TK6-', '-', df3$group)
df3$group = gsub('\\-\\S+microM\\-', '-', df3$group)
df3$group = gsub('ATAD5_MSH2', 'MSH2_ATAD5', df3$group)
df3$group = gsub('FANCD2_MSH2', 'MSH2_FANCD2', df3$group)

df3$group = factor(df3$group, levels=c('MSH2-nodrug', 'MSH2-TMZ',
									'MSH2_ALKBH3-nodrug', 'MSH2_ALKBH3-TMZ',
									'MSH2_ALKBH2-nodrug', 'MSH2_ALKBH2-TMZ',
									'MSH2_XRCC1-nodrug', 'MSH2_XRCC1-TMZ',
									'MSH2_MPG-nodrug', 'MSH2_MPG-TMZ',
									'MSH2_RAD18-nodrug', 'MSH2_RAD18-TMZ',
									'MSH2_REV1-nodrug', 'MSH2_REV1-TMZ',
									'MSH2_REV7-nodrug', 'MSH2_REV7-TMZ',
									'MSH2_FANCD2-nodrug', 'MSH2_FANCD2-TMZ',
									'MSH2_p53-nodrug', 'MSH2_p53-TMZ',
									'MSH2_ATAD5-nodrug', 'MSH2_ATAD5-TMZ'))


facet.name = c('MSH2-nodrug', 'MSH2-TMZ',
'MSH2_ALKBH3-nodrug', 'MSH2_ALKBH3-TMZ',
'MSH2_ALKBH2-nodrug', 'MSH2_ALKBH2-TMZ',
'MSH2_XRCC1-nodrug', 'MSH2_XRCC1-TMZ',
'MSH2_MPG-nodrug', 'MSH2_MPG-TMZ',
'MSH2_RAD18-nodrug', 'MSH2_RAD18-TMZ',
'MSH2_REV1-nodrug', 'MSH2_REV1-TMZ',
'MSH2_REV7-nodrug', 'MSH2_REV7-TMZ',
'MSH2_FANCD2-nodrug', 'MSH2_FANCD2-TMZ',
'MSH2_p53-nodrug', 'MSH2_p53-TMZ',
'MSH2_ATAD5-nodrug', 'MSH2_ATAD5-TMZ')
facet.label = toupper(facet.name)
facet.label = gsub('\\-', ' ', facet.label)
facet.label = gsub('NODRUG', 'No Drug', facet.label)
names(facet.label) = facet.name

colors = c('#01bbed', '#121212', '#db2d2a', '#ccc8c9', '#a2cd61', '#ecc6c4')

pval.df = df2[df2$variable == 'transcribed_versus_untranscribed_p_value',]

pval.df$group = gsub('-TK6-', '-', pval.df$group)
pval.df$group = gsub('\\-\\S+microM\\-', '-', pval.df$group)
pval.df$group = gsub('ATAD5_MSH2', 'MSH2_ATAD5', pval.df$group)
pval.df$group = gsub('FANCD2_MSH2', 'MSH2_FANCD2', pval.df$group)

pval.df$symbol = ifelse(pval.df$value > 0.1, 'NS',
					ifelse(pval.df$value > 0.05, '#',
					ifelse(pval.df$value > 0.01, '*',
					ifelse(pval.df$value > 0.001, '**', '***'))))

pval.df2 = pval.df[pval.df$value <= '0.1',]

pval.df3 = merge(pval.df2, df3, by=c('group', 'type'))
pval.df3$size = ifelse(pval.df3$value.x > 0.05, 5, 10)

pval.df3$group = factor(pval.df3$group, levels=c('MSH2-nodrug', 'MSH2-TMZ',
									'MSH2_ALKBH3-nodrug', 'MSH2_ALKBH3-TMZ',
									'MSH2_ALKBH2-nodrug', 'MSH2_ALKBH2-TMZ',
									'MSH2_XRCC1-nodrug', 'MSH2_XRCC1-TMZ',
									'MSH2_MPG-nodrug', 'MSH2_MPG-TMZ',
									'MSH2_RAD18-nodrug', 'MSH2_RAD18-TMZ',
									'MSH2_REV1-nodrug', 'MSH2_REV1-TMZ',
									'MSH2_REV7-nodrug', 'MSH2_REV7-TMZ',
									'MSH2_FANCD2-nodrug', 'MSH2_FANCD2-TMZ',
									'MSH2_p53-nodrug', 'MSH2_p53-TMZ',
									'MSH2_ATAD5-nodrug', 'MSH2_ATAD5-TMZ'))
colnames(pval.df3)[7] = 'value'

output.dir2 = paste0(output.dir, '/sbs.transcribed.untranscribed.pval.txt')
write.table(x=pval.df3, file=output.dir2, quote=F, sep='\t', row.names=F, col.names=T)

pdf(file=file.path(output.dir, 'sbs.odds.ratio.transcribed.untranscribed.pdf'), width=20, height=15)
ggplot(df3, aes(x=value, y=type, color=type)) +
	geom_point(size=4, shape=19) +
	geom_text(data=pval.df3, aes(x=value, y=type, label=symbol, size=size), color='black', nudge_y=0.2) +
	scale_size(range=c(6, 10), guide='none') +
	geom_vline(xintercept=1, linetype='longdash') +
	scale_colour_manual(values=colors) +
	facet_wrap(group~., ncol=4, labeller=as_labeller(facet.label)) +
	xlab(NULL) + ylab(NULL) +
	theme_light(base_size=10, base_family='sans') +
	theme(strip.background=element_rect(colour='white', fill='white')) +
	theme(strip.text=element_text(colour='black', size=15, face='bold')) +
	theme(axis.text=element_text(colour='black')) +
	guides(color=guide_legend(title='Mutation\nType', size=15)) +
	scale_x_continuous(trans='log2')
dev.off()

df3 = df2[df2$variable %in% c('Transcribed_real_count', 'UnTranscribed_real_count', 'Transcribed_mean_sims_count', 'UnTranscribed_mean_sims_count'),]

df3$group = gsub('-TK6-', '-', df3$group)
df3$group = gsub('\\-\\S+microM\\-', '-', df3$group)
df3$group = gsub('ATAD5_MSH2', 'MSH2_ATAD5', df3$group)
df3$group = gsub('FANCD2_MSH2', 'MSH2_FANCD2', df3$group)

df4 = df3
df4$sim = ifelse(grepl('sim', df4$variable), 'Simulated', 'Real') 

df4$group = factor(df4$group, levels=c('MSH2-nodrug', 'MSH2-TMZ',
									'MSH2_ALKBH3-nodrug', 'MSH2_ALKBH3-TMZ',
									'MSH2_ALKBH2-nodrug', 'MSH2_ALKBH2-TMZ',
									'MSH2_XRCC1-nodrug', 'MSH2_XRCC1-TMZ',
									'MSH2_MPG-nodrug', 'MSH2_MPG-TMZ',
									'MSH2_RAD18-nodrug', 'MSH2_RAD18-TMZ',
									'MSH2_REV1-nodrug', 'MSH2_REV1-TMZ',
									'MSH2_REV7-nodrug', 'MSH2_REV7-TMZ',
									'MSH2_FANCD2-nodrug', 'MSH2_FANCD2-TMZ',
									'MSH2_p53-nodrug', 'MSH2_p53-TMZ',
									'MSH2_ATAD5-nodrug', 'MSH2_ATAD5-TMZ'))

facet.name = c('MSH2-nodrug', 'MSH2-TMZ',
'MSH2_ALKBH3-nodrug', 'MSH2_ALKBH3-TMZ',
'MSH2_ALKBH2-nodrug', 'MSH2_ALKBH2-TMZ',
'MSH2_XRCC1-nodrug', 'MSH2_XRCC1-TMZ',
'MSH2_MPG-nodrug', 'MSH2_MPG-TMZ',
'MSH2_RAD18-nodrug', 'MSH2_RAD18-TMZ',
'MSH2_REV1-nodrug', 'MSH2_REV1-TMZ',
'MSH2_REV7-nodrug', 'MSH2_REV7-TMZ',
'MSH2_FANCD2-nodrug', 'MSH2_FANCD2-TMZ',
'MSH2_p53-nodrug', 'MSH2_p53-TMZ',
'MSH2_ATAD5-nodrug', 'MSH2_ATAD5-TMZ')
facet.label = toupper(facet.name)
facet.label = gsub('\\-', ' ', facet.label)
facet.label = gsub('NODRUG', 'No Drug', facet.label)
names(facet.label) = facet.name

colors = c('#01bbed', '#db2d2a', '#01bbed', '#db2d2a')

df5 = aggregate(value ~ type + group, data=df3, mean)

df5$value = ifelse(grepl('TMZ', df5$group), df5$value + df5$value*0.01, df5$value + df5$value*0.1)
colnames(df5)[3] = 'yvalue'

pval.df3 = merge(pval.df2, df5, by=c('group', 'type'))
pval.df3$size = ifelse(pval.df3$value > 0.05, 5, 10)
pval.df3$sim = 'Simulated'
pval.df3$variable = 'Transcribed_mean_sims_count'

pval.df3$group = factor(pval.df3$group, levels=c('MSH2-nodrug', 'MSH2-TMZ',
									'MSH2_ALKBH3-nodrug', 'MSH2_ALKBH3-TMZ',
									'MSH2_ALKBH2-nodrug', 'MSH2_ALKBH2-TMZ',
									'MSH2_XRCC1-nodrug', 'MSH2_XRCC1-TMZ',
									'MSH2_MPG-nodrug', 'MSH2_MPG-TMZ',
									'MSH2_RAD18-nodrug', 'MSH2_RAD18-TMZ',
									'MSH2_REV1-nodrug', 'MSH2_REV1-TMZ',
									'MSH2_REV7-nodrug', 'MSH2_REV7-TMZ',
									'MSH2_FANCD2-nodrug', 'MSH2_FANCD2-TMZ',
									'MSH2_p53-nodrug', 'MSH2_p53-TMZ',
									'MSH2_ATAD5-nodrug', 'MSH2_ATAD5-TMZ'))

pdf(file=file.path(output.dir, 'sbs.count.transcribed.untranscribed.pdf'), width=20, height=15)
ggplot(df4, aes(x=type, y=value, fill=variable, pattern=sim)) +
	geom_col_pattern(position='dodge', colour='black', pattern_fill='black', pattern_density=0.1, pattern_key_scale_factor=0.1) +
	geom_text(data=pval.df3, aes(x=type, y=yvalue, label=symbol, size=size), color='black', nudge_y=0.2) +
	scale_size(range=c(6, 10), guide='none') +
	scale_fill_manual(values=colors) +
	scale_pattern_manual(values=c(Simulated='stripe', Real='none')) +
	facet_wrap(group~., ncol=4, labeller=as_labeller(facet.label), scales='free_y') +
	xlab(NULL) + ylab(NULL) +
	theme_light(base_size=10, base_family='sans') +
	theme(strip.background=element_rect(colour='white', fill='white')) +
	theme(strip.text=element_text(colour='black', size=15, face='bold')) +
	theme(axis.text=element_text(colour='black')) +
	guides(pattern=guide_legend(override.aes=list(fill='white')),
			fill=guide_legend(override.aes=list(pattern='none'))) +
	labs(pattern='Data\nType', fill='Mutation\nType')
dev.off()

df[df$group == 'atad5.tmz' & df$type == 'T>C',]
df[df$group == 'atad5.nodrug',]

df[df$transcribed_versus_untranscribed_p_value <= '0.1',]

   type Transcribed_real_count UnTranscribed_real_count
1:  C>T                  44052                    41932
2:  C>A                     24                        9
   Transcribed_mean_sims_count UnTranscribed_mean_sims_count
1:                     41287.7                       40486.9
2:                        22.2                          24.9
   transcribed_versus_untranscribed_p_value real_ratio sims_ratio odds_ratio
1:                              0.002376202   1.050558  1.0197792   1.030182
2:                              0.037470668   2.666667  0.8915663   2.990991
            group
1: msh2.atad5.tmz
2:      wt.nodrug
