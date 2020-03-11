window = 1000000
fn_DSS_adeno_benign
grdss  = dss_met_bp
samples2use  = setdiff(sample_ids_wgbs,samples_hypermut)

rowsannovarfull = annovar_full$Chr %in% chromosomes & annovar_full$sample %in% samples2use
mutgr = makeGRangesFromDataFrame(annovar_full[rowsannovarfull,])

rows2use = rows_manta_pass & data_manta_slim$sample_id %in% samples2use
dfmanta1 = data_manta_slim[rows2use,c('chr1','pos1','pos1')]
colnames(dfmanta1) = c('chr','start','end')
dfmanta2 = data_manta_slim[rows2use,c('chr2','pos2','pos2')]
colnames(dfmanta2) = c('chr','start','end')
rowsmissing = is.na(dfmanta2$start)
dfmanta2[rowsmissing,'chr'] = dfmanta1[rowsmissing,'chr']
dfmanta2[rowsmissing,'start'] = dfmanta1[rowsmissing,'start']
dfmanta2[rowsmissing,'end'] = dfmanta1[rowsmissing,'end']
grmanta = makeGRangesFromDataFrame(rbind.data.frame(dfmanta1,dfmanta2))

grdss = removetc(grdss)
mutgr = removetc(mutgr)
grmanta = removetc(grmanta)

df = data.frame()
for(i in 1:length(chrsizes)) {
    chr = names(chrsizes)[i]
    chrsize = chrsizes[i]
    numbins = ceiling(chrsize/window)
    dfchr = data.frame(chr=chr,start=(0:(numbins-1))*window)
    dfchr$end = dfchr$start+window-1
    df = rbind.data.frame(df,dfchr)
}

gr = makeGRangesFromDataFrame(df,keep.extra.columns=T)

exclude = getpctoverlap(gr,c(centromeres,telomeres))
rownum2exclude = exclude[exclude$x>0.5,'Category']
print(paste('exclude',length(rownum2exclude)))


df$mut = countOverlaps(gr,mutgr)/length(samples2use)
df$sv = countOverlaps(gr,grmanta)
df$cnbreaks = countOverlaps(gr,copycat_gr)
cntab = getpctoverlap(gr,copycat_gr)
df$cn = 0
df[cntab$Category,'cn'] = cntab$x
df$cn = (df$cn/100)/window

df$mutlog = log10(df$mut+1)
df$svlog = log10(df$sv+1)
genetab = getpctoverlap(gr,tracks[['genebody']])
df$genes = 0
df[genetab$Category,'genes'] = genetab$x

df$score = 0
overlaps = findOverlaps(gr,grdss)
#windows = aggregate(grdss[subjectHits(overlaps)]$diff.Methy, by=list(Category=queryHits(overlaps)), FUN=sum)
windows = aggregate(grdss[subjectHits(overlaps)]$areaStat, by=list(Category=queryHits(overlaps)), FUN=sum)
df[windows$Category,'score'] = windows$x

df = df[-rownum2exclude,]

print(cor.test(df$score,df$cnbreaks,method='spearman'))
print(cor.test(df$score,df$genes,method='spearman'))
print(cor.test(df$score,df$sv,method='spearman'))
corobj = cor.test(df$score,df$mutlog,method='spearman')
print(corobj)

ggplot(df, aes(x=mutlog, y=score)) + geom_point(alpha=0.4) + 
    theme_classic() +
    ggtitle(paste('p =',corobj$p.value,'rho =',corobj$estimate)) +
    xlab('Log10(Mut#)') + 
    ylab('Differential Hypo-methylation')+
    geom_hline(yintercept=0)
ggsave(fn_figure6b)

# trackgene = tracks[['genebody']][tracks[['genebody']]$gene %in% importantgenes]
# overlaps = findOverlaps(gr,trackgene)
# genes = aggregate(trackgene[subjectHits(overlaps)]$gene, by=list(Category=queryHits(overlaps)), FUN=paste, collapse='|')
# df[genes$Category,'genes'] = genes$x
# df = df[!is.na(df$score),]
# df = df[order(df$score),]
# write.table(df,'C:/Users/Admin/Desktop/dss.xls',sep='\t',col.names=T,row.names=F)
# 
