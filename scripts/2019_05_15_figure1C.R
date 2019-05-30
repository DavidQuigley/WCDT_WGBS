#source('/Volumes/datasets_1/human_sequence_prostate_WGBS/reproduce/scripts/2019_05_15_prepare_environment.r')

#Plots figure 1C

samples = samples_cmp
groups = rep('Non-CMP',length(sample_ids_wgbs))
names(groups) = sample_ids_wgbs
groups[sample_ids_wgbs %in% samples] = 'CMP'

hmrcpgi = c()
hmrnoncpgi = c()
hmrtsg = c()
hmrog = c()
cpgrhmrcount = sum( countOverlaps( tracks[['HMRseg']], tracks[['cpg_shelf']] )>0 )
noncpgrhmrcount = sum( countOverlaps( tracks[['HMRseg']], tracks[['cpg_shelf']] )==0 )
stopifnot(cpgrhmrcount+noncpgrhmrcount==length(tracks[['HMRseg']]))

tsgrhmrcount = sum(countOverlaps(tracks[['HMRseg']],gr_tsg)>0)
ogcpgrhmrcount = sum(countOverlaps(tracks[['HMRseg']],gr_oncogenes)==0)

for(i in 1:length(hmr)) {
    print(names(hmr)[i])
    hmri = hmr[[i]]
    hmri = hmri[countOverlaps(hmri,tracks[['HMRseg']])>0]
    hmrcpgi[i] = (cpgrhmrcount-sum(countOverlaps(hmri,tracks[['cpg_shelf']])>0))/cpgrhmrcount
    hmrnoncpgi[i] = (noncpgrhmrcount-sum(countOverlaps(hmri,tracks[['cpg_shelf']])==0))/noncpgrhmrcount
    hmrtsg[i] = tsgrhmrcount-sum(countOverlaps(hmri,gr_tsg)>0)
    hmrog[i] = ogcpgrhmrcount-sum(countOverlaps(hmri,gr_oncogenes)>0)
}
names(hmrcpgi) = toupper(names(hmr))
names(hmrnoncpgi) = toupper(names(hmr))
names(hmrtsg) = toupper(names(hmr))
names(hmrog) = toupper(names(hmr))

df = data.frame(groups)
df$hmrcpgi = hmrcpgi[sample_ids_wgbs]
df$hmrnoncpgi = hmrnoncpgi[sample_ids_wgbs]

df$cn = cellsummary[sample_ids_wgbs,'percent.copy.altered']
df$hmrtsg = hmrtsg[sample_ids_wgbs]
df$hmrog = hmrog[sample_ids_wgbs]

df$mut = cellsummary[sample_ids_wgbs,'mutation_count']
df$ploidy = cellsummary[sample_ids_wgbs,'estimated.ploidy']
df$purhist = tumorpurity[sample_ids_wgbs,'Tumor.Purity.Histo']
df$purcomp = tumorpurity[sample_ids_wgbs,'Tumor.Purity.Comp']

#Use this line of code to make stratified boxplots
#df = df[df$purcomp>50,]

plist = list()
namelist = c('CpG island / shore / shelf methylation','CpG open-sea methylation',
              '% Copy Altered',
              'HMR TSG promoter','HMR Oncogene promoter',
              'Mutation count','Ploidy','Purity_hist','Purity_comp')

df = df[,c(1,2,3,9,10)]
namelist = namelist[c(1,2,8,9)]

colors = rep('black',dim(df)[1])
colors[ df$groups=='CMP' ] = 'mediumturquoise'
df$groups = factor( df$groups, levels=c('Non-CMP','CMP') )

for(i in 2:length(df)) {
    colname = colnames(df)[i]
    test = wilcox.test(df[,colname]~df$groups)
    plist[[i-1]] = ggplot(df,aes_string(x='groups', y=colname))+
        xlab('')+ylab(namelist[i-1])+ggtitle(paste('Wilcox p:',signif(test$p.value,4)))+
        geom_boxplot(outlier.shape=NA)+
        geom_jitter(size=1,position=position_jitter(height=0,width=.2),color=colors)+
        theme_classic()
}
ggarrange(plotlist=plist,ncol=length(plist),nrow=1)
#ggsave(fileout,width=8,height=6,useDingbats=F)



#sink('C:/Users/Admin/Desktop/box_cmp_sink.txt')
#print(summary(lm(df$hmrcpgi~df$purcomp+df$groups)))
#print(summary(lm(df$hmrnoncpgi~df$purcomp+df$groups)))
#print(summary(lm(df$hmrcpgi~df$purhist+df$groups)))
#print(summary(lm(df$hmrnoncpgi~df$purhist+df$groups)))
#sink()






