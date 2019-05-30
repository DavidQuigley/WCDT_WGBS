#source('/Volumes/datasets_1/human_sequence_prostate_WGBS/reproduce/scripts/2019_05_15_prepare_environment.r')

# Plots Figure 1A

colstf <- c('CRPC_ARE','LOC_ARE','chip_FOXA1','chip_HOXB13','chip_ERG')
pcttotal=FALSE


output <- data.frame()
for( i in 1 : length(hmr) ) {
    print(i)
    sample <- toupper(names(hmr)[i])
    hmri <- hmr[[i]]
    
    totalhmrs <- length(hmri)
    output[sample,'Sample'] <- sample
    output[sample,'hmr'] <- totalhmrs
    
    divideby <- 1
    if(pcttotal) {
        divideby <- totalhmrs
    }
    in_tf <- rep(F,totalhmrs)
    for(j in 1:length(colstf)) {
        coltf <- colstf[j]
        in_tf <- in_tf | countOverlaps(hmri,tracks[[coltf]],ignore.strand=T)
    }
    in_genebody = countOverlaps(hmri,tracks[['genebody']],ignore.strand=T)>0
    in_promoter = countOverlaps(hmri,tracks[['promoter']],ignore.strand=T)>0
    in_h3k27ac = countOverlaps(hmri,tracks[['pca100_H3K27ac']],ignore.strand=T)>0
    in_h3k27me3 = countOverlaps(hmri,tracks[['pca100_H3K27me3']],ignore.strand=T)>0
    
    output[sample,'Promoter'] <- sum(in_promoter) / divideby
    output[sample,'Genebody'] <- sum(in_genebody & !in_promoter) / divideby
    output[sample,'TFBS_not_promoter'] <- sum(in_tf & !in_genebody & !in_promoter) / divideby
    output[sample,'H3K27ac_not_promoter_tfbs'] <- sum(in_h3k27ac & !in_tf & !in_genebody & !in_promoter) / divideby
    output[sample,'H3K27me3_not_others'] <- sum(in_h3k27me3 & !in_h3k27ac & !in_tf & !in_genebody & !in_promoter) / divideby
    
    promoter_oncogenes = countOverlaps(hmri,gr_oncogenes, ignore.strand=T)
    output[sample,'oncogene_promoter'] = sum(promoter_oncogenes>0)/divideby
    promoter_tsg = countOverlaps(hmri, gr_tsg, ignore.strand=T)
    output[sample,'TS_promoter'] = sum(promoter_tsg>0)/divideby
    
    in_cpg_island = countOverlaps(hmri,tracks[['cpg_island']],ignore.strand=T)>0
    in_cpg_shore = countOverlaps(hmri,tracks[['cpg_shore']],ignore.strand=T)>0
    in_cpg_shelf = countOverlaps(hmri,tracks[['cpg_shelf']],ignore.strand=T)>0
    
    output[sample,'CpG_island'] = sum(in_cpg_island)/divideby
    output[sample,'CpG_shore'] = sum(in_cpg_shore & !in_cpg_island)/divideby
    output[sample,'CpG_shelf'] = sum(in_cpg_shelf & !in_cpg_shore & !in_cpg_island)/divideby
}

#output <- output[rownames(tumorpurity)[tumorpurity[samples_wgbs,'Tumor.Purity.Comp']>80],]
output <- output[rownames(output) %in% sample_ids_wgbs,]

lvl <- output[order(output$hmr),'Sample']
output$Sample <- factor(output$Sample,levels=lvl)

if(pcttotal) {
    countsother <- 1
} else {
    countsother <- output$hmr
}
countsother <- output$hmr-(output$Promoter+output$Genebody+output$TFBS_not_promoter+
                               output$H3K27me3_not_others+output$H3K27ac_not_promoter_tfbs)

df1 <- rbind.data.frame(
    cbind.data.frame(sample=output$Sample,count=countsother,group='Other'),
    cbind.data.frame(sample=output$Sample,count=output$Promoter,group='Promoter'),
    cbind.data.frame(sample=output$Sample,count=output$Genebody,group='Genebody'),
    cbind.data.frame(sample=output$Sample,count=output$TFBS_not_promoter,group='TFBS'),
    cbind.data.frame(sample=output$Sample,count=output$H3K27ac_not_promoter_tfbs,group='H3K27ac'),
    cbind.data.frame(sample=output$Sample,count=output$H3K27me3_not_others,group='H3K27me3'))
df1$group <- factor(df1$group,levels=c('Other','H3K27me3','H3K27ac','TFBS','Genebody','Promoter'))

df3 <- rbind.data.frame(
    cbind.data.frame(sample=output$Sample,count=output$CpG_island,group='CpG_island'),
    cbind.data.frame(sample=output$Sample,count=output$CpG_shelf,group='CpG_shelf'),
    cbind.data.frame(sample=output$Sample,count=output$CpG_shore,group='CpG_shore'))
df3$group <- factor(df3$group,levels=c('CpG_shelf','CpG_shore','CpG_island'))

df5 <- cellsummary[lvl,]
df5$mutation_count <- df5$mutation_count / (HG38_SIZE / 1000000)
df5$sample <- factor(rownames(df5),levels=lvl)
df5$namecn <- '%Altered'
df5$namemut <- 'Mut/Mb'

df6slice <- manta_counts[lvl,]
df6slice$sample <- factor(rownames(df6slice),levels=lvl)
df6 <- rbind.data.frame(
    cbind.data.frame(sample=df6slice$sample,count=df6slice$BND,group='BND'),
    cbind.data.frame(sample=df6slice$sample,count=df6slice$DEL,group='DEL'),
    cbind.data.frame(sample=df6slice$sample,count=df6slice$DUP,group='DUP'),
    cbind.data.frame(sample=df6slice$sample,count=df6slice$INS,group='INS'),
    cbind.data.frame(sample=df6slice$sample,count=df6slice$INV,group='INV'))

cols <- brewer.pal(n=5, name='Set1')
plist <- list()
plist[[1]] <- ggplot(df1,aes(x=sample,y=count,fill=group))+theme_classic()+geom_bar(stat='identity')+
    theme(legend.position = "none", 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    ylab('HMR#')+
    scale_fill_manual(values=c('darkgray',cols))
plist[[2]] <- ggplot(df3,aes(x=sample,y=count,fill=group))+theme_classic()+geom_bar(stat='identity')+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")+
    scale_fill_brewer(palette='Dark2')+ylab('CpG#')
plist[[3]] <- ggplot(df5,aes(x=sample,y=percent.copy.altered,fill=namecn))+theme_classic()+geom_bar(stat='identity')+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")+
    scale_fill_manual(values='black')+ylab('% CN Altered')
plist[[4]] <- ggplot(df5,aes(x=sample,y=mutation_count,fill=namemut))+theme_classic()+geom_bar(stat='identity')+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position = "none")+
    scale_fill_manual(values='black')+ylab('Mut/Mb')+scale_y_log10() 

plist[[5]] <- ggplot(df6,aes(x=sample,y=count,fill=group))+theme_classic()+geom_bar(stat='identity')+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
          axis.ticks.x=element_blank(),
          legend.position = "none")+ 
    scale_fill_brewer(palette='Set2')+ylab('# SV')

ggarrange(plotlist=plist,
          ncol=1,
          nrow=length(plist),
          heights=c(3,1,1,1,2) ,
          align="v")

