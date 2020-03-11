#Plots figure 1C, 1D
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
df3 <- rbind.data.frame(
    cbind.data.frame(sample=output$Sample,count=output$CpG_island,group='CpG_island'),
    cbind.data.frame(sample=output$Sample,count=output$CpG_shelf,group='CpG_shelf'),
    cbind.data.frame(sample=output$Sample,count=output$CpG_shore,group='CpG_shore'))
df3$group <- factor(df3$group,levels=c('CpG_shelf','CpG_shore','CpG_island'))

ssi = df3
ssi$CMP = rep("Non-CMP", dim(ssi)[1] )
ssi$CMP[ match.idx( samples_cmp, ssi$sample, allow.multiple.B = TRUE)$idx.B ] = "CMP"

n_in_shelf = c()
n_outside_shelf = c()
for(i in 1:length(hmr)) {
    hmri = hmr[[i]]
    n_in_shelf[i] = sum( countOverlaps(hmri,tracks[['cpg_shelf']])>0 )
    n_outside_shelf[i] = sum( countOverlaps(hmri,tracks[['cpg_shelf']])==0 )
}

df=data.frame( sample_id = names(hmr), 
               n_in_shelf, 
               n_outside_shelf, 
               sum=n_in_shelf+n_outside_shelf)
df = df[match.idx( df$sample_id, sample_ids_wgbs )$idx.A,]
df$CMP = rep("Non-CMP", 100)
df$CMP[ match.idx( df$sample_id, samples_cmp)$idx.A] = "CMP"
df$groups = factor( df$CMP, levels=c('Non-CMP','CMP') )
df$colors = rep("black", dim(df)[1])
df$colors[df$CMP=="CMP"] = "cornflowerblue"
df = df[order(df$sum),]

ssi$count = ssi$count / 1000
df$n_outside_shelf = df$n_outside_shelf / 1000

shoreshelf = ssi[ssi$group=="CpG_shelf" | ssi$group=="CpG_shore",]
shoreshelf$groups = factor( shoreshelf$CMP, levels=c('Non-CMP','CMP') )
shoreshelf$colors = rep("black", dim(shoreshelf)[1])
shoreshelf$colors[shoreshelf$CMP=="CMP"] = "cornflowerblue"

plist[[1]] = ggplot(shoreshelf,aes_string(x='groups', y="count") )+
    xlab('') + ylab( '' ) +
    geom_boxplot( outlier.shape=NA) +
    geom_jitter( size=1, position=position_jitter(height=0,width=.2), colour=shoreshelf$colors) +
    theme_classic()
ggarrange(plotlist=plist,ncol=length(plist),nrow=1)
ggsave(fn_figure1c,width=2, height=3,useDingbats=FALSE)


plist[[1]] = ggplot(df,aes_string(x='groups', y="n_outside_shelf") )+
    xlab('') + ylab( '' ) +
    geom_boxplot( outlier.shape=NA) +
    geom_jitter( size=1, position=position_jitter(height=0,width=.2), color=df$colors) +
    theme_classic()
ggarrange(plotlist=plist,ncol=length(plist),nrow=1)
ggsave(fn_figure1d,width=2, height=3,useDingbats=FALSE)



wilcox.test( ssi$count[ ssi$group=="CpG_island" & !ssi$CMP], 
             ssi$count[ssi$group=="CpG_island" & ssi$CMP] )
# W = 770, p-value = 0.4666
wilcox.test(ssi$count[ (ssi$group=="CpG_shore"|ssi$group=="CpG_shelf") & !ssi$CMP], 
            ssi$count[( ssi$group=="CpG_shore" | ssi$group=="CpG_shelf") & ssi$CMP])
# W = 6154.5, p-value = 9.928e-16
wilcox.test( c(df$n_outside_shelf[!df$CMP]  ), 
             c(df$n_outside_shelf[df$CMP]  ))
# W = 1707, p-value = 1.661e-12



