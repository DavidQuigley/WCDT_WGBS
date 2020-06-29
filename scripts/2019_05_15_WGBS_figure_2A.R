#Loop through HMRs per sample and assign them to groups
output <- data.frame()
samples2use <- sample_ids_wgbs
valleydf <- read.delim(file=fn_valley,sep='\t',header=F,stringsAsFactors=F)
colnames(valleydf) <- c('chr','start','end','sample')
valley <- makeGRangesFromDataFrame(valleydf,keep.extra.columns=T)

for(i in 1:length(samples2use)) {
    sample <- samples2use[i]
    valleyi <- valley[valley$sample==sample]
    
    valleywidth <- sum(width(valleyi))
    output[sample,'Sample'] <- sample
    output[sample,'width'] <- valleywidth
    output[sample,'count'] <- length(valleyi)
    
    print(paste(i, output[sample,'count']))
    
    in_h3k4me3 <- countOverlaps(tracks[['pca100_H3K4me3']],valleyi,ignore.strand=T)>0
    in_h3k27me3 <- countOverlaps(tracks[['pca100_H3K27me3']],valleyi,ignore.strand=T)>0
    in_cpg_island <- countOverlaps(tracks[['cpg_island']],valleyi,ignore.strand=T)>0
    
    count_h3k4me3_valley <- sum(in_h3k4me3)
    count_h3k27me3_valley <- sum(in_h3k27me3)
    count_h3k4me3_novalley <- sum(!in_h3k4me3)
    count_h3k27me3_novalley <- sum(!in_h3k27me3)
    
    tab <- rbind(c(count_h3k4me3_valley,count_h3k4me3_novalley),
                 c(count_h3k27me3_valley,count_h3k27me3_novalley))
    ft <- fisher.test(tab)
    
    output[sample,'OR'] <- log2(ft$estimate)
    output[sample,'p'] <- ft$p.value
}

output <- output[rownames(output) %in% sample_ids_wgbs,]
output$fdr <- p.adjust(output$p,method='fdr')
output$logfdr <- -log10(output$fdr)

lvl <- output[order(output$count),'Sample']
#lvl <- output[order(output$width),'Sample']
output$Sample <- factor(output$Sample,levels=lvl)

#Set up data frames for plotting
df1 <- cbind.data.frame(sample=output$Sample,count=output$OR,SS=output$fdr<=0.05)
df3 <- cbind.data.frame(sample=output$Sample,count=output$width/HG38_SIZE*100)
df4 <- cbind.data.frame(sample=output$Sample,count=output$count)

cols <- brewer.pal(n=4, name='Dark2')
plist <- list()
plist[[1]] <- ggplot(df1,aes(x=sample,y=count))+theme_classic()+geom_bar(stat='identity')+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ylab('Log2 Odds Ratio')+
    scale_fill_manual(values=c('darkgray',cols[1]))
plist[[3]] <- ggplot(df3,aes(x=sample,y=count))+theme_classic()+geom_bar(stat='identity',fill=cols[2])+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    ylab('DMV %Genome')
plist[[2]] <- ggplot(df4,aes(x=sample,y=count))+theme_classic()+geom_bar(stat='identity',fill=cols[3])+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    ylab('DMV Count')

ggarrange(plotlist=plist,ncol=1,nrow=length(plist),heights=c(2,1,1))
ggsave(fn_figure2a,width=15,height=12)



