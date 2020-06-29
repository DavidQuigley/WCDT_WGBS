genename <- 'MYC' #plots figure 5D
hmrsegs <- tracks[['HMRseg']]

rowensembl <- ensembl2sym$name==genename
ensemblid <- rownames(ensembl2sym)[rowensembl]
stopifnot(sum(rowensembl)==1)
chr <- ensembl2sym[rowensembl,'chr']
gene_start <- ensembl2sym[rowensembl,'start']
gene_end <- ensembl2sym[rowensembl,'end']

samples2use <- sample_ids_wgbs

start <- gene_start-10000
end <- gene_end+160000
gr <- GRanges(chr,IRanges(start,end))
expr <- as.numeric(tpm[ensemblid,samples2use])
covs <- data.frame(expr)
rownames(covs) <- samples2use

hmrgene <- hmrsegs[countOverlaps(hmrsegs,gr)>0]
covs$add2 <- cnregion(chr,gene_start,gene_end)[samples2use]
add2name <- 'CN'
print(cor.test(as.numeric(covs$add2),as.numeric(expr),method='spearman'))
         
#Make linear models and compare
for(i in 1:length(hmrgene)) {
    print(i)
    hmr_start <- start(hmrgene)[i]
    hmr_end <- end(hmrgene)[i]
    methylmean <- meanmethyl(ensemblid,dir_gene_output_mcrpc,hmr_start,hmr_end)[samples2use]
    covs[,i+2+as.numeric(genename=='ERG')] <- methylmean
}

lm1 <- lm(expr~add2,covs)

print(summary(lm1))
lm2 <- lm(expr~.,covs[,-1])
pval <- anova(lm1,lm2)[2,'Pr(>F)']
print(summary(lm1)$adj.r.squared)
print(summary(lm2)$adj.r.squared)

pred1 <- predict(lm1,interval='confidence')
pred2 <- predict(lm2,interval='confidence')

covs$lm1 <- pred1[,1]
covs$lm1up <- pred1[,2]
covs$lm1lw <- pred1[,3]
covs$lm2 <- pred2[,1]
covs$lm2up <- pred2[,2]
covs$lm2lw <- pred2[,3]
covs$x <- as.numeric(factor(samples2use,levels=samples2use[order(expr)]))

#Make the plots
dfplot <- rbind.data.frame(
    cbind.data.frame(x=covs$x,y=covs$expr,group=paste('Actual',genename,'TPM')),
    cbind.data.frame(x=covs$x,y=covs$lm1,group=paste('LM predictions with',add2name)),
    cbind.data.frame(x=covs$x,y=covs$lm2,group='LM predictions with Methylation')
)

p <- ggplot(dfplot,aes(x=x,y=y,color=group))+
    geom_smooth(method='loess',se=F)+
    geom_point()+
    xlab(paste('ANOVA P:',signif(pval,4)))+ylab(paste(genename,'TPM'))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme_classic()+
    scale_color_manual(values=c("darkred", "darkgreen", "darkblue"))
ggsave(fn_figure5d,p,width=7,height=7,useDingbats=F)




