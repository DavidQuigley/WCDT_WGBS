n_hmr = rep(0, 100)
for(i in 1:length(sample_ids_wgbs)){
    n_hmr[i] = length( hmr[[ sample_ids_wgbs[i]  ]] )
}

m = match.idx( sample_ids_wgbs, samples_cmp )
is_cmp = rep(FALSE, 100)
is_cmp[m$idx.A]=TRUE

gg = rep("non-CMP", 100)
gg[is_cmp] = "CMP"
df = data.frame( n_hmr/1000, groups=gg )
colors = rep('black',dim(df)[1])
colors[ df$groups=='CMP' ] = 'cornflowerblue'
df$groups = factor( df$groups, levels=c('non-CMP','CMP') )
colname = colnames(df)[1]
test = wilcox.test(df[,colname]~df$groups)
plist = ggplot(df,aes_string(x='groups', y=colname))+
    xlab('') + ylab( '' ) +
    geom_boxplot( outlier.shape=NA) +
    geom_jitter( size=1, position=position_jitter(height=0,width=.2), color=colors) +
    theme_classic() +
    scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20), expand = c(0, 0))

#ggarrange(plotlist=plist,ncol=length(plist),nrow=1)
ggsave(fn_figure1c,width=2,height=3,useDingbats=FALSE)

wilcox.test(n_hmr~is_cmp)
# W = 1686, p-value = 5.756e-12