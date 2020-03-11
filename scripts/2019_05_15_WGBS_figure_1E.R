n_hmr = rep(0, 100)
for(i in 1:length(sample_ids_wgbs)){
    n_hmr[i] = length( hmr[[ sample_ids_wgbs[i]  ]] )
}

stopifnot( sum(sample_ids_wgbs==dimnames(mat_hmr_recurrent_full)[[2]]) == 100 )
m = match.idx( dimnames(mat_hmr_recurrent_full)[[2]], samples_cmp )
is_cmp = rep(FALSE, 100)
is_cmp[m$idx.A]=TRUE

gg = rep("Non-CMP", 100)
gg[is_cmp] = "CMP"
df = data.frame( n_hmr/1000, groups=gg )
colors = rep('black',dim(df)[1])
colors[ df$groups=='CMP' ] = 'cornflowerblue'
df$groups = factor( df$groups, levels=c('Non-CMP','CMP') )
colname = colnames(df)[1]
test = wilcox.test(df[,colname]~df$groups)
plist[[1]] = ggplot(df,aes_string(x='groups', y=colname))+
    xlab('') + ylab( '' ) +
    geom_boxplot( outlier.shape=NA) +
    geom_jitter( size=1, position=position_jitter(height=0,width=.2), color=colors) +
    theme_classic()

ggarrange(plotlist=plist,ncol=length(plist),nrow=1)
ggsave(fn_figure1e,width=2,height=3,useDingbats=F)


wilcox.test(n_hmr~is_cmp)
# W = 1686, p-value = 5.756e-12