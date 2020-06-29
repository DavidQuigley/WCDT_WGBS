df_hmr_recurrent_full = read.delim(fn_hmr_recurrent, sep='\t', header=T, check.names=F)
df_hmr_recurrent = df_hmr_recurrent_full
samples_plot = intersect(sample_ids_wgbs,colnames(df_hmr_recurrent))
mat_hmr_recurrent_full = as.matrix(df_hmr_recurrent_full[,samples_plot])
mat_hmr_recurrent_full[is.nan(mat_hmr_recurrent_full)] = NA


stopifnot( sum(sample_ids_wgbs==dimnames(mat_hmr_recurrent_full)[[2]]) == 100 )
m = match.idx( dimnames(mat_hmr_recurrent_full)[[2]], samples_cmp )
is_cmp = rep(FALSE, 100)
is_cmp[m$idx.A]=TRUE

meth_medians = colMedians( mat_hmr_recurrent_full, na.rm=TRUE)

fn_figure1d_test='~/mnt/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/figure_1d_test.pdf'
postscript(fn_figure1d_test, height=3,width=2)
par(mar=c(3,3,2,1))
boxplot(meth_medians~is_cmp, names=c("", ""), las=1, 
        ylab="", ylim=c(0,60))
lines(c(1,2), c(94,94))
lines(c(1,1), c(94,93))
lines(c(2,2), c(94,93))
text( 1.5, 97, "*", font=2, cex=2)
dev.off()

wilcox.test(meth_medians~is_cmp)
# W = 41, p-value = 1.09e-11

gg = rep("non-CMP", 100)
gg[is_cmp] = "CMP"
df = data.frame( median=meth_medians, groups=gg )
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
    scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10), expand = c(0, 0))

#ggarrange(plotlist=plist,ncol=length(plist),nrow=1)
ggsave(fn_figure1d,width=2,height=3 ,useDingbats=F)
