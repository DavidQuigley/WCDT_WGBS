df_mcrpc = read.delim(fn_recurrent_hmr_mcrpc,sep='\t',header=TRUE,check.names=FALSE)
df_upitt = read.delim(fn_recurrent_hmr_upitt,sep='\t',header=TRUE,check.names=FALSE)

stopifnot(identical(df_mcrpc$chr,df_upitt$chr))
stopifnot(identical(df_mcrpc$start,df_upitt$start))
stopifnot(identical(df_mcrpc$end,df_upitt$end))

df = cbind.data.frame(df_mcrpc,df_upitt[,-1:-3])

samplesnocmp = setdiff(sample_ids_wgbs,samples_cmp)
samples2use = toupper( 
               c(samples_upitt_benign_prostate,
                samples_upitt_localized_tumor,
                samples_cmp,
                samplesnocmp) )
group = rep('mCRPC',length(samples2use))
names(group) = samples2use
group[samples_tscnc] = 'tSCNC'
group[toupper(samples_upitt_benign_prostate)] = 'Benign'
group[toupper(samples_upitt_localized_tumor)] = 'Primary'
group[samples_cmp] = 'CMP'
group[samples_braf] = 'BRAF'
group[samples_idh1] = 'IDH1'
group[samples_tet2] = 'TET2'
group[samples_dnmt3b] = 'DNMT3B'

M = t( df[ , match.idx( samples2use, colnames(df))$idx.B ] )
M = M[ ,colSums(is.na(M))==0]
seed = 2
numpca = 10
perp = 30
iter = 1000

titlestr = paste('seed:',seed,', PCA:',numpca,', perplexity:',perp,', iter:',iter)

set.seed(seed)
tsne_obj = Rtsne(M,normalize=FALSE,verbose=TRUE,num_threads=7,
                  initial_dims=numpca,
                  perplexity=perp,
                  theta=0,#0 for exact, 1 for approximations
                  max_iter=iter)
tsnedf = as.data.frame(tsne_obj$Y)
colnames(tsnedf) = c('t1','t2')

labeltxt = rep('',dim(tsnedf)[1])
loctxt = metastasis_locations[rownames(M)]
names(loctxt) = rownames(M)
loctxt[toupper(samples_upitt_benign_prostate)] = 'Prostate'
loctxt[toupper(samples_upitt_localized_tumor)] = 'Prostate'

p = ggplot(data=tsnedf, aes(x=t1,y=t2,color=group,label=labeltxt)) +
    geom_point(aes(shape=factor(loctxt)),alpha=1,size=2) +
    geom_text(colour="black",alpha=1,size=2,nudge_y=0.1) +
    theme_classic() + 
    scale_shape_discrete(name='Site') +
    scale_color_manual(values=c(brewer.pal(8,'Dark2'),'black'),name='') +
    ggtitle(paste('seed:',seed,', PCA:',numpca,', perplexity:',perp,', iter:',iter))

ggsave(fn_figure_tsne_a, p, width=9, height=7)

