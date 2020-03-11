source('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/scripts/2019_05_15_prepare_environment.r')

plot_expr_meth = function( symbol, meth_str, mc, bed ){
  ensembl_id = dimnames(ensembl2sym)[[1]][ ensembl2sym$name==symbol]
    idx_tpm = match.idx( sample_ids_wgbs, dimnames(tpm)[[2]] )$idx.B
  xp = log10(1+as.numeric( tpm[ensembl_id,idx_tpm] ) )
  meth = mc[which(bed$str==meth_str),]
  plot( xp, meth, xlab=paste(symbol, "expression"), ylab=paste("% meth ",meth_str) )
  idx_purity = match.idx( sample_ids_wgbs, dimnames(tumorpurity)[[1]] )$idx.B
  purity=tumorpurity$Tumor.Purity.Comp[idx_purity]
  plot( purity, meth )
  cor_xp = round( cor.test( xp, meth, method="spearman")$estimate, 2) 
  cor_pure = round( cor.test( purity, meth, method="spearman")$estimate, 2)
  print( paste( "cor XP:", cor_xp, "cor purity:", cor_pure ))
}

hclust2 = function(x, method="average") {
    hclust(x, method=method)
}
dist2 = function(x) {
    as.dist(1-cor(t(x), method="spearman"))
}

generate_rsidebar=function( M ){
    cgi = rep('white', dim(M)[1] )
    t(matrix(cgi))   
}

generate_csidebar = function( M, samples_plot ){
    toplot = data.frame(samples_plot)
    toplot$BRCA2 = allele_effect('BRCA2')$alleles[samples_plot,'bi']
    toplot$ERG = allele_effect('ERG')$alleles[samples_plot,'activating_sv']
    toplot$ETV1 = allele_effect('ETV1')$alleles[samples_plot,'activating_sv']
    toplot$ETV4 = allele_effect('ETV4')$alleles[samples_plot,'activating_sv']
    toplot$ETV5 = allele_effect('ETV5')$alleles[samples_plot,'activating_sv']
    toplot$CHD1 = allele_effect('CHD1')$alleles[samples_plot,'n_alleles_inactivated']==2
    toplot$SPOP = allele_effect('SPOP')$alleles[samples_plot,'inactivating_missense']
    toplot$RB1 = allele_effect('RB1')$alleles[samples_plot,'n_alleles_inactivated']==2
    toplot$PTEN = allele_effect('PTEN')$alleles[samples_plot,'n_alleles_inactivated']==2
    toplot$MYC = allele_effect('MYC')$alleles[samples_plot,'CNA_amp']
    toplot$TP53 = allele_effect('TP53')$alleles[samples_plot,'n_alleles_inactivated']==2
    CMP = samples_plot %in% samples_cmp
    ETS = toplot$ERG | toplot$ETV1 | toplot$ETV4 | toplot$ETV5
    BRCA2 = toplot$BRCA2 
    CHD1_SPOP = toplot$CHD1 | toplot$SPOP
    PTEN = toplot$PTEN
    RB1 = toplot$RB1
    TP53 = toplot$TP53
    MYC = toplot$MYC
    Mutation = samples_plot %in% c(samples_braf, samples_idh1, samples_tet2, samples_dnmt3b)
    Sites = metastasis_locations[samples_plot]
    colside = rep('white',dim(mat_hmr_recurrent)[2])
    names(colside) = colnames(mat_hmr_recurrent)
    stopifnot(length(colside)==length(samples_plot))
    tSCNC = colside
    tSCNC[samples_tscnc] = 'gray30'
    stopifnot(length(tSCNC)==length(samples_plot))
    
    ETS = colside
    ETS[toplot$ERG] = 'gray30'
    ETS[toplot$ETV1] = 'gray30'
    ETS[toplot$ETV4] = 'gray30'
    ETS[toplot$ETV5] = 'gray30'
    stopifnot(length(ETS)==length(samples_plot))
    
    MYC = colside
    MYC[toplot$MYC] = 'gray30'
    stopifnot(length(MYC)==length(samples_plot))
    
    SPOP = colside
    SPOP[toplot$SPOP] = 'gray30'
    stopifnot(length(SPOP)==length(samples_plot))
    
    BRCA2 = colside
    BRCA2[toplot$BRCA2] = 'gray30'
    stopifnot(length(BRCA2)==length(samples_plot))
    
    RB1 = colside
    RB1[toplot$RB1] = 'gray30'
    stopifnot(length(RB1)==length(samples_plot))
    
    PTEN = colside
    PTEN[toplot$PTEN] = 'gray30'
    stopifnot(length(PTEN)==length(samples_plot))
    
    TP53 = colside
    TP53[toplot$TP53] = 'gray30'
    stopifnot(length(TP53)==length(samples_plot))
    
    CMP = colside
    CMP[samples_cmp] = 'lightsteelblue'
    CMP[samples_braf] = 'orange'
    CMP[samples_idh1] = 'green'
    CMP[samples_tet2] = 'purple'
    CMP[samples_dnmt3b] = 'deeppink'
    #CMP['DTB-252-BL'] = 'black'
    stopifnot(length(CMP)==length(samples_plot))
    Site = colside
    Site[intersect(names(metastasis_locations)[metastasis_locations=='Bone'],samples_plot)] = wes_palette('Chevalier1')[1] #darkgreen
    Site[intersect(names(metastasis_locations)[metastasis_locations=='Lymph_node'],samples_plot)] = wes_palette('Chevalier1')[2] #gold
    Site[intersect(names(metastasis_locations)[metastasis_locations=='Liver'],samples_plot)] = wes_palette('Chevalier1')[4] #brown
    Site[intersect(names(metastasis_locations)[metastasis_locations=='Other'],samples_plot)] = wes_palette('Chevalier1')[3] #gray
    
    purity = colside
    purity[tumorpurity[colnames(mat_hmr_recurrent),'Tumor.Purity.Histo']<50] = 'black'
    stopifnot(length(purity)==length(samples_plot))
    
    csidebar = cbind(Site,CMP,tSCNC,TP53,RB1,PTEN,MYC,ETS,BRCA2)
    csidebar
}


read_bed = function(fn_bed){
    bed = read.table(fn_bed, header=FALSE, sep='\t', stringsAsFactors=FALSE)
    bed = bed[,1:4]
    names(bed) = c("chrom", "start", "end", "score")
    bed
}
read_methylation_counts = function( fn_bed, feature_type, dir_output ){
    bed = read_bed(fn_bed)
    
    hsh_bed = hsh_from_vectors( paste(bed$chrom,bed$start,sep=":"), 1:dim(bed)[1] )
    M = matrix(NA, nrow=dim(bed)[1], ncol=100 )
    chroms = c(as.character(1:22),'X','Y')
    for(i in 1:length( sample_ids_wgbs )){
        print( paste( sample_ids_wgbs[i], i, "of 100" ) )
        for( cc in 1:length( chroms )){
            fn = paste( dir_output, '/', sample_ids_wgbs[i], '_methylseekr_chr', chroms[cc], '_counts.txt', sep='', collapse='')
            cnt = read.table( fn, stringsAsFactors=FALSE, header=TRUE, sep='\t' )
            cnt = cnt[cnt$feature_type==feature_type,]
            idx = hsh_get( hsh_bed, paste( cnt$chrom, cnt$start, sep=":"))
            M[idx, i] = cnt$median
        }
    }
    dimnames(M)[[2]] = sample_ids_wgbs
    data.matrix(M)
}

heatmap_decile = function( M, quantile_val, main="" ){
    sds = apply(M[,samples_plot],1,sd)
    rows_top_decile = sds >= quantile( sds, quantile_val, na.rm=T )
    rows_top_decile[ is.na(rows_top_decile) ] = FALSE
    M = M[rows_top_decile,]
    
    h=heatmap.3(M, 
                hclustfun=hclust, 
                distfun=dist,
                na.rm=TRUE, 
                scale='none', 
                dendrogram="column", margins=c(4,4),
                Rowv=TRUE, Colv=TRUE, 
                ColSideColors=generate_csidebar(M, samples_plot), 
                RowSideColors=generate_rsidebar(M), 
                symbreaks=FALSE, 
                key=TRUE, 
                symkey=FALSE,
                density.info="none", trace="none",
                labCol=colnames(M), 
                labRow=rep("",(dim(M)[1])), 
                main=main, cexRow=1, cexCol=1, 
                col=colorRampPalette(c('blue','white','red'))(n = 1000),
                ColSideColorsSize=8, RowSideColorsSize=1)
}

df_hmr_recurrent_full = read.delim(fn_hmr_recurrent, sep='\t', header=T, check.names=F)
df_hmr_recurrent = df_hmr_recurrent_full

samples_plot = intersect(sample_ids_wgbs,colnames(df_hmr_recurrent))
mat_hmr_recurrent_full = as.matrix(df_hmr_recurrent_full[,samples_plot])
mat_hmr_recurrent_full[is.nan(mat_hmr_recurrent_full)] = NA
mat_hmr_recurrent_full[is.infinite(mat_hmr_recurrent_full)] = NA

# restrict to the top 10% of variability in HMR by standard deviation
sds = apply(mat_hmr_recurrent_full,1,sd)
rows_top_decile = sds >= quantile( sds, 0.9, na.rm=T )
rows_top_decile[ is.na(rows_top_decile) ] = FALSE
mat_hmr_recurrent = mat_hmr_recurrent_full[rows_top_decile,]

n_hmr = rep(0, 100)
for(i in 1:length(sample_ids_wgbs)){
  n_hmr[i] = length( hmr[[ sample_ids_wgbs[i]  ]] )
}
stopifnot( sum(sample_ids_wgbs==dimnames(mat_hmr_recurrent_full)[[2]]) == 100 )


# to answer:
# 1) identify unique HMR regions
 hmr_tumor = GRanges()
 for(i in 1:100){
   hmr_tumor = c(hmr_tumor, hmr[[ sample_ids_wgbs[i] ]] )
 }
 hmr_compact = reduce(hmr_tumor)

# 2) write to bed
fn_hmr_reduced='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/hmrs_reduced.bed'
bed_hmrs = data.frame(
   chrom = seqnames(hmr_compact),
   start=start(hmr_compact),
   end=end(hmr_compact),
   score = rep(1, length(hmr_compact)),
   stringsAsFactors = FALSE)
# write.table( bed_hmrs, fn_hmr_reduced,
#              row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)

# 3) extract median hmr from all loci in bed file

# FN_BED=/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/hmrs_reduced.bed
# FN_META=/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/methyl_per_HMR/metadata/meth_metadata.txt
# DIR_OUT=/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/methyl_per_HMR/processed
# python /data1/marlowe/marlowe/generate_cluster_workflow_from_metadata.py \
# -w quantify_feature_methylation \
# -i ${FN_META} \
# -v "sample_id,SAMPLE_ID;chrom,CHROM;dir_in,DIR_IN;fn_in,FN_IN" \
# --constants "FEATURE_TYPES,bed;DIR_OUT,${DIR_OUT};FNS_FEATURE_BEDS,${FN_BED};QUANTIFY_BY,cpg;FEATURE_TYPE,window;OUTPUT_TYPE,stats" \
# --json_file quant.json \
# --mem_per_core 2 \
# --mem_monolithic 2 \
# --n_cores 1 \
# -r quant \
# -d ${DIR_OUT}
# python /data1/marlowe/marlowe/process_workflow.py -i ${DIR_OUT}/quant.json

# load results by sample

dir_output='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/methyl_per_HMR/processed'
#m_hmr_mp = read_methylation_counts( fn_hmr_reduced, "bed", dir_output )
load('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/methyl_per_HMR/processed/m_hmr_mp.Rdata')    

# for Question 8, 9

fn_hmr_pearson='~/mnt/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/hmrseg_gene_pearson_spearman.txt'
hmr_pear = read.table( fn_hmr_pearson, header=TRUE, sep='\t', stringsAsFactors = FALSE)
iterm = -1*log10( as.numeric(hmr_pear$p_interact_cmp) )
hmr_pear[which(iterm>10),]

# source('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/scripts/calculate_meth_expr_correlation.R')

fn_corh3k='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_h3k27ac_1mb.txt'
fn_corhmr= '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_hmr_1mb.txt'
fn_bed_h3k27ac='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/H3K27ac_sum_hg38.bed'
fn_h3k27ac_medians='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/methyl_per_regulatory/processed/h3k27ac_methylation_medians.txt'
fn_hmr_medians='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/methyl_per_regulatory/processed/hmr_methylation_medians.txt'
fn_bed_hmr = '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/hmr_methylation_medians.bed'

#fn_bed_h3k27ac='/notebook/human_sequence_prostate_WGBS/reproduce/WCDT_WGBS/secondary_data/recurrent_HMRs/H3K27ac_sum_hg38.bed'
#fn_h3k27ac_medians='/notebook/human_sequence_prostate_WGBS/reproduce/WCDT_WGBS/secondary_data/recurrent_HMRs/h3k27ac_methylation_medians.txt'
#fn_corh3k='/notebook/human_sequence_prostate_WGBS/reproduce/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_h3k27ac_1mb.txt'
#fn_corhmr= '/notebook/human_sequence_prostate_WGBS/reproduce/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_hmr_1mb.txt'


# read in methylation values and correlations 
cor_h3k = read.table(fn_corh3k, stringsAsFactors = FALSE, sep='\t', header=TRUE)
cor_h3k$fdr = p.adjust(cor_h3k$pval, "BH")
cor_hmr = read.table(fn_corhmr, stringsAsFactors = FALSE, sep='\t', header=TRUE)
cor_hmr$fdr = p.adjust(cor_hmr$pval, "BH")
names(cor_hmr)[ which( names(cor_hmr)=="chr" )] = "chrom"

bed_h3k27ac = read.table(fn_bed_h3k27ac, stringsAsFactors = FALSE, header=FALSE, sep='\t')
names(bed_h3k27ac) = c("chrom", "start", "end", "score")
bed_h3k27ac$str = paste(bed_h3k27ac$chrom,':',bed_h3k27ac$start,'-',bed_h3k27ac$end,sep='')
mc_h3k27ac = data.matrix( read.table( fn_h3k27ac_medians, stringsAsFactors = FALSE, header=TRUE,sep='\t', check.names = FALSE) )

bed_hmr = read.table(fn_bed_hmr, stringsAsFactors = FALSE, header=TRUE, sep='\t')
names(bed_hmr) = c("chrom", "start", "end")
bed_hmr$str = paste(bed_hmr$chrom,':',bed_hmr$start,'-',bed_hmr$end,sep='')
mc_hmr = data.matrix( read.table( fn_hmr_medians, stringsAsFactors = FALSE, header=TRUE,sep='\t', check.names = FALSE) )

# Annotate correlations with the association between the enhancer methylation value and tumor purity
idx_purity = match.idx( sample_ids_wgbs, dimnames(tumorpurity)[[1]] )$idx.B
purity=tumorpurity$Tumor.Purity.Comp[idx_purity]
cor_purity = rep(NA, dim(bed_h3k27ac)[1])
has_meth = rowSums( !is.na(mc_h3k27ac) ) >10
for(i in 1:length(cor_purity)){
  if( i %% 10000 == 0  ){
    print( paste(i, length(cor_purity)))
  }
  if(has_meth[i]){
    cor_purity[i] = round( cor.test( purity, mc_h3k27ac[i,], method="spearman")$estimate, 2)  
  }
}
hsh_enh2purity = hsh_from_vectors( paste( bed_h3k27ac$chrom, ":", bed_h3k27ac$start, "-", bed_h3k27ac$end, sep=''), cor_purity )
cor_h3k$cor_purity = hsh_get( hsh_enh2purity, paste(cor_h3k$chrom,':',cor_h3k$start,'-',cor_h3k$end,sep='') )
rm(hsh_enh2purity)

cor_purity = rep(NA, dim(bed_hmr)[1])
has_meth = rowSums( !is.na(mc_hmr) ) >10
for(i in 1:length(cor_purity)){
  if( i %% 10000 == 0  ){
    print( paste(i, length(cor_purity)))
  }
  if(has_meth[i]){
    cor_purity[i] = round( cor.test( purity, mc_hmr[i,], method="spearman")$estimate, 2)  
  }
}
hsh_enh2purity = hsh_from_vectors( paste( bed_hmr$chrom, ":", bed_hmr$start, "-", bed_hmr$end, sep=''), cor_purity )
cor_hmr$cor_purity = hsh_get( hsh_enh2purity, paste(cor_hmr$chrom,':',cor_hmr$start,'-',cor_hmr$end,sep='') )
rm(hsh_enh2purity)


ensembl_left = ensembl2sym$start - PROMOTER_RANGE
ensembl_right = ensembl2sym$start + PROMOTER_RANGE
#ensembl_left = ensembl2sym$start
#ensembl_right = ensembl2sym$end
#ensembl_left[ensembl2sym$strand=="+"] = ensembl_left[ensembl2sym$strand=="+"] - PROMOTER_RANGE
#ensembl_right[ensembl2sym$strand=="-"] = ensembl_right[ensembl2sym$strand=="-"] + PROMOTER_RANGE
ensembl2left = hsh_from_vectors( dimnames(ensembl2sym)[[1]], ensembl_left)
ensembl2right = hsh_from_vectors( dimnames(ensembl2sym)[[1]], ensembl_right)

cor_h3k$in_promoter = rep(NA, dim(cor_h3k)[1])
cor_left = hsh_get( ensembl2left, cor_h3k$ensembl_id )
cor_right = hsh_get( ensembl2right, cor_h3k$ensembl_id )
cor_h3k$in_promoter[ cor_h3k$start > cor_right | 
                       cor_h3k$end   < cor_left ] = FALSE
cor_h3k$in_promoter[ cor_h3k$start > cor_left & cor_h3k$start < cor_right ] = TRUE
cor_h3k$in_promoter[ cor_h3k$end   > cor_left & cor_h3k$end   < cor_right ] = TRUE
table(cor_h3k$in_promoter)

cor_hmr$in_promoter = rep(NA, dim(cor_hmr)[1])
cor_left = hsh_get( ensembl2left, cor_hmr$ensembl_id )
cor_right = hsh_get( ensembl2right, cor_hmr$ensembl_id )
cor_hmr$in_promoter[ cor_hmr$start > cor_right | 
                       cor_hmr$end < cor_left ] = FALSE
cor_hmr$in_promoter[ cor_hmr$start > cor_left & cor_hmr$start < cor_right ] = TRUE
cor_hmr$in_promoter[ cor_hmr$end   > cor_left & cor_hmr$end   < cor_right ] = TRUE

rm(ensembl2left)
rm(ensembl2right)



# filters
FDR5_H3K = cor_h3k$fdr<0.05
FDR5_H3K[ is.na(FDR5_H3K)] = FALSE
NEG_H3K = cor_h3k$cor < (-0.5)
NEG_H3K[ is.na(NEG_H3K) ] = FALSE
PROMOTER_H3K = cor_h3k$in_promoter

FDR5_HMR = cor_hmr$fdr<0.05
FDR5_HMR[ is.na(FDR5_HMR)] = FALSE
NEG_HMR = cor_hmr$cor < (-0.5)
NEG_HMR[ is.na(NEG_HMR) ] = FALSE
PROMOTER_HMR = cor_hmr$in_promoter

tss = ensembl2sym$start
tss[ ensembl2sym$strand=="-" ] = ensembl2sym$end[ ensembl2sym$strand=="-" ]
ensembl2tss = hsh_from_vectors( rownames(ensembl2sym), tss )
cor_h3k$tss = hsh_get( ensembl2tss, cor_h3k$ensembl_id)
cor_hmr$tss = hsh_get( ensembl2tss, cor_hmr$ensembl_id)

# what is the distance from the TSS for each significant methylation association?
enhancer_mid =  rowMeans( cor_h3k[, c("start", "end")] )
cor_h3k$tss_distance = abs( cor_h3k$tss - enhancer_mid )

enhancer_mid =  rowMeans( cor_hmr[, c("start", "end")] )
cor_hmr$tss_distance = abs( cor_hmr$tss - enhancer_mid )
rm(ensembl2tss)



# ********************************************************************************
# Question 1
# As previously noted, the definition of the ‘CMP’ group in figure 1B is not convincing 
# and appears to be a highly subjective demarcation based on the presence of mutations 
# rather than a distinct hypermethylation group – the flanking nodes indeed appear to 
# be hypermethylated compared to the so-called CMP grouping.
# ********************************************************************************

pdf('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/response/Q1_dendrogram_full.pdf',
    height=10,width=8)
h=heatmap.3(mat_hmr_recurrent, 
          hclustfun=hclust, 
          distfun=dist,
          na.rm=TRUE, 
          scale='none', 
          dendrogram="column", margins=c(4,4),
          Rowv=TRUE, Colv=TRUE, 
          ColSideColors=generate_csidebar(mat_hmr_recurrent, samples_plot), 
          RowSideColors=generate_rsidebar(mat_hmr_recurrent), 
          symbreaks=FALSE, 
          key=TRUE, 
          symkey=FALSE,
          density.info="none", trace="none",
          labCol=colnames(mat_hmr_recurrent), 
          labRow=rep("",(dim(mat_hmr_recurrent)[1])), 
          main="", cexRow=1, cexCol=1, 
          col=colorRampPalette(c('blue','white','red'))(n = 1000),
          ColSideColorsSize=8, RowSideColorsSize=1)
dev.off()


dendro_CMP = h$colInd[ 79:100]
dendro_tSCNC = h$colInd[ 74:78 ]
dendro_3rd = h$colInd[ 67:73 ]
dendro_4th = h$colInd[ 40:66 ]
dendro_left = h$colInd[ 1:39 ]

cluster_assignment = rep(0, 100)
cluster_assignment[dendro_CMP] = 4
cluster_assignment[dendro_tSCNC] = 3
cluster_assignment[dendro_3rd] = 2
cluster_assignment[dendro_4th] = 1
cluster_median = colMedians( mat_hmr_recurrent_full, na.rm=TRUE )

wilcox.test( cluster_median[dendro_CMP], cluster_median[dendro_3rd])
# W = 128, p-value = 0.007695
# alternative hypothesis: true location shift is not equal to 0

wilcox.test( cluster_median[dendro_CMP], cluster_median[dendro_tSCNC])
# W = 96, p-value = 0.008027
# alternative hypothesis: true location shift is not equal to 0

wilcox.test( cluster_median[dendro_CMP], cluster_median[dendro_4th])
# W = 593, p-value = 8.048e-14
# alternative hypothesis: true location shift is not equal to 0

wilcox.test( cluster_median[dendro_CMP], cluster_median[dendro_left])
# W = 858, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

pdf('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/response/Q1_dendrogram_boxplots.pdf',
    height=8,width=10)
par(mar=c(3,5,1,1))
boxplot( cluster_median~cluster_assignment, las=1, ylim=c(15,60), ylab="% methylated at HMR", names=c("C1", "C2", "C3", "C4", "CMP") )

lines(c(1,4.8), c(60,60))
lines(c(1,1), c(59,60))
lines(c(4.8,4.8), c(59,60))
text( mean(c(1,4.8)), 61, "*", font=2, cex=2)

lines(c(2,4.8), c(58,58))
lines(c(2,2), c(57,58))
lines(c(4.8,4.8), c(57,58))
text( mean(c(2,4.8)), 59, "*", font=2, cex=2)

lines(c(3,4.8), c(56,56))
lines(c(3,3), c(55,56))
lines(c(4.8,4.8), c(55,56))
text( mean(c(3,4.8)), 57, "*", font=2, cex=2)

lines(c(4,4.8), c(54,54))
lines(c(4,4), c(53,54))
lines(c(4.8,4.8), c(53,54))
text( mean(c(4,4.8)), 55, "*", font=2, cex=2)
dev.off()

# ********************************************************************************
# QUESTION 2
# The classification of CMP is dependent on the definition of HMRs and thus there 
# needs to be a better description of HMRs, such as the distribution of fraction 
# methylation value within these regions per sample or at least the median across 
# samples like for PMDs in Fig S8. 
# ********************************************************************************

meds = colMedians( mat_hmr_recurrent_full, na.rm=TRUE )

m = match.idx( dimnames(mat_hmr_recurrent_full)[[2]], samples_cmp )
is_cmp = rep(FALSE, 100)
is_cmp[m$idx.A]=TRUE
stopifnot( length( setdiff( sample_ids_wgbs[is_cmp], samples_cmp)) == 0 )
cols = rep("white", 100)
cols[ is_cmp[ order(meds) ]] = "lightblue"

pdf('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/response/Q2_boxplot.pdf',
    height=4,width=12)
layout(matrix(1,1,1))
par(mar=c(9,5,2,1))
boxplot( mat_hmr_recurrent_full[,order( meds ) ], 
         las=2, ylab="% methylated at all HMRs",
         col=cols, 
         outline=FALSE,
         cex.axis=0.5, 
         xaxs="i")
dev.off()


# ********************************************************************************
# QUESTION 3 
# Given that the number of HMRs range from ~20K to ~90K, are CMP samples associated
# with lower numbers of HMRs generally – related to previous comments? 
# ********************************************************************************

m = match.idx( dimnames(mat_hmr_recurrent_full)[[2]], samples_cmp )
is_cmp = rep(FALSE, 100)
is_cmp[m$idx.A]=TRUE
pdf('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/response/Q3_boxplot.pdf',
    height=5,width=5)
par(mar=c(3,6,2,1))
boxplot(n_hmr/1000~is_cmp, names=c("non-CMP", "CMP"), las=1, 
        ylab="number of HMR per sample (thousands)", ylim=c(0,100))
lines(c(1,2), c(94,94))
lines(c(1,1), c(94,93))
lines(c(2,2), c(94,93))
text( 1.5, 97, "*", font=2, cex=2)
dev.off()

wilcox.test(n_hmr/1000~is_cmp)
# W = 1686, p-value = 5.756e-12


# ********************************************************************************
# QUESTION 4
# ********************************************************************************

# Manually updated figure 1 to include CMP status as a row

# ********************************************************************************
# QUESTION 5
# ********************************************************************************

# rewrote figures 1CDEF, code in main figure sections
    
# ********************************************************************************
# QUESTION 6
# A more direct approach would be to plot the mean fractional methylation across 
# all CGIs etc.. in the CMP vs. non-CMP sample – or the same for the top 10% 
# variable HMRs.
# ********************************************************************************

sds = apply(mat_hmr_recurrent_full,1,sd)
rows_top_decile = sds >= quantile( sds, 0.9, na.rm=T )
rows_top_decile[ is.na(rows_top_decile) ] = FALSE

meds = colMedians( mat_hmr_recurrent_full[rows_top_decile,], na.rm=TRUE )
pdf('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/response/Q6_boxplot.pdf',
    height=5,width=4)
boxplot( meds ~ is_cmp, 
         names=c("not CMP", "CMP"), las=1, ylab="% methylated at recurrent HMR" )
dev.off()

# ********************************************************************************
# QUESTION 7
# In this light, what is the deltaFM between the ‘methylated’ and unmethylated HMRs? 
# How do these plots compare to the flanking 3 nodes from figure 1B. 
# ********************************************************************************

# create binary matrix n_row x 100
# fill with overlaps for a given sample's HMR vs. hmr_compact
m_has_hmr = matrix(FALSE, nrow=dim(bed_hmrs)[1], ncol=100)
for(i in 1:100){
    h = hmr[[i]]   
    hits = findOverlaps(h, hmr_compact)
    m_has_hmr[to(hits), i] = TRUE
}
m_hmr_mp_present = m_hmr_mp
m_hmr_mp_absent = m_hmr_mp
m_hmr_mp_present[!m_has_hmr] = NA
m_hmr_mp_absent[m_has_hmr] = NA
mu_with_hmr = rowMeans(m_hmr_mp_present, na.rm=TRUE )
mu_without_hmr = rowMeans(m_hmr_mp_absent, na.rm=TRUE )

pdf('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/response/Q7_scatter_HRM_present_absent.pdf',
    height=10,width=8)
layout(matrix(1,1,1))
par(mar=c(5,5,2,1))
plot( 100*mu_with_hmr, 100*mu_without_hmr, pch=19, cex=0.1, col="#33333333", las=1, 
      xlab="rHMR is present",
      ylab="rHMR is absent",
      main="Mean % methylation at rHMR loci", xaxs="i", yaxs="i")
dev.off()

# calculate difference in means for each

summary( mu_without_hmr - mu_with_hmr)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -0.691   0.280   0.433   0.423   0.586   0.998   11388 

#The average methylation difference between a sample with a HMR and a sample without a HMR is
# a 48% reduction in methylation.


# ********************************************************************************
# QUESTION 8
# As previously noted, methylation gains associated with IDH/TET2 mutations are typically associated with 
# regulatory features, most notably enhancers. Analysis for Fig3 (eHMR) is a missed opportunity not to look 
# at HMR (or at least eHMR) in the context of promoter and non-promoter and their association with gene 
# expression. While there are some handpicked genes at end of paper this needs to be done in an unbiased way. 
# Given the focus on CMP, it would be useful in this context then to compare CMP and non-CMP at these eHMR 
# regions. For example are there differences between CMP and non-CMP at AR regulatory regions and does 
# that correlate with expression?
# ********************************************************************************

# how many significant methylation associations?
sum(FDR5_H3K, na.rm=TRUE)
# [1] 59395
sum(FDR5_HMR, na.rm=TRUE)
# [1] 58436

# how many transcripts have a significant methylation association?
length(unique(cor_h3k$ensembl_id[ which(FDR5_H3K & cor_h3k$cor<0) ] ) )
# [1] 10,412
length(unique(cor_hmr$ensembl_id[ which(FDR5_HMR & cor_hmr$cor<0) ] ) )
# [1] 11,928

# combined overall 
length( which(FDR5_H3K & cor_h3k$cor<0) ) + length( which(FDR5_HMR & cor_hmr$cor<0) )
# [1] 71163



# ------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------

pdf('~/mnt/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/response/Q8_correlation_plots.pdf',
    height=8,width=10)
layout(matrix(1:8,4,2),heights=c(1,3,1,3))
idx=which(PROMOTER_H3K & FDR5_H3K & cor_h3k$cor<0 & cor_h3k$tss_distance<=5000 )
par(mar=c(0,5,2,2))
h=hist(cor_h3k$tss_distance[idx], xlim=c(0,5000), breaks=200, col="grey",
     xaxs="i", axes=FALSE,ylab="Frequency", xlab="", 
     main="Correlated H3k27ac loci, within promoter", freq = TRUE, 
     ylim=c(0,60))
axis(2,at=c(0,30,60), labels = c(0, 30,60),las=1)
par(mar=c(4,5,0,2))
plot( cor_h3k$tss_distance[idx], cor_h3k$cor[idx], 
      main="",
      pch=19,col="#333333ee", cex=0.2, las=1, 
      ylab="correlation coefficient", xlab="distance from locus center to TSS",
      xlim=c(0,5000), ylim=c(-0.8, -0.3), 
      xaxs="i", yaxs="i", axes=FALSE)
axis(2, at=seq(from=-0.8,to=-0.4, by=0.1), las=1 )
axis(1, at=seq(from=0, to=5000,by=1000), labels = paste( seq(from=0,to=5,by=1), "kb") )
box()


#layout(matrix(1:4,2,2),heights=c(1,3))
idx=which(PROMOTER_HMR & FDR5_HMR & cor_hmr$cor<0 & cor_hmr$tss_distance<=5000 )
par(mar=c(0,5,2,2))
h=hist(cor_hmr$tss_distance[idx], xlim=c(0,5000), col="grey", breaks=200,
       xaxs="i", axes=FALSE,ylab="Frequency", xlab="", 
       main="Correlated hypomethylated loci, within promoter", ylim=c(0,60))
axis(2,at=c(0,30,60), labels = c(0, 30,60),las=1)

par(mar=c(4,5,0,2))
plot( cor_hmr$tss_distance[idx], cor_hmr$cor[idx], 
      main="",
      pch=19,col="#333333ee", cex=0.2, las=1, 
      ylab="correlation coefficient", xlab="distance from locus center to TSS",
      xlim=c(0,5000), ylim=c(-0.8, -0.3), 
      xaxs="i", yaxs="i", axes=FALSE)
axis(2, at=seq(from=-0.8,to=-0.3, by=0.1), las=1 )
axis(1, at=seq(from=0, to=5000,by=1000), labels = paste( seq(from=0,to=5,by=1), "kb") )
box()

idx=which(!PROMOTER_H3K & FDR5_H3K & cor_h3k$cor<0 )
par(mar=c(0,5,2,2))
h=hist(cor_h3k$tss_distance[idx], xlim=c(0,900000), breaks=200, col="grey",
       xaxs="i", axes=FALSE,ylab="Frequency", xlab="", 
       main="Correlated H3K27ac loci, outside promoter")
axis(2,at=c(0,1500,3000), labels = c(0, 1500,3000),las=1)

par(mar=c(4,5,0,2))
plot( cor_h3k$tss_distance[idx], cor_h3k$cor[idx], 
      main="",
      pch=19,col="#333333ee", cex=0.2, las=1, 
      ylab="correlation coefficient", xlab="distance from locus center to TSS",
      xlim=c(0,900000), ylim=c(-0.8, -0.3), 
      xaxs="i", yaxs="i", axes=FALSE)
axis(2, at=seq(from=-0.8,to=-0.4, by=0.1), las=1 )
axis(1, at=seq(from=0, to=900000,by=100000), labels = paste( seq(from=0,to=900,by=100), "kb") )
box()

idx=which(!PROMOTER_HMR & FDR5_HMR & cor_hmr$cor<0 )
par(mar=c(0,5,2,2))
h=hist(cor_hmr$tss_distance[idx], xlim=c(0,900000), breaks=200, col="grey",
       xaxs="i", axes=FALSE,ylab="Frequency", xlab="", 
       main="Correlated hypomethylated loci, outside promoter",
       ylim=c(0,4000))
axis(2,at=c(0,2000,4000), labels = c(0, 2000,4000),las=1)

par(mar=c(4,5,0,2))
plot( cor_hmr$tss_distance[idx], cor_hmr$cor[idx], 
      main="",
      pch=19,col="#333333ee", cex=0.2, las=1, 
      ylab="correlation coefficient", xlab="distance from locus center to TSS",
      xlim=c(0,900000), ylim=c(-0.8, -0.3), 
      xaxs="i", yaxs="i", axes=FALSE)
axis(2, at=seq(from=-0.8,to=-0.4, by=0.1), las=1 )
axis(1, at=seq(from=0, to=900000,by=100000), labels = paste( seq(from=0,to=900,by=100), "kb") )
box()
dev.off()

# ------------------------------------------------------------
# TABLES
# ------------------------------------------------------------

idx=which(FDR5_H3K & cor_h3k$cor<0 )
h1 = cor_h3k[idx,c("ensembl_id","symbol","cor","pval","logp","type","chrom","start","end","fdr","tss_distance") ]

idx=which(FDR5_HMR & !is.na( cor_hmr$cor) & cor_hmr$cor<0 )
h2 = cor_hmr[idx, c("ensembl_id","symbol","cor","pval","logp","type","chrom","start","end","fdr","tss_distance") ]
comb = rbind(h1, h2)
comb$cor = round( comb$cor,3)
comb$pval = signif( comb$pval,3)
comb$fdr = round( comb$fdr,3)
comb = comb[order(comb$chrom, comb$start),]
write.table( comb, fn_table_ehmr, row.names=FALSE, sep='\t', quote=FALSE)


# ##############################################################################
# Question 9
# Given the focus on CMP, it would be useful in this context then to compare 
# CMP and non-CMP at these eHMR regions. For example are there differences 
# between CMP and non-CMP at AR regulatory regions and does that correlate with 
# expression?
# ##############################################################################


idx = FDR5_HMR & cor_hmr$cor<0
strs = paste( cor_hmr$chrom[idx],":",cor_hmr$start[idx],"-",cor_hmr$end[idx],sep='')
cors = cor_hmr$cor[idx]
distances = cor_hmr$tss_distance[idx]
m=match.idx( bed_hmr$str, strs, allow.multiple.B = TRUE )
mu_cmp = rowMeans( mc_hmr[,is_cmp], na.rm=TRUE)[m$idx.A]
mu_not = rowMeans( mc_hmr[,!is_cmp], na.rm=TRUE)[m$idx.A]
dm = mu_cmp - mu_not 
layout(matrix(1,1,1))
par(mar=c(4,4,3,1))
cols = rep("#33333333", length(dm))
cols[ distances < 1500 ] = "blue"
layout(matrix(1:2,2,1))

plot( dm[distances<=1500], cors[distances<=1500], pch=19, cex=0.2, las=1,
      xlab="mean difference in methylation when CMP",
      ylab="Correlation between locus methylation and expression", col="#33333333",
      main="rHMR within promoters (<1500 bp to TSS)")

plot( dm[distances>1500], cors[distances>1500], pch=19, cex=0.2, las=1,
      xlab="mean difference in methylation when CMP",
      ylab="Correlation between locus methylation and expression", col="#33333333",
      main="rHMR outside of promoters")
cor.test( dm, cors, method="spearman")

