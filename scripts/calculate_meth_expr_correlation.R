#qmod -d all.q@node1
library(doParallel)

source('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/scripts/2019_05_15_prepare_environment.r')

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

ensembl_ids = rownames(tpm)
chroms = as.character( ensembl2sym$chr )
starts = ensembl2sym$start
ends = ensembl2sym$end
strands = as.character( ensembl2sym$strand )

fn_bed_h3k27ac = '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/chipseq/pca100/H3K27ac_sum_hg38.bed'
#fn_bed_promoter = '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/metadata/hg38_promoters.bed'
#fn_bed_hmr = '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/hmrs_reduced.bed'
#dir_output_hmr='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/methyl_per_HMR/processed'

dir_output='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/methyl_per_regulatory/processed'
#fn_mc_h3k27ac = paste(dir_output,"h3k27ac_methylation_medians.txt", sep="/")
#fn_mc_hmr = paste(dir_output,"hmr_methylation_medians.txt", sep="/")


mc_h3k27ac = read_methylation_counts( fn_bed_h3k27ac, "H3K27AC", dir_output )
#mc_promoter = read_methylation_counts( fn_bed_promoter, "PROMOTER", dir_output )
#mc_HMR = read_methylation_counts( fn_bed_hmr, "bed", dir_output_hmr)
#bed_promoter = read_bed(fn_bed_promoter)
bed_h3k27ac = read_bed(fn_bed_h3k27ac )
#bed_hmr = read_bed(fn_bed_hmr)

outputfile <- '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_h3k27ac_1mb.txt'

# Set up cores for parallel processing
#cl <- makeCluster( 60, outfile="")
#registerDoParallel(cl)
samples_cmp = c(
    "DTB-018-BL", "DTB-021-BL", "DTB-022-BL",  "DTB-023-BL",
    "DTB-024-PRO", "DTB-037-BL","DTB-042-BL","DTB-074-BL",
    "DTB-089-BL", "DTB-094-BL", "DTB-102-PRO", "DTB-112-BL",
    "DTB-124-BL", "DTB-137-PRO", "DTB-141-BL", "DTB-151-BL",
    "DTB-188-BL", "DTB-190-BL", "DTB-194-PRO", "DTB-202-BL", 
    "DTB-204-BL", "DTB-260-BL")

fn_hmr_recurrent='/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/hmr_recurrent.txt'
df_hmr_recurrent_full = read.delim(fn_hmr_recurrent, sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors = FALSE)
df_hmr_recurrent = df_hmr_recurrent_full

samples_plot = intersect(sample_ids_wgbs,colnames(df_hmr_recurrent))
mat_hmr_recurrent_full = as.matrix(df_hmr_recurrent_full[,samples_plot])
mat_hmr_recurrent_full[is.nan(mat_hmr_recurrent_full)] = NA
mat_hmr_recurrent_full[is.infinite(mat_hmr_recurrent_full)] = NA


for(i in 1:10000 ){
    ensembl_id = ensembl_ids[i]
    idx_symbol = which( rownames(ensembl2sym) == ensembl_id)
    symbol = ensembl2sym$name[ which( rownames(ensembl2sym)==ensembl_id ) ][1]
    print( paste(i, symbol, ensembl_id))
    
    cur_chrom = chroms[idx_symbol]
    cur_start = starts[idx_symbol]
    cur_end = ends[idx_symbol]
    cur_strand = strands[ idx_symbol ]
    idx_samples = match.idx(sample_ids_wgbs, dimnames(tpm)[[2]])$idx.B
    xp = log10(1+as.numeric( tpm[ensembl_id,idx_samples] ))
    if( cur_strand=="+"){
        idx_h3 = which(bed_h3k27ac$chrom==cur_chrom &
                           bed_h3k27ac$start > (cur_start - 900000) &
                           bed_h3k27ac$end < (cur_end + 100000))
    }else{
        idx_h3 = which(bed_h3k27ac$chrom==cur_chrom &
                           bed_h3k27ac$start > (cur_start - 100000) &
                           bed_h3k27ac$end < (cur_end + 900000))
    }
    cors = rep(NA, length(idx_h3))
    pvals = rep(NA, length(idx_h3))
    for(i in 1:length(idx_h3)){
        mp = mc_h3k27ac[ idx_h3[i],]
        has_both = !is.na( xp ) & !is.na( mp )
        if(sum(has_both)>2){
            cc=cor.test( xp, mp, method="spearman" )
            cors[i]=cc$estimate
            pvals[i]=cc$p.value
        }
    }
    if( sum( !is.na(cors))>0 ){
        res=data.frame( 
            ensembl_id=rep(ensembl_id,length(cors)),
            symbol=rep(symbol,length(cors)),
            cor=cors,
            pval=pvals,
            logp = round(-1*log10(pvals),3),
            type=rep( "H3K27AC", length(cors)),
            stringsAsFactors=FALSE
        )
        cbind(res, bed_h3k27ac[idx_h3,])
    }
}
stopCluster(cl)

write.table(res_all,outputfile,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
write.table( mc_h3k27ac, fn_mc_h3k27ac,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')





res_hr = c()

#for( i in 1:length(ensembl_ids) ){
for(i in 40001:58381 ){
    ensembl_id = ensembl_ids[i]
    idx_symbol = which( rownames(ensembl2sym) == ensembl_id)
    symbol = ensembl2sym$name[ which( rownames(ensembl2sym)==ensembl_id ) ][1]
    print( paste(i, length(ensembl_ids), symbol, ensembl_id))
    
    cur_chrom = chroms[idx_symbol]
    cur_start = starts[idx_symbol]
    cur_end = ends[idx_symbol]
    cur_strand = strands[ idx_symbol ]
    idx_samples = match.idx(sample_ids_wgbs, dimnames(tpm)[[2]])$idx.B
    xp = log10(1+as.numeric( tpm[ensembl_id,idx_samples] ))
    if( cur_strand=="+"){
        idx_hmr = which(df_hmr_recurrent$chr==cur_chrom &
                        df_hmr_recurrent$start > (cur_start - 900000) &
                        df_hmr_recurrent$end < (cur_end + 100000))
    }else{
        idx_hmr = which(df_hmr_recurrent$chr==cur_chrom &
                        df_hmr_recurrent$start > (cur_start - 100000) &
                        df_hmr_recurrent$end < (cur_end + 900000))
    }
    cors = rep(NA, length(idx_hmr))
    pvals = rep(NA, length(idx_hmr))
    for(i in 1:length(idx_hmr)){
        mp = mat_hmr_recurrent_full[ idx_hmr[i],]
        has_both = !is.na( xp ) & !is.na( mp )
        if(sum(has_both)>2){
            cc=cor.test( xp, mp, method="spearman" )
            cors[i]=cc$estimate
            pvals[i]=cc$p.value
        }
    }
    if( sum( !is.na(cors))>0 ){
        res=data.frame( 
            ensembl_id=rep(ensembl_id,length(cors)),
            symbol=rep(symbol,length(cors)),
            cor=cors,
            pval=pvals,
            logp = round(-1*log10(pvals),3),
            type=rep( "HMR", length(cors)),
            stringsAsFactors=FALSE
        )
        res_hr = rbind( res_hr, cbind(res, df_hmr_recurrent[idx_hmr,1:3]), stringsAsFactors=FALSE)
    }
}
outputfile_hmr1 <- '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_hmr_1mb_1.txt'
outputfile_hmr2 <- '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_hmr_1mb_2.txt'
outputfile_hmr3 <- '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_hmr_1mb_3.txt'
outputfile_hmr4 <- '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_hmr_1mb_4.txt'
outputfile_hmr5 <- '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/cor_meth_expr_at_hmr_1mb_5.txt'

write.table(res_hr,outputfile_hmr5,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
write.table( mat_hmr_recurrent_full, fn_mc_hmr,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')

r1=read.table(outputfile_hmr1,header=TRUE,sep='\t', stringsAsFactors=FALSE)
r2=read.table(outputfile_hmr2,header=TRUE,sep='\t', stringsAsFactors=FALSE)
r3=read.table(outputfile_hmr3,header=TRUE,sep='\t', stringsAsFactors=FALSE)
r4=read.table(outputfile_hmr4,header=TRUE,sep='\t', stringsAsFactors=FALSE)
r5=read.table(outputfile_hmr5,header=TRUE,sep='\t', stringsAsFactors=FALSE)

res_all = rbind(r1, r2, r3, r4, r5, stringsAsFactors=FALSE)
write.table(res_hr,outputfile_hmr,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')


#qmod -e all.q@node1
