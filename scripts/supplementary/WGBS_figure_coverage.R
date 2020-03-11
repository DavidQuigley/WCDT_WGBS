#available.genomes()
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(forcats)
library(ggplot2)
library(plyr)
library(BSgenome)
library(rtracklayer)
library(DNAcopy)
library("BSgenome.Hsapiens.UCSC.hg38")
library(MethylSeekR)
library(parallel)

sLengths=seqlengths(Hsapiens)

match.idx=function(A, B, allow.multiple.B=F){
    # return dataframe of indices into A and B restricted to perfect matches
    # between A and B, where idx.A[i] == idx.B[i] for each i in matched pairs
    if( allow.multiple.B ){
        idx.B = which(B %in% A)
        idx.A = match(B[idx.B], A)
    }
    else{
        in.both = intersect(A,B)
        idx.A = match(in.both, A)
        idx.B = match(in.both, B)
    }
    C= data.frame(idx.A, idx.B)
    if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
        stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
    C
}

BASE_DIR = "/data1/projects/WCDT_WGBS_2019/WCDT_WGBS"
make_path = function(BASE_DIR, s){
    paste(BASE_DIR, s, sep='/' )
}
fn_cov_sam = make_path(BASE_DIR, 'metadata/coverage/cpg_depth_by_sample_both_strands.txt' )
fn_cov_chr = make_path(BASE_DIR, 'metadata/coverage/cpg_coverage_per_chromosome_both_strands.txt' )

chrom_lengths = read.table(make_path(BASE_DIR, 'metadata/genome_builds/HG38_chromosome_lengths.txt'),
                           sep='\t')
                           
c_sam = read.table(fn_cov_sam, sep='\t',header=TRUE, stringsAsFactors=FALSE)
c_sam = c_sam[c_sam$chrom !='NC_000024.10',]
c_chr = read.table(fn_cov_chr, sep='\t',header=TRUE, stringsAsFactors=FALSE)

n0 = ddply( c_sam, .(sample_id), summarize, sum(n_0) )
n10 = ddply( c_sam, .(sample_id), summarize, sum(n_10) )
n20 = ddply( c_sam, .(sample_id), summarize, sum(n_20) )
n30 = ddply( c_sam, .(sample_id), summarize, sum(n_30) )
names(n0) = c("sample_id", "total")
names(n10) = c("sample_id", "total")
names(n20) = c("sample_id", "total")
names(n30) = c("sample_id", "total")

layout(matrix(1,1,1))
boxplot( n0$total, n10$total, n20$total, n30$total,
         boxwex=0.25, las=1, col="lightblue",
         main="CpG Coverage",
         names=c("CpGs tested", "≥ 10x", "≥ 20x", "≥ 30x"))

summary(n10$total)
sd(n10$total)
# mean 26.3M +/- 1.59M


barplot( round(n10$total / n0$total,2)*100, ylim=c(0,100), col="lightblue", las=1,
         ylab="Percent assayed at >= 10x")

write.table( 
    data.frame( sample_id = n0$sample_id, 
            percent_10x=round(n10$total / n0$total,2), stringsAsFactors=FALSE),
    make_path(BASE_DIR, 'metadata/coverage/percent_10x_per_sample.txt'),
              sep='\t', quote=FALSE, row.names=FALSE)


#--------------------

round(rowSums( c_chr[,92:101]) / rowSums(c_chr[2:101]), 2)
sum( c_chr[,92:101]) / sum( c_chr[2:101])
sum( c_chr[,96:101]) / sum( c_chr[2:101])

sum( c_chr[1:23, 96:101]) / sum( c_chr[1:23, 2:101])

sum( c_chr[1:23, 101]) / sum( c_chr[1:23, 2:101])

#----------------------------------
# CNA WGBS vs WGS
#----------------------------------

# Load binned copy number data across all WGS tumors
# -------------------------------------------------
fn_binned_CN = "/data1/projects/WCDT_WGS_2018/build_2018_04_15/2018_04_15_matrix_binned_weighted_CN_copycat.txt"
cnbin = read.table(fn_binned_CN,
                   header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                   check.names = FALSE)


# calculate cnbin (3Mb windows) for all chromosomes
# -------------------------------------------------

#How did I calculate CNA_meth/processed?
# I called quantify_feature_methylation on bins.
# qfm does not combine adjacent nucleotides correctly.

sample_ids_wgs = c("DTB-003-BL", "DTB-005-BL", "DTB-008-BL", "DTB-009-BL", "DTB-011-BL", "DTB-018-BL", "DTB-019-PRO", "DTB-020-BL", "DTB-021-BL", 
                   "DTB-022-BL", "DTB-023-BL", "DTB-024-PRO", "DTB-034-BL", "DTB-035-BL", "DTB-036-BL", "DTB-037-BL", "DTB-040-BL", "DTB-042-BL", 
                   "DTB-053-BL", "DTB-055-PRO", "DTB-059-BL", "DTB-060-BL", "DTB-061-BL", "DTB-063-BL", "DTB-064-BL", "DTB-067-PRO", "DTB-069-BL", 
                   "DTB-071-BL", "DTB-074-BL", "DTB-077-PRO", "DTB-080-BL", "DTB-083-BL", "DTB-085-BL", "DTB-089-BL", "DTB-090-PRO", "DTB-091-BL", 
                   "DTB-092-BL", "DTB-094-BL", "DTB-097-PRO", "DTB-098-PRO2", "DTB-100-BL", "DTB-101-BL", "DTB-102-PRO", "DTB-104-BL", "DTB-111-PRO", 
                   "DTB-112-BL", "DTB-119-PRO", "DTB-121-BL", "DTB-124-BL", "DTB-126-BL", "DTB-127-PRO", "DTB-128-BL", "DTB-129-BL", "DTB-132-BL", 
                   "DTB-135-PRO", "DTB-137-PRO", "DTB-138-BL", "DTB-140-BL", "DTB-141-BL", "DTB-143-BL", "DTB-146-BL", "DTB-149-BL", "DTB-151-BL", 
                   "DTB-156-BL", "DTB-159-BL", "DTB-165-PRO", "DTB-167-PRO", "DTB-170-BL", "DTB-172-BL", "DTB-173-BL", "DTB-175-BL", "DTB-176-BL", 
                   "DTB-183-BL", "DTB-186-BL", "DTB-187-BL", "DTB-188-BL", "DTB-190-BL", "DTB-193-BL", "DTB-194-PRO", "DTB-201-PRO", "DTB-202-BL", 
                   "DTB-204-BL", "DTB-205-BL", "DTB-206-BL", "DTB-210-BL", "DTB-213-BL", "DTB-214-BL", "DTB-216-PRO", "DTB-220-BL", "DTB-222-BL", 
                   "DTB-223-BL", "DTB-232-PRO", "DTB-234-BL", "DTB-251-BL", "DTB-252-BL", "DTB-255-BL", "DTB-258-BL", "DTB-260-BL", "DTB-261-BL", 
                   "DTB-265-PRO", "DTB-266-BL")
sample_ids_wgbs = setdiff(sample_ids_wgs, "DTB-193-BL")

cnbin_wgbs = cnbin
cnbin_wgbs[,3:103] = NA
dir_cnbin = '/data1/projects/WCDT_WGBS_2019/CNA_meth/processed'
chroms = unique(cnbin_wgbs$chrom)
keys = paste( cnbin_wgbs$chrom, cnbin_wgbs$bin_start)
for( i in 1:100){
    sample_id = sample_ids_wgbs[i]
    print(sample_id)
    idx = which( dimnames(cnbin_wgbs)[[2]] == sample_id )
    for( j in 1:24){
        chrom = chroms[j]
        fn=paste(dir_cnbin,"/",sample_id,"_methylseekr_",chrom,"_counts.txt",sep='')
        data=read.table( fn, header=TRUE, stringsAsFactor=FALSE)
        m = match.idx( keys, paste( chrom, data$start))
        cnbin_wgbs[m$idx.A,idx] = data$median_coverage[m$idx.B]
    }
}

# convert to vector
y_wgbs = c()
y_wgs = c()
for(i in 1:length(sample_ids_wgbs)){
    idx_cnbin= which(dimnames(cnbin)[[2]] == sample_ids_wgs[i] )
    y_wgs = c(y_wgs, cnbin[,idx_cnbin] )
    idx_wgbs = which( dimnames(cnbin_wgbs)[[2]] == sample_ids_wgbs[i] )
    y_wgbs = c(y_wgbs, cnbin_wgbs[,idx_wgbs] )
}

fn_wgs_wgbs_scatter = make_path(BASE_DIR, 'figures/supplementary_figure_WGBS_vs_WGS_CNA.pdf' )
pdf(fn_wgs_wgbs_scatter, width = 10, height = 10)

layout(matrix(1,1,1))
par(mar=c(4,4,2,1))
plot( y_wgs, y_wgbs, col="#0000ff33", cex=0.5, pch=19,
      xlab="WGS copy number", ylab="WGBS depth", xlim=c(0,10), ylim=c(0,200))
cor.test( y_wgs, y_wgbs, method="spearman")
abline( lm( y_wgbs~y_wgs ))
dev.off()


#-------------

fn_wgs_wgbs = make_path(BASE_DIR, 'figures/supplementary_figure_CNA.pdf' )

pdf(fn_wgs_wgbs, width = 14, height = 4)
samples = sample_ids_wgs[c(11,12,13,23)]
layout(matrix(1:8,nrow=2,ncol=4,byrow=FALSE))
idx_chrX = which( cnbin$chrom=="chrX" )

for(i in 1:4){
    par(mar=c(1,3,2,1))
    idx_cnbin = which(dimnames(cnbin)[[2]] == samples[i] )
    plot( cnbin[idx_chrX,idx_cnbin], ylim=c(0,6),las=1, pch=19, cex=0.5)
    text(0.5,.75,"DNA",adj=0, font=2)
    text(0.5,.3,samples[i], adj=0, font=2)
    par(mar=c(2,3,2,1))
    idx_wgbs = which(dimnames(cnbin_wgbs)[[2]] == samples[i] )
    plot( cnbin_wgbs[idx_chrX,idx_wgbs],  ylim=c(0,60),las=1, pch=19, cex=0.5)
    text(0.5,7,"WGBS",adj=0, font=2)
    text(0.5,3,samples[i],adj=0, font=2)
}
dev.off()

