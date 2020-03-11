#########
# GLOBAL VARIABLES
#######
window_szs <- 1000
xscale <- 'norm'

# Relative widths of gene body vs upstream and downstream regions in methylation to RNA expression correlation plot
norm_genesz <- 1000
norm_upstreamsz <- 1000
norm_downstreamsz <- norm_upstreamsz # Required to be the same as norm_upstreamsz for now

##################
# Helper functions
###################

# TODO: decide whether we want to use filename coordinates or actual methylation min/max coordinates as start and end
extract_metadata <- function(fname) {
  retval <- list()
  retval$ensembl_name <- strsplit(fname, '_')[[1]][1]
  retval$hgnc_name <- strsplit(fname, '_')[[1]][2]
  retval$chr <- strsplit(fname, '_')[[1]][3]
  retval$start <- as.numeric(strsplit(fname, '_')[[1]][4])
  retval$end <- as.numeric(strsplit(fname, '_')[[1]][5])
  retval$strand <- substr(strsplit(fname, '_')[[1]][6],1,1)
  return(retval)
}

#
extract_metadata_are <- function(fname) {
  retval <- list()
  retval$tier <- strsplit(fname, '_')[[1]][1]
  retval$chr <- strsplit(fname, '_')[[1]][3]
  retval$start <- as.numeric(strsplit(fname, '_')[[1]][4])
  retval$end <- as.numeric(strsplit(fname, '_')[[1]][5])
  return(retval)
}

# normalize coordinates
normalize_coords <- function(data, gene_start, gene_end, val=NA) {
  if(is.na(val)) {
    cur_vals <- data$X
  } else {
    cur_vals <- val
  }
  # scale gene body to be standard width (i.e. norm_genesz)
  retval <- (cur_vals - gene_start) / (gene_end - gene_start) * norm_genesz
  # scale upstream region to be standard width (i.e. norm_upstreamsz)
  if(sum(retval < 0) > 0) {
    retval[retval < 0] <- retval[retval < 0] / (0 - min(retval)) * norm_upstreamsz
  }
  # scale downstream region to be standard width (i.e. norm_downstreamsz)
  if(sum(retval > norm_genesz) > 0) {
    downstream_coords_norm <- retval[retval > norm_genesz] - norm_genesz
    downstream_coords_norm <- downstream_coords_norm / max(downstream_coords_norm) * norm_downstreamsz
    retval[retval > norm_genesz] <- downstream_coords_norm + norm_genesz
  }
  
  #retval <- (cur_vals - min(data$X)) / (max(data$X)- min(data$X))*10000
  return(retval)
}

########################
# Experimental analysis
########################
library(ggplot2)

curdate <- '2019-08-19'
#setwd('/Users/Will/Desktop/med_school_homebase/feng_lab/wgbs/data/gene_selected_output_cor_dec15/cor/')
setwd('/media/will/Elements/wgbs_2019-08-10/cor_gene_output/')
savedest <- sprintf('/home/will/Desktop/wgbs/revision.1/plots/%s/', curdate)

filenames <- list.files(".", pattern="*.txt", full.names=FALSE)
all_data <- list()

cur_filename <- filenames[1]

# Store methyl data with normalized xcoordinates for each gene
all_methyl_norm <- list()
counter <- 1 
#Broken filenames to exclude
idx.excl <- c(13532, 13933, 20234, 21805, 21949, 21977, 22087, 22991, 23754, 24079, 24610, 25062, 26278, 26503, 26753, 27522, 28455, 28529,
             29655, 30211, 30336, 30458, 32576, 33306, 33502, 40764, 43128, 45092, 47550, 49539)
for(i in counter:length(filenames)) {
  if(i %in% idx.excl) next
  cur_filename <- filenames[i]
  if(counter %% 100 == 0) print(counter)
  cur_meth_tab <- read.delim(file = cur_filename, sep='\t', header=T)
  head(cur_meth_tab)
  if(nrow(cur_meth_tab) == 0) {
    print(sprintf('Skipping file %s', cur_filename))
    counter <- counter + 1
    next
  }
  
  cur_metadata <- extract_metadata(cur_filename) 
  
  
  # Ensure that we're only using coordinates that are within 10kb of gene start and gene end (e.g. for AR)
  cur_meth_tab <- cur_meth_tab[which(cur_meth_tab$X >= cur_metadata$start - 10000 & cur_meth_tab$X <= cur_metadata$end + 10000),]
  
  # TODO: allow oncogene vs tumor suppressor designation
  
  # adjusted X coordinates for track that is alignable across genes
  cur_meth_tab$coordAdj <- normalize_coords(cur_meth_tab, gene_start=cur_metadata$start, gene_end=cur_metadata$end)
  head(cur_meth_tab)
  
  #all_data[[cur_metadata$hgnc_name]] <- cur_meth_tab
  cur_meth_tab$spearmanRho <- cur_meth_tab$spearman_rho #spearman correlation value
  cur_meth_tab$chr <- cur_metadata$chr
  cur_meth_tab$pos <- cur_meth_tab$X
  all_methyl_norm[[cur_filename]] <- cur_meth_tab[,c('chr', 'pos', 'coordAdj', 'spearmanRho', 'spearman_p')]
  counter <- counter+1
}

# Align normalized negative sense gene coordinates with normalized positive sense gene coordinates
for(curfile in names(all_methyl_norm)) {
  cur_metadata <- extract_metadata(curfile) 
  if(cur_metadata$strand == '-') {
    curdata <- all_methyl_norm[[curfile]]
    curdata$coordAdj <- curdata$coordAdj * -1 + norm_genesz
    all_methyl_norm[[curfile]] <- curdata
  }
}
head(all_methyl_norm[[cur_filename]])
save(all_methyl_norm, file='/media/will/Elements/wgbs_2019-08-10/normCoordAllGenesCor-10kb.RData')

load('/media/will/Elements/wgbs_2019-08-10/normCoordAllGenesCor-10kb.RData')

######################
# Load additional data
###################

# Get Ensembl gene names in order that they occur in all_methyl_norm list
gene_names <- sapply(names(all_methyl_norm), function(x) {
  strsplit(x, '_')[[1]][1]
})

# Load methyl cor segments and identify genes with statistically significant EMRs
methyl_cor_segments <- read.delim(file='/media/will/Elements/wgbs_2019-08-10/hmrseg_gene.txt', header=T, sep='\t', stringsAsFactors=F)
methyl_cor_segments_sig <- methyl_cor_segments[which(p.adjust(methyl_cor_segments$avg_p, method='fdr') < 0.01),]

# Positively correlated with expression
methyl_cor_segments_pos_sig <- methyl_cor_segments_sig[which(methyl_cor_segments_sig$avg_cor > 0),]
genes_pos_sig <- unique(methyl_cor_segments_pos_sig$gene_id)
pos_sig_idx <- which(gene_names %in% genes_pos_sig)

# Negatively correlated with expression
methyl_cor_segments_neg_sig <- methyl_cor_segments_sig[which(methyl_cor_segments_sig$avg_cor < 0),]
genes_neg_sig <- unique(methyl_cor_segments_neg_sig$gene_id)
neg_sig_idx <- which(gene_names %in% genes_neg_sig)

length(pos_sig_idx)
length(neg_sig_idx)

all_metadata <- lapply(names(all_methyl_norm), extract_metadata)

# Fit 1 spline for all Cpgs concatenated (1 spline for all genes) - positively correlated
direction <- 'pos'
all_cpgs_list <- list()
counter <- 1
library(GenomicRanges)
for(j in pos_sig_idx) {
  if(counter %% 1000 == 0) print(counter)
  curtab <- all_methyl_norm[[j]]
  curtab <- curtab[complete.cases(curtab),]
  
  #Order in terms of increasing genomic position (upstream -> downstream)
  curtab <- curtab[order(curtab$coordAdj, decreasing = F),] 
  curgene_metadata <- all_metadata[[j]]
  
  if(direction == 'pos') {
    # Only use positive correlations
    curtab <- curtab[which(curtab$spearmanRho >= 0 & curtab$spearman_p < 0.05),]
  } else if(direction == 'neg') {
    # Only use positive correlations
    curtab <- curtab[which(curtab$spearmanRho <= 0 & curtab$spearman_p < 0.05),]
  }
  
  all_cpgs_list[[counter]] <- curtab
  counter <- counter + 1
}

library(data.table)
all_cpgs_tab <- rbindlist(all_cpgs_list)
all_cpgs_tab <- all_cpgs_tab[order(all_cpgs_tab$coordAdj, decreasing=F),]
cur.methyl.smooth <- smooth.spline(all_cpgs_tab[,c('coordAdj', 'spearmanRho')], spar=0.9)
cur_x <- seq(-1000,2000, by=0.1)
cur_y <- predict(cur.methyl.smooth, cur_x)$y
pdf(sprintf('/home/will/Desktop/wgbs/revision.1/plots/%sgenes_spline_cor_aggregate_spar9.pdf', direction), width=5, height=5)
plot(cur_x, cur_y, pch=19, xlab='', ylab='Spearman\'s rho',  xaxt='n', bty='n', ylim=c(0.32,0.38), col='black')
lines(cur_x, cur_y)
abline(v=0, col='green', lty=2)
abline(v=1000, col='red', lty=2)
dev.off()


# Fit 1 spline for all Cpgs concatenated (1 spline for all genes) - negatively correlated
all_cpgs_list_neg <- list()
counter <- 1
direction <- 'neg'
library(GenomicRanges)
for(j in neg_sig_idx) {
  if(counter %% 1000 == 0) print(counter)
  curtab <- all_methyl_norm[[j]]
  curtab <- curtab[complete.cases(curtab),]
  
  #Order in terms of increasing genomic position (upstream -> downstream)
  curtab <- curtab[order(curtab$coordAdj, decreasing = F),] 
  curgene_metadata <- all_metadata[[j]]
  
  if(direction == 'pos') {
    # Only use positive correlations
    curtab <- curtab[which(curtab$spearmanRho >= 0 & curtab$spearman_p < 0.05),]
  } else if(direction == 'neg') {
    # Only use positive correlations
    curtab <- curtab[which(curtab$spearmanRho <= 0 & curtab$spearman_p < 0.05),]
  }
  
  all_cpgs_list_neg[[counter]] <- curtab
  counter <- counter + 1
}

all_cpgs_tab <- rbindlist(all_cpgs_list_neg)
all_cpgs_tab <- all_cpgs_tab[order(all_cpgs_tab$coordAdj, decreasing=F),]
cur.methyl.smooth <- smooth.spline(all_cpgs_tab[,c('coordAdj', 'spearmanRho')], spar=0.9)
cur_x <- seq(-1000,2000, by=0.1)
cur_y <- predict(cur.methyl.smooth, cur_x)$y
pdf(sprintf('/home/will/Desktop/wgbs/revision.1/plots/%sgenes_spline_cor_aggregate_spar9.pdf', direction), width=5, height=5)
plot(cur_x, cur_y, pch=19, xlab='', ylab='Spearman\'s rho',  xaxt='n', bty='n', ylim=c(-0.36,-0.25), col='black')
abline(v=0, col='green', lty=2)
abline(v=1000, col='red', lty=2)
dev.off()