###################
# Helper functions
####################

# normalize gene coordinates to positions relative to TSS and gene end
normalize_coord = function(curpos, gene.start, gene.end) {
    if(curpos < gene.start) {
        # curpos (HMR segment boundary) is upstream of gene start
        curpos.norm = (curpos-gene.start)/upstream_window*upstream_sz
    } else if(curpos >= gene.start & curpos <= gene.end) {
        gene.len = gene.end-gene.start
        curpos.norm = (curpos-gene.start)/gene.len*gb_sz
    } else if(curpos > gene.end) {
        curpos.norm = gb_sz + (curpos-gene.end)/downstream_window*downstream_sz
    } else {
        stop('Something is wrong')
    }
    return(curpos.norm)
}

##############
# Parameters
#############
upstream_sz = 1000 #plot upstream region with width = 1000 normalized units
gb_sz = 1000 # plot gene body with width = 1000 normalized units
downstream_sz = 1000 # plot downstream region with width = 1000 normalized units
upstream_window = 10000 #upstream window is 10kb (unnormalized) units
downstream_window = 10000 #downstream window is 10kb unnormalized units

# Note: for the current implementation, upstream_sz and downstream_sz need to be the same

##################
# Load data
##################
hmr.tab = read.delim(file=fn_hmrseg_gene, sep='\t', header=T, stringsAsFactors = F)
ensembl2sym = read.delim(file=fn_ensembl2sym, sep='\t', header=T, stringsAsFactors=F)
expressed.genes = read.delim(file=fn_expressed_genes, sep='\t', header=F, stringsAsFactors=F)
expressed.genes = as.character(t(expressed.genes))

# Ensure that none of the HMRs are in expressed genes (should be true because George prefiltered those out)
stopifnot(length(which(!hmr.tab$gene_id %in% expressed.genes)) == 0)

# Rename columns to be more specific
colnames(hmr.tab)[match(c('chr', 'start', 'end'), colnames(hmr.tab))] = c('hmr.chr', 'hmr.start', 'hmr.end')

# Get map from hmr rows to gene metadata rows
hmr2metadata_map = match(hmr.tab$gene_id, ensembl2sym$ensembl_ID)
gene.metadata  = ensembl2sym[hmr2metadata_map,c('chr', 'start', 'end', 'strand')]
colnames(gene.metadata) = c('gene.chr', 'gene.start', 'gene.end', 'strand')

# Attach gene metadata to each HMR
hmr.tab.annot = cbind.data.frame(hmr.tab, gene.metadata)
# Ensure all annotations are correct (i.e. all HMRs assigned to a gene are within 10kb window of gene body)
stopifnot(identical(hmr.tab.annot$hmr.chr, hmr.tab.annot$gene.chr))
stopifnot(hmr.tab.annot$gene.start - hmr.tab.annot$hmr.end <= 10000)
stopifnot(hmr.tab.annot$hmr.start - hmr.tab.annot$gene.end <= 10000)


###################################
# Analysis and plotting code below
###################################

# Get normalized coordinates for each HMR segment
coords.norm = t(apply(hmr.tab.annot, 1, function(currow) {
    start.coord.norm = normalize_coord(as.numeric(currow['hmr.start']), as.numeric(currow['gene.start']), as.numeric(currow['gene.end']))
    end.coord.norm = normalize_coord(as.numeric(currow['hmr.end']), as.numeric(currow['gene.start']), as.numeric(currow['gene.end']))
    return(c(start.coord.norm, end.coord.norm))
}))
colnames(coords.norm) = c('hmr.start.norm', 'hmr.end.norm')

stopifnot(nrow(coords.norm) == nrow(hmr.tab.annot))
print(min(coords.norm[,'hmr.start.norm'])) #should be roughly -1000
print(max(coords.norm[,'hmr.end.norm'])) # should be roughly 2000

# Truncate segment for plotting purposes
coords.norm[coords.norm[,'hmr.start.norm'] < -upstream_sz, 'hmr.start.norm'] = -upstream_sz
coords.norm[coords.norm[,'hmr.end.norm'] > gb_sz + downstream_sz, 'hmr.end.norm'] = gb_sz + downstream_sz

print(min(coords.norm[,'hmr.start.norm'])) #should be >= -1000
print(max(coords.norm[,'hmr.end.norm'])) # should be <= 2000

hmr.tab.annot.final = cbind.data.frame(hmr.tab.annot, coords.norm)

# Take into account strandedness (i.e. the fact that gene end < gene start for these genes)
negsense.idx = which(hmr.tab.annot.final$strand == '-')
hmr.tab.annot.final[negsense.idx,c('hmr.start.norm', 'hmr.end.norm')] = (-1*hmr.tab.annot.final[negsense.idx,c('hmr.end.norm', 'hmr.start.norm')])+gb_sz 

# Keep only HMRs with significant FDRs (i.e., significantly associated with expression)
hmr.tab.annot.final.sig = hmr.tab.annot.final[which(hmr.tab.annot.final$fdr < 0.05),]

# Separate into positive and negative EMRs
hmr.tab.annot.final.sig.pos = hmr.tab.annot.final.sig[hmr.tab.annot.final.sig$avg_cor > 0,]
hmr.tab.annot.final.sig.neg = hmr.tab.annot.final.sig[hmr.tab.annot.final.sig$avg_cor < 0,]


# Show the degree of correlation in negative HMRs
hmr.tab.annot.final.sig.neg.pts.list = apply(hmr.tab.annot.final.sig.neg, 1, function(currow) {
    x = floor(as.numeric(currow['hmr.start.norm'])):ceiling(as.numeric(currow['hmr.end.norm']))
    y = as.numeric(currow['avg_cor'])
    return(cbind.data.frame(x,y))
})

hmr.tab.annot.final.sig.neg.pts = rbindlist(hmr.tab.annot.final.sig.neg.pts.list)

# Generate point estimate for mean correlation at every 1/1000th of the gene body and at every10 bases of upstream / downstream regions
x_norm = unique(hmr.tab.annot.final.sig.neg.pts$x)
x_norm = x_norm[order(x_norm)]
y_norm = rep(NA, length(x_norm))
for(i in 1:length(x_norm)) {
    curx = x_norm[i]
    y_norm[i] = mean(hmr.tab.annot.final.sig.neg.pts[hmr.tab.annot.final.sig.neg.pts$x == curx,y])
}
plot(x_norm, y_norm, pch=19)
data_neg = data.frame(xcoords=x_norm, cor_avg=y_norm)

# Generate smoothed plot for final figure
curplot = ggplot(data_neg, aes(x=xcoords, y=cor_avg)) +
    coord_cartesian(ylim=c(-0.42,-0.32)) +
    #geom_point(size=1, colour='black') +
    geom_smooth(col='black')+
    labs(x="Relative genomic position", y = "Spearman's rho") +
    theme_classic()+
    theme(axis.text.x=element_blank(),
          axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank())

curplot = curplot + geom_vline(xintercept=0, colour='green', linetype = "dashed") 
curplot = curplot + geom_vline(xintercept=gb_sz, colour='red', linetype = "dashed") 
ggsave(filename=fn_figure_neg_ehmr_cor_allgenes, plot=curplot, width=8, height=5, device='pdf')

#################################################################
# Plot density of negative eHMR segments relative to gene start and end
#########################################################

# Negative eHMRs
curhist = hist(hmr.tab.annot.final.sig.neg.pts$x, 3000, xlab='Relative genomic position', ylab='Density', main='')
xcoords = curhist$mids
ycoords = curhist$counts
data_neg_density = cbind.data.frame(pos=xcoords, density=ycoords)

# # Histogram (line plot)
# curplot = ggplot(data_neg_density, aes(x=pos, y=density)) +
#   geom_line() +
#   labs(x="Relative genomic position", y = "Density") +
#   theme_classic() 

# Histogram (shaded in)
curplot = ggplot(hmr.tab.annot.final.sig.neg.pts, aes(x=x)) + 
    geom_histogram(binwidth=1, color='black', fill='black') +
    labs(x="Relative genomic position", y = "Density") +
    ylim(0,6000) +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank())

curplot = curplot + geom_vline(xintercept=0, colour='green', linetype = "dashed") 
curplot = curplot + geom_vline(xintercept=gb_sz, colour='red', linetype = "dashed") 
ggsave(filename=fn_figure_neg_ehmr_cordensity_allgenes, plot=curplot, width=8, height=5, device='pdf')


#####################################
# Make similar plots for positive eHMRs
####################################
# Show the degree of correlation in negative HMRs
hmr.tab.annot.final.sig.pos.pts.list = apply(hmr.tab.annot.final.sig.pos, 1, function(currow) {
    x = floor(as.numeric(currow['hmr.start.norm'])):ceiling(as.numeric(currow['hmr.end.norm']))
    y = as.numeric(currow['avg_cor'])
    return(cbind.data.frame(x,y))
})
library(data.table)
hmr.tab.annot.final.sig.pos.pts = rbindlist(hmr.tab.annot.final.sig.pos.pts.list)

# Generate point estimate for mean correlation at every 1/1000 of the gene body and at every 10 bases of upstream / downstream regions
x_norm = unique(hmr.tab.annot.final.sig.pos.pts$x)
x_norm = x_norm[order(x_norm)]
y_norm = rep(NA, length(x_norm))
for(i in 1:length(x_norm)) {
    curx = x_norm[i]
    y_norm[i] = mean(hmr.tab.annot.final.sig.pos.pts[hmr.tab.annot.final.sig.pos.pts$x == curx,y])
}
#plot(x_norm, y_norm, pch=19)
data_pos = data.frame(xcoords=x_norm, cor_avg=y_norm)

# Generate smoothed plot for final figure
curplot = ggplot(data_pos, aes(x=xcoords, y=cor_avg)) +
    #geom_point(size=1, colour='black') +
    geom_smooth(col='black',se=F, data=data_pos)+
    coord_cartesian(ylim = c(0.34,0.38)) +
    labs(x="Relative genomic position", y = "Spearman's rho") +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank())

curplot = curplot + geom_vline(xintercept=0, colour='green', linetype = "dashed") 
curplot = curplot + geom_vline(xintercept=gb_sz, colour='red', linetype = "dashed") 
curplot
ggsave(filename=fn_figure_pos_ehmr_cor_allgenes, plot=curplot, width=8, height=5, device='pdf')


#################################################################
# Plot density of positive eHMR segments relative to gene start and end
#########################################################
curhist = hist(hmr.tab.annot.final.sig.pos.pts$x, 3000, xlab='Relative genomic position', ylab='Density', main='', ylim=c(0,5500))
xcoords = curhist$mids
ycoords = curhist$counts
data_pos_density = cbind.data.frame(pos=xcoords, density=ycoords)

# #Histogram (line plot)
# curplot = ggplot(data_pos_density, aes(x=pos, y=density)) +
#   #geom_point(size=1, colour='black') +
#   geom_line() +
#   labs(x="Relative genomic position", y = "Density") +
#   theme_classic() 

# Histogram (shaded in)
curplot = ggplot(hmr.tab.annot.final.sig.pos.pts, aes(x=x)) + 
    geom_histogram(binwidth=1, color='black', fill='black') +
    labs(x="Relative genomic position", y = "Density") +
    ylim(0,6000) +
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.line.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank())


curplot = curplot + geom_vline(xintercept=0, colour='green', linetype = "dashed") 
curplot = curplot + geom_vline(xintercept=gb_sz, colour='red', linetype = "dashed") 
curplot


ggsave(filename=fn_figure_pos_ehmr_cordensity_allgenes, plot=curplot, width=8, height=5, device='pdf')

