##################
# Helper function
###################
extract_metadata = function(fname) {
    retval = list()
    retval$ensembl_name = strsplit(fname, '_')[[1]][1]
    retval$hgnc_name = strsplit(fname, '_')[[1]][2]
    retval$chr = strsplit(fname, '_')[[1]][3]
    
    retval$start = strsplit(fname, '_')[[1]][4]
    retval$end = strsplit(fname, '_')[[1]][5]
    retval$strand = substr(strsplit(fname, '_')[[1]][6],1,1)
    return(retval)
}

#################################
# Paths and sample annotations
#################################
incl.nls = TRUE
region.opts = c('full', 'upstream', 'promoter', 'body', 'downstream')
curdate = '2019-08-21'
tscnc_pids = c('DTB.003.BL','DTB.036.BL', 'DTB.040.BL', 'DTB.135.PRO', 'DTB.205.BL')

# Get patient identifiers
benign_pids = c('Bis159AT', 'Bis165AT', 'Bis171AT', 'Bis49AT')
localized_pids = c('Bis158T', 'Bis159T', 'Bis165T', 'Bis171T', 'Bis49T')

###################
# Perform analysis
###################
# Set directory to be where gene-level correlation files are
# Zoomed out AR file

# Methylation data for subset of important genes
filenames = 'ENSG00000169083.16_AR_chrX_67544032_67730619_+.txt'
fn_wide_normal = paste(dir_normal_wide_slice, filenames[1], sep='/')
fn_wide_mcrpc = paste(dir_mcrpc_wide_slice, filenames[1], sep='/')

cur_region = region.opts[1]

###################
# GENERATE PLOT
###################

methylation.tab = read.delim(file=fn_wide_mcrpc, sep='\t', header=T, stringsAsFactors = F)
colnames(methylation.tab)[1] = 'POS'
curgene.data = extract_metadata(filenames[1])

# Plot normals alongside methylation
if(incl.nls) {
    methyl.tab.upitt = read.delim(file=fn_wide_normal, sep='\t', header=T, stringsAsFactors = F)
    colnames(methyl.tab.upitt) = sapply(colnames(methyl.tab.upitt), function(x) {
        strsplit(x, '_')[[1]][1]
    })
    colnames(methyl.tab.upitt)[1] = 'POS'
    
    benign_idx = which(colnames(methyl.tab.upitt) %in% benign_pids)
    localized_idx = which(colnames(methyl.tab.upitt) %in% localized_pids)
}



# Special upstream/downstream window for AR
if(curgene.data$hgnc_name %in% c('AR', 'ROBO2')) {
    methylation.tab = methylation.tab[methylation.tab$POS >= as.numeric(curgene.data$start) - 2000000 & methylation.tab$POS <= as.numeric(curgene.data$end) + 800000,]
    if(incl.nls) {
        methyl.tab.upitt = methyl.tab.upitt[methyl.tab.upitt$POS >= as.numeric(curgene.data$start) - 2000000 & methyl.tab.upitt$POS <= as.numeric(curgene.data$end) + 800000,]
    }
}

stopifnot(curgene.data$strand %in% c('+', '-'))
if(curgene.data$strand == '+') {
    gene_start = as.numeric(curgene.data$start)
    gene_end = as.numeric(curgene.data$end)
} else {
    gene_start = as.numeric(curgene.data$end) #
    gene_end = as.numeric(curgene.data$start)
}

minpos.x = min(methylation.tab$POS)
maxpos.x = max(methylation.tab$POS)

# Get mean methylation value for all tumors at each genomic coordinate
# Divide into t-SCNC  and non-tscnc
tscnc_idx = match(tscnc_pids, colnames(methylation.tab))
tscnc_idx = tscnc_idx[!is.na(tscnc_idx)]
adeno_idx  = setdiff(2:101, tscnc_idx)
group1.means = apply(methylation.tab[,tscnc_idx], 1, mean, na.rm=T)
group2.means = apply(methylation.tab[,adeno_idx], 1, mean, na.rm=T)
methyl.mean.tab.g1 = cbind.data.frame(POS=methylation.tab$POS, methyl.pct=group1.means)
methyl.mean.tab.g2 = cbind.data.frame(POS=methylation.tab$POS, methyl.pct=group2.means)
# include benign prostate and localized PCa
if(incl.nls) {
    group3.means = apply(methyl.tab.upitt[,localized_idx], 1, mean, na.rm=T)
    methyl.mean.tab.g3 = cbind.data.frame(POS=methyl.tab.upitt$POS, methyl.pct=group3.means)
    group4.means = apply(methyl.tab.upitt[,benign_idx], 1, mean, na.rm=T)
    methyl.mean.tab.g4 = cbind.data.frame(POS=methyl.tab.upitt$POS, methyl.pct=group4.means)
}

linecols = colorRampPalette(brewer.pal(9,'Blues'))(12)
cur.cmap = linecols[c(5,7,9,11)]
names(cur.cmap) = c("Benign Prostate", "Localized PCa", "mCRPC Adeno", 't-SCNC')

# Draw splines
cur.methylplot = ggplot(methyl.mean.tab.g2, aes(x=POS, y=methyl.pct, color='grey')) +
    labs(x="", y = "", color = "Methylation value") +
    ylim(0, 100) +
    labs(x="Genomic coordinate", y = "Methylation percent (%)") +
    # Draw smoothed fit on average
    geom_smooth(data=methyl.mean.tab.g1, method='loess', span=0.1, se=F, n=80, fill='grey20', aes(colour="t-SCNC")) + #n= window size (default 80), span = jaggedness
    geom_smooth(data=methyl.mean.tab.g2,  method='loess', span=0.1, se=F, n=80, fill='grey20', aes(colour="mCRPC Adeno")) +
    geom_smooth(data=methyl.mean.tab.g3, method='loess', span=0.1, se=F, n=80, fill='grey20', aes(colour="Localized PCa")) + 
    geom_smooth(data=methyl.mean.tab.g4, method='loess', span=0.1, se=F, n=80, fill='grey20',  aes(colour="Benign Prostate")) +
    scale_colour_manual(name="Sample type",breaks=c("Benign Prostate", "Localized PCa", "mCRPC Adeno", 't-SCNC'),  values=cur.cmap) +
    theme_classic() +
    theme(axis.title.y=element_text(size=8.5), axis.text.y=element_text(size=9.5), legend.position="bottom")

# Mark promoter 
if(cur_region %in% c('promoter', 'body', 'upstream', 'full')) {
    cur.methylplot = cur.methylplot + geom_vline(xintercept=gene_start, colour='green', linetype = "dashed") 
}

# Mark terminator
if(cur_region %in% c('body', 'downstream', 'full')) {
    cur.methylplot = cur.methylplot + geom_vline(xintercept=gene_end, colour='red', linetype = "dashed") 
}
# Plot AR enhancer (coordinates from Viswanathan et al, Suppl. Table 7)
# hg19 coordinates reported were 66115000-66130000; hg38 coordinates from Liftover = 66895158-66910158
if(curgene.data$hgnc_name == 'AR' && cur_region %in% c('body', 'upstream', 'full')) {
    cur.methylplot = cur.methylplot + geom_vline(xintercept=66895158, colour='black', linetype = "dashed") 
    cur.methylplot = cur.methylplot + geom_vline(xintercept=66910158, colour='black', linetype = "dashed") 
}

# Save plot
ggsave(filename=fn_figure2b, plot=cur.methylplot, width=10, height=5)



