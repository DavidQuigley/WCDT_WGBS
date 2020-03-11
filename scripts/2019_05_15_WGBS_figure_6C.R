###################
# HELPER FUNCTIONS
###################
getDSSVals <- function(curmuts.gr, curdss.gr, window.sz) {
    # Identify mapping from each mutwindow segment to overlapping dss segments
    curwindows_overlap_hits <- findOverlaps(curmuts.gr, curdss.gr, ignore.strand=T)
    hits.tab <- as.data.frame(curwindows_overlap_hits)
    
    max(queryHits(curwindows_overlap_hits))
    max(subjectHits(curwindows_overlap_hits))
    
    
    dssmat <- matrix(0,nrow=length(curmuts.gr), ncol=window.sz)
    counter <- 1
    
    # Compute score as the sum of all HMR widths in each 5kb window
    for(j in unique(hits.tab$queryHits)) {
        if(counter %% 100 == 0) print(counter)
        curwin <- curmuts.gr[j,]@ranges
        curwin@start
        
        currows <- hits.tab[hits.tab$queryHits == j,'subjectHits']
        cur.dss.segs <- data.frame(curdss.gr[currows,]) 
        cur.dss.segs[,c('start', 'end')] <- cur.dss.segs[,c('start', 'end')] - curwin@start + 1
        # Truncate range to be constrained to window size
        cur.dss.segs$start[cur.dss.segs$start <= 0] <- 1
        cur.dss.segs$end[cur.dss.segs$end >= window.sz] <- window.sz
        for(k in seq_len(nrow(cur.dss.segs))) {
            dssmat[j, cur.dss.segs$start[k]:cur.dss.segs$end[k]] <- cur.dss.segs$diff.Methy[k]
        }
        
        #cur.dss.segs <- intersect(curmuts.gr[j,], curdss.gr)
        #data.frame(cur.dss.segs@ranges)- curwin@start + 1
        counter <- counter+1
    }
    return(dssmat)
}


getWindowedGR <- function(curmuts.gr, window.sz) {
    curmuts.df <- as.data.frame(curmuts.gr)
    #stopifnot(max(curmuts.df$width) <= 2500)
    
    # get middle of range for mutation / indel / segment
    midpt <- (curmuts.df$start + curmuts.df$end)/2
    curmuts.df$end <- midpt + window.sz/2
    curmuts.df$start <- midpt - (window.sz/2-1)
    curmuts.window.df <- makeGRangesFromDataFrame(curmuts.df)
    return(curmuts.window.df)
}


getSubsampleGR <- function(cur.gr, subsample.sz) {
    set.seed(1)
    if(length(cur.gr) <= subsample.sz) {
        return(cur.gr)
    }
    return(cur.gr[sample(1:length(cur.gr), subsample.sz),])
}

##############
# LOAD DATA
###############
#setwd('/Users/Will/Desktop/med_school_homebase/feng_lab/wgbs/scripts/from_george')
#source('/Users/Will/Desktop/med_school_homebase/feng_lab/wgbs/scripts/from_george/loadranges_will.R')

##############################
# Slice out areas around TFBS
##############################
window.sz <- 20000

curtest.tfbs.ar.gr <- tracks[['CRPC_ARE']]
curtest.tfbs.foxa1.gr <- tracks[['chip_FOXA1']]
curtest.tfbs.hoxb13.gr <- tracks[['chip_HOXB13']]
curtest.tfbs.erg.gr <- tracks[['chip_ERG']]

curtest.tfbs.ar.window.gr <- getWindowedGR(curtest.tfbs.ar.gr, window.sz)
curtest.tfbs.foxa1.window.gr <- getWindowedGR(curtest.tfbs.foxa1.gr, window.sz)
curtest.tfbs.hoxb13.window.gr <- getWindowedGR(curtest.tfbs.hoxb13.gr, window.sz)
curtest.tfbs.erg.window.gr <- getWindowedGR(curtest.tfbs.erg.gr, window.sz)

# Get DSS values at each mutation window
dssmat.tfbs.ar <- getDSSVals(curtest.tfbs.ar.window.gr, dss_metbp, window.sz)
dssmat.tfbs.foxa1 <- getDSSVals(curtest.tfbs.foxa1.window.gr, dss_metbp, window.sz)
dssmat.tfbs.hoxb13 <- getDSSVals(curtest.tfbs.hoxb13.window.gr, dss_metbp, window.sz)
dssmat.tfbs.erg <- getDSSVals(curtest.tfbs.erg.window.gr, dss_metbp, window.sz)

################
# Plot spline
################
# Save spline for each TFBS
dssmat.final <- cbind.data.frame(Position=1:window.sz - window.sz/2, DiffMethyl=colMeans(dssmat.tfbs.foxa1))
dssmat.final[,2] <- dssmat.final[,2] * 100 #convert decimal to percentage
dssmat.final.foxa1 <- dssmat.final
save(dssmat.final.foxa1, file=fn_spline_foxa1)

dssmat.final <- cbind.data.frame(Position=1:window.sz - window.sz/2, DiffMethyl=colMeans(dssmat.tfbs.hoxb13))
dssmat.final[,2] <- dssmat.final[,2] * 100 #convert decimal to percentage
dssmat.final.hoxb13 <- dssmat.final
save(dssmat.final.hoxb13, file=fn_spline_hoxb13)

dssmat.final <- cbind.data.frame(Position=1:window.sz - window.sz/2, DiffMethyl=colMeans(dssmat.tfbs.ar))
dssmat.final[,2] <- dssmat.final[,2] * 100 #convert decimal to percentage
dssmat.final.ar <- dssmat.final
save(dssmat.final.ar, file=fn_spline_ar)

dssmat.final <- cbind.data.frame(Position=1:window.sz - window.sz/2, DiffMethyl=colMeans(dssmat.tfbs.erg))
dssmat.final[,2] <- dssmat.final[,2] * 100 #convert decimal to percentage
dssmat.final.erg <- dssmat.final
save(dssmat.final.erg, file=fn_spline_erg)


# Plot FOXA1
cur.ddsplot <- ggplot(dssmat.final.foxa1, aes(x=Position, y=DiffMethyl)) +
    labs(x="Relative genomic position", y = "Mean differential methylation (%)") +
    geom_line(col='dodgerblue1', size=1.5) +
    geom_vline(xintercept=0, colour='black', linetype=2) +
    theme_classic() +
    theme(axis.title.y=element_text(size=8.5), axis.text.y=element_text(size=9.5))
cur.ddsplot
ggsave(filename=fn_figure6d_foxa1, plot=cur.ddsplot, width=4,h=3)

# Plot HOXB13
cur.ddsplot <- ggplot(dssmat.final.hoxb13, aes(x=Position, y=DiffMethyl)) +
    labs(x="Relative genomic position", y = "Mean differential methylation (%)") +
    geom_line(col='dodgerblue1', size=1.5) +
    geom_vline(xintercept=0, colour='black', linetype=2) +
    theme_classic() +
    theme(axis.title.y=element_text(size=8.5), axis.text.y=element_text(size=9.5))
cur.ddsplot
ggsave(filename=fn_figure6d_hoxb13, plot=cur.ddsplot, width=4,h=3)

# Plot AR
cur.ddsplot <- ggplot(dssmat.final.ar, aes(x=Position, y=DiffMethyl)) +
    labs(x="Relative genomic position", y = "Mean differential methylation (%)") +
    geom_line(col='dodgerblue1', size=1.5) +
    geom_vline(xintercept=0, colour='black', linetype=2) +
    theme_classic() +
    theme(axis.title.y=element_text(size=8.5), axis.text.y=element_text(size=9.5))
cur.ddsplot
ggsave(filename=fn_figure6d_ar, plot=cur.ddsplot, width=4,h=3)


# Plot ERG
cur.ddsplot <- ggplot(dssmat.final.erg, aes(x=Position, y=DiffMethyl)) +
    labs(x="Relative genomic position", y = "Mean differential methylation (%)") +
    geom_line(col='dodgerblue1', size=1.5) +
    geom_vline(xintercept=0, colour='black', linetype=2) +
    theme_classic() +
    theme(axis.title.y=element_text(size=8.5), axis.text.y=element_text(size=9.5))
cur.ddsplot
ggsave(filename=fn_figure6d_erg, plot=cur.ddsplot, width=4,h=3)


################################
# H3K27ac (requires more memory)
################################
curtest.tfbs.gr <- tracks[['pca100_H3K27ac']]
curtest.tfbs.gr <- curtest.tfbs.gr[1:length(tracks[['pca100_H3K27ac']]),]
curtest.tfbs.window.gr <- getWindowedGR(curtest.tfbs.gr, window.sz)
# Get DSS values at each mutation window
dssmat.curtest.tfbs <- getDSSVals(curtest.tfbs.window.gr, dss_metbp, window.sz)
dssmat.final <- cbind.data.frame(Position=1:window.sz - window.sz/2, DiffMethyl=colMeans(dssmat.curtest.tfbs))
dssmat.final[,2] <- dssmat.final[,2] * 100 #convert decimal to percentage
dssmat.final.h3k27ac <- dssmat.final
save(dssmat.final.h3k27ac, file=fn_spline_h3k27ac)

library(ggplot2)
cur.ddsplot <- ggplot(dssmat.final.h3k27ac, aes(x=Position, y=DiffMethyl)) +
    labs(x="Relative genomic position", y = "Mean differential methylation (%)") +
    geom_line(col='dodgerblue1', size=1.5) +
    geom_vline(xintercept=0, colour='black', linetype=2) +
    theme_classic() +
    theme(axis.title.y=element_text(size=8.5), axis.text.y=element_text(size=9.5))
cur.ddsplot
ggsave(filename=sprintf('/Users/Will/Desktop/med_school_homebase/feng_lab/wgbs/results/global.analyses/dss.tfbs/2019-apr12/H3K27ac.dss.pdf'), plot=cur.ddsplot, width=4,h=3)


