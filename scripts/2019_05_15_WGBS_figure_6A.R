###################
# Helper functions
####################

# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='', curwidth=10) {
    scale = (length(lut)-1)/(max-min)
    
    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', main=title, ylab = '')
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min
        rect(0,y,curwidth,y+1/scale, col=lut[i], border=lut[i])
    }
    mtext(text = "Differential methylation (%)",
          side = 2, #side 2 = left
          line = 2.5)
}


get_pid <- function(fname) {
    fname <- as.character(fname)
    return(paste(strsplit(fname,'-')[[1]][1], strsplit(fname,'-')[[1]][2], sep='-'))
}

get_sname <- function(fname) {
    fname <- as.character(fname)
    return(strsplit(fname,'_')[[1]][1])
}

# Helper function for plotting
ctscolor <- function(vect,min_range=NA,max_range=NA,colors=c('white',rep('black',5)),scale=1) {
    numbers <- round(vect*scale)
    if(is.na(min_range)) {
        min_range <- min(numbers, na.rm=T)
    }
    if(is.na(max_range)) {
        max_range <- max(numbers, na.rm=T)
    }
    number_colors <- colorRampPalette(colors)(n=(max_range-min_range)+1)
    names(number_colors) <- paste(min_range:max_range)
    colors_out <- number_colors[paste(numbers)]
    stopifnot(length(colors_out)==length(vect))
    return(colors_out)
}

## Plot ideogram with 2 tracks + horizontal barplot showing mutation frequency
plotideogram_2col_mutdensity <- function(filename, background=NULL, left=NULL, mid=NULL, right=NULL, 
                                         barwidth=10000, genes.highlight=NULL, mut.density=NULL, 
                                         flank.width=0.25, left.margin = 0.1) {
    ymax <- 250000000/barwidth
    
    numcol <- 2
    methyl.width <- (1 - flank.width - left.margin) / numcol
    
    if(strsplit(filename, '\\.')[[1]][length(strsplit(fname,'\\.')[[1]])] == 'pdf') {
        pdf(filename, height=6, width=12)
    } else {
        png(filename, height=6, width=12, res=600, units='in')
    }
    layout(matrix(1,1,1))
    par(mar=c(2,4,1,1))
    chr_width <- 1
    
    plot( -100, -100, xlim=c(1, 25.8), ylim=c(0, ymax), axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
    for(i in 1:24 ) {
        rect( i + left.margin, 0, i+chr_width - flank.width, 1 + (hg38[i,'length'] / barwidth), col='white' ) # Note that the 0.2 is to make the width narrower
    }
    axis(2, at=seq(from=0, to=ymax, by=ymax/25), labels= seq(from=0, to=ymax, by=ymax/25) / 100, las=1)
    axis(1, at=(2:25)-0.5, labels = c(1:22, 'X', 'Y'),las=1)
    
    if(!is.null(left)) {
        #stopifnot(sum(colnames(left)!=c('chr','start','end','value'))==0)
        left$value <- (left$value/max(left$value)) * methyl.width
    }
    if(!is.null(mid)) {
        #stopifnot(sum(colnames(mid)!=c('chr','start','end','value'))==0)
        mid$value <- (mid$value/max(mid$value)) * methyl.width
    }
    if(!is.null(right)) {
        #stopifnot(sum(colnames(right)!=c('chr','start','end','value'))==0)
        right$value <- (right$value/max(right$value)) * methyl.width
    }
    if(!is.null(mut.density)) {
        mut.density$mutd <- (mut.density$mutd/max(mut.density$mutd)) * flank.width * 1.5
        mut.density$mutd[mut.density$mutd >= flank.width * 1.5] <- flank.width * 1.5 # Truncate outliers
    }
    
    for(i in 1:24) {
        chri <- paste(i)
        if(i==23) {chri <- 'X'}
        else if(i==24) {chri <- 'Y'}
        print(paste('chromosome',chri))
        
        if(!is.null(background)) {
            #stopifnot(sum(colnames(background)!=c('chr','start','end','value'))==0)
            backgroundcolors <- ctscolor(background$value, 0, max(unlist(background)))
            rowsbackground <- background$chr==chri
            print(paste('Number of background bins =', sum(rowsbackground)))
            backgroundi <- background[rowsbackground,]
            backgroundcolorsi <- backgroundcolors[rowsbackground]
            for(j in 1:length(backgroundcolorsi)) {
                rowindex <- ceiling(backgroundstarti[j,'start'] / barwidth)
                rect( i+0.05, rowindex, i + 0.95, rowindex+1, col=backgroundcolorsi[j], border=NA)
            }
        }
        
        # Centromeric region
        #centro_start <- hg38[i,'centromerStart']
        #centro_end <- hg38[i,'centromerEnd']
        
        centro_start <- centromere.coords[centromere.coords$chrom == paste('chr', chri, sep=''), 'start']
        centro_end <- centromere.coords[centromere.coords$chrom == paste('chr', chri, sep=''), 'end']
        print(centro_end)
        print(centro_start)
        print(barwidth)
        # Ensure that centromere is at least some visible width
        if(centro_end - centro_start <= 200 * barwidth) {
            centro_start <- centro_end - 200 * barwidth
        }
        
        bins_centromere <- floor( round( centro_start / barwidth, 0 ) ) :
            ceiling( round( centro_end / barwidth, 0 ) )
        top_cent <- max(bins_centromere)
        bottom_cent <- min(bins_centromere)
        
        
        if(!is.null(left) && !is.null(mid) && !is.null(mut.density)) {
            rowsleft <- left$chr==chri
            rowsmid <- mid$chr==chri
            rowsright <- right$chr==chri
            print(paste('Number of left bins =', sum(rowsleft)))
            print(paste('Number of mid bins =', sum(rowsmid)))
            #print(paste('Number of right bins =', sum(rowsright)))
            lefti <- left[rowsleft,]
            midi <- mid[rowsmid,]
            #righti <- right[rowsright,]
            genesi <- genes.highlight[genes.highlight$chr == chri,]
            mdi <- mut.density[mut.density$chr == chri,]
            
            #Don't draw mutation density in centromeric regions (not reliable)
            mdi[mdi$start < centro_start & mdi$end > centro_start,'end'] <- centro_start 
            mdi[mdi$start >= centro_start & mdi$end <= centro_end,'mutd'] <- 0
            mdi[mdi$start <= centro_end & mdi$end > centro_end,'start'] <- centro_end
            
            for(j in seq_len(nrow(lefti))){
                rowindex_start <- ceiling(lefti[j,'start'] / barwidth)
                # ensure not to draw anything in telomeric regions
                rowindex_start <- max(rowindex_start, ceiling(telomere.coords.5prime[telomere.coords.5prime$chrom == paste('chr', chri, sep=''), 'chromEnd']/barwidth))
                
                rowindex_end <- max(ceiling(lefti[j,'end'] / barwidth), rowindex_start+1)
                # ensure not to draw anything in telomeric regions
                rowindex_end <- min(rowindex_end, floor(telomere.coords.3prime[telomere.coords.3prime$chrom == paste('chr', chri, sep=''), 'chromStart']/barwidth))
                if(rowindex_end > rowindex_start) {
                    rect( i + left.margin, rowindex_start, i + left.margin + abs(lefti[j,'value']), rowindex_end, col=lefti[j,'color'], border=NA)
                }
            }
            
            for(j in seq_len(nrow(midi))){
                rowindex_start <- ceiling(midi[j,'start'] / barwidth)
                # ensure not to draw anything in telomeric regions
                rowindex_start <- max(rowindex_start, ceiling(telomere.coords.5prime[telomere.coords.5prime$chrom == paste('chr', chri, sep=''), 'chromEnd']/barwidth))
                
                rowindex_end <- max(ceiling(midi[j,'end'] / barwidth), rowindex_start+1)
                # ensure not to draw anything in telomeric regions
                rowindex_end <- min(rowindex_end, floor(telomere.coords.3prime[telomere.coords.3prime$chrom == paste('chr', chri, sep=''), 'chromStart']/barwidth))
                if(rowindex_end > rowindex_start) {
                    rect( i+ left.margin + methyl.width, rowindex_start, i + left.margin + methyl.width + abs(midi[j,'value']), rowindex_end, col=midi[j,'color'], border=NA)
                }
            }
            # Plot mutation density
            for(j in seq_len(nrow(mdi))) {
                rowindex_start <- ceiling(mdi[j,'start']/barwidth)
                # ensure not to draw anything in telomeric regions
                rowindex_start <- max(rowindex_start, ceiling(telomere.coords.5prime[telomere.coords.5prime$chrom == paste('chr', chri, sep=''), 'chromEnd']/barwidth))
                
                rowindex_end <- ceiling(mdi[j,'end']/barwidth)
                # ensure not to draw anything in telomeric regions
                rowindex_end <- min(rowindex_end, floor(telomere.coords.3prime[telomere.coords.3prime$chrom == paste('chr', chri, sep=''), 'chromStart']/barwidth))
                if(rowindex_end > rowindex_start) {
                    rect( i+ left.margin + 2*methyl.width, rowindex_start, i +left.margin + 2*methyl.width + mdi[j, 'mutd'], rowindex_end, col='lavenderblush', lwd=0.5)
                }
            }
        }
        
        # White out centromeric region
        rect(i+0.15,bottom_cent+15,i+0.85,top_cent-15, col="white", border=NA)
        
        #polygon( c( i, i+1, mean(c(i, i+1))),
        polygon( c( i+0.38, i+0.62, mean(c(i, i+1))),       # Make centromere triangle  
                 c( top_cent, top_cent, mean(c(top_cent, bottom_cent)) ), col="gold")
        #polygon( c( i, i+1, mean(c(i, i+1))),
        polygon( c( i+0.38, i+0.62, mean(c(i, i+1))),
                 c( bottom_cent, bottom_cent, mean(c(top_cent, bottom_cent)) ), col="gold")
        lines( c(i+left.margin,i+1-flank.width), c(top_cent, top_cent), col="black")
        lines( c(i+left.margin,i+1-flank.width), c(bottom_cent, bottom_cent), col="black")
    }
    dev.off()
}


####################
# Load dependencies
####################
library(rCGH) # Load hg38 chromosomal info
head(hg38)

# Load centromere and telomere info

centromere.coords = read.delim(file=fn_centromeres, sep='\t', header=TRUE,stringsAsFactors=FALSE)
telomere.coords = read.delim(file=fn_telomeres, sep='\t', header=TRUE,stringsAsFactors=FALSE)
telomere.coords.5prime = telomere.coords[seq(from=1,to=nrow(telomere.coords), by=2),]
telomere.coords.3prime = telomere.coords[seq(from=2,to=nrow(telomere.coords), by=2),]

###############
# LOAD DATA
##############

# Load CNA data
filenames_dss = c( fn_DSS_localized_benign, fn_DSS_adeno_localized)

dss_data <- list()
dss_data_gr <- list()
# TODO: loop through to load DSS data
for(fname in filenames_dss) {
    cur_data <- read.delim(fname, header=T, sep='\t', stringsAsFactors = F)
    cur_data <- cur_data[,c(1,2,3,8)]
    colnames(cur_data) <- c('Chr', 'Start', 'End', 'value')
    cur_data$Chr <- gsub('chr', '', cur_data$Chr)
    #cur_data <- as.matrix(cur_data)
    cur_data[,2] <- as.numeric(cur_data[,2])
    cur_data[,3] <- as.numeric(cur_data[,3])
    cur_data[,4] <- as.numeric(cur_data[,4])
    cur_data$Chr[cur_data$Chr == 23] <- 'X'
    cur_data$Chr[cur_data$Chr == 24] <- 'Y'
    dss_data[[fname]] <- cur_data
    
    cur_data$Chr <- paste('chr', cur_data$Chr, sep='')
    cur_data.gr <- makeGRangesFromDataFrame(cur_data, keep.extra.columns = T)
    dss_data_gr[[fname]] <- cur_data.gr
}

# remove samples not included in wgbs manuscript
samples.annovar <- unique(annovar_full$sample)
samples.rm <- c("DTB-053-RP", "DTB-265-BL", "DTB-193-BL") # Remove the samples not included in WGBS
samples.rm <- c(samples.rm, c('DTB-083-BL','DTB-126-BL')) # Remove hypermutated samples
samples.annovar.final <- samples.annovar[-which(samples.annovar %in% samples.rm)]
annovar_full <- annovar_full[-which(annovar_full$sample %in% samples.rm),]
length(unique(annovar_full$sample))
allmuts_gr <- makeGRangesFromDataFrame(annovar_full)

# Load blacklisted, poorly-mapped regions
blacklist.regions <- read.delim(file=fn_gap, sep='\t', header=T)
blacklist.regions <- blacklist.regions[,1:3]
colnames(blacklist.regions)[1:3] <- c('chrom', 'start', 'end')
blacklist.regions <- rbind.data.frame(blacklist.regions, centromere.coords) #exclude centromeres
colnames(blacklist.regions)[1:3] <- c('Chr', 'Start', 'End')
blacklist.regions.gr <- makeGRangesFromDataFrame(blacklist.regions)

# Exclude mutations called in these regions
allmuts_gr <- allmuts_gr[-queryHits(findOverlaps(allmuts_gr, blacklist.regions.gr, type='any')),]
# Exclude any DSS segment overlapping one of these regions
for(i in 1:length(dss_data_gr)) {
    curgr.data <- dss_data_gr[[i]]
    curgr.data <- curgr.data[-queryHits(findOverlaps(curgr.data, blacklist.regions.gr, type='any')),]
    dss_data_gr[[i]] <- curgr.data
}
for(curname in names(dss_data)) {
    curgr <- dss_data_gr[[curname]]
    curdf <- as.data.frame(curgr, stringsAsFactors=F)
    curdf <- curdf[,c(1,2,3,6)]
    colnames(curdf) <- c('chr', 'start', 'end', 'value')
    curdf[,'chr'] <- gsub('chr', '', curdf[,'chr'])
    dss_data[[curname]] <- curdf
}

#################
# GENERATE PLOT
#################
# Plot DSS data as ideogram
dss.threshold <- 0 # Only plot DSS segments with a diffMethyl >= this threshold

left <- dss_data[[1]]
mid <- dss_data[[2]]

# Get color intensities of DSS segments
left$value[abs(left$value) < dss.threshold] <- 0
mid$value[abs(mid$value) < dss.threshold] <- 0
left$intensity <- left$value
mid$intensity <- mid$value


blues.gradient <- colorRampPalette(c("white", "dodgerblue1"))(50)
reds.gradient <- colorRampPalette(c("white", "red"))(50)

# Colors corresponding to degree of differential methylation
left$color <- sapply(left$intensity, function(x) {
    curval <- round(abs(x)*100)
    if(curval == 0) curval <- 1
    if(curval > 50) curval <- 50
    if(x >= 0) reds.gradient[curval]
    else blues.gradient[curval]
})

mid$color <- sapply(mid$intensity, function(x) {
    curval <- round(abs(x)*100)
    if(curval == 0) curval <- 1
    if(curval > 50) curval <- 50
    if(x >= 0) reds.gradient[curval]
    else blues.gradient[curval]
})

# Set toggle as -1 (hypomethylated), 0 (not significant), or 1 (hypermethylated)
left$value = 0
left$value[dss_data[[1]]$value > dss.threshold] = 1
left$value[dss_data[[1]]$value < -dss.threshold] = -1

mid$value <- 0
mid$value[dss_data[[2]]$value > dss.threshold] = 1
mid$value[dss_data[[2]]$value < -dss.threshold] = -1

left <- left[left$value != 0,]
mid <- mid[mid$value != 0,]


colors.blue <- sapply(seq(0.8,0,by=-0.001), function(x) {
    adjustcolor( "dodgerblue1", alpha.f = x)
})
colors.red <- sapply(seq(0,0.8,by=0.001), function(x) {
    adjustcolor( "red", alpha.f = x)
})

pdf(file=fn_figure6a_colorbar, w=3,h=4)
color.bar(colorRampPalette(c("dodgerblue1", "white", "red"))(200), -50, 50, nticks=11, curwidth=1)
dev.off()

# Generate histogram bins tallying number of genes in each 1MB segment of the genome
chrlist <- list()
density_binsz <- 1000000
for(i in 1:nrow(hg38)) {
    currow <- hg38[i,]
    chrlen <- currow$length
    chrbreaks <- seq(from=1,to=chrlen, by=density_binsz)
    binstart <- chrbreaks
    binend <- chrbreaks + density_binsz - 1
    binend <- sapply(binend, min, chrlen)
    
    # Reformat chromosome names
    curchr <- currow$chrom
    if(curchr == 23) curchr <- 'X'
    if(curchr == 24) curchr <- 'Y'
    curchr <- paste('chr', curchr, sep='')
    
    chrtab <- cbind.data.frame(chr=curchr,start=binstart,end=binend)
    chrlist[[i]] <- chrtab
}

chrtab_full <- rbindlist(chrlist)
chrbins_full_gr <- makeGRangesFromDataFrame(chrtab_full)

#######################################
# Generate mutation density histograms
######################################
mutd_bins <- countOverlaps(chrbins_full_gr, allmuts_gr)
# Ensure counts correspond to number of genes in each bin (sanity check)
stopifnot(length(mutd_bins) == nrow(chrtab_full))
chrtab_full$mut.counts <- mutd_bins
chrtab_full$mutd <- chrtab_full$mut.counts / (chrtab_full$end - chrtab_full$start + 1)

## Try plotting mutational frequency as color column
head(chrtab_full)
chrtab_full$mutd.norm <- chrtab_full$mutd / max(chrtab_full$mutd)
chrtab_full$value <- 1
chrtab_full$color <- sapply(chrtab_full$mutd.norm, function(x) {
    adjustcolor( "lightsalmon1", alpha.f = abs(x) * 2)
})

#######################################################
# Plot 2-col ideogram with mutational frequency as bars
bar.ht <- 0.3 # height of histogram bars showing mutational frequency
chrtab_full$mutd[chrtab_full$mutd > 0.001] <- 0.001 #Truncate bins with average per-sample mut density of > 10 mutations / MB
chrtab_full$chr <- gsub('chr', '', chrtab_full$chr)

# Generally want to save as png because pdf has weird saturation issues in terms of plot colors
plotideogram_2col_mutdensity(fn_figure6a, background=NULL, left=left, mid=mid, barwidth=10000, genes.highlight=NULL, mut.density = chrtab_full, flank.width = bar.ht, left.margin = bar.ht)
