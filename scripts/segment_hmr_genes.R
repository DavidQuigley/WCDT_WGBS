#qmod -d all.q@node1
library(doParallel)
#source('helper.R')
#source('helper_peaks.R')

window <- 100
maxwidth <- 10000
cutoff <- 5
within <- 100

filepath <- '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/gene_output/gene_output/'
outputfile <- '/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/recurrent_HMRs/hmrseg_gene_pearson_spearman.txt'
source('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/scripts/2019_05_15_prepare_environment.r')
# Set up cores for parallel processing
cl <- makeCluster( 60, outfile="")
registerDoParallel(cl)

filenames <- list.files(filepath, pattern='.txt')

samples_cmp = c(
    "DTB-018-BL", "DTB-021-BL", "DTB-022-BL",  "DTB-023-BL",
    "DTB-024-PRO", "DTB-037-BL","DTB-042-BL","DTB-074-BL",
    "DTB-089-BL", "DTB-094-BL", "DTB-102-PRO", "DTB-112-BL",
    "DTB-124-BL", "DTB-137-PRO", "DTB-141-BL", "DTB-151-BL",
    "DTB-188-BL", "DTB-190-BL", "DTB-194-PRO", "DTB-202-BL", 
    "DTB-204-BL", "DTB-260-BL")

segments <- foreach( i = 1:length(filenames), .combine=rbind) %dopar% {
	library(stringr)
	library(GenomicRanges)
	
	split <- strsplit( filenames[i], '_', fixed=TRUE)[[1]]
	ensemblid <- split[1]
	rest <- paste( split[-1], collapse='_')
	split <- strsplit(rest, '_chr', fixed=TRUE)[[1]]
	name <- split[1]
	chr <- paste('chr', strsplit(split[2],'_', fixed=TRUE )[[1]][1],sep='')
	stopifnot(chr %in% chromosomes)
	
	filenamefull <- paste(filepath, filenames[i],sep='')
	methyl <- read.delim(filenamefull,sep='',row.names=1,header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
	colnames(methyl) <- toupper(str_split_fixed(colnames(methyl),'\\.',2)[,1])
	stopifnot(dim(methyl)[1]>0)
	is_cmp = factor( colnames(methyl) %in% samples_cmp )
	methylpos <- as.numeric(rownames(methyl))
	start <- floor(min(methylpos)/window)*window
	end <- ceiling(max(methylpos)/window)*window
	stopifnot(start<end)
	
	hmrbin <- data.frame(bins=1:((end-start)/window))
	hmrbin$chr <- chr
	hmrbin$start <- start+((hmrbin$bins-1)*window)
	hmrbin$end <- start+(hmrbin$bins*window)
	bins <- makeGRangesFromDataFrame(hmrbin)
	for(j in 1:length(hmr)) {
		sample <- names(hmr)[[j]]
		rowshmrj <- width(hmr[[j]]) <= maxwidth
		hmrbin[,sample] <- countOverlaps(bins,hmr[[j]][rowshmrj])>0
	}
	hmrbin$sum <- rowSums(hmrbin[,(1:length(hmr))+(dim(hmrbin)[2]-length(hmr))])
	hmrbin <- hmrbin[,c('chr','start','end','sum')]
	colnames(hmrbin) <- c('chrom','start','end','value')
	peaks <- findpeaks(hmrbin,cutoff,minwidth=window,within=within,minsize=cutoff)
	if( !is.null(peaks) ) {
	    tryCatch(
	      {
    		numpeaks <- dim(peaks)[1]
	    	#print( paste( 'found', numpeaks, 'peaks in', filenamefull))
		    peaks$chr <- chr
    		peaks$gene_id <- ensemblid
	    	peaks$gene <- name
		
    		expr <- as.numeric(tpm[ensemblid,colnames(methyl)])
	    	for(j in 1:numpeaks) {
		    	methylrows <- methylpos >= peaks[j,'start'] & methylpos <= peaks[j,'end']
			    if( sum(methylrows) > 1) {
				    methylslice <- methyl[methylrows,]
    				avgmethyl <- colMeans(methylslice,na.rm=T)
	    			both_present = !is.na( avgmethyl ) & !is.na( expr ) & !is.infinite( avgmethyl ) & !is.infinite( expr )
		    		if( sum( both_present ) >1 ){
    		    		cor <- cor.test(avgmethyl,expr,method='spearman')
	    		    	peaks[j,'avg_cor'] <- cor$estimate
		    		    peaks[j,'avg_p'] <- cor$p.value
    			    	cor <- cor.test(avgmethyl,expr)
	    			    peaks[j,'avg_cor_pearson'] <- cor$estimate
    	    			peaks[j,'avg_p_pearson'] <- cor$p.value
    	    			m1 = aov( expr~avgmethyl+is_cmp)
    	    			m2 = aov( expr~avgmethyl*is_cmp)
    	    			peaks[j, 'p_interact_cmp'] = anova(m1, m2)[2,6]
    		    	}
    			}
	    	}
	      
		    if(!is.null(peaks$avg_cor) & !is.null(peaks$avg_p)) {
			    peaks <- peaks[,c('gene_id','gene','chr','start','end','max','avg_cor','avg_p','avg_cor_pearson','avg_p_pearson','p_interact_cmp')]
    			peaks
	    	}
	      }
	    , error=function(c){print( paste( "ERROR in ", filenamefull ) )} )
	}
}
stopCluster(cl)

segments$fdr <- p.adjust(segments$avg_p, method='fdr')
segments$fdr_pearson <- p.adjust(segments$avg_p_pearson, method='fdr')

segments$avg_cor = signif( as.numeric(segments$avg_cor), 3)
segments$avg_p = signif( as.numeric(segments$avg_p), 3)
segments$avg_cor_pearson = signif( as.numeric(segments$avg_cor_pearson), 3)
segments$avg_p_pearson = signif( as.numeric(segments$avg_p_pearson), 3)

print(sum(segments$fdr<=0.05,na.rm=T)/dim(segments)[1])

write.table(segments,outputfile,row.names=F,col.names=T,quote=F,sep='\t')

#qmod -e all.q@node1
