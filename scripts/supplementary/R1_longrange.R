library(gplots)

filepath <- 'C:/Users/Admin/Desktop/longrange/'

range <- 2
mincor <- 0.3
minmaxcor <- 0.3
mincor2 <- 0.1

maxna <- 0.2
filename <- '../methylseq/hmrseg/cgi.tsv'
cgi <- read.delim(filename,sep='\t',header=T,check.names=F)
cgiorder <- order(match(cgi$chr,chromosomes),cgi$start)
cgi <- cgi[cgiorder,]
rownames(cgi) <- 1:length(cgiorder)

rows2keep <- rowSums(is.na(cgi[,samples_wgbs]))<=maxna
cgi <- cgi[rows2keep,]
cgi$pos <- (cgi$start+cgi$end)/2



plot_hclust_cor <- function(matrix, filename, use_pdf=F) {
	color <- colorRampPalette(c("blue","white","red"))(n = 1000)

	if(use_pdf) {
		filename <- paste(filename, ".pdf", sep="")
		pdf(file=filename, onefile=FALSE)
	} else {
		filename <- paste(filename, ".png", sep="")
		png(filename, width=1000, height=1000)
	}

	heatmap.2(matrix, col=color, 
		Rowv=T, Colv=T,
		rowsep=0:dim(matrix)[1], colsep=0:dim(matrix)[2],
		#sepcolor="black", sepwidth=c(0.001,0.001),
		margins = c(20, 20), cexRow=1, cexCol=1,
		symm=T,symkey=T,symbreaks=T,
		trace='none', na.rm=TRUE)
	dev.off()
}


findpeaks <- function(values, threshold=mincor, minwidth=(range+range+1), within=range, minht=minmaxcor) {
	results <- data.frame()
	inpeak <- F
	currentpeak <- 0
	
	if(max(values,na.rm=T)<threshold) {
		return(NULL)
	}
	
	for(i in 1:length(values)) {
		stopifnot(!is.na(currentpeak))
		val <- values[i]
		if(!is.na(val) & val>=threshold) {
			if(!inpeak) { #new peak?
				inpeak <- T
				if(currentpeak==0) { #first peak
					currentpeak <- currentpeak+1
					results[currentpeak,'start'] <- i
					results[currentpeak,'max'] <- val
				} else if(i-results[currentpeak,'end']>within) { #new peak!
					currentpeak <- currentpeak+1
					results[currentpeak,'start'] <- i
					results[currentpeak,'max'] <- val
				} else {
					#keep current peak, and continue
					results[currentpeak,'end'] <- NA
				}
			} else { 
				#continue peak
				results[currentpeak,'max'] <- max(results[currentpeak,'max'],val)
			}
			if(i==length(values)) { #stop if last one
				results[currentpeak,'end'] <- i
				inpeak <- F
			}
		} else {
			if(inpeak) { #done with peak
				results[currentpeak,'end'] <- i-1
				inpeak <- F
			} else {
				#do nothing and keep going
			}
		}
		#print(results)
	}
	stopifnot(sum(is.na(results$start))==0)
	stopifnot(sum(is.na(results$end))==0)
	stopifnot(sum(is.na(results$max))==0)
	results$width <- results$end-results$start
	results <- results[results$width>=minwidth & results$max>=minht,]
	return(results)
}

results <- data.frame()
for(chr in chromosomes) {
	for(i in 1:2) {
		if(i==1) {
			rowscgi <- cgi$chr==chr & cgi$end < centromeresdf[chr,'start']
			rowsgenes <- ensembl2sym$chr==chr & ensembl2sym$end < centromeresdf[chr,'start']
			fullname <- paste(chr,'p',sep='_')
		} else {
			rowscgi <- cgi$chr==chr & cgi$start > centromeresdf[chr,'end']
			rowsgenes <- ensembl2sym$chr==chr & ensembl2sym$start > centromeresdf[chr,'end']
			fullname <- paste(chr,'q',sep='_')
		}
		print(fullname)
		
		if(sum(rowscgi)>5 && sum(rowsgenes)>5) {
			print('CGI cor')
			cgi_chr <- cgi[rowscgi,c('start','end','pos',samples_wgbs)]
			rownames(cgi_chr) <- paste('CGI',rownames(cgi_chr),sep='')
			cgi_chr <- cgi_chr[order(cgi_chr$pos),]
			tocorcgi <- t(as.matrix(cgi_chr[,samples_wgbs]))
			cormatcgi <- cor(tocorcgi,use='complete.obs',method='spearman')
			meancorcgi <- rep(NA,dim(cormatcgi)[1])
			for(k in (1+range):(dim(cormatcgi)[1])-range) {
				if(k%%1000==0) print(k)
				prevcor <- as.numeric(cormatcgi[k,(k+1):(k+range)])
				nextcor <- as.numeric(cormatcgi[(k+1):(k+range),k])
				meancorcgi[k] <- mean(c(prevcor,nextcor),na.rm=T)
			}
			cgi_peaks <- findpeaks(meancorcgi)
			cgi_peaks$start <- cgi_chr[cgi_peaks$start,'start']
			cgi_peaks$end <- cgi_chr[cgi_peaks$end,'end']
			
			print('TPM cor')
			genes_chr <- intersect(rownames(ensembl2sym)[rowsgenes],expressed)
			tpm_chr <- tpm[genes_chr,samples_wgbs]
			tpm_chr$start <- ensembl2sym[genes_chr,'start']
			tpm_chr$end <- ensembl2sym[genes_chr,'end']
			tpm_chr$pos <- (tpm_chr$start+tpm_chr$end)/2
			tpm_chr <- tpm_chr[order(tpm_chr$pos),]
			tocortpm <- t(as.matrix(tpm_chr[,samples_wgbs]))
			cormattpm <- cor(tocortpm,use='complete.obs',method='spearman')
			meancortpm <- rep(NA,dim(cormattpm)[1])
			for(k in (1+range):(dim(cormattpm)[1])-range) {
				if(k%%1000==0) print(k)
				prevcor <- as.numeric(cormattpm[k,(k+1):(k+range)])
				nextcor <- as.numeric(cormattpm[(k+1):(k+range),k])
				meancortpm[k] <- mean(c(prevcor,nextcor),na.rm=T)
			}
			tpm_peaks <- findpeaks(meancortpm)
			tpm_peaks$start <- tpm_chr[tpm_peaks$start,'start']
			tpm_peaks$end <- tpm_chr[tpm_peaks$end,'end']
			
			if(dim(cgi_peaks)[1]>0 && dim(tpm_peaks)[1]>0) {
				for(k in 1:dim(tpm_peaks)[1]) {
					tpm_start <- tpm_peaks[k,'start']
					tpm_end <- tpm_peaks[k,'end']
					rowsoverlap <- overlap(cgi_peaks$start,cgi_peaks$end,tpm_start,tpm_end)
					stopifnot(sum(is.na(rowsoverlap))==0)
					if(sum(rowsoverlap)>0) {
						print(k)
						print(tpm_peaks[k,])
						print(cgi_peaks[rowsoverlap,])
						
						chri <- chr
						starti_tpm <- tpm_peaks[k,'start']
						endi_tpm <- tpm_peaks[k,'end']
						rowsgene <- ensembl2sym$chr==chri & 
							overlap(ensembl2sym$start,ensembl2sym$end,starti_tpm,endi_tpm) &
							rownames(ensembl2sym) %in% expressed
						genes2cor <- rownames(ensembl2sym)[rowsgene]
						genecor <- tpm[genes2cor,samples_wgbs]
						rownames(genecor) <-  make.unique(ensembl2sym[rownames(genecor),'name'])
						CNi <- cnregion(chri,starti_tpm,endi_tpm)[samples_wgbs]
						tocor <- rbind(genecor,CN=CNi)
						
						rowscytoband <- cytobandsdf$chr==chri & overlap(cytobandsdf$start,cytobandsdf$end,starti_tpm,endi_tpm)
						cbi <- cytobandsdf[rowscytoband,'band']
						
						rownums <- (1:length(rowsoverlap))[rowsoverlap]
						found <- F
						for(rownum in rownums) {
							starti_cgi<- cgi_peaks[rownum,'start']
							endi_cgi <- cgi_peaks[rownum,'end']
							rowscgi <- cgi$chr==chri & overlap(cgi$start,cgi$end,starti_cgi,endi_cgi)
							cgicor <- cgi[rowscgi,samples_wgbs]
							rownames(cgicor) <- paste('CGI',rownames(cgicor),sep='')
							
							corcgitpm <- cor(t(genecor),t(cgicor),use='complete.obs',method='spearman')
							diag(corcgitpm) <- NA
							meancor <- mean(corcgitpm,na.rm=T)
							if(abs(meancor) >= mincor2) {
								found <- T
								tocor <- rbind(tocor,cgicor)
								resultsrow <- dim(results)[1]+1
								results[resultsrow,'chr'] <- chr
								results[resultsrow,'start'] <- min(starti_tpm,starti_cgi)
								results[resultsrow,'end'] <- max(endi_tpm,endi_cgi)
								results[resultsrow,'cytoband'] <- paste(cbi,collapse='|')
								results[resultsrow,'genes'] <- paste(ensembl2sym[genes2cor,'name'],collapse='|')
								results[resultsrow,'CGI'] <- paste(rownames(cgicor),collapse='|')
								if(meancor>=0) {
									results[resultsrow,'correlation'] <- 'positive'
								} else {
									results[resultsrow,'correlation'] <- 'negative'
								}
								results[resultsrow,'rho'] <- meancor
							}
						}
						
						if(found) {
							tocor <- t(tocor)
							cormat <- cor(tocor,use='complete.obs',method='spearman')
							
							print('plotting')
							filename <- paste(filepath,chri,'_',starti_tpm,'_',endi_tpm,sep='')
							plot_hclust_cor(cormat,filename)
						}
					}
				}
			}
		}
	}
}

filename <- paste(filepath,'longrange.xls',sep='')
write.table(results,filename,sep='\t',col.names=T,row.names=F)



