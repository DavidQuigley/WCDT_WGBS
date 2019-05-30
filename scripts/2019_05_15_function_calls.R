#######################################################################################################
# Function calls
#######################################################################################################
##################
### BEGIN HASH ###
##################

load.matrix = function( fn ){
    read.table(file=fn, sep="\t",row.names=1, check.names=F, header=T, na.strings=c('NA', '-99999'), stringsAsFactors=F)
}


get.split.col = function(v, string, col=0, last=F, first=F){
    if( last & first )
        stop("Cannot request both last and first column")
    if( col==0 & !last & !first)
        stop("Must request either a column by index, first, or last")
    
    for(i in 1:length(v)){
        x = strsplit( v[i], string, fixed=T)[[1]]
        if(last){
            v[i] = x[length(x)]
        }
        else if(first){
            v[i] = x[1]
        }
        else{
            v[i] = x[col]
        }
    }
    v
}



hsh_new = function(){
    new.env(hash=TRUE, parent=emptyenv()) 
}

hsh_in = function(H, key){
    exists(key, H)
}

hsh_get = function( H, key, na.if.not.found=F ){
    if( length(key)==1 ){
        if( na.if.not.found ){
            if( exists(key, H) )
                get(key, H)
            else
                NA
        }
        else{
            get(key, H)
        }
    }
    else{
        results = rep(0, length(key) )
        if( !na.if.not.found ){
            for(i in 1:length(key) ){
                if( exists(key[i], H) ){
                    results[i] = get(key[i], H )
                }
                else{
                    results[i] = NA
                }
            }
        }
        else{
            for(i in 1:length(key) ){
                results[i] = get(key[i], H )
            }
        }
        results
    }
}

hsh_set = function( H, key, value ){
    assign(key, value, envir=H)
}

hsh_keys = function( H ){
    return(sort(ls(H)))
}

hsh_keys_values = function( H ){
    keys = ls(H)
    values = hsh_get(H, keys)
    data.frame( keys, values, stringsAsFactors=F)
}

hsh_from_vectors = function( v1, v2=NULL ){
    # Create a hash from vectors v1, v2 with keys from v1 and values from v2
    # if v2 is null, set it to 1:length(v1)
    if( is.null(v2) )
        v2 = 1:length(v1)
    if( length(v1) != length(v2) ){
        stop("Length of v1 != length of v2")
    }
    H = hsh_new()
    for( i in 1:length(v1) ){
        hsh_set(H, v1[i], v2[i] )
    }
    H
}


hsh_unique_values = function( hsh, keys=NULL ){
    # given hsh where values are vectors, iterate through keys and 
    # identify the set of unique values
    all_values = hsh_new()
    if(is.null(keys)){
        keys = hsh_keys(hsh)
    }
    for(i in 1:length(keys)){
        if( hsh_in(hsh, keys[i] ) ){
            values = hsh_get( hsh, keys[i])
            for( j in 1:length(values)){
                hsh_set( all_values, values[j], 1)
            }    
        }
    }
    sort(hsh_keys(all_values))
}

match.idx = function(A, B, allow.multiple.B=F){
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

allele_effect=function( symbol, do_plot=FALSE, axis_override=NULL, order=NULL ){
    # missense mutations all come from curated data 
    # external references:
    #sample_ids gene_locs 
    #matrix_SV matrix_tpm matrix_germline matrix_CNA_int_ploidy
    #curated_fs curated_missense
    
    has_somatic_inactivation = rep(FALSE,N_SAMPLES)
    has_activating_missense = rep(FALSE, length(sample_ids))
    has_inactivating_missense = rep(FALSE, length(sample_ids))
    has_activating_sv = rep(FALSE, length(sample_ids))
    has_inactivating_sv = rep(FALSE, length(sample_ids))
    has_inactive = matrix_inactive[symbol,] != '.'
    has_germline = matrix_germline[symbol,] != '.'
    has_loh = as.character(matrix_CNA[symbol,]) == "LOH"
    
    chrom = gene_locs$chrom[ which( rownames(gene_locs)==symbol ) ]
    
    # copy number
    cna = as.numeric( matrix_CNA_int_ploidy[symbol,] )
    threshold_single =  LOSS_SINGLE_NONSEX   
    threshold_double =  LOSS_DOUBLE_NONSEX   
    threshold_amp = GAIN_NONSEX   
    cn.hi.threshold = GAIN_NONSEX
    if( chrom=="chrX" | chrom=="chrY" ){
        threshold_single =  LOSS_SEX
        threshold_double =  LOSS_SEX   
        threshold_amp = GAIN_SEX
        cn.hi.threshold = GAIN_SEX
    }
    
    # expression
    xp =  as.numeric( matrix_tpm[symbol,] )
    val_mad_all=mad(xp, na.rm=TRUE)
    val_median_all = median(xp, na.rm=TRUE)
    
    # mutations
    if( symbol %in% unique(curated_fs$symbol) ){
        samples_inactive = curated_fs$sample_id[curated_fs$gene==symbol]
        has_somatic_inactivation[ match.idx( sample_ids, samples_inactive)$idx.A ] = TRUE
    }
    
    if( symbol %in% unique( curated_missense$gene ) ){
        samples_active = curated_missense$sample_id[curated_missense$gene==symbol & curated_missense$consequence=="activate"] 
        has_activating_missense[ match.idx( sample_ids, samples_active)$idx.A ] = TRUE
        samples_inactive = curated_missense$sample_id[curated_missense$gene==symbol & curated_missense$consequence!="activate"] 
        has_inactivating_missense[ match.idx( sample_ids, samples_inactive)$idx.A ] = TRUE
        matrix_missense[symbol, !has_activating_missense & !has_inactivating_missense] = '.'
    }
    
    # sv
    idx = which( curated_sv$threeprime==symbol & curated_sv$consequence=="break" ) 
    if(length(idx)>0){
        for(i in 1:length(idx)){
            sample_id = curated_sv$sample[idx][i]
            idx_sample_id = which(sample_ids==sample_id)
            has_inactivating_sv[idx_sample_id] = TRUE
            matrix_SV[ symbol, idx_sample_id ] = curated_sv$mechanism[idx[i]]
        }
    }
    idx = which( curated_sv$threeprime==symbol & curated_sv$consequence=="activate" ) 
    if(length(idx)>0){
        for(i in 1:length(idx)){
            sample_id = curated_sv$sample[idx][i]
            idx_sample_id = which(sample_ids==sample_id)
            has_activating_sv[idx_sample_id] = TRUE
            matrix_SV[ symbol, idx_sample_id ] = curated_sv$mechanism[idx[i]]
        }
    }
    
    df = data.frame(
        activating_missense = has_activating_missense,
        inactivating_missense = has_inactivating_missense,
        nonsense = has_inactive,
        activating_sv = has_activating_sv,
        inactivating_sv = has_inactivating_sv,
        inactivating_germline = has_germline,
        CNA_1 = cna <= threshold_single,
        CNA_2 = cna <= threshold_double,
        CNA_amp = cna >= threshold_amp,
        LOH = has_loh,
        xp=xp,
        activating_cna = cna >= cn.hi.threshold & !is.na(xp) & xp >= 100
    )
    df$activated = df$activating_missense | df$activating_cna | df$activating_sv
    
    alleles_gone = rowSums( df[,c("inactivating_missense", "LOH", "nonsense", 
                                  "inactivating_sv", "inactivating_germline", 
                                  "CNA_1")])
    has_biallelic =    df$CNA_2 | alleles_gone > 1
    has_monoallelic = !df$CNA_2 & alleles_gone == 1
    
    # special case for CDK12 with two somatic mutations
    if( symbol=="CDK12" ){
        has_biallelic[ which(sample_ids=="DTB-214-BL")]=TRUE
        has_monoallelic[ which(sample_ids=="DTB-214-BL")]=FALSE
    }
    
    n_mono = sum( has_monoallelic )
    n_bi = sum( has_biallelic )
    n_amp = sum( df$CNA_amp )
    biallelic_because_sv = has_biallelic & 
        !df$CNA_2 & df$inactivating_sv  & 
        rowSums( df[,c("inactivating_missense", "LOH", "nonsense", 
                       "inactivating_germline", "CNA_1")] ) == 1
    df$bi = has_biallelic
    df$mono = has_monoallelic
    df$bi_from_sv = biallelic_because_sv
    n_bi_from_sv = sum( biallelic_because_sv )
    
    n_alleles = rep(0, N_SAMPLES)
    n_alleles[has_biallelic] = 2
    n_alleles[has_monoallelic] = 1
    logxp = log(xp+1)
    if(do_plot){
        if( is.null(axis_override)){
            boxplot( logxp[n_alleles==2], logxp[n_alleles==1], logxp[n_alleles==0], 
                     las=1, box.wex=0.25, 
                     cex.axis=0.75, lwd=0.5,
                     col=c("darkgrey", "lightgrey", "white"),
                     names=c("","",""))
        }else{
            boxplot( logxp[n_alleles==2], logxp[n_alleles==1], logxp[n_alleles==0], 
                     las=1, box.wex=0.25, 
                     cex.axis=0.75, lwd=0.5, axes=FALSE,
                     col=c("darkgrey", "lightgrey", "white"),
                     names=c("","",""), ylim=c( min(axis_override), max(axis_override) ))
            axis( 4, axis_override, las=1, cex.axis=0.75)
            axis( 1, c(1,2,3), labels = c("","","") )
            box()
        }
    }
    df$n_alleles_inactivated = n_alleles
    cor_estimate = NA
    cor_pval = NA
    if( sum( !is.na(xp) ) > 50 & ( n_bi > 2 | n_mono > 2) ){
        cc = cor.test( xp, n_alleles )
        cor_estimate = cc$estimate
        cor_pval = cc$p.value
    }
    if( !is.null(order)){
        n_bi = n_bi[order]
        n_mono = n_mono[order]
        n_bi_from_sv = n_bi_from_sv[order]
        n_amp=n_amp[order]
        cor_estimate = cor_estimate[order]
        cor_pval = cor_pval[order]
        df = df[order,]
    }
    list( 
        n_bi = n_bi, n_mono = n_mono, n_bi_from_sv = n_bi_from_sv, n_amp=n_amp,
        cor_alleles_xp = signif( as.numeric(cor_estimate), 3),
        pval_alleles_xp = signif( cor_pval, 3 ),
        alleles = df
    )
}


plot_gene_methylation_across_samples = function(genename){
    
    #How many of the top eHMRs do you want to plot? 1 for all figures in paper
    numbest <- 1
    
    scalerange = TRUE #show raw methyl or scale to min/max
    logexpr = TRUE
    groupalts = FALSE #group the DNA altered T/F groups separately
    use_sc = FALSE #mark the tSCNC samples? may need for revisions
    
    isgof = FALSE
    hmrsegs = TRUEracks[['HMRseggene']]
    
    rowensembl <- ensembl2sym$name==genename
    ensemblid <- rownames(ensembl2sym)[rowensembl]
    stopifnot(sum(rowensembl)==1)
    chr <- ensembl2sym[rowensembl,'chr']
    gene_start <- ensembl2sym[rowensembl,'start']
    gene_end <- ensembl2sym[rowensembl,'end']
    
    samples2use <- samples_wgbs
    if(genename=='MYC') {
        start <- gene_start-5000
        end <- gene_end+120000
    } else if(genename=='AR') {
        start <- gene_start-850000
        end <- gene_end+100000
    } else if(genename=='ERG') {
        start <- gene_start+3100000
        end <- gene_end+2930000
        samples2use <- rownames(wcdt_genes)[wcdt_genes$TMPRSS2_ERG==3]
        #promoter/genebody eHMRs
    } else if(ensembl2sym[ensemblid,'strand']=='+') {
        start <- gene_start-1500
        end <- gene_end
    } else {
        start <- gene_start
        end <- gene_end+1500
    }
    gr <- GRanges(chr,IRanges(start,end))
    
    expr <- as.numeric(tpm[ensemblid,samples2use])
    if(logexpr) {
        expr <- log2(expr+1)
    }
    
    expr <- expr / max(expr)
    covs <- data.frame(expr)
    rownames(covs) <- samples2use
    
    gof <- c(
        'activating_missense',
        'activating_sv',
        'CNA_amp')
    lof <- c(
        'inactivating_missense',
        'nonsense',
        'inactivating_germline',
        'inactivating_sv',
        'CNA_1',
        'CNA_2')
    grouplvls <- c(gof,lof,'HMR','TPM','tSCNC')
    colors <- c(gofcol,lofcol,
                'green','darkgray','black')
    
    #Set up data frames for the plots
    df <- data.frame()
    
    all <- allele_effect(genename)$alleles
    
    df1 <- data.frame(
        x=samples2use,
        y=3,
        ht=0.95,
        fill=rep(NA,length(samples2use)),
        stringsAsFactors=F)
    df2 <- data.frame(
        x=samples2use,
        y=4,
        ht=0.95,
        fill=rep(NA,length(samples2use)),
        stringsAsFactors=F)
    df3 <- data.frame()
    dfexpr <- data.frame(
        x=samples2use,
        y=2-expr,
        ht=expr*2,
        fill=rep('TPM',length(samples2use)))
    if(use_sc) {
        dfexpr$fill <- paste(dfexpr$fill)
        dfexpr[samples2use %in% tscnc,'fill'] <- 'tSCNC'
    }
    
    rownames(df1) <- samples2use
    rownames(df2) <- samples2use
    rownames(dfexpr) <- samples2use
    
    promotermethyl <- hmrsegs[countOverlaps(hmrsegs,gr)>0 & hmrsegs$gene_id==ensemblid]
    roworder <- rev(order(abs(promotermethyl$avg_cor)))
    emr_best <- promotermethyl[roworder]
    emr_best <- emr_best[!is.na(emr_best$avg_cor)]
    print(range(emr_best$avg_cor))
    print(emr_best)
    
    colramp <- colorRampPalette(c('blue','white','red'))(n = 100)
    grouplvls <- c(grouplvls,paste(1:100))
    colors <- c(colors,colramp)
    names(colors) <- grouplvls
    
    if(length(emr_best)>=numbest) {
        #Loop through all recurrent HMRs and add them to model
        for(i in 1:length(emr_best)) {
            emr_start <- start(emr_best)[i]
            emr_end <- end(emr_best)[i]
            methylmean <- meanmethyl(ensemblid, dir_gene_output_mcrpc, emr_start, emr_end)[samples2use]
            covs[,i+1] <- methylmean
            #Only plot the best numbest recurrent HMRs
            if(i <= numbest) {
                df3i <- data.frame(
                    x=samples2use,
                    y=i+4,
                    ht=0.95,
                    fill=rep(NA,length(samples2use)),
                    stringsAsFactors=F)
                rownames(df3i) <- samples2use
                if(scalerange) {
                    meanmethylround <- round((methylmean/max(methylmean,na.rm=T))*100)
                    print(range(methylmean,na.rm=T))
                } else {
                    meanmethylround <- round(methylmean)
                }
                meanmethylround[meanmethylround==0] <- 1
                df3i$fill <- paste(meanmethylround)
                df3 <- rbind.data.frame(df3,df3i)
            }
        }
    } else {
        print('No EMR')
    }
    
    #Add DNA alterations to plot dataframes
    if(isgof) {
        covs$n_alleles <- 2
        for(j in 1:length(gof)) {
            namegof <- gof[j]
            samplesalt <- samples2use[all[samples2use,namegof]]
            for(samplealt in samplesalt) {
                covs[samplealt,'n_alleles'] <- 3
                count <- is.na(df1[samplealt,'fill']) + is.na(df2[samplealt,'fill'])
                if(count==2) {
                    df1[samplealt,'fill'] <- namegof
                } else if(count==1) {
                    df2[samplealt,'fill'] <- namegof
                }
            }
        }
    } else {
        covs$n_alleles <- 2-all[samples2use,'n_alleles_inactivated']
        for(j in 1:length(lof)) {
            namelof <- lof[j]
            samplesalt <- samples2use[all[samples2use,namelof]]
            for(samplealt in samplesalt) {
                if(namelof=='CNA_2') {
                    df1[samplealt,'fill'] <- namelof
                    df2[samplealt,'fill'] <- namelof
                } else {
                    count <- is.na(df1[samplealt,'fill']) + is.na(df2[samplealt,'fill'])
                    if(count==2) {
                        df1[samplealt,'fill'] <- namelof
                    } else if(count==1) {
                        df2[samplealt,'fill'] <- namelof
                    }
                }
            }
        }
    }
    
    if(groupalts) {
        collvls <- samples2use[order(is.na(df1$fill),expr)]
    } else {
        collvls <- samples2use[order(expr)]
    }
    
    df <- rbind.data.frame(df1,df2,df3,dfexpr)
    df$x = FALSEactor(df$x,levels=collvls)
    df$fill = FALSEactor(df$fill,levels=grouplvls)
    rowsna <- is.na(df$fill)
    df <- df[!rowsna,]
    
    #Compare linear models
    rows2keep <- rowSums(is.na(covs))==0
    
    brks <- c(1,3:(4+numbest))
    rowlabs <- c('TPM','DNA1','DNA2',paste('Methyl',1:numbest,sep=''))
    lm1 <- lm(expr~n_alleles,covs[rows2keep,])
    lm2 <- lm(expr~.,covs[rows2keep,])
    pval <- anova(lm1,lm2)[2,'Pr(>F)']
    print(summary(lm1)$adj.r.squared)
    print(summary(lm2)$adj.r.squared)
    
    colors <- colors[names(colors) %in% df$fill]
    
    if(logexpr) {
        rowlabs[1] <- 'Log2(TPM)+1'
    }
    
    p <- ggplot(df,aes(x=x,y=y,width=0.5,height=ht))+geom_tile(aes(fill=fill))+
        scale_fill_manual(values=colors,breaks=c(gof,lof))+theme_classic()+
        #theme(axis.text.x=element_text(angle=90, hjust=1))+
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
        xlab(paste('ANOVA P:',signif(pval,4)))+
        scale_y_reverse(name=genename,minor_breaks=c(2.05,4.5),breaks=brks,labels=rowlabs)+
        theme(panel.grid.minor=element_line(colour='black',size=0.5))
    #ggsave(fileout,p,width=10,height=4)
    
    covs$n_alleles = FALSEactor(covs$n_alleles)
    p <- ggplot(covs,aes(x=n_alleles,y=V2))+geom_boxplot()+theme_classic()
    #ggsave(filebox,p)
}

#gets the mean methylation for a region given the path to the gene-level slice files
meanmethyl = function(id, genepath, start, end) {
    stopifnot(length(id)==1)
    file_methyl = list.files(genepath, pattern=id)
    stopifnot(length(file_methyl)==1)
    filenamefull = paste(genepath,file_methyl,sep='')
    methyl = read.delim(filenamefull,sep='',row.names=1,header=T,stringsAsFactors=F,check.names=F)
    colnames(methyl) = toupper(str_split_fixed(colnames(methyl),'\\.',2)[,1])
    stopifnot(dim(methyl)[1]>0)
    methylpos = as.numeric(rownames(methyl))
    methylrows = methylpos >= start & methylpos <= end
    means = colMeans(methyl[methylrows,],na.rm=T)
    return(means)
}

#Overlap between two pairs of vectors
overlappairs = function(x,y,a,b) {
    stopifnot(sum(x>=y)==0)
    stopifnot(sum(a>=b)==0)
    return((x>=a & x<=b) | (y>=a & y<=b) | (x<a & y>b))
}


# Helper function to extract information from gene level filename
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

#Converts chromosome number to name
convertchr = function(chrin) {
    temp = gsub('chr','',chrin,fixed=T)
    temp = gsub('X','23',temp,fixed=T)
    temp = gsub('Y','24',temp,fixed=T)
    return(as.numeric(temp))
}

#Calculate the 1D position given chr/pos
calcpos = function(chrin,pos1,pos2) {
    stopifnot(pos2>=pos1)
    chrin = convertchr(chrin)
    pos = (pos1+pos2)/2
    for(i in 1:length(chrsizes)) {
        rows2add = chrin > i
        pos[rows2add] = pos[rows2add] + chrsizes[i]
    }
    return(pos)
}

#Calculate the chr/pos given the 1D position
calcposrev = function(pos1,pos2) {
    stopifnot(pos2>=pos1)
    initsize = 0
    for(i in 1:length(chrsizes)) {
        chrsize = as.numeric(chrsizes[i])
        if(pos1>initsize & pos1<=initsize+chrsize) { #found chr
            if(pos2>initsize & pos2<=initsize+chrsize) { #both in chr
                return(c(i,pos1-initsize,pos2-initsize))
            } else { #split across chr, not OK
                return(NULL)
            }
        }
        initsize = initsize+chrsize
    }
    return(NULL)
}

heatmap.3 = function(x,
                     Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                     distfun = dist,
                     hclustfun = hclust,
                     dendrogram = c("both","row", "column", "none"),
                     symm = FALSE,
                     scale = c("none","row", "column"),
                     na.rm = TRUE,
                     revC = identical(Colv,"Rowv"),
                     add.expr,
                     breaks,
                     symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                     col = "heat.colors",
                     colsep,
                     rowsep,
                     sepcolor = "white",
                     sepwidth = c(0.05, 0.05),
                     cellnote,
                     notecex = 1,
                     notecol = "cyan",
                     na.color = par("bg"),
                     trace = c("none", "column","row", "both"),
                     tracecol = "cyan",
                     hline = median(breaks),
                     vline = median(breaks),
                     linecol = tracecol,
                     margins = c(5,5),
                     ColSideColors,
                     RowSideColors,
                     side.height.fraction=0.3,
                     cexRow = 0.2 + 1/log10(nr),
                     cexCol = 0.2 + 1/log10(nc),
                     labRow = NULL,
                     labCol = NULL,
                     key = TRUE,
                     keysize = 1.5,
                     density.info = c("none", "histogram", "density"),
                     denscol = tracecol,
                     symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                     densadj = 0.25,
                     main = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     lmat = NULL,
                     lhei = NULL,
                     lwid = NULL,
                     ColSideColorsSize = 1,
                     RowSideColorsSize = 1,
                     KeyValueName="Value",...){
    
    invalid = function (x) {
        if (missing(x) || is.null(x) || length(x) == 0)
            return(TRUE)
        if (is.list(x))
            return(all(sapply(x, invalid)))
        else if (is.vector(x))
            return(all(is.na(x)))
        else return(FALSE)
    }
    
    x = as.matrix(x)
    scale01 = function(x, low = min(x), high = max(x)) {
        x = (x - low)/(high - low)
        x
    }
    retval = list()
    scale = if (symm && missing(scale))
        "none"
    else match.arg(scale)
    dendrogram = match.arg(dendrogram)
    trace = match.arg(trace)
    density.info = match.arg(density.info)
    if (length(col) == 1 && is.character(col))
        col = get(col, mode = "function")
    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
                "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv))
        Rowv = FALSE
    if (is.null(Colv) || is.na(Colv))
        Colv = FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv = FALSE
    di = dim(x)
    if (length(di) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")
    nr = di[1]
    nc = di[2]
    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote))
        cellnote = matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                     c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram = "column"
            else dedrogram = "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                     c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram = "row"
            else dendrogram = "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                    dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr = Rowv
        rowInd = order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr = hclustfun(distfun(x))
        ddr = as.dendrogram(hcr)
        ddr = reorder(ddr, Rowv)
        rowInd = order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv = rowMeans(x, na.rm = na.rm)
        hcr = hclustfun(distfun(x))
        ddr = as.dendrogram(hcr)
        ddr = reorder(ddr, Rowv)
        rowInd = order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd = nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc = Colv
        colInd = order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc = ddr
            colInd = order.dendrogram(ddc)
        }
        else colInd = rowInd
    }
    else if (is.integer(Colv)) {
        hcc = hclustfun(distfun(if (symm)
            x
            else t(x)))
        ddc = as.dendrogram(hcc)
        ddc = reorder(ddc, Colv)
        colInd = order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv = colMeans(x, na.rm = na.rm)
        hcc = hclustfun(distfun(if (symm)
            x
            else t(x)))
        ddc = as.dendrogram(hcc)
        ddc = reorder(ddc, Colv)
        colInd = order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd = 1:nc
    }
    retval$rowInd = rowInd
    retval$colInd = colInd
    retval$call = match.call()
    x = x[rowInd, colInd]
    x.unscaled = x
    cellnote = cellnote[rowInd, colInd]
    if (is.null(labRow))
        labRow = if (is.null(rownames(x)))
            (1:nr)[rowInd]
    else rownames(x)
    else labRow = labRow[rowInd]
    if (is.null(labCol))
        labCol = if (is.null(colnames(x)))
            (1:nc)[colInd]
    else colnames(x)
    else labCol = labCol[colInd]
    if (scale == "row") {
        retval$rowMeans = rm = rowMeans(x, na.rm = na.rm)
        x = sweep(x, 1, rm)
        retval$rowSDs = sx = apply(x, 1, sd, na.rm = na.rm)
        x = sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans = rm = colMeans(x, na.rm = na.rm)
        x = sweep(x, 2, rm)
        retval$colSDs = sx = apply(x, 2, sd, na.rm = na.rm)
        x = sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks = 16
        else breaks = length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks)
            breaks = seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                         length = breaks)
        else {
            extreme = max(abs(x), na.rm = TRUE)
            breaks = seq(-extreme, extreme, length = breaks)
        }
    }
    nbr = length(breaks)
    ncol = length(breaks) - 1
    if (class(col) == "function")
        col = col(ncol)
    min.breaks = min(breaks)
    max.breaks = max(breaks)
    x[x < min.breaks] = min.breaks
    x[x > max.breaks] = max.breaks
    if (missing(lhei) || is.null(lhei))
        lhei = c(keysize, 4)
    if (missing(lwid) || is.null(lwid))
        lwid = c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat = rbind(4:3, 2:1)
        
        if (!missing(ColSideColors)) {
            #if (!is.matrix(ColSideColors))
            #stop("'ColSideColors' must be a matrix")
            if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of nrow(x) rows")
            lmat = rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei = c(lhei[1], 0.2, lhei[2])
            lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
        }
        
        if (!missing(RowSideColors)) {
            #if (!is.matrix(RowSideColors))
            #stop("'RowSideColors' must be a matrix")
            if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of ncol(x) columns")
            lmat = cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid = c(lwid[1], 0.2, lwid[2])
            lwid = c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
        }
        lmat[is.na(lmat)] = 0
    }
    
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op = par(no.readonly = TRUE)
    on.exit(par(op))
    
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    
    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
            par(mar = c(margins[1], 0, 0, 0.5))
            image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[,rowInd, drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
            if (length(rownames(RowSideColors)) > 0) {
                axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), cex.axis=2, las = 2, tick = FALSE)
            }
        }
    }
    
    if (!missing(ColSideColors)) {
        
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[colInd, , drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            image(csc, col = as.vector(csc.colors), axes = FALSE)
            if (length(colnames(ColSideColors)) > 0) {
                axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), cex.axis=2, las = 2, tick = FALSE)
            }
        }
    }
    
    par(mar = c(margins[1], 0, 0, margins[2]))
    x = t(x)
    cellnote = t(cellnote)
    if (revC) {
        iy = nr:1
        if (exists("ddr"))
            ddr = rev(ddr)
        x = x[, iy]
        cellnote = cellnote[, iy]
    }
    else iy = 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    retval$carpet = x
    if (exists("ddr"))
        retval$rowDendrogram = ddr
    if (exists("ddc"))
        retval$colDendrogram = ddc
    retval$breaks = breaks
    retval$col = col
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat = ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
              col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
         cex.axis = cexCol)
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
         cex.axis = cexRow)
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    min.scale = min(breaks)
    max.scale = max(breaks)
    x.scaled = scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline = vline
        vline.vals = scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol,
                       lty = 2)
            }
            xv = rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv = c(xv[1], xv)
            yv = 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline = hline
        hline.vals = scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv = rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv = rev(c(yv[1], yv))
            xv = length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
             col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main))
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks = breaks
        if (symkey) {
            max.raw = max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw = -max.raw
            tmpbreaks[1] = -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] = max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw = min(x, na.rm = TRUE)
            max.raw = max(x, na.rm = TRUE)
        }
        
        z = seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
              xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv = pretty(breaks)
        xv = scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens = density(x, adjust = densadj, na.rm = TRUE)
            omit = dens$x < min(breaks) | dens$x > max(breaks)
            dens$x = dens$x[-omit]
            dens$y = dens$y[-omit]
            dens$x = scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                  lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h = hist(x, plot = FALSE, breaks = breaks)
            hx = scale01(breaks, min.raw, max.raw)
            hy = c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                  col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable = data.frame(low = retval$breaks[-length(retval$breaks)],
                                   high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

convertClinSamplenames = function(samples) {
    samples[samples=='DTB-77Pro'] ='DTB-077Pro'
    samples = gsub('PRO', '-PRO', toupper(samples))
    samples = gsub('PRO2', 'PRO2', toupper(samples))
    notpro = !grepl('PRO', samples)
    samples[notpro] = paste(samples[notpro],'BL', sep='-')
    return(samples)
}

## slices out rows of a dataframe which match str2find in a specific colname
slicerows = function(dataframe, colname, str2find) {
    rows = toupper(dataframe[,colname])==str2find
    rows[is.na(rows)] = F
    if(sum(rows)>0) {
        return(dataframe[rows,])
    } else {
        return(NULL)
    }
}

getnum = function(txt) {
    return(as.numeric(gsub("[^0-9]", "", txt)))
}

getexpr = function(gene, sample, dataset) {
    if(is.na(gene) || is.na(sample)) {return(c(NA,NA))}
    if(gene %in% rownames(dataset) && sample %in% colnames(dataset)) {
        expr = dataset[gene, sample]
        pctile = (sum(dataset[gene,]<expr)/(dim(dataset)[2])) * 100
        return(c(expr,round(pctile)))
    } else {
        return(c(NA,NA))
    }
}

getexprs = function(gene, samples, dataset) {
    stopifnot(gene %in% rownames(dataset))
    expr = rep(NA, length(samples))
    
    exprows = samples %in% colnames(dataset)
    samples2find = samples[exprows]
    expr[exprows] = as.numeric(dataset[gene, samples2find])
    
    pctile = expr
    for(i in 1:length(expr)) {
        pctile[i] = (sum(dataset[gene,]<expr[i])/(dim(dataset)[2])) * 100
    }
    return(cbind(expr,round(pctile)))
}

inrange = function(x, y, z) {
    #y = start coord
    #z = end coord 
    stopifnot(y<=z)
    return(x>=y & x<=z)
}

overlap = function(x1, x2, y, z) {
    #y = start coord
    #z = end coord 
    stopifnot(length(x1)==length(x2))
    if(length(x1)==1 && !is.na(x1) && !is.na(x2) && x1>x2) {
        x3 = x2
        x2 = x1
        x1 = x3
    } else if(length(x1)>1) {
        rowswrong = x1>x2
        rowswrong[is.na(rowswrong)] = F
        x3 = x2[rowswrong]
        x2[rowswrong] = x1[rowswrong]
        x1[rowswrong] = x3
    }
    stopifnot(y<=z)
    x1in = x1>=y & x1<=z
    x2in = x2>=y & x2<=z
    xencompasses = x1<y & x2>z
    return(x1in | x2in | xencompasses)
}

splitgenes = function(genes, sep='|',max=10) {
    split = strsplit(genes,sep,fixed=T)
    splitn = sapply(split,length)
    splitfinal = split[splitn <= max]
    genesfinal = unique(unlist(splitfinal))
    genesfinal = genesfinal[!is.na(genesfinal) & genesfinal!='' & genesfinal!='NA']
    return(genesfinal)
}

gene2samples = function(gene2find, dataset, colsgene, colsample) {
    stopifnot(length(colsgene)>=1)
    rows = rep(F,dim(dataset)[1])
    for(i in 1:length(colsgene)) {
        rows = rows | dataset[,colsgene[i]]==gene2find |
            grepl(paste('^',gene2find,'\\|',sep=''),dataset[,colsgene[i]]) |
            grepl(paste('\\|',gene2find,'$',sep=''),dataset[,colsgene[i]]) |
            grepl(paste('|',gene2find,'|',sep=''),dataset[,colsgene[i]],fixed=T)
    }
    
    results = dataset[rows,colsample]
    results = results[!is.na(results)]
    results = results[results!='']
    return(results)
}

ensemblstrip = function(ensembldot) {
    return(str_split_fixed(ensembldot, "\\.", 2)[,1])
}

# Copycat functions
#------------------------------------------------------------------------------------------------------------------

countoverlap = function(x1,y1,x2,y2) {
    stopifnot(x1<=y1)
    stopifnot(x2<=y2)
    return(max(0,min(y1,y2)-max(x1,x2)+1))
}

cnregion = function(chr,start,end) {
    rowsbed = copycat_bed$chrom==chr & overlap(copycat_bed$start,copycat_bed$end,start,end)
    stopifnot(sum(rowsbed)>0)
    cc = copycat_bed[rowsbed,]
    genewidth = end-start+1
    output = rep(0,length(unique(cc$sample_id)))
    names(output) = unique(cc$sample_id)
    for(i in 1:dim(cc)[1]) {
        samplei = cc[i,'sample_id']
        width = countoverlap(cc[i,'start'],cc[i,'end'],start,end)
        output[samplei] = output[samplei] + (cc[i,'cn']*width/genewidth)
    }
    return(output)
}

# Annovar functions
#------------------------------------------------------------------------------------------------------------------

annotate_annovar = function(annovardata) {
    stopifnot(sum(grepl(',',annovardata$Alt,fixed=T))==0)
    annovardata$muttype = annovardata$Func.ensGene
    rows_exon = grepl('exonic', annovardata$Func.ensGene, fixed=T)
    annovardata[rows_exon,'muttype'] = 'Exonic'
    rows_intron = grepl('intronic', annovardata$Func.ensGene, fixed=T)
    annovardata[rows_intron,'muttype'] = 'Intronic'
    rows_splice = grepl('splicing', annovardata$Func.ensGene, fixed=T)
    annovardata[rows_splice,'muttype'] = 'Splice'
    
    annovardata[annovardata$Func.ensGene=='UTR3','muttype'] = 'UTR'
    annovardata[annovardata$Func.ensGene=='UTR5','muttype'] = 'UTR'
    annovardata[annovardata$Func.ensGene=='UTR5;UTR3','muttype'] = 'UTR'
    annovardata[annovardata$Func.ensGene=='downstream','muttype'] = 'Down1k'
    annovardata[annovardata$Func.ensGene=='upstream','muttype'] = 'Up1k'
    annovardata[annovardata$Func.ensGene=='upstream;downstream','muttype'] = 'Up1k'
    annovardata[annovardata$ExonicFunc.ensGene=='synonymous SNV','muttype'] = 'Silent'
    annovardata[annovardata$ExonicFunc.ensGene=='nonsynonymous SNV','muttype'] = 'Missense'
    annovardata[annovardata$ExonicFunc.ensGene=='frameshift substitution','muttype'] = 'Frameshift'
    annovardata[annovardata$ExonicFunc.ensGene=='stopgain','muttype'] = 'Stopgain'
    annovardata[annovardata$ExonicFunc.ensGene=='nonframeshift substitution','muttype'] = 'Indel'
    annovardata[annovardata$ExonicFunc.ensGene=='stoploss','muttype'] = 'Stoploss'
    
    annovardata$funcmut = F
    annovardata$funcmut[annovardata$muttype == 'Frameshift' |
                            annovardata$muttype == 'Splice' |
                            annovardata$muttype == 'Stopgain' |
                            annovardata$muttype == 'Stoploss' |
                            annovardata$SIFT_pred=='D' | annovardata$Polyphen2_HVAR_pred=='D' |
                            grepl('Pathogenic',annovardata$CLINSIG,fixed=T)] = T
    cols2keep = c(1:5,7:10,50:55)
    return(annovardata[,cols2keep])
}


# RNAseq functions
#------------------------------------------------------------------------------------------------------------------

modify_RNAseq_names = function(exprdf) {
    colnames(exprdf) = toupper(gsub('-D-T-RNA', '', colnames(exprdf), fixed=T))
    colnames(exprdf) = gsub('-T-RNA', '', colnames(exprdf), fixed=T)
    samples2remove = NULL
    for(i in 1:dim(exprdf)[1]) {
        samplename = colnames(exprdf)[i]
        if(grepl('-NS',samplename,fixed=T)) {
            nameclean = gsub('-NS', '', samplename)
            if(nameclean %in% colnames(exprdf)) {
                samples2remove = c(samples2remove,nameclean)
            }
        }
    }
    colnames(exprdf)[colnames(exprdf)=='DTB-095-BL'] = 'DTB-095-1-BL'
    colnames(exprdf)[colnames(exprdf)=='DTB-095-BL2-BONE'] = 'DTB-095-2-BL'
    exprdf = exprdf[,!(colnames(exprdf) %in% samples2remove)]
    colnames(exprdf) = gsub('-NS', '', colnames(exprdf), fixed=T)
    return(exprdf)
}


# MethylseqR functions

extract_methylseqR_metadata = function(fname) {
    # Extract MethylSeqR metadata from filename
    retval = list()
    retval$sname = strsplit(fname, '_')[[1]][1]
    retval$chr = strsplit(fname, '_')[[1]][3]
    return(retval)
}
