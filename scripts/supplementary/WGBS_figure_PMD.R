# python /opt/software/marlowe/pipeline_scripts/generate_cluster_json_quantify_methylation.py \
#   -i /home/dquigley/notebook/human_sequence_prostate_WGBS/methylseekr/feature_beds/BEDs_for_PMD_tiling_10K.txt \
#   -t tile \
#   -j quantify_localized.json \
#   -d localized \
#   -o /data1/queuescripts/quantify_methylation

#python /opt/software/marlowe/marlowe/process_workflow.py -r  \
#  -i /data1/queuescripts/quantify_methylation/quantify_localized.json > /data1/queuescripts/quantify_methylation/quantify_localized.sh

source('/notebook/code/src/R/quantitative_genetics.R')

plot.color.grid = function(M, block.height=20, block.width=10, space.X=3, space.Y=10, 
         cex.x=1, cex.y=1, border=TRUE, show.axis.X=TRUE, show.axis.Y=TRUE, font.Y=1){
    
    n.rows = dim(M)[1]
    n.cols = dim(M)[2]
    total.width =  ( block.width*n.cols) + ( (n.cols-1) * space.X)
    total.height = ( block.height * n.rows ) + ( (n.rows-1) * space.Y)
    plot(0,0,col="white", xlim=c(0,total.width), ylim=c(0,total.height), 
         axes=F, xlab="", ylab="", bg="azure2", xaxs="i", yaxs="i")
    xlab_locs = rep(0, n.cols)
    ylab_locs = rep(0, n.rows)
    cur.y = total.height
    for(rr in 1:n.rows){
        for(cc in 1:n.cols){  
            this.x.left = (cc-1)*block.width + (cc-1)*space.X
            this.x.right = this.x.left+block.width
            xlab_locs[cc] = (this.x.right + this.x.left)/2
            #xlab_locs[cc] = this.x.right - (block.width)
            if( border )
                rect( this.x.left , cur.y - block.height, this.x.right, cur.y, 
                      col=M[rr,cc], border="azure2")
            else
                rect( this.x.left , cur.y - block.height, this.x.right, cur.y, 
                      col=M[rr,cc], border=NA)
        }
        ylab_locs[rr] = cur.y - (0.5*block.height)
        cur.y = cur.y - block.height - space.Y
    }
    if( show.axis.X ){
        axis(1, at=xlab_locs, labels=dimnames(M)[[2]], las=2, cex.axis=cex.x, tick=FALSE, padj=1 )
    }
    if( show.axis.Y ){   
        axis(2, at=ylab_locs, labels=rownames(M), las=2, cex.axis=cex.y, tick=FALSE, hadj=1, font=font.Y)
    }
}

color_scale = function( V, color_map, color_bounds=NA, color_NA=NA ){
    if( is.na( color_bounds[1] ) ){
        Vmax = max(V, na.rm=TRUE)
        Vmin = min(V, na.rm=TRUE)
    }else{
        Vmin = color_bounds[1]
        Vmax = color_bounds[2]
    }
    increment = (Vmax-Vmin) / (length(color_map)-1)
    lookup = seq(from=Vmin, to=Vmax, by=increment )
    if( sum(V>Vmax, na.rm=TRUE)>0 ){
        warning("clipped colors that exceed Vmax")   
    }
    V[V>Vmax] = Vmax-0.001
    if( is.vector(V) ){
        out = rep("", length(V))
        for(i in 1:length(V)){
            out[i] = color_map[ which( lookup>V[i] )[1] ]
        }
    }else if( is.matrix(V) ){
        out = matrix("", nrow=dim(V)[1], ncol=dim(V)[2], 
                     dimnames=list( dimnames(V)[[1]], dimnames(V)[[2]]))
        for(rr in 1:dim(V)[1]){
            for(cc in 1:dim(V)[2]){
                col=color_map[ which( lookup>=V[rr,cc] )[1]]
                out[rr,cc] = col 
            }
        }
    }
    if( !is.na( color_NA) )
        out[is.na(out)] = color_NA
    out
}

pmd_to_means=function(pmd, chrom,  ncol, binwidth=250000){
    M = matrix(NA, nrow=3, ncol=ncol)
    n = ceiling( chrom_lengths$V2[ which(dimnames(chrom_lengths)[[1]] == chrom ) ] / binwidth )
    for(i in 1:n ){
        binend = i * binwidth
        binstart = binend - binwidth
        idx = which( pmd$chrom==chrom & pmd$pos>=binstart & pmd$pos<=binend )
        M[1,i] = mean( pmd$mu_n[ idx ], na.rm = TRUE)
        M[2,i] = mean( pmd$mu_loc[ idx ], na.rm = TRUE)
        M[3,i] = mean( pmd$mu_t[ idx ], na.rm = TRUE)
    }
    M
}

medians_to_colors=function( M, chrom, col.low, col.med, col.high, binwidth=250000 ){
    ncol = dim(M)[2]
    n = ceiling( chrom_lengths$V2[ which(dimnames(chrom_lengths)[[1]] == chrom ) ] / binwidth )
    cmap = grDevices::colorRampPalette( colors=c(col.low, col.med, col.high) )(50)
    cs=color_scale(M, cmap, color_bounds = c(0,1), color_NA = "lightgrey")
    cs[,(n+1):ncol] = "white"
    cs
}

methyl_matrix = function(D){
    sid = sort(unique(D$sample_id))
    nrow = length(sid)
    ncol = length(D$chrom[D$sample_id==sid[1]])
    M = matrix(NA, nrow=nrow, ncol=ncol)
    dimnames(M)[[1]] = sid
    sample_id_list = sort(unique(D$sample_id))
    for(i in 1:nrow){
        M[i,] = D$mean[ D$sample_id==sample_id_list[i] ]
    }
    M
}

add_rowsum = function(D){
    M = methyl_matrix(D)
    vals = rep(NA, dim(D)[1])
    sums = rowSums(M, na.rm=TRUE)
    for( i in 1:length( sums )){ 
        vals[ D$sample_id == dimnames(M)[[1]][i] ] = sums[i]
    }
    D$order = vals
    D
}

add_cluster = function(D){
    M = methyl_matrix(D)
    h = hclust(dist(M))
    vals = rep(NA, dim(D)[1])
    for( i in 1:length( h$order )){ 
        vals[ D$sample_id == dimnames(M)[[1]][i] ] = h$order[i]
    }
    D$order = vals
    D
}


create_pmd_plot = function(D, chrom, order="sum"){
    DX=D
    DX$order = order(D$sample_id)
    if(order=="sum"){
        DX = add_rowsum(D)
    }else if( order=="cluster"){
        DX = add_cluster(DX)
    }
    ggplot( DX, 
            aes(x=start, y=fct_reorder(sample_id, order, .desc=TRUE), fill=mean) ) + 
        geom_tile() + 
        #scale_fill_gradient2(low = "#ffffff", high = "#253494", mid = "#41b6c4", 
        #                     midpoint = 0.5, limit = c(0,1)) +
        scale_fill_gradient2(low = "#253494", high = "yellow", mid = "lightgrey", 
                             midpoint = 0.5, limit = c(0,1)) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_text(size=14, face="bold")) +
        xlab("position") + ylab("sample") +
        ggtitle(chrom)
}

create_label_plot = function( S, sample_order ){
    df = data.frame( x=rep(1, length(S)),
                     y=seq(from=length(S), to=1, by=-1),
                     z=factor( S ))
    df = df[order(sample_order),]
    ggplot(df, aes(x, y)) +
        geom_tile(aes(fill = z), colour = "white") +
        scale_x_continuous(expand = c(0.0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.title = element_text(size=14, face="bold"),
              legend.position="none",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
        ggtitle("")
}

create_combined_pmd_plot = function(chrom, xlims=NULL, order_by="sum", DN=NULL, DT=NULL){
    if(is.null(DN)){
        DN = read.pmd.chromosome(chrom, "localized", xlim=xlims)
    }else{
        DN = DN[DN$start >= xlims[1] & DN$stop <= xlims[2],]     
    }
    DN$source = rep("normal", dim(DN)[1])
    DN = add_rowsum(DN)
    DN$order = -1*DN$order
    if( is.null(DT)){
        DT=read.pmd.chromosome(chrom, "tumor", xlim=xlims)
    }else{
        DT = DT[DT$start >= xlims[1] & DT$stop <= xlims[2],]     
    }
    DT$source = rep("tumor", dim(DT)[1])
    DTO=DT
    DTO$order = order(DT$sample_id)
    if(order_by=="sum"){
        DTO = add_rowsum(DTO)
    }else if( order_by=="cluster"){
        DTO = add_cluster(DTO)
    }
    D = rbind(DN, DTO)
    nt=ggplot( D, 
               aes(x=start, 
                   y=fct_reorder(sample_id, order, .desc=TRUE), 
                   fill=mean) ) + 
        geom_tile() + 
        scale_x_continuous(expand = c(0.01, 0)) +
        scale_fill_gradient2(low = "#253494", high = "yellow", mid = "lightgrey", 
                             midpoint = 0.5, limit = c(0,1)) +
        theme(axis.title.y=element_blank(),
              axis.text.y = element_text(color = "grey20", size = 6 ),
              axis.ticks.y=element_blank(),
              plot.title = element_text(size=14, face="bold"),
              legend.position="none",
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank() ) +
        xlab("position") + ylab("sample") +
        ggtitle(chrom)
    
    sids = sort( unique( DT$sample_id ))
    m = match.idx( sids, sample_ids_wgs )
    N_normal = length(unique(DN$sample_id))
    N_tumor = length(unique(DTO$sample_id))
    ords = rep(0, length(sids))
    for(i in 1:length(sids)){
        ords[i] = min(DTO$order[DTO$sample_id==sids[i] ], na.rm=TRUE)
    }
    sids = c( unique(DN$sample_id), sids )
    ords = c( -1 * (1:length( unique(DN$sample_id) )), ords )
    sample_type = c(rep( "normal", N_normal), rep( "tumor", N_tumor) )
    ETSs = rep( "neg", N_tumor )
    ETSs[ matrix_samples$has_ETS[m$idx.A] ] = "pos"
    ETSs = c( rep("normal", N_normal), ETSs )
    ETSs = ETSs[ order(ords)] 
    
    ids_tSCNC_WGS = c("DTB-205-BL","DTB-003-BL","DTB-036-BL", 
                      "DTB-135-PRO","DTB-040-BL")
    sids = sort( unique( DT$sample_id ))
    m = match.idx( sids, ids_tSCNC_WGS)
    neuros = rep( "adeno", N_tumor )
    neuros[ m$idx.A ] = "t-SCNC"
    neuros = c( rep("normal", N_normal), neuros )
    neuros = neuros[ order(ords)] 
    
    m = match.idx( sids, samples_cmp )
    CMPs = rep("adeno", N_tumor)
    CMPs[ m$idx.A ] = "CMP"
    CMPs = c( rep("normal", N_normal), CMPs)
    CMPs = CMPs[ order( ords ) ]
    
    label_sample_type = create_label_plot( sample_type, ords )
    label_ETS = create_label_plot( ETSs, ords )
    label_tscnc = create_label_plot( neuros, ords )
    label_CMP = create_label_plot( CMPs, ords )
    list(data=nt, 
         label_sample_type=label_sample_type, 
         label_ETS=label_ETS, 
         label_tscnc = label_tscnc,
         label_CMP = label_CMP)
    
}

read.pmd.chromosome = function(chrom, sample_type="tumor", xlim=NULL){
    D = data.frame()
    if( sample_type=="tumor" ){
        sample_set = sample_ids_wgbs
    }else if( sample_type == "normal" ){
        sample_set = samples_normal_filecase
    }else if( sample_type == "localized" ){
        sample_set = setdiff( samples_upitt_localized_tumor, "Bis158T")
    }else{
        stop("unrecognized sample type")   
    }
    for(sample_id in sample_set){
        fn=paste(dir_pmd, '/', sample_id,"_methylseekr_",chrom,"_counts.txt",sep="" )
        d = read.table( fn, header=TRUE,stringsAsFactors = FALSE )    
        if(dim(D)[1]==0){
            D = d
        }else{
            D = rbind(D,d)
        }
    }
    if( !is.null(xlim)){
        D = D[D$start>=xlim[1] & D$stop <=xlim[2],]   
    }
    D
}

##########################################################################################


for( chrom in c( paste("chr",1:22,sep=""), "chrX", "chrY") ) {
    print(chrom)
    pmd_full_t = read.pmd.chromosome(chrom, sample_type="tumor", xlim=NULL)
    pmd_full_n = read.pmd.chromosome(chrom, sample_type="normal", xlim=NULL)
    pmd_full_loc = read.pmd.chromosome(chrom, sample_type="localized", xlim=NULL)
    print("Loaded")
    starts = sort(unique(pmd_full_t$start))
    pmd_full_t = pmd_full_t[order(pmd_full_t$start, pmd_full_t$sample_id),]
    pmd_full_n = pmd_full_n[order(pmd_full_n$start, pmd_full_n$sample_id),]
    pmd_full_loc = pmd_full_loc[order(pmd_full_loc$start, pmd_full_loc$sample_id),]
     
    mus_n = rep(NA, length(starts))
    mus_t = rep(NA, length(starts))
    mus_loc = rep(NA, length(starts))
    ns_t = length(sample_ids_wgbs)
    ns_l = length(samples_upitt_localized_tumor)
    ns_b = length(samples_normal_filecase)
    for(i in 1:length(starts)){
        if( i %% 5000 == 0 ){
            print(paste(i,"of",length(starts)))   
        }
        mus_t[i] = mean( pmd_full_t$mean[ ((ns_t*i)-(ns_t-1)) : (ns_t*i) ], na.rm=TRUE)
        mus_loc[i] = mean( pmd_full_loc$mean[ (( ns_l*i)-(ns_l-1)) : (ns_l*i) ], na.rm=TRUE)
        mus_n[i] = mean( pmd_full_n$mean[ (( ns_b*i)-(ns_b-1)) : (ns_l*i) ], na.rm=TRUE)
    }
    if( chrom == "chr1" ){
        res = data.frame( 
            chrom = rep(chrom, length(mus_n)),
            pos = starts,
            mu_n = mus_n,
            mu_t = mus_t,
            mu_loc = mus_loc,
            stringsAsFactors = FALSE)
    }else{
        res = rbind( res, 
                     data.frame( 
                         chrom = rep(chrom, length(mus_n)),
                         pos = starts,
                         mu_n = mus_n,
                         mu_t = mus_t,
                         mu_loc = mus_loc,
                         stringsAsFactors = FALSE)
        )
    }
}

pmd_mu_summary = res

write.table(pmd_mu_summary, fn_pmd_mu_summary,
            sep='\t', quote=FALSE, row.names=FALSE)

# Plot chromosome 16 
# ----------------------------------------


ntl = create_combined_pmd_plot( chrom="chr16", xlim=c(0,35000000), order="sum")
plot( ntl$data )
ggsave(filename = fn_figure_pmd_chr16, height=4, width=8)

# calculate means
LM = list()
for( chrom in dimnames(chrom_lengths)[[1]]){
    print(chrom)
    LM = c(LM, list( pmd_to_means(pmd_mu_summary, chrom, ncol=1000) ) )
}
# convert to colors
LC = list()
for(i in 1:24){
    chrom=dimnames(chrom_lengths)[[1]][ i ]
    LC = c(LC, list( medians_to_colors(LM[[i]], chrom, "#2c7fb8", "#7fcdbb", "#edf8b1") ) )
}

pdf( fn_figure_pmd_chrom_summary, height=5, width=10)
layout(matrix(1:24,24,1))
par(mar=c(0.4,1,0,1))
for(i in 1:24){
    if(i==1){
        par(mar=c(0.4, 1, 0, 1))
    }else{
        par(mar=c(0.4, 1, 0, 1 ) )   
    }
    plot.color.grid(LC[[i]], block.height = 10, block.width=1, space.X = 0, space.Y=2, 
                show.axis.X = FALSE, show.axis.Y=FALSE,border = FALSE)
}
dev.off()

