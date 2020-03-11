chroms = c( paste("chr",1:22,sep=""), "chrX", "chrY" )
beds = data.frame()

for(i in 1:24){
    locs = seq(from=1, to=chrom_lengths$V2[i] - 100000, by=100000) 
    bb = data.frame( 
        chrom = rep(chroms[i], length(locs)),
        start = locs,
        end = locs+100000,
        stringsAsFactors = FALSE)
    if(i == 1 ){
        beds = bb
    }else{
        beds = rbind(beds, bb)   
    }
}
write.table(beds, 
            '/notebook/human_sequence_prostate_WGBS/drafts/response/BEDs_for_PMD_tiling_100K.txt', 
            sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)


beds = data.frame()
for(i in 1:24){
    locs = seq(from=1, to=chrom_lengths$V2[i] - 10000, by=10000) 
    bb = data.frame( 
        chrom = rep(chroms[i], length(locs)),
        start = locs,
        end = locs+10000,
        stringsAsFactors = FALSE)
    if(i == 1 ){
        beds = bb
    }else{
        beds = rbind(beds, bb)   
    }
}
write.table(beds, 
            '/notebook/human_sequence_prostate_WGBS/drafts/response/BEDs_for_PMD_tiling_10K.txt', 
            sep='\t', quote=FALSE, row.names = FALSE, col.names = FALSE)

