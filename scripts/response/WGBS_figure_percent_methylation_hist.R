sample_ids_wgs = c("DTB-003-BL", "DTB-005-BL", "DTB-008-BL", "DTB-009-BL", "DTB-011-BL", "DTB-018-BL", "DTB-019-PRO", "DTB-020-BL", "DTB-021-BL", 
                   "DTB-022-BL", "DTB-023-BL", "DTB-024-PRO", "DTB-034-BL", "DTB-035-BL", "DTB-036-BL", "DTB-037-BL", "DTB-040-BL", "DTB-042-BL", 
                   "DTB-053-BL", "DTB-055-PRO", "DTB-059-BL", "DTB-060-BL", "DTB-061-BL", "DTB-063-BL", "DTB-064-BL", "DTB-067-PRO", "DTB-069-BL", 
                   "DTB-071-BL", "DTB-074-BL", "DTB-077-PRO", "DTB-080-BL", "DTB-083-BL", "DTB-085-BL", "DTB-089-BL", "DTB-090-PRO", "DTB-091-BL", 
                   "DTB-092-BL", "DTB-094-BL", "DTB-097-PRO", "DTB-098-PRO2", "DTB-100-BL", "DTB-101-BL", "DTB-102-PRO", "DTB-104-BL", "DTB-111-PRO", 
                   "DTB-112-BL", "DTB-119-PRO", "DTB-121-BL", "DTB-124-BL", "DTB-126-BL", "DTB-127-PRO", "DTB-128-BL", "DTB-129-BL", "DTB-132-BL", 
                   "DTB-135-PRO", "DTB-137-PRO", "DTB-138-BL", "DTB-140-BL", "DTB-141-BL", "DTB-143-BL", "DTB-146-BL", "DTB-149-BL", "DTB-151-BL", 
                   "DTB-156-BL", "DTB-159-BL", "DTB-165-PRO", "DTB-167-PRO", "DTB-170-BL", "DTB-172-BL", "DTB-173-BL", "DTB-175-BL", "DTB-176-BL", 
                   "DTB-183-BL", "DTB-186-BL", "DTB-187-BL", "DTB-188-BL", "DTB-190-BL", "DTB-193-BL", "DTB-194-PRO", "DTB-201-PRO", "DTB-202-BL", 
                   "DTB-204-BL", "DTB-205-BL", "DTB-206-BL", "DTB-210-BL", "DTB-213-BL", "DTB-214-BL", "DTB-216-PRO", "DTB-220-BL", "DTB-222-BL", 
                   "DTB-223-BL", "DTB-232-PRO", "DTB-234-BL", "DTB-251-BL", "DTB-252-BL", "DTB-255-BL", "DTB-258-BL", "DTB-260-BL", "DTB-261-BL", 
                   "DTB-265-PRO", "DTB-266-BL")
sample_ids_wgbs = setdiff(sample_ids_wgs, "DTB-193-BL")

read_meth = function( chrom, sample){
    fn = paste( '/data1/projects/WCDT_WGBS_2019/methylseekr/input_data/tumors/',sample,'_methylseekr_',chrom,'.txt',sep='')
    M = read.table(fn,
                   sep='\t', stringsAsFactors=FALSE)
    names(M) = c("chrom", "pos", "nT", "nM")
    M
}


load_methylation_genomewide = function(sample){
    chroms = c( paste("chr",1:22,sep=''), "chrX", "chrY")
    for( chrom in chroms ){
        print(chrom)
        M=read_meth( chrom, sample )
        if( chrom=="chr1"){
            Mall=M
        }else{
            Mall = rbind(Mall, M)   
        }
    }
    Mall
}

sample = "DTB-213-BL"
M = load_methylation_genomewide(sample)
percents = M$nM/M$nT

hist( percents, main=sample,col="lightblue", border='blue', freq=FALSE, breaks=30,
      xlab="% methylated")

pdf('/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/DTB-213-BL.pdf',
    height=5, width=5)
hist( percents, main=sample,col="lightblue", border='blue', freq=FALSE, breaks=40,
      xlab="% methylated")
box()
dev.off()