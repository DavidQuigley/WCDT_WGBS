use2=FALSE
rangeup = 850000
rangedown = 100000
gene2find = 'AR'

rowensembl = ensembl2sym$name==gene2find
ensemblid = rownames(ensembl2sym)[rowensembl]
stopifnot(sum(rowensembl)==1)

chr = ensembl2sym[rowensembl,'chr']
gene_start = ensembl2sym[rowensembl,'start']
gene_end = ensembl2sym[rowensembl,'end']

start = round(gene_start-rangeup,-3)
end = round(gene_end+rangedown,-3)
stopifnot(start<end)

grange = GRanges(chr, IRanges(start, end))
ch = import.chain(fn_lift_38_to_19)
grangehg19 = liftOver(grange,ch)[[1]]

name = paste(gene2find,'loci',sep='_')
print(paste(name, ' - ', chr, ':', start, '-', end, sep=''))
numplots = 11
window = 100
maxwidth = 10000
cutoff = 5
within = 100

minscore = 2
cutgroups = 4
 
########################################################################
######################### LOAD IN VARIOUS DATA #########################
########################################################################

bedgenes = ensembl2sym[,c('chr','start','end','name','type','strand')]
bedgenes$strand = as.numeric(bedgenes$strand=='+')
bedgenes[bedgenes$strand==0,'strand'] = -1
colnames(bedgenes) = c('chrom','start','stop','gene','score','strand')

exons = read.delim(file=fn_gencode_exons, header=TRUE, stringsAsFactors = F) 
colnames(exons)[1] = 'exon_ID'
bedgene = exons[,c(5,6,7,3,1,4)]
bedgene$strand = as.numeric(bedgene$strand=='+')
bedgene[bedgene$strand==0,'strand'] = -1

## CpG islands
cgi = as.data.frame(tracks[['cpg_island']])[1:3]
colnames(cgi)[1] = 'chr'
cgi$row = 4

## CHIPseq tracks
chrXlist = list()
dir_hg19 = dir_ar_enhancer_bw
files_hg19 = list.files(dir_hg19,pattern='.bw')

coords = start:end
h3k27ac_t = rep(0,length(coords))
h3k27ac_m = rep(0,length(coords))
hoxb13 = rep(0,length(coords))
foxa1 = rep(0,length(coords))

for(i in 1:(length(files_hg19)-2)) {
    filei = paste(dir_hg19,files_hg19[i],sep='/')
    chrx = import.bw(filei,as='NumericList',
                     selection=BigWigSelection(ranges=grangehg19))[[1]]
    stopifnot(length(coords)==length(chrx))
    if(i %in% 1:6) {
        h3k27ac_t = h3k27ac_t + chrx
    } else if(i %in% 7:10) {
        h3k27ac_m = h3k27ac_m + chrx
    } else if(i %in% 11:14) {
        hoxb13 = hoxb13 + chrx
    } else if(i %in% 15:18) {
        foxa1 = foxa1 + chrx
    }
}
h3k27ac_t = h3k27ac_t/6
h3k27ac_m = h3k27ac_m/4
hoxb13 = hoxb13/4
foxa1 = foxa1/4

h3k27ac_t = cbind.data.frame(chr,coords,coords,h3k27ac_t)
colnames(h3k27ac_t) = c('chrom','start','end','value')
h3k27ac_m = cbind.data.frame(chr,coords,coords,h3k27ac_m)
colnames(h3k27ac_m) = c('chrom','start','end','value')
hoxb13 = cbind.data.frame(chr,coords,coords,hoxb13)
colnames(hoxb13) = c('chrom','start','end','value')
foxa1 = cbind.data.frame(chr,coords,coords,foxa1)
colnames(foxa1) = c('chrom','start','end','value')

erg = read.delim(fn_chipseq_erg_hg38, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(erg) = c('chrom','start','end','value')

ergvcap = read.delim(fn_chipseq_erg_vcap, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(ergvcap) = c('chrom','start','end','value')

## MANTA tracks
# print('Sliding windows for SVs')
# fusionrows = (data_manta_slim$chr1==chr & inrange(data_manta_slim$pos1,start,end)) |
#     (data_manta_slim$chr2==chr & inrange(data_manta_slim$pos2,start,end)) |
#     (data_manta_slim$chr1==data_manta_slim$chr2 & data_manta_slim$chr1==chr &
#          overlap(data_manta_slim$pos1,data_manta_slim$pos2,start,end)) &
#     data_manta_slim$sample_id %in% sample_ids_wgbs
# fusionrows = fusionrows & rows_manta_pass
# fusionrows[is.na(fusionrows)] = FALSE
# 
# cols2keep = c('chr1','pos1','chr2','pos2','svtype','sample_id')
# plot_sv = data_manta_slim[fusionrows,cols2keep]
# rowsdel = plot_sv$svtype=='MantaDEL'
# rowsinv = plot_sv$svtype=='MantaINV'
# rowsins = plot_sv$svtype=='MantaINS'
# rowsbnd = plot_sv$svtype=='MantaBND'
# rowsdup = plot_sv$svtype=='MantaDUP'
# bed_mantadel = svwindow(plot_sv[rowsdel,],chr,start,end,endsonly=FALSE, extra=0, window=1000)
# bed_mantainv = svwindow(plot_sv[rowsinv,],chr,start,end, extra=0, window=1000)
# bed_mantabnd = svwindow(plot_sv[rowsbnd,],chr,start,end, extra=0, window=1000)
# bed_mantadup = svwindow(plot_sv[rowsdup,],chr,start,end,endsonly=FALSE, extra=0, window=1000)
# bed_mantains = svwindow(plot_sv[rowsins,],chr,start,end, extra=0, window=1000)

c('chr1','pos1','chr2','pos2','svtype','sample_id')
dups = list_sv_m[list_sv_m$svtype=="TANDEM",]
dups = dups[, c('chrom_start','pos_start','chrom_end', 'pos_end', 'svtype', 'sample_id')] 
names(dups) = c('chr1','pos1','chr2','pos2','svtype','sample_id')
bed_mantadup = svwindow(dups,chr,start,end,endsonly=FALSE, extra=0, window=1000)

## CN tracks
rowsbedcn = copycat_bed$sample_id %in% sample_ids_wgbs & copycat_bed$chr==chr
rowsbedgain = rowsbedcn & copycat_bed$call=='GAIN'
rowsbedloss = rowsbedcn & copycat_bed$call=='LOSS'
bed_gain = slidingwindow(copycat_bed[rowsbedgain,-4:-6],chr,start,end, extra=0, window=1000)
bed_loss = slidingwindow(copycat_bed[rowsbedloss,-4:-6],chr,start,end, extra=0, window=1000)

## Mutation tracks
rowsannovar = annovar_full$sample %in% sample_ids_wgbs &
    overlap(annovar_full$Start,annovar_full$End,start,end) &
    annovar_full$Chr==chr
bed_mut_temp = annovar_full[rowsannovar,c(1:3)]
colnames(bed_mut_temp) = c('chrom','start','end')
bed_mut_temp$value = 1
bed_mut = slidingwindow(bed_mut_temp,chr,start,end,FALSE, extra=0, window=1000)

# Load UMRs
print('LOAD UMRS')
hmrbin = data.frame(bins=1:((end-start)/window))
hmrbin$chr = chr
hmrbin$start = start+((hmrbin$bins-1)*window)
hmrbin$end = start+(hmrbin$bins*window)
bins = makeGRangesFromDataFrame(hmrbin)
for(i in 1:length(hmr)) {
    sample = names(hmr)[[i]]
    if(sample %in% sample_ids_wgbs) {
        rowshmri = width(hmr[[i]]) <= maxwidth
        hmrbin[,sample] = countOverlaps(bins,hmr[[i]][rowshmri])>0
    }
}
hmrbin$sum = rowSums(hmrbin[,(1:length(sample_ids_wgbs))+(dim(hmrbin)[2]-length(sample_ids_wgbs))])
hmrbin = hmrbin[,c('chr','start','end','sum')]
colnames(hmrbin) = c('chrom','start','end','value')
peaks = findpeaks(hmrbin,cutoff,minwidth=window,within=within,minsize=cutoff)
peaks$chr = chr
numpeaks = dim(peaks)[1]
print(paste('peaks',numpeaks))

file_methyl = list.files(dir_gene_output_mcrpc, pattern=ensemblid)
stopifnot(length(file_methyl)==1)
filenamefull = paste(dir_gene_output_mcrpc,file_methyl,sep='/')
methyl = read.delim(filenamefull,sep='',row.names=1,header=TRUE,stringsAsFactors=FALSE,check.names=F)
colnames(methyl) = toupper(str_split_fixed(colnames(methyl),'\\.',2)[,1])
stopifnot(dim(methyl)[1]>0)
methylpos = as.numeric(rownames(methyl))
expr = as.numeric(tpm[ensemblid,colnames(methyl)])
expr = log2(expr+1)
names(expr) = colnames(methyl)

filenames = list.files(dir_methylseqR, pattern=paste(chr,'UMRsLMRs.txt',sep='_'), full.names=F)
all_hmrs = list()
for(curfile in filenames) {
    nametab = strsplit(curfile,'/',fixed=T)[[1]]
    samplename = strsplit(nametab[length(nametab)],'_',fixed=T)[[1]][1]
    curfile = paste(dir_methylseqR,curfile,sep='/')
    curtab = read.delim(file=curfile,sep='\t',header=TRUE,stringsAsFactors=F)
    curtab$sample = toupper(samplename)
    all_hmrs[[samplename]] = curtab
}
hmrdf = rbindlist(all_hmrs)
rowshmr = overlap(hmrdf$start,hmrdf$end,start,end) & hmrdf$sample %in% 
    c(sample_ids_wgbs,samples_normal,samples_upitt_benign_prostate,toupper(samples_upitt_localized_tumor))

# Load EMRs
if(use2) {
    if(gene2find=='ERG') {
        samplesgroup = sample_ids_wgbs[allele_effect('ERG')$alleles[colnames(methyl),'activating_sv']]
    } else if (gene2find=='RB1') {
        samplesgroup = sample_ids_wgbs[allele_effect('RB1')$alleles[colnames(methyl),'n_alleles_inactivated']==2]
    } else if (gene2find=='TP53') {
        samplesgroup = sample_ids_wgbs[allele_effect('TP53')$alleles[colnames(methyl),'n_alleles_inactivated']==2]
    }
    samplesnotgroup = setdiff(colnames(methyl),samplesgroup)
    print(range(expr[samplesgroup]))
    print(range(expr[samplesnotgroup]))
    for(i in 1:numpeaks) {
        print(i)
        starti = peaks[i,'start']
        endi = peaks[i,'end']
        methylrows = methylpos >= starti & methylpos <= endi
        methylslice = methyl[methylrows,]
        means1 = colMeans(methylslice[,samplesgroup],na.rm=T)
        corobj1 = cor.test(means1,expr[samplesgroup],method='spearman')
        peaks[i,'avg_p1'] = corobj1$p.value
        peaks[i,'avg_cor1'] = corobj1$estimate
        means2 = colMeans(methylslice[,samplesnotgroup],na.rm=T)
        corobj2 = cor.test(means2,expr[samplesnotgroup],method='spearman')
        peaks[i,'avg_p2'] = corobj2$p.value
        peaks[i,'avg_cor2'] = corobj2$estimate
    }
    peakscols1 = c('chr','start','end','avg_cor1')
    peakscols2 = c('chr','start','end','avg_cor2')
    peaks$fdr = p.adjust(peaks$avg_p1,method='fdr')
    peaks$rowsig = peaks$avg_cor1<0 & peaks$fdr<0.05
}else{
    for(i in 1:numpeaks) {
        print(i)
        starti = peaks[i,'start']
        endi = peaks[i,'end']
        methylrows = methylpos >= starti & methylpos <= endi
        if(sum(methylrows)>0) {
            methylslice = methyl[methylrows,]
            means = colMeans(methylslice,na.rm=T)
            corobj = cor.test(means,expr,method='spearman')
            #corobj = cor.test(means[samples_noargain],expr[samples_noargain],method='spearman')
            peaks[i,'avg_p'] = corobj$p.value
            peaks[i,'avg_cor'] = corobj$estimate
        }
    }
    peakscols = c('chr','start','end','avg_cor')
    peaks$fdr = p.adjust(peaks$avg_p,method='fdr')
    peaks$rowsig = peaks$avg_cor<0 & peaks$fdr<0.05
}
print(peaks)

# Load AREs
crpc_are = read.delim(fn_crpc_are, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=TRUE)
colnames(crpc_are) = c('chrom','start','end','score')
crpc_are$score = crpc_are$score/5

loc_are = read.delim(fn_loc_are, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(loc_are) = c('chrom','start','end','score')
loc_are$score = loc_are$score/3

# Load PC 100 chipseq
pc100_are = read.delim(fn_pc100_are, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(pc100_are) = c('chrom','start','end','score')

pc100_h3k4me3 = read.delim(fn_pc100_h3k4me3, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(pc100_h3k4me3) = c('chrom','start','end','score')

pc100_h3k27ac = read.delim(fn_pc100_h3k27ac, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(pc100_h3k27ac) = c('chrom','start','end','score')

# Load AR ERG Chiapet
print('Loading chiapet')
ch = import.chain(fn_lift_19_to_38)

chiapet_ar = readchiapet(fn_chiapet_ar)
chiapet_ar_bind = read.delim(fn_chia_ar_bind, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(chiapet_ar_bind) = c('chrom','start','end','score')

chiapet_erg = readchiapet(fn_chiapet_erg)
chiapet_erg_bind = read.delim(fn_chiapet_erg_bind, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=T)
colnames(chiapet_erg_bind) = c('chrom','start','end','score')

# Load DMRs
#tvnslice = data.frame(tvn[countOverlaps(tvn,grange)>0])
#tvnsliceall = tvnslice[order(tvnslice$start),c(1:3,10)]
locslice = data.frame(dss_loc[countOverlaps(dss_loc,grange)>0])
locslice = locslice[order(locslice$start),c(1:3,10)]
metslice = data.frame(dss_met[countOverlaps(dss_met,grange)>0])
metslice = metslice[order(metslice$start),c(1:3,10)]

## Load HiC
hic_lncap = readhic(fn_hic_lncap)
hic_pc3 = readhic(fn_hic_pc3)
hic_prec = readhic(fn_hic_prec)


#######################################################################
############################# PLOT TRACKS #############################
#######################################################################

pdf(file=fn_figure4c, onefile=TRUE, width=12, height=10)
matrows = c(1,1,rep(2,10),3,3,4:numplots)
layout(matrix(matrows, length(matrows), 1, byrow=TRUE))
par(mai=c(0.05,1,0.05,0.5))

## Plot genes
par(mai=c(0.4,1,0.05,0))
exonrows = bedgene$chr==chr & overlap(bedgene$start, bedgene$end, start, end)
if(sum(exonrows)>0) {
    plotGenes(bedgene[exonrows,],chr,start,end)
    labelgenome(chr,start,end,n=20,scale="Mb")
} else {
    plot.new()
}

## Plot HMR
par(mai=c(0.05,1,0.05,0))
hmrslice = hmrdf[rowshmr,c(1:3,9)]
hmrslice$start = pmax(hmrslice$start,start+1)
hmrslice$end = pmin(hmrslice$end,end-1)
colnames(hmrslice) = c('chrom','start','end','sample')
hmrslice$score = 1
hmrslice$strand = 1

if(use2) {
    hmrslice$color = 'thistle'
    rowsgroup = hmrslice$sample %in% samplesgroup
    hmrslice[rowsgroup,'color'] = 'black'
} else {
    hmrslice$color = 'black'
}
rows_loc = hmrslice$sample %in% samples_upitt_localized_tumor
hmrslice[rows_loc,'color'] = 'gray50'
rows_ben = hmrslice$sample %in% samples_upitt_benign_prostate
hmrslice[rows_ben,'color'] = 'gray80'
rows_nt = hmrslice$sample %in% samples_normal
hmrslice[rows_nt,'color'] = 'slategray'


if(use2) {
    roworder = as.numeric(factor(hmrslice$sample,
                                  levels=c(samples_normal,samples_upitt_benign_prostate,
                                           samples_upitt_localized_tumor,
                                           samplesgroup,samplesnotgroup)))
} else {
    roworder = as.numeric(factor(hmrslice$sample,
                                  levels=c(samples_normal,samples_upitt_benign_prostate,
                                           samples_upitt_localized_tumor,sample_ids_wgbs)))
}
rows2keep = !is.na(roworder)
hmrslice = hmrslice[rows2keep,]
roworder = roworder[rows2keep]
plotBed(hmrslice,chr,start,end,row='supplied',rownumber=roworder,color=hmrslice$color)
axis(side=2,las=2,tcl=.2)
mtext('HMR',side=2,line=4,cex=0.5)

## Plot EMR
emr_range = c(-0.8,0.6)
if(use2) {
    cols = c('black','thistle')
    plotBedgraph(peaks[,peakscols2],chr,start,end,color=cols[2],range=emr_range,transparency=0.6)
    plotBedgraph(peaks[,peakscols1],chr,start,end,color=cols[1],range=emr_range,transparency=0.6,overlay=T)
    legend("topright",inset=0.025,legend=c('altered','wt'),
           fill=opaque(cols),border=cols,cex=0.5)
} else {
    plotBedgraph(peaks[,peakscols],chr,start,end,color='black',range=emr_range)
}
axis(side=2,las=2,tcl=.2)
mtext('Rho',side=2,line=4,cex=0.5)
abline(h=0, col='black')

## Plot chipseq
cols = c('red','gray50')
if(chr=='chrX' && length(grangehg19)==1) {
    plotBedgraph(h3k27ac_m,chr,start,end,color=cols[1])
    plotBedgraph(h3k27ac_t,chr,start,end,color=cols[2],overlay=T)
    axis(side=2,las=2,tcl=.2)
    mtext('H3K27ac',side=2,line=4,cex=0.5)
    legend("topright",inset=0.025,legend=c('mCRPC','Primary'),
           fill=opaque(cols),border=cols,cex=0.5)
} else {
    plotBedgraph(cres,chr,start,end,color='red')
    axis(side=2,las=2,tcl=.2)
    mtext('H3K27ac',side=2,line=4,cex=0.5)
}
abline(h=0, col='black')

plotBedgraph(hoxb13,chr,start,end,color='darkgreen')
axis(side=2,las=2,tcl=.2)
mtext('HOXB13',side=2,line=4,cex=0.5)
abline(h=0, col='black')

plotBedgraph(foxa1,chr,start,end,color='darkgreen')
axis(side=2,las=2,tcl=.2)
mtext('FOXA1',side=2,line=4,cex=0.5)
abline(h=0, col='black')

cols = c('darkgreen','gray50')
plotBedgraph(crpc_are,chr,start,end,color=cols[1])
plotBedgraph(loc_are,chr,start,end,color=cols[2],overlay=TRUE)
axis(side=2,las=2,tcl=.2)
mtext("ARE",side=2,line=4,cex=0.5)
abline(h=0, col='black')
legend("topright",inset=0.025,legend=c('mCRPC','Primary'),
       fill=opaque(cols),border=cols,cex=0.5)

cols = c('darkgreen','black')
plotBedgraph(ergvcap,chr,start,end,color=cols[1])
plotBedgraph(erg,chr,start,end,color=cols[2],overlay=TRUE)
axis(side=2,las=2,tcl=.2)
mtext('ERG',side=2,line=4,cex=0.5)
abline(h=0, col='black')
legend("topright",inset=0.025,legend=c('VCaP','Metastasis'),
       fill=opaque(cols),border=cols,cex=0.5)

plotBedpe(chiapet_ar,chr,start,end,heights=chiapet_ar$score,plottype='loops',
              color='navy',flip=F)
axis(side=2,las=2,tcl=.2)
mtext('AR',side=2,line=4,cex=0.5)
abline(h=0, col='black')


par(mai=c(0.05,1,0,0))

plotBedpe(chiapet_erg,chr,start,end,heights=chiapet_erg$score,plottype='loops',
          color='royalblue3',flip=T)

axis(side=2,las=2,tcl=.2)
mtext('ERG',side=2,line=4,cex=0.5)
abline(h=0, col='black')


plotBedgraph(bed_mantadup,chr,start,end,color='plum4')
axis(side=2,las=2,tcl=.2)
mtext("Tandem Dup",side=2,line=4,cex=0.5)
dev.off()