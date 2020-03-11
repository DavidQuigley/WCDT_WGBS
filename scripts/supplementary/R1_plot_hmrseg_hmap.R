library(gplots)
library(wesanderson)
library(fpc)

#Plots fig 1B

filename <- 'C:/Users/Admin/Desktop/datasets/wgbs/rhmr.tsv'
df1 <- read.delim(filename,sep='\t',header=T,check.names=F)

filename <- 'C:/Users/Admin/Desktop/datasets/wgbs/rhmr_nt.tsv'
df2 <- read.delim(filename,sep='\t',header=T,check.names=F)

filename <- 'C:/Users/Admin/Desktop/datasets/wgbs/rhmr_upitt.tsv'
df3 <- read.delim(filename,sep='\t',header=T,check.names=F)

stopifnot(identical(df1$chr,df2$chr))
stopifnot(identical(df1$chr,df3$chr))
stopifnot(identical(df1$start,df2$start))
stopifnot(identical(df1$start,df3$start))
stopifnot(identical(df1$end,df2$end))
stopifnot(identical(df1$end,df3$end))

df <- cbind.data.frame(df1,df2[,-1:-3],df3[,-1:-3])

#samples2use <- c(samples_wgbs,samples_normal,upitt_benign_prostate,upitt_localized_tumor)
#samples2use <- c(upitt_benign_prostate,upitt_localized_tumor,samples_wgbs)
samples2use <- samples_wgbs
samples2use <- samples2use[samples2use %in% colnames(df)]

#sds <- apply(df[,samples2use],1,sd)
sds <- apply(df1[,intersect(colnames(df1),samples2use)],1,sd)
rows2keep <- sds>=quantile(sds,0.9,na.rm=T)
rows2keep[is.na(rows2keep)] <- F
df <- df[rows2keep,]

print('converting to matrix')
rowsreg <- rep(F,dim(df)[1])
gr <- makeGRangesFromDataFrame(df)
colsreg <- c('CRPC_ARE','LOC_ARE','chip_FOXA1','chip_HOXB13',
	'chip_ERG','chip_ERG_VCaP','promoter','genebody',
	'pca100_H3K27ac','pca100_H3K27me3')
for(i in 1:length(colsreg)) {
	rowsreg[countOverlaps(gr,tracks[[colsreg[i]]])>0] <- T
}
#rowsreg <- rep(T,dim(df)[1])
mat <- as.matrix(df[rowsreg,samples2use])
print(dim(mat))

rsidebar <- t(matrix(rep('white',dim(mat)[1])))

toplot <- data.frame(samples2use)
toplot$ERG <- allele_effect('ERG')$alleles[samples2use,'activating_sv']
toplot$ETV1 <- allele_effect('ETV1')$alleles[samples2use,'activating_sv']
toplot$ETV4 <- allele_effect('ETV4')$alleles[samples2use,'activating_sv']
toplot$ETV5 <- allele_effect('ETV5')$alleles[samples2use,'activating_sv']
toplot$CHD1 <- allele_effect('CHD1')$alleles[samples2use,'n_alleles_inactivated']==2
toplot$SPOP <- allele_effect('SPOP')$alleles[samples2use,'inactivating_missense']
toplot$RB1 <- allele_effect('RB1')$alleles[samples2use,'n_alleles_inactivated']==2
toplot$PTEN <- allele_effect('PTEN')$alleles[samples2use,'n_alleles_inactivated']==2
toplot$MYC <- allele_effect('MYC')$alleles[samples2use,'CNA_amp']
toplot$TP53 <- allele_effect('TP53')$alleles[samples2use,'n_alleles_inactivated']==2

CMP <- samples2use %in% samples_cmp
ETS <- toplot$ERG|toplot$ETV1|toplot$ETV4|toplot$ETV5
CHD1_SPOP <- toplot$CHD1 | toplot$SPOP
PTEN <- toplot$PTEN
RB1 <- toplot$RB1
TP53 <- toplot$TP53
MYC <- toplot$MYC
Mutation <- samples2use %in% c(samples_braf,samples_idh1,samples_tet2,samples_dnmt3b)
Sites <- locations[samples2use]
print(fisher.test(CMP,ETS))
print(fisher.test(CMP,TP53))
print(fisher.test(CMP,MYC))
print(fisher.test(CMP,Mutation))
print(fisher.test(CMP,Sites))

colside <- rep('white',dim(mat)[2])
names(colside) <- colnames(mat)
stopifnot(length(colside)==length(samples2use))

Sample <- colside
Sample[samples2use] <- wes_palette('BottleRocket2')[5] #gray
Sample[tscnc] <- wes_palette('BottleRocket2')[1] #yellow
#Sample[samples_normal] <- wes_palette('BottleRocket2')[3] #slate
#Sample[upitt_benign_prostate] <- wes_palette('BottleRocket2')[2] #red
#Sample[upitt_localized_tumor] <- wes_palette('BottleRocket2')[4] #darkgreen
stopifnot(length(Sample)==length(samples2use))

tSCNC <- colside
tSCNC[tscnc] <- 'gray30'
stopifnot(length(tSCNC)==length(samples2use))

ETS <- colside
ETS[toplot$ERG] <- 'gray30'
ETS[toplot$ETV1] <- 'gray30'
ETS[toplot$ETV4] <- 'gray30'
ETS[toplot$ETV5] <- 'gray30'
stopifnot(length(ETS)==length(samples2use))

MYC <- colside
MYC[toplot$MYC] <- 'gray30'
stopifnot(length(MYC)==length(samples2use))

SPOP <- colside
SPOP[toplot$SPOP] <- 'gray30'
stopifnot(length(SPOP)==length(samples2use))

RB1 <- colside
RB1[toplot$RB1] <- 'gray30'
stopifnot(length(RB1)==length(samples2use))

PTEN <- colside
PTEN[toplot$PTEN] <- 'gray30'
stopifnot(length(PTEN)==length(samples2use))

TP53 <- colside
TP53[toplot$TP53] <- 'gray30'
stopifnot(length(TP53)==length(samples2use))

CMP <- colside
CMP[samples_cmp] <- 'mediumturquoise'#'lightsteelblue'
CMP[samples_braf] <- 'orange'
CMP[samples_idh1] <- 'green'
CMP[samples_tet2] <- 'purple'
CMP[samples_dnmt3b] <- 'deeppink'
#CMP['DTB-252-BL'] <- 'black'
stopifnot(length(CMP)==length(samples2use))

Site <- colside
Site[intersect(names(locations)[locations=='Bone'],samples2use)] <- wes_palette('Chevalier1')[1] #darkgreen
Site[intersect(names(locations)[locations=='Lymph_node'],samples2use)] <- wes_palette('Chevalier1')[2] #gold
Site[intersect(names(locations)[locations=='Liver'],samples2use)] <- wes_palette('Chevalier1')[4] #brown
Site[intersect(names(locations)[locations=='Other'],samples2use)] <- wes_palette('Chevalier1')[3] #gray

purity <- colside
purcolname <- 'Tumor.Purity.Histo'
#purcolname <- 'Tumor.Purity.Comp'
tp <- tumorpurity[colnames(mat),purcolname]
names(tp) <- colnames(mat)
#tp[samples_normal] <- 0
#tp[upitt_benign_prostate] <- 0
#tp[upitt_localized_tumor] <- 80
tpcol <- colorRampPalette(c('gray100','black'))(n = 101)
purity <- tpcol[tp+1]
stopifnot(length(purity)==length(samples2use))

samples_noncmp <- setdiff(samples_wgbs,samples_cmp)
sampleorder <- c(#upitt_benign_prostate,
	#upitt_localized_tumor[c(1,3,4,5,2)],
	samples_cmp[order(tp[samples_cmp])],
	samples_noncmp[order(tp[samples_noncmp])])

csidebar <- cbind(purity,Site,tSCNC,CMP,TP53,RB1,PTEN,MYC,ETS)
color <- colorRampPalette(c('blue','white','red'))(n = 1000)

print('Plotting heatmap')
fileout <- 'C:/Users/Admin/Desktop/rhmr_hmap.pdf'
source('heatmap3.R')

hclust2 <- function(x, method="ward.D") {
	hclust(x, method=method)
}
dist2 <- function(x) {
	as.dist(1-cor(t(x), method="spearman"))
}

pdf(file=fileout,height=10,width=15)
main_title <- ""
heatmap.3(mat[,sampleorder], hclustfun=hclust, distfun=dist,
	na.rm=T, scale='none', dendrogram="none", margins=c(4,4),
	Rowv=T, Colv=F, ColSideColors=csidebar[sampleorder,], 
	RowSideColors=rsidebar, symbreaks=F, key=T, symkey=F,
	density.info="none", trace="none",labCol=colnames(mat)[sampleorder],
	labRow=rep("",(dim(mat)[1])), 
	main=main_title, cexRow=1, cexCol=0.8, col=color,
	ColSideColorsSize=8, RowSideColorsSize=1)
dev.off()

numk <- 12
cb <- clusterboot(t(mat),B=100,bootmethod='boot',clustermethod=hclustCBI,k=numk,method='complete',seed=7)
for(j in 1:numk) {
	samplesk <- names(cb$partition)[cb$partition==j]
	if(length(intersect(samplesk,samples_cmp))==length(samples_cmp)) {
		print('CMP')
		print(cb$bootmean[j])
	}
}






