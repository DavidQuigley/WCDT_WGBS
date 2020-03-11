#Analysis of partially methylated domains (PMDs) to confirm previous findings in other cancer types
#Martin Sjostrom
#2019-09-20

# Prepare environment -----------------------------------------------------

source("/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/scripts/2019_05_15_prepare_environment.r")

# Libraries ---------------------------------------------------------------
library(pryr)

# Directories and filenames -----------------------------------------------
DIR_METHYLSEEKR <- "/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/MethylSeekR/" 
DIR_PMD_MEANS <- "/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/secondary_data/partially_methylated_domains/output/sample_PMDs/"

OUTDIR <- "/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/PMD"
OUTDIR_BED <- "/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/PMD/bed"
OUTDIR_GEX_IN_PMD <- "/data1/projects/WCDT_WGBS_2019/WCDT_WGBS/figures/PMD/gex_in_pmd"


#fn_tumor_purity <- "F:/Shared VM external/WGBS/build_2019_04_05_v2/metadata/WGS/WGS_tumor_purity_AF_012318.csv"
#fn_exons <- "F:/Shared VM external/WGBS/build_2019_04_05_v2/metadata/gencode28_exons.txt"
#fn_genes <- "F:/Shared VM external/WGBS/build_2019_04_05_v2/metadata/ensembl2sym.txt"
#fn_WGS_data <- "F:/Shared VM external/WGBS/build_2019_04_05_v2/metadata/WGS/2018_04_25_build_image.RData"

# Constants and gene regions ----------------------------------------------
#load exons and genes
exons <- read.table(fn_gencode_exons)
exons.gr <- makeGRangesFromDataFrame(exons, keep.extra.columns = T)

coding.genome <- reduce(exons.gr, ignore.strand = T) #this removes directionality
coding.genome <- coding.genome[!grepl("chrM",seqnames(coding.genome))] #remove M

stopifnot(all(rownames(ensembl2sym)==rownames(tpm)))
genes <- ensembl2sym[rownames(ensembl2sym) %in% rownames(tpm),] #keep genes for which we have gene expression data
genes.gr <- makeGRangesFromDataFrame(genes,keep.extra.columns = T)
genes.gr.protein <- genes.gr[genes.gr$type %in% gencode_protein]

#colors
col_ben_loc_mcrpc <- c("darkolivegreen2","goldenrod2","indianred2")
col_no_pmd_pmd <- c("deepskyblue1","brown1")
col_noCMP_CMP <- c("palegreen2","tomato2")

# Sample annotations ------------------------------------------------------
#keep upitt and two normals as mixed case to match filenames, 
samples_upitt_benign_prostate = c('Bis49AT','Bis159AT','Bis165AT','Bis171AT')
samples_upitt_localized_tumor <- c('Bis49T','Bis158T','Bis159T','Bis165T','Bis171T')
samples_normal[samples_normal %in% "NT-DTB-077-PRO"] <- "NT-DTB-077-Pro"
samples_normal[samples_normal %in% "NT-DTB-194-PRO"] <- "NT-DTB-194-Pro"
included_samples <- c(sample_ids_wgbs, samples_upitt_localized_tumor, samples_upitt_benign_prostate, samples_normal)

# Load and organize data --------------------------------------------------
# Load and PMD information as from original methylseekr run
fileNames_methylseekr <- list.files(DIR_METHYLSEEKR,"*PMD.rds")
sampleNames <- unique(gsub("_methylseekr.+$","",fileNames_methylseekr))
stopifnot(all(included_samples %in% sampleNames) & all(sampleNames %in% included_samples))

setwd(DIR_METHYLSEEKR)
all.samples <- list()

for(i in 1:length(sampleNames)){
  sampleName <- sampleNames[i]
  filesToRead <- fileNames_methylseekr[grep(sampleName,fileNames_methylseekr)]
  allFiles <- list()
  
  for(j in 1:length(filesToRead)){
    file_name <- filesToRead[j]
    newFile <- readRDS(file_name)
    allFiles[[j]] <- newFile
  }
  
  allFiles <- do.call(c,unlist(allFiles))
  all.samples[[i]] <- allFiles
  names(all.samples)[i] <- sampleName
}

#export a bedfile per sample for extraction of mean methylation levels
exportBed <- function(gr,sampleName){
  export(object = gr[gr$type=="PMD"], con = paste0(sampleName,".PMD.bed"),format = "bed")
}
# setwd(OUTDIR_BED)
# for(i in 1:length(sampleNames)){
#   sampleName <- sampleNames[i]
#   gr <- all.samples[[sampleName]]
#   exportBed(gr,sampleName)
# }

#Load mean methylation in PMDs
#This was extracted using the bed files created above
fileNames <- list.files(DIR_PMD_MEANS,"*.txt")
sampleNames <- unique(gsub("_methylseekr.+$","",fileNames))

length(sampleNames) #one is missing

pmdMeth <- list()

setwd(DIR_PMD_MEANS)
for(i in 1:length(sampleNames)){
  sampleName <- sampleNames[i]
  print(paste(i, length(sampleNames)))
  pmdData <- data.frame(matrix(NA,ncol = 11, nrow = 0))
  colnames(pmdData) <- c("sample_id","chrom","start","stop","score","name","feature_type","median","mean","median_coverage","n_cpg")
  
  for(chrom in chromosomes){
    newData <- as.data.frame(fread(file = paste0(sampleName,"_methylseekr_",chrom,"_counts.txt")))
    pmdData <- rbind(pmdData,newData)
  }
  pmdMeth[[i]] <- pmdData
  names(pmdMeth)[i] <- sampleName 
}

# Calculate % genome PMD and mean methylation -----------------------------
# Potential further considerations: 
# -filter away centromeric (and telomeric) regions
# -remove CGIs (and shores) and promoters from calculation of mean PMD methylation levles ("weighted methylation")

#Calculate percentage of genmoe with PMDs for each sample
calcPercPMD <- function(sampleName){
  percPMD <- sum(width(sampleName)[sampleName$type=="PMD"])/sum(width(sampleName))
  return(percPMD)
}

percPMD <- unlist(lapply(all.samples,calcPercPMD))

#Calculate mean level of methylation in PMDs (the mean of all PMD means)
meanPMD <- unlist(lapply(pmdMeth, FUN = function(x){mean(x[,"mean"],na.rm=T)}))

# Overview and correlation plots ------------------------------------------

#### Histogram/boxplot %PMD (breast cancer has 10-50% with mean of 32%)
#### Boxplot PMD level in benign, localized, mCRPC, not our adjacent normal
fig.perPMD %<a-% {
  group <- ifelse(names(percPMD) %in% sample_ids_wgbs,3,
                  ifelse(names(percPMD) %in% samples_normal,NA,
                         ifelse(names(percPMD) %in% samples_upitt_localized_tumor,2,
                                ifelse(names(percPMD) %in% samples_upitt_benign_prostate,1,NA))))
  boxplot(percPMD~group,
          names = c("Benign prostate", "Localized PCa","mCRPC"), 
          xlab = "",
          ylab = "Prop. of genome with PMD",
          col = col_ben_loc_mcrpc,
          main = "Fraction of genome in PMDs")
}

fig.meanPMDMeth %<a-% {
  group <- ifelse(names(meanPMD) %in% sample_ids_wgbs,3,
                  ifelse(names(meanPMD) %in% samples_normal,NA,
                         ifelse(names(meanPMD) %in% samples_upitt_localized_tumor,2,
                                ifelse(names(meanPMD) %in% samples_upitt_benign_prostate,1,NA))))
  boxplot(meanPMD~group,
          names = c("Benign prostate", "Localized PCa", "mCRPC"), 
          xlab = "",
          ylab = "Mean PMD methylation",
          col = col_ben_loc_mcrpc,
          main = "Mean PMD methylation per disease stage")
  ps <- pairwise.wilcox.test(meanPMD,group, p.adjust.method = 'none')
  text(2,0.55,labels = paste0("Benign - localized: p = ",signif(ps$p.value["2","1"],2),"\n",
                              "Localized - mCRPC: p = ",signif(ps$p.value["3","2"],2),"\n",
                              "Benign - mCRPC: p = ",signif(ps$p.value["3","1"],2)))
}


#### Correlate PMD level with %TC content
percPMD.mcrpc <- percPMD[names(percPMD) %in% sample_ids_wgbs]

fig.PMDvsTCcont %<a-%{
  percPMD.mcrpc <- percPMD[names(percPMD) %in% sample_ids_wgbs]

  s.cor <- signif(cor(tumorpurity[names(percPMD.mcrpc),"Tumor.Purity.Histo"],percPMD.mcrpc,method = "spearman"),2)
  pval <- signif(cor.test(tumorpurity[names(percPMD.mcrpc),"Tumor.Purity.Histo"],percPMD.mcrpc,method = "spearman")$p.value,2)
  lmod <- summary(lm(percPMD.mcrpc~tumorpurity[names(percPMD.mcrpc),"Tumor.Purity.Histo"]))
  #pval <- signif(lmod$coefficients[2,4],2) #this corresponds to cor.test with method="pearson"?
  interc <- lmod$coefficients[1,1]
  slope <- lmod$coefficients[2,1]
  
  plot(tumorpurity[names(percPMD.mcrpc),"Tumor.Purity.Histo"],percPMD.mcrpc,
       ylab = "Prop. genome with PMD",
       xlab = "Tumor purity, histology assessment",
       main = "Tumor purity and PMDs")
  #abline(a = interc,b = slope,col=2,lty=2)
  text(37,0.23,labels = paste0("P-value = ",pval,"\n Spearman's Rho = ",s.cor))
}

meanPMD.mcrpc <- meanPMD[names(meanPMD) %in% sample_ids_wgbs]

fig.PMDMeanvsTCcont %<a-%{
  meanPMD.mcrpc <- meanPMD[names(meanPMD) %in% sample_ids_wgbs]
  
  s.cor <- signif(cor(tumorpurity[names(meanPMD.mcrpc),"Tumor.Purity.Histo"],meanPMD.mcrpc,method = "spearman"),2)
  pval <- signif(cor.test(tumorpurity[names(meanPMD.mcrpc),"Tumor.Purity.Histo"],meanPMD.mcrpc,method = "spearman")$p.value,2)
  lmod <- summary(lm(meanPMD.mcrpc~tumorpurity[names(meanPMD.mcrpc),"Tumor.Purity.Histo"]))
  #pval <- signif(lmod$coefficients[2,4],2) #this corresponds to cor.test with method="pearson"?
  interc <- lmod$coefficients[1,1]
  slope <- lmod$coefficients[2,1]
  
  plot(tumorpurity[names(meanPMD.mcrpc),"Tumor.Purity.Histo"],meanPMD.mcrpc,
       ylab = "Mean PMD methylation",
       xlab = "Tumor purity, histological assessment",
       main = "Tumor purity and PMDs")
  #abline(a = interc,b = slope,col=2,lty=2)
  text(50,0.50,labels = paste0("P-value = ",pval,"\n Spearman's Rho = ",s.cor))
}

#### Correlate PMD % genome with TMB
fig.PMDvsMuts %<a-%{
  tmb <- matrix_mutcount[names(percPMD.mcrpc),]
  tmb <- tmb[!rownames(tmb) %in% samples_hypermut,] #remove 2 outlier samples
  percPMD.mcrpc.tmb <- percPMD.mcrpc[rownames(tmb)]
  stopifnot(all(names(percPMD.mcrpc.tmb)==rownames(tmb)))
  
  s.cor <- signif(cor(percPMD.mcrpc.tmb,tmb$mutation_count,method = "spearman"),2)
  pval <- signif(cor.test(percPMD.mcrpc.tmb,tmb$mutation_count,method = "spearman")$p.value,2)
  lmod <- summary(lm(tmb$mutation_count~percPMD.mcrpc.tmb))
  #pval <- signif(lmod$coefficients[2,4],2)
  interc <- lmod$coefficients[1,1]
  slope <- lmod$coefficients[2,1]
  
  plot(percPMD.mcrpc.tmb,tmb$mutation_count,
       xlab = "Prop. genome with PMD",
       ylab = "Number of mutations",
       main = "Number of mutations and prop. genome with PMDs")
  #abline(a = interc,b = slope,col=2,lty=2)
  text(0.27,30000,labels = paste0("P-value = ",pval,"\n Spearman's Rho = ",s.cor))
}

#### Correlate PMD methylation level with TMB
fig.PMDMeanvsMuts %<a-%{
  tmb <- matrix_mutcount[names(meanPMD.mcrpc),]
  tmb <- tmb[!rownames(tmb) %in% samples_hypermut,] #remove 2 outlier samples
  meanPMD.mcrpc.tmb <- meanPMD.mcrpc[rownames(tmb)]
  stopifnot(all(names(meanPMD.mcrpc.tmb)==rownames(tmb)))
  
  s.cor <- signif(cor(meanPMD.mcrpc.tmb,tmb$mutation_count,method = "spearman"),2)
  pval <- signif(cor.test(meanPMD.mcrpc.tmb,tmb$mutation_count,method = "spearman")$p.value,2)
  lmod <- summary(lm(tmb$mutation_count~meanPMD.mcrpc.tmb))
  #pval <- signif(lmod$coefficients[2,4],2)
  interc <- lmod$coefficients[1,1]
  slope <- lmod$coefficients[2,1]
  
  plot(meanPMD.mcrpc.tmb,tmb$mutation_count,
       xlab = "Mean PMD methylation",
       ylab = "Number of mutations",
       main = "Number of mutations and PMD methylation level")
  #abline(a = interc,b = slope,col=2,lty=2)
  text(0.55,32000,labels = paste0("P-value = ",pval,"\n Spearman's Rho = ",s.cor))
}

#### Correlate PMD % genome with %CN altered
fig.PMDvsCNA %<a-%{
  pcna <- percent_CNA[names(percPMD.mcrpc)]
  
  s.cor <- signif(cor(percPMD.mcrpc,pcna,method = "spearman"),2)
  pval <- signif(cor.test(percPMD.mcrpc,pcna,method = "spearman")$p.value,2)
  lmod <- summary(lm(pcna~percPMD.mcrpc))
  #pval <- signif(lmod$coefficients[2,4],2)
  interc <- lmod$coefficients[1,1]
  slope <- lmod$coefficients[2,1]
  
  plot(percPMD.mcrpc,pcna,
       xlab = "Mean PMD methylation",
       ylab = "% CNA",
       main = "% CN altered and prop. genome with PMDs")
  #abline(a = interc,b = slope,col=2,lty=2)
  text(0.27,0.43,labels = paste0("P-value = ",pval,"\n Spearman's Rho = ",s.cor))
}

#### Correlate PMD mean level with %CN altered
fig.PMDMeanvsCNA %<a-%{
  pcna <- percent_CNA[names(meanPMD.mcrpc)]
  
  s.cor <- signif(cor(meanPMD.mcrpc,pcna,method = "spearman"),2)
  pval <- signif(cor.test(meanPMD.mcrpc,pcna,method = "spearman")$p.value,2)
  lmod <- summary(lm(pcna~meanPMD.mcrpc))
  #pval <- signif(lmod$coefficients[2,4],2)
  interc <- lmod$coefficients[1,1]
  slope <- lmod$coefficients[2,1]
  
  plot(meanPMD.mcrpc,pcna,
       xlab = "Mean PMD methylation",
       ylab = "% CNA",
       main = "% CN altered and PMD methylation level")
  #abline(a = interc,b = slope,col=2,lty=2)
  text(0.55,0.45,labels = paste0("P-value = ",pval,"\n Spearman's Rho = ",s.cor))
}

#### Compare PMD levels with different phenotypes, e.g. COMP, mutations etc 
fig.PMDperCMP %<a-%{
  group <- (names(percPMD.mcrpc) %in% samples_cmp)
  boxplot(percPMD.mcrpc~group,
          names=c("No CMP","CMP"),
          ylab = "Prop. of genome with PMD",
          xlab = "CMP status",
          main = "Prop. genome with PMD per CMP status",
          col = col_noCMP_CMP)
}

fig.PMDMeanperCMP %<a-%{
  group <- (names(meanPMD.mcrpc) %in% samples_cmp)
  pval <- signif(wilcox.test(meanPMD.mcrpc~group)$p.value,2)
  boxplot(meanPMD.mcrpc~group,
          names=c("No CMP","CMP"),
          ylab = "Mean PMD methylation",
          xlab = "CMP status",
          main = "PMD methylation level per CMP status",
          col = col_noCMP_CMP)
  text(2,0.53,labels = paste0("P-value = ", pval))
}

# Gene expression analysis ------------------------------------------------
# Correlation of PMD methylation levels and gene expression
# Perform GSEA/pathway enrichment analysis for top genes (proliferation?) externally

sub <- intersect(names(meanPMD),colnames(tpm)) #subset samples with both PMD data and gex data, #99 
meanPMD.sub <- meanPMD[sub]
gex <- tpm[,sub]
all(rownames(gex)==rownames(ensembl2sym))
stopifnot(all(names(meanPMD.sub)==colnames(gex)))

gex <- gex[rownames(gex) %in% rownames(ensembl2sym)[ensembl2sym$type %in% gencode_protein],] #restrict to coding genes

genes.cor <- data.frame(matrix(NA,nrow=nrow(gex),ncol=3))
colnames(genes.cor) <- c("gene","cor","pval")
  
for(i in 1:nrow(gex)){
  if(i %% 1000==0){
    print(i) 
  }
  ensembl_transcript <- rownames(gex)[i]
  gene <- ensembl2sym[ensembl_transcript,"name"]
  pval <- cor.test(meanPMD.sub, as.numeric(gex[i,]), method = "spearman", use="complete.obs")$p.value
  cor <- cor(meanPMD.sub, as.numeric(gex[i,]), method = "spearman",use="complete.obs")
  rownames(genes.cor)[i] <- ensembl_transcript
  genes.cor[i,"gene"] <- gene
  genes.cor[i,"cor"] <- cor
  genes.cor[i,"pval"] <- pval
}  

genes.cor[,2] <- as.numeric(genes.cor[,2])
genes.cor[,3] <- as.numeric(genes.cor[,3])  

#How to handel NAs? Add 0 and 1 to allow the list to be run with GSEA etc. 
genes.cor[,2][is.na(genes.cor[,2])] <- 0
genes.cor[,3][is.na(genes.cor[,3])] <- 1

#export table for external analyses. 
# write.table(x = genes.cor[order(genes.cor$cor,decreasing = F),1:2], file = "~/Shared_VM_external/WGBS/pmd/output/cor_gex_PMDmeth.txt",sep = "\t",row.names = F, quote = F)

# Properties of PMDs ------------------------------------------------------
# Check if mutations are more common in PMDs vs not PMDs
# Overlap mutations with GRanges PMD object, summarize TMB in PMD/noPMD for each sample  

mutations <- read.table(file = fn_muts ,header=T)  
mutations$start <- mutations$end <- mutations$pos
muts <- makeGRangesFromDataFrame(mutations,keep.extra.columns = T)
muts.all <- split(muts,muts$sample_id)

#calc number of mutations in pmd vs non pmd, and convert to TMB
calcTMBenr <- function(sampleName, mutData = muts.all, PMDdata = all.samples){
  
  pmd.gr <- all.samples[[sampleName]]
  muts.gr <- muts.all[[sampleName]]
  
  yesPMD.gr <- pmd.gr[pmd.gr$type=="PMD"]
  noPMD.gr <- pmd.gr[pmd.gr$type=="notPMD"]
  
  sizePMD <- sum(width(yesPMD.gr))
  sizeNoPMD <- sum(width(noPMD.gr))
  
  overlapsPMD <- sum(overlapsAny(muts.gr,yesPMD.gr))
  overlapsNoPMD <- sum(overlapsAny(muts.gr,noPMD.gr))
  
  tmbPMD <- (overlapsPMD/sizePMD)*1000000  
  tmbNoPMD <- (overlapsNoPMD/sizeNoPMD)*1000000  
  
  return(c(tmbNoPMD,tmbPMD))
}  

TMBenr <- data.frame(matrix(data = NA,nrow = 1,ncol = 2))
colnames(TMBenr) <- c("tmbNoPMD","tmbPMD")

for(i in 1:length(sample_ids_wgbs)){
  sampleName <- sample_ids_wgbs[i]
  enr <- calcTMBenr(sampleName)
  TMBenr[i,] <- enr
  rownames(TMBenr)[i] <- sampleName
}

fig.TMBinPMD %<a-%{
  TMBenr.nooutliers <- TMBenr[! rownames(TMBenr) %in% samples_hypermut,] #remove hypermutated
  
  pval <- signif(wilcox.test(TMBenr.nooutliers[,1],TMBenr.nooutliers[,2])$p.value,2)
  
  boxplot((TMBenr.nooutliers[,1]),(TMBenr.nooutliers[,2]),
          names=c("Outside PMD","Inside PMD"),
          xlab="Region",
          ylab="Mutations per Mb",
          col = col_no_pmd_pmd,
          main = "TMB in and outside of PMDs")
  text(1.5,11,paste0("P-value = ",pval))
}

fig.enrTMB %<a-%{
  TMBenr.nooutliers <- TMBenr[! rownames(TMBenr) %in% samples_hypermut,] #remove hypermutated
  
  boxplot(TMBenr.nooutliers[,2]/TMBenr.nooutliers[,1],
          main = "Enrichent of mutations inside PMDs",
          ylab = "TMB in PMD / TMB outside PMD",
          col = col_no_pmd_pmd[2])
}

#### Check if important regions fall into PMDs (e.g. gene density)
coding.density <- data.frame(matrix(NA, ncol = 4, nrow = length(sample_ids_wgbs)))
colnames(coding.density) <- c("sample","frac_PMD","frac_noPMD","ratio")

#restrict this analysis to mCRPC samples
for(i in 1:length(sample_ids_wgbs)){
  sampleName <- sample_ids_wgbs[i]
  sample.gr <- all.samples[[sampleName]]
  
  length_PMD <- sum(width(sample.gr[sample.gr$type=="PMD"]))
  length_notPMD <- sum(width(sample.gr[sample.gr$type=="notPMD"]))
  
  coding_in_PMD <- sum(width(intersect(coding.genome,sample.gr[sample.gr$type=="PMD"]))) 
  coding_not_in_PMD <- sum(width(intersect(coding.genome,sample.gr[sample.gr$type=="notPMD"])))
  
  frac_PMD <- coding_in_PMD/length_PMD
  frac_noPMD <- coding_not_in_PMD/length_notPMD
  
  ratio <- frac_noPMD/frac_PMD
  
  coding.density[i,] <- c(sampleName,frac_PMD,frac_noPMD,ratio)
}

coding.density[,2] <- as.numeric(coding.density[,2])
coding.density[,3] <- as.numeric(coding.density[,3])
coding.density[,4] <- as.numeric(coding.density[,4])

fig.fractionCoding %<a-%{
  pval <- signif(wilcox.test(coding.density[,2],coding.density[,3])$p.value,2)
  boxplot(coding.density[,2],coding.density[,3],
          names=c("In PMD","Outside PMD"),
          ylab = "Fraction coding genome",
          col = col_no_pmd_pmd,
          main = "Fraction coding genome in PMDs")
  text(1,0.07,labels = paste0("P-value = ",pval))
}

fig.fractionCodingEnr %<a-%{
  boxplot(coding.density[,4],
          names="",
          ylab = "Fraction coding genome outside/inside PMDs",
          col = col_no_pmd_pmd[[2]],
          main = "Enrichment of coding genome in PMDs")
}

#superset PMDs, this makes up ~91% of the genome and is not very informative
#could potentially replicate breast analysis with number of samples with PMD in a region of the genome
all.gr <- GRangesList(unlist(all.samples))
all.pmd.gr <- GRangesList(lapply(all.gr, FUN = function(x){x[x$type=="PMD"]}))
supersetPMD <- reduce(c(unlist(all.pmd.gr)))

#### Plot gene expression of genes in and outside of PMDs
ratioGex <- vector() #save ratios for plotting enrichment
setwd(OUTDIR_GEX_IN_PMD)

new.gex <- gex + 0.01 #avoid having zero

for(i in 1:length(sub)){
  
  sampleName <- sub[i]
  
  genes.gr.sample <- genes.gr.protein
  genes.gr.sample$geneInPMD <- overlapsAny(genes.gr.sample, all.samples[[sampleName]][all.samples[[sampleName]]$type=="PMD"])

  stopifnot(all(names(genes.gr.sample)==rownames(new.gex)))
  
  ratio <- mean(new.gex[genes.gr.sample$geneInPMD==FALSE,sampleName])/mean(new.gex[genes.gr.sample$geneInPMD==TRUE,sampleName])
  ratioGex[i] <- ratio
  names(ratioGex)[i] <- sampleName
  
  # png(paste0(sampleName,".png"))
  # boxplot(log2(new.gex[,sampleName])~genes.gr.sample$geneInPMD,
  #        main = sampleName,
  #        ylab = "Gene expression (log2(TPM))",
  #        xlab = "Overlap of gene with PMD",
  #        names = c("Outside PMD", "Inside PMD"),
  #        col = col_no_pmd_pmd)
  # dev.off()
}

fig.fractionGexIncrease %<a-%{
  boxplot(ratioGex,
          ylab = "Mean gex outside PMD / Mean gex inside PMD",
          main = "Increase in gene expression outide of PMDs",
          col = col_no_pmd_pmd[2])
}

sampleName <- sampleNames[52] #random example
fig.fractionGexIncreaseExample %<a-%{
  pval <- wilcox.test(log2(new.gex[,sampleName])~genes.gr.sample$geneInPMD)$p.value
  boxplot(log2(new.gex[,sampleName])~genes.gr.sample$geneInPMD,
          main = sampleName,
          ylab = "Gene expression (log2(TPM))",
          xlab = "Overlap of gene with PMD",
          names = c("Outside PMD", "Inside PMD"),
          col = col_no_pmd_pmd)
  text(1.5,8,labels = paste0("P-value = ",pval))
}

# Final figures -----------------------------------------------------------
setwd(OUTDIR)
png(filename = "PMD_frac_and_meth.png", width = 12, height = 7, units = "in", res = 300)
par(mfrow=c(1,2))
fig.perPMD
fig.meanPMDMeth
dev.off()

png(filename = "PMD_tmb_cna_tcc.png", width = 12, height = 5, units = "in", res = 300)
par(mfrow=c(1,3))
fig.PMDMeanvsTCcont
fig.PMDMeanvsMuts
fig.PMDMeanvsCNA
dev.off()

png(filename = "PMD_tmb_enrich.png", width = 10, height = 7, units = "in", res = 300)
par(mfrow=c(1,2))
fig.TMBinPMD
fig.enrTMB
dev.off()

png(filename = "PMD_gene_density.png", width = 10, height = 7, units = "in", res = 300)
par(mfrow=c(1,2))
fig.fractionCoding
fig.fractionCodingEnr
dev.off()

png(filename = "PMD_gex_in_pmd.png", width = 10, height = 7, units = "in", res = 300)
par(mfrow=c(1,2))
fig.fractionGexIncreaseExample
fig.fractionGexIncrease
dev.off()

png(filename = "PMD_in_CMP.png", width = 5, height = 7, units = "in", res = 300)
fig.PMDMeanperCMP
dev.off()
