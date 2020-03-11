library(ggplot2)
library(grid)
library(wesanderson)
source('wgbs/helper.R')

#Plots figure 2

fileout <- 'C:/Users/Admin/Desktop/genes.pdf'
filepath <- 'E:/wgbs/gene_output/'

#1=NT,2=prostate
ctrlgroup <- 2 #use benign prostate as control group
if(ctrlgroup==1) {
	filepathctrl <- 'E:/wgbs/gene_output_NT/'
	samples_ctrl <- samples_normal
} else {
	filepathctrl <- 'E:/wgbs/gene_output_upitt/'
	samples_ctrl <- upitt_benign_prostate
	#samples_ctrl <- upitt_localized_tumor
}

tosort <- data.frame(samples_wgbs)
tosort$ERG <- !allele_effect('ERG')$alleles[samples_wgbs,'activating_sv']
tosort$ETV1 <- !allele_effect('ETV1')$alleles[samples_wgbs,'activating_sv']
tosort$ETV4 <- !allele_effect('ETV4')$alleles[samples_wgbs,'activating_sv']
tosort$ETV5 <- !allele_effect('ETV5')$alleles[samples_wgbs,'activating_sv']
tosort$CHD1 <- !allele_effect('CHD1')$alleles[samples_wgbs,'n_alleles_inactivated']==2
tosort$SPOP <- !allele_effect('SPOP')$alleles[samples_wgbs,'inactivating_missense']
tosort$CDK12 <- !allele_effect('CDK12')$alleles[samples_wgbs,'n_alleles_inactivated']==2
tosort$RB1 <- !allele_effect('RB1')$alleles[samples_wgbs,'n_alleles_inactivated']==2
tosort$PTEN <- !allele_effect('PTEN')$alleles[samples_wgbs,'n_alleles_inactivated']==2
tosort$MYC <- !allele_effect('MYC')$alleles[samples_wgbs,'CNA_amp']
tosort$TP53 <- !allele_effect('TP53')$alleles[samples_wgbs,'n_alleles_inactivated']==2

print(wilcox.test(hmrlengths[samples_wgbs]~tosort$TP53))
print(wilcox.test(hmrlengths[samples_wgbs]~tosort$MYC))

map <- c('Mutation','SV','CN Gain','Mutation','Mutation','Mutation','SV','CN Loss','CN Loss')
names(map) <- c(gof,lof)
site <- c('Bone','Lymph_node','Liver','Other')

height <- 0.9

groups <- c(5,15,25)
numgroups <- length(groups)*2
groupnames <- c(rev(paste('-',groups,'%',sep='')),paste('+',groups,'%',sep=''))

pvals <- c()

#What genes to use, and if DNA alterations are gain or loss of function
genes2use <- rbind.data.frame(
	cbind.data.frame(gene='ERG',gof=T),
	cbind.data.frame(gene='ETV1',gof=T),
	cbind.data.frame(gene='ETV4',gof=T),
	#cbind.data.frame(gene='ETV5',gof=T),
	cbind.data.frame(gene='AR',gof=T),
	cbind.data.frame(gene='KLK3',gof=T),
	cbind.data.frame(gene='NKX3-1',gof=F),
	cbind.data.frame(gene='FOLH1',gof=T),
	#cbind.data.frame(gene='FOXA1',gof=T),
	#cbind.data.frame(gene='NCOR1',gof=F),
	#cbind.data.frame(gene='NCOR2',gof=F),
	#cbind.data.frame(gene='ASXL2',gof=F),
	#cbind.data.frame(gene='ZBTB16',gof=F),
	cbind.data.frame(gene='SCHLAP1',gof=T),
	#cbind.data.frame(gene='PCAT14',gof=T),
	#cbind.data.frame(gene='PCAT1',gof=T),
	#cbind.data.frame(gene='MALAT1',gof=T),
	#cbind.data.frame(gene='NEAT1',gof=T),
	#cbind.data.frame(gene='PCGEM1',gof=T), #PCAT9
	#cbind.data.frame(gene='PRNCR1',gof=T),
	#cbind.data.frame(gene='CTBP1-AS',gof=T), #PCAT10
	#cbind.data.frame(gene='HOTAIR',gof=T),
	#cbind.data.frame(gene='H19',gof=T),
	#cbind.data.frame(gene='GAS5',gof=T),
	cbind.data.frame(gene='PIK3CA',gof=T),
	cbind.data.frame(gene='MYC',gof=T),
	#cbind.data.frame(gene='AXL',gof=T),
	cbind.data.frame(gene='SPOP',gof=F),
	cbind.data.frame(gene='CHD1',gof=F),
	cbind.data.frame(gene='RB1',gof=F),
	#cbind.data.frame(gene='CDKN2A',gof=F),
	#cbind.data.frame(gene='APC',gof=F),
	#cbind.data.frame(gene='CCND1',gof=T),
	#cbind.data.frame(gene='BRAF',gof=T),
	#cbind.data.frame(gene='HRAS',gof=T),
	cbind.data.frame(gene='TP53',gof=F),
	cbind.data.frame(gene='BRCA1',gof=F),
	cbind.data.frame(gene='BRCA2',gof=F),
	#cbind.data.frame(gene='ZFHX3',gof=F),
	#cbind.data.frame(gene='GSTP1',gof=F)
	#cbind.data.frame(gene='ATM',gof=F),
	#cbind.data.frame(gene='CDK12',gof=F),
	#cbind.data.frame(gene='ERCC2',gof=F),
	#cbind.data.frame(gene='PRKDC',gof=F),
	#cbind.data.frame(gene='MLH1',gof=F),
	#cbind.data.frame(gene='MSH2',gof=F),
	#cbind.data.frame(gene='MSH6',gof=F)
	#cbind.data.frame(gene='DNMT1',gof=T),
	#cbind.data.frame(gene='DNMT3A',gof=T),
	#cbind.data.frame(gene='DNMT3B',gof=T),
	#cbind.data.frame(gene='TET1',gof=F),
	#cbind.data.frame(gene='TET2',gof=F),
	#cbind.data.frame(gene='TET3',gof=F),
	#cbind.data.frame(gene='CTLA4',gof=F),
	#cbind.data.frame(gene='CD274',gof=F),
	#cbind.data.frame(gene='PDCD1LG2',gof=F),
	#cbind.data.frame(gene='IDH1',gof=F),
	#cbind.data.frame(gene='PTEN',gof=F),
	#cbind.data.frame(gene='AKT1',gof=T),
	#cbind.data.frame(gene='CTNNB1',gof=T),
	#cbind.data.frame(gene='ZNRF3',gof=F),
	#cbind.data.frame(gene='KMT2C',gof=F),
	#cbind.data.frame(gene='KMT2D',gof=F),
	#cbind.data.frame(gene='KDM6A',gof=F),
	#cbind.data.frame(gene='HDAC4',gof=F),
	#cbind.data.frame(gene='MED12',gof=T),
	#cbind.data.frame(gene='GNAS',gof=F),
	cbind.data.frame(gene='GSTP1',gof=F),
	cbind.data.frame(gene='WT1',gof=F),
	cbind.data.frame(gene='GRASP',gof=F),
	#cbind.data.frame(gene='RHCG',gof=F),#no EMR
	cbind.data.frame(gene='PHOX2A',gof=F),
	cbind.data.frame(gene='POU3F3',gof=F),
	cbind.data.frame(gene='VAX1',gof=F),
	cbind.data.frame(gene='TMEM106A',gof=F),
	cbind.data.frame(gene='CYP27A1',gof=F),
	cbind.data.frame(gene='SIX6',gof=F),
	cbind.data.frame(gene='IRX1',gof=F),
	#cbind.data.frame(gene='CCDC8',gof=F),#no EMR
	cbind.data.frame(gene='LBX1',gof=F),
	cbind.data.frame(gene='HOXA3',gof=F),
	cbind.data.frame(gene='NKX2-5',gof=F),
	cbind.data.frame(gene='NKX2-1',gof=F),
	cbind.data.frame(gene='BARHL2',gof=F),
	cbind.data.frame(gene='HOXD8',gof=F),
	cbind.data.frame(gene='FOXE3',gof=F),
	cbind.data.frame(gene='LHX2',gof=F),
	cbind.data.frame(gene='KIT',gof=F),
	cbind.data.frame(gene='HHEX',gof=F),
	cbind.data.frame(gene='EN2',gof=F),
	cbind.data.frame(gene='FBLN1',gof=F),
	#cbind.data.frame(gene='HOXD4',gof=F),#no EMR
	cbind.data.frame(gene='NKX2-2',gof=F),
	#cbind.data.frame(gene='GP5',gof=F),#no EMR
	cbind.data.frame(gene='WNT2',gof=F),
	cbind.data.frame(gene='BDNF',gof=F),
	#cbind.data.frame(gene='TBX3',gof=F),#hypo-methylated
	#cbind.data.frame(gene='DLX1',gof=F),#hypo-methylated
	#cbind.data.frame(gene='ESR1',gof=F),#hypo-methylated
	cbind.data.frame(gene='EFS',gof=F),
	cbind.data.frame(gene='HOXC11',gof=F),
	cbind.data.frame(gene='LAMP5',gof=F),
	cbind.data.frame(gene='TBX15',gof=F),
	cbind.data.frame(gene='LHX9',gof=F),
	#cbind.data.frame(gene='RUNX3',gof=F),#no EMR
	#cbind.data.frame(gene='RARB',gof=F),#no EMR
	cbind.data.frame(gene='AOX1',gof=F),
	cbind.data.frame(gene='ZNF154',gof=F),
	cbind.data.frame(gene='APC',gof=F),
	cbind.data.frame(gene='CDKN2A',gof=F),
	cbind.data.frame(gene='HOXD3',gof=F),
	cbind.data.frame(gene='PTGS2',gof=F),
	cbind.data.frame(gene='WT1',gof=F)
)

genes2use$gene <- paste(genes2use$gene)

filename <- '../methylseq/hmrseg/summary_expr.txt'
sumexpr <- read.table(filename,sep='\t',stringsAsFactors=F,header=T,row.names=1)
sumexpr$anovafdr <- p.adjust(sumexpr$anovap,method='fdr')
rowsfdr <- sumexpr$anovafdr<=0.05
rowsfdr[is.na(rowsfdr)] <- F
genesfdr <- sumexpr[rowsfdr,'name']

numgenes <- dim(genes2use)[1]

collvls <- paste(tosort[order(tosort$ERG,tosort$ETV1,tosort$ETV4,tosort$ETV5,tosort$CHD1,tosort$SPOP,
	tosort$CDK12,tosort$RB1,tosort$PTEN,tosort$MYC),'samples_wgbs'])

df <- data.frame(
	x=samples_wgbs,
	y=rep(-1,length(samples_wgbs)),
	ht=height,
	fill=rep(NA,length(samples_wgbs)),
	stringsAsFactors=F)
rownames(df) <- samples_wgbs
df[intersect(names(locations)[locations=='Bone'],samples_wgbs),'fill'] <- site[1]
df[intersect(names(locations)[locations=='Lymph_node'],samples_wgbs),'fill'] <- site[2]
df[intersect(names(locations)[locations=='Liver'],samples_wgbs),'fill'] <- site[3]
df[intersect(names(locations)[locations=='Other'],samples_wgbs),'fill'] <- site[4]

#Loop through genes and set up plotting data frames
for(i in 1:numgenes) {
	genei <- paste(genes2use[i,'gene'])
	print(genei)
	
	gofi <- genes2use[i,'gof']
	all <- allele_effect(genei)$alleles
	
	df1 <- data.frame(
		x=samples_wgbs,
		y=rep((i*4)-3,length(samples_wgbs)),
		ht=height,
		fill=rep(NA,length(samples_wgbs)),
		stringsAsFactors=F)
	df2 <- data.frame(
		x=samples_wgbs,
		y=rep((i*4)-2,length(samples_wgbs)),
		ht=height,
		fill=rep(NA,length(samples_wgbs)),
		stringsAsFactors=F)
	df3 <- data.frame(
		x=samples_wgbs,
		y=rep((i*4)-1,length(samples_wgbs)),
		ht=height,
		fill=rep(NA,length(samples_wgbs)),
		stringsAsFactors=F)
	rownames(df1) <- samples_wgbs
	rownames(df2) <- samples_wgbs
	rownames(df3) <- samples_wgbs
	
	rowensembl <- ensembl2sym$name==genei
	ensemblid <- rownames(ensembl2sym)[rowensembl]
	stopifnot(sum(rowensembl)==1)
	chr <- ensembl2sym[rowensembl,'chr']
	gene_start <- ensembl2sym[rowensembl,'start']
	gene_end <- ensembl2sym[rowensembl,'end']
	
	#Just plot the promoter
	if(ensembl2sym[ensemblid,'strand']=='+') {
		start <- gene_start-1500
		end <- gene_start+1500
		#end <- gene_end
	} else {
		#start <- gene_start
		start <- gene_end-1500
		end <- gene_end+1500
	}
	gr <- GRanges(chr,IRanges(start,end))
	
	#get recurrent HMRs
	emr_gene <- tracks[['HMRseggene']][tracks[['HMRseggene']]$gene==genei]
	emr_gene <- emr_gene[countOverlaps(emr_gene,gr)>0 & !is.na(emr_gene$avg_cor)]
	print(emr_gene)
	if(length(emr_gene)>0) {
		#Get the best eHMR
		#print(emr_gene)
		bestcor <- min(emr_gene$avg_p)
		maxrow <- emr_gene$avg_p==bestcor
		emr_best <- emr_gene[maxrow]
		stopifnot(length(emr_best)==1)
		
		#Get the methylation in the eHMR
		emr_start <- start(emr_best)
		emr_end <- end(emr_best)
		methylmean <- meanmethyl(ensemblid,filepath,emr_start,emr_end)[samples_wgbs]
		methylmean_nt <- meanmethyl(ensemblid,filepathctrl,emr_start,emr_end)
		names(methylmean_nt) <- gsub('_R1_BISMARK_BT2_PE','',names(methylmean_nt),fixed=T)
		methylmean_nt <- methylmean_nt[samples_ctrl]
		
		#Compare against the control
		if(ctrlgroup==1) {
			methyldiff <- rep(NA,length(methylmean))
			for(loci in unique(locations)) {
				samplestosub <- names(methylmean) %in% names(locations)[locations==loci]
				if(loci %in% locs_normal) {
					tosub <- mean(methylmean_nt[names(locs_normal)[locs_normal==loci]],na.rm=T)
				} else {
					tosub <- mean(methylmean_nt,na.rm=T)
				}
				methyldiff[samplestosub] <- methylmean[samplestosub] - tosub
			}
		} else {
			methyldiff <- methylmean-(mean(methylmean_nt,na.rm=T))
		}
		methylmeanadeno <- methylmean[samples_wgbs]
		print(mean(methylmeanadeno))
		print(mean(methylmean_nt))
		print(wilcox.test(methylmeanadeno,methylmean_nt))
		pvals[genei] <- wilcox.test(methylmeanadeno,methylmean_nt)$p.value
		#if(genei == 'ERG') {
		#	sampleserg <- samples_wgbs[allele_effect('ERG')$alleles[samples_wgbs,'activating_sv']]
		#	methylmeanerg <- methylmean[sampleserg]
		#	print(wilcox.test(methylmeanerg,methylmean_nt))
		#}
		
		for(cutoff in groups) {
			rowspos <- methyldiff>=cutoff
			rowspos[is.na(rowspos)] <- F
			df3[rowspos,'fill'] <- paste('+',cutoff,'%',sep='')
			rowsneg <- methyldiff<=(-cutoff)
			rowsneg[is.na(rowsneg)] <- F
			df3[rowsneg,'fill'] <- paste('-',cutoff,'%',sep='')
		}
	} else {
		print(paste('No EMR for',genei))
	}
	
	#Add new AR enhancer amplification
	if(genei=='AR') { #Add all eHMR enhancers too
		chrx_1 <- cnregion('chrX',66744900,66745900)[samples_wgbs]
		df1[names(chrx_1)[chrx_1 >=GAIN_SEX],'fill'] <- 'CNA_amp'
		chrx_2 <- cnregion('chrX',66900400,66902700)[samples_wgbs]
		df1[names(chrx_2)[chrx_2 >=GAIN_SEX],'fill'] <- 'CNA_amp'
		chrx_3 <- cnregion('chrX',66906200,66908000)[samples_wgbs]
		df1[names(chrx_3)[chrx_3 >=GAIN_SEX],'fill'] <- 'CNA_amp'
		chrx_4 <- cnregion('chrX',67043300,67045400)[samples_wgbs]
		df1[names(chrx_4)[chrx_4 >=GAIN_SEX],'fill'] <- 'CNA_amp'
		chrx_5 <- cnregion('chrX',67361700,67362400)[samples_wgbs]
		df1[names(chrx_5)[chrx_5 >=GAIN_SEX],'fill'] <- 'CNA_amp'
		chrx_6 <- cnregion('chrX',67543100,67552000)[samples_wgbs]
		df1[names(chrx_6)[chrx_6 >=GAIN_SEX],'fill'] <- 'CNA_amp'
		chrx_7 <- cnregion('chrX',67799500,67800400)[samples_wgbs]
		df1[names(chrx_7)[chrx_7 >=GAIN_SEX],'fill'] <- 'CNA_amp'
	}
	
	#Add DNA alterations
	if(gofi) {
		for(j in 1:length(gof)) {
			namegof <- gof[j]
			samplesalt <- samples_wgbs[all[samples_wgbs,namegof]]
			if(genei=='AR' && namegof=='CNA_amp') {
				df1[samplealt,'fill'] <- namegof
			} else {
				for(samplealt in samplesalt) {
					count <- is.na(df1[samplealt,'fill']) + is.na(df2[samplealt,'fill'])
					if(count==2) {
						df1[samplealt,'fill'] <- namegof
					} else if(count==1) {
						df2[samplealt,'fill'] <- namegof
					}
				}
			}
		}
	} else {
		for(j in 1:length(lof)) {
			namelof <- lof[j]
			samplesalt <- samples_wgbs[all[samples_wgbs,namelof]]
			for(samplealt in samplesalt) {
				if(namelof=='CNA_2' || (namelof=='inactivating_missense' && genei=='CDK12')) {
					df1[samplealt,'fill'] <- namelof
					df2[samplealt,'fill'] <- namelof
				} else {
					count <- is.na(df1[samplealt,'fill']) + is.na(df2[samplealt,'fill'])
					if(count==2) {
						df1[samplealt,'fill'] <- namelof
					} else if(count==1) {
						df2[samplealt,'fill'] <- namelof
					} else {
						print(paste(genei,samplealt,namelof))
					}
				}
			}
		}
	}
	df <- rbind.data.frame(df,df1,df2,df3)
}
df$x <- factor(df$x,levels=collvls)
rowsna <- is.na(df$fill)
df <- df[!rowsna,]

rowsalt <- df$fill %in% c(gof,lof)
df[rowsalt,'fill'] <- map[df[rowsalt,'fill']]
df$fill <- factor(df$fill,levels=c(site,' ',unique(map),'  ',groupnames))

reds <- colorRampPalette(c('white','red'))(n=10)
blues <- colorRampPalette(c('white','blue'))(n=10)
colpal <- c(blues[10],blues[7],blues[5],reds[5],reds[7],reds[10])

xpos <- -10.5
xcex <- 0.7
p <- ggplot(df,aes(x,y,width=0.5,height=ht))+geom_tile(aes(fill=fill))+
	scale_fill_manual(drop=F,name="",
		values=c(
		'Bone'=wes_palette('Chevalier1')[1], #darkgreen
		'Lymph_node'=wes_palette('Chevalier1')[2], #gold
		'Liver'=wes_palette('Chevalier1')[4], #brown
		'Other'=wes_palette('Chevalier1')[3], #gray
		" "="white",
		"Mutation"="purple",
		"SV"="orange",
		"CN Gain"="green",
		"CN Loss"="cyan4",
		"  "="white",
		"-25%"=colpal[1],
		"-15%"=colpal[2],
		"-5%"=colpal[3],
		"+5%"=colpal[4],
		"+15%"=colpal[5],
		"+25%"=colpal[6]))+
	guides(fill=F)+
	theme_classic()+xlab('Samples')+
	theme(plot.margin=unit(c(0.2,3.3,0,0.2),"cm"))+
	theme(axis.text.x=element_text(angle=90, hjust=1,size=5))+
	scale_y_reverse(name='',minor_breaks=(0:numgenes)*4,
		breaks=(c(1,(1:numgenes)*4)-2),labels=c('Biopsy Site',paste(genes2use$gene)))+
	coord_cartesian(ylim=c(205,2))+
	theme(panel.grid.minor=element_line(colour='black',size=0.5))+
	annotation_custom(grob=textGrob(expression(underline('  ETS  ')),
		hjust='center',rot=90,gp=gpar(cex=xcex)),
		ymin=-6,ymax=-6,xmin=xpos,xmax=xpos)+
	annotation_custom(grob=textGrob(expression(underline('     AR     ')),
		hjust='center',rot=90,gp=gpar(cex=xcex)),
		ymin=-20,ymax=-20,xmin=xpos,xmax=xpos)+
	#annotation_custom(grob=textGrob(expression(underline('  lncRNA  ')),
	#	hjust='center',rot=90,gp=gpar(cex=xcex)),
	#	ymin=-36,ymax=-36,xmin=xpos,xmax=xpos)+
	annotation_custom(grob=textGrob(expression(underline('Oncogene')),
		hjust='center',rot=90,gp=gpar(cex=xcex)),
		ymin=-34,ymax=-34,xmin=xpos,xmax=xpos)+
	annotation_custom(grob=textGrob(expression(underline('Tumor suppressor')),
		hjust='center',rot=90,gp=gpar(cex=xcex)),
		ymin=-52,ymax=-52,xmin=xpos,xmax=xpos)+
	annotation_custom(grob=textGrob(expression(underline('                                                       Hyper-methylated genes in PCa from the literature                                                      ')),
		hjust='center',rot=90,gp=gpar(cex=xcex)),
		ymin=-138,ymax=-138,xmin=xpos,xmax=xpos)+
	annotation_custom(grob=linesGrob(),xmin=101,xmax=103,ymin=-61,ymax=-57)+
	annotation_custom(grob=linesGrob(),xmin=101,xmax=103,ymin=-62,ymax=-60)+
	annotation_custom(grob=linesGrob(),xmin=101,xmax=103,ymin=-63,ymax=-63)+
	annotation_custom(grob=textGrob('DNA alteration 1',
		hjust=0,gp=gpar(cex=0.7)),
		ymin=-57,ymax=-57,xmin=103.5,xmax=103.5)+
	annotation_custom(grob=textGrob('DNA alteration 2',
		hjust=0,gp=gpar(cex=0.7)),
		ymin=-60,ymax=-60,xmin=103.5,xmax=103.5)+
	annotation_custom(grob=textGrob('Methylation',
		hjust=0,gp=gpar(cex=0.7)),
		ymin=-63,ymax=-63,xmin=103.5,xmax=103.5)

#Add custom legends
legx <- 105
legtx <- 107
width <- 5
fontsize <- 0.8

ystart <- 3
p <- p +
	annotation_custom(grob=linesGrob(gp=gpar(col=wes_palette('Chevalier1')[1],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-5,ymax=ystart-9)+
	annotation_custom(grob=linesGrob(gp=gpar(col=wes_palette('Chevalier1')[2],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-10,ymax=ystart-14)+
	annotation_custom(grob=linesGrob(gp=gpar(col=wes_palette('Chevalier1')[4],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-15,ymax=ystart-19)+
	annotation_custom(grob=linesGrob(gp=gpar(col=wes_palette('Chevalier1')[3],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-20,ymax=ystart-24)+
	annotation_custom(grob=textGrob('Bone',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-7,ymax=ystart-7,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('LN',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-12,ymax=ystart-12,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('Liver',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-17,ymax=ystart-17,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('Other',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-22,ymax=ystart-22,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('Biopsy site',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-2,ymax=ystart-2,xmin=legx-2,xmax=legx-2)

ystart <- -28
p <- p +
	annotation_custom(grob=linesGrob(gp=gpar(col='purple',lwd=width)),xmin=legx,xmax=legx,ymin=ystart-5,ymax=ystart-9)+
	annotation_custom(grob=linesGrob(gp=gpar(col='orange',lwd=width)),xmin=legx,xmax=legx,ymin=ystart-10,ymax=ystart-14)+
	annotation_custom(grob=linesGrob(gp=gpar(col='green',lwd=width)),xmin=legx,xmax=legx,ymin=ystart-15,ymax=ystart-19)+
	annotation_custom(grob=linesGrob(gp=gpar(col='cyan4',lwd=width)),xmin=legx,xmax=legx,ymin=ystart-20,ymax=ystart-24)+
	annotation_custom(grob=textGrob('Mutation',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-7,ymax=ystart-7,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('SV',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-12,ymax=ystart-12,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('CN Gain',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-17,ymax=ystart-17,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('CN Loss',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-22,ymax=ystart-22,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('DNA alteration',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-2,ymax=ystart-2,xmin=legx-2,xmax=legx-2)

ystart <- -80
p <- p +
	annotation_custom(grob=linesGrob(gp=gpar(col=colpal[1],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-5,ymax=ystart-9)+
	annotation_custom(grob=linesGrob(gp=gpar(col=colpal[2],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-10,ymax=ystart-14)+
	annotation_custom(grob=linesGrob(gp=gpar(col=colpal[3],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-15,ymax=ystart-19)+
	annotation_custom(grob=linesGrob(gp=gpar(col=colpal[4],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-20,ymax=ystart-24)+
	annotation_custom(grob=linesGrob(gp=gpar(col=colpal[5],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-25,ymax=ystart-29)+
	annotation_custom(grob=linesGrob(gp=gpar(col=colpal[6],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-30,ymax=ystart-34)+
	annotation_custom(grob=textGrob(expression(""<="-25%"),hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-7,ymax=ystart-7,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob(expression(""<="-15%"),hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-12,ymax=ystart-12,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob(expression(""<="-5%"),hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-17,ymax=ystart-17,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob(expression("">="+5%"),hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-22,ymax=ystart-22,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob(expression("">="+15%"),hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-27,ymax=ystart-27,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob(expression("">="+25%"),hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-32,ymax=ystart-32,xmin=legtx,xmax=legtx)+
	annotation_custom(grob=textGrob('Differential',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart+6,ymax=ystart+6,xmin=legx-2,xmax=legx-2)+
	annotation_custom(grob=textGrob('Promoter',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart+2,ymax=ystart+2,xmin=legx-2,xmax=legx-2)+
	annotation_custom(grob=textGrob('Methylation',hjust=0,gp=gpar(cex=fontsize)),ymin=ystart-2,ymax=ystart-2,xmin=legx-2,xmax=legx-2)

g = ggplotGrob(p)
g$layout$clip[g$layout$name=="panel"] <- "off"
g$layout$clip[g$layout$name=="legend"] <- "off"
grid.draw(g)

ggsave(fileout,g,width=10,height=10)
dev.off()

print(p.adjust(pvals,method='fdr'))
