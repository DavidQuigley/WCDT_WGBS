#source('/Volumes/datasets_1/human_sequence_prostate_WGBS/reproduce/scripts/2019_05_15_prepare_environment.r')

#Plots figure 2

#1=NT,2=prostate
ctrlgroup = 2 #use benign prostate as control group
if(ctrlgroup==1){
    filepathctrl = dir_gene_output_NT 
    samples_ctrl = samples_normal
}else{
    filepathctrl = dir_gene_output_upitt
    samples_ctrl = samples_upitt_benign_prostate
}

tosort = data.frame(sample_ids_wgbs)
tosort$ERG = !allele_effect('ERG')$alleles[sample_ids_wgbs,'activating_sv']
tosort$ETV1 = !allele_effect('ETV1')$alleles[sample_ids_wgbs,'activating_sv']
tosort$ETV4 = !allele_effect('ETV4')$alleles[sample_ids_wgbs,'activating_sv']
tosort$ETV5 = !allele_effect('ETV5')$alleles[sample_ids_wgbs,'activating_sv']
tosort$CHD1 = !allele_effect('CHD1')$alleles[sample_ids_wgbs,'n_alleles_inactivated']==2
tosort$SPOP = !allele_effect('SPOP')$alleles[sample_ids_wgbs,'inactivating_missense']
#tosort$CDK12 = !allele_effect('CDK12')$alleles[sample_ids_wgbs,'n_alleles_inactivated']==2
tosort$RB1 = !allele_effect('RB1')$alleles[sample_ids_wgbs,'n_alleles_inactivated']==2
tosort$PTEN = !allele_effect('PTEN')$alleles[sample_ids_wgbs,'n_alleles_inactivated']==2
tosort$MYC = !allele_effect('MYC')$alleles[sample_ids_wgbs,'CNA_amp']
tosort$TP53 = !allele_effect('TP53')$alleles[sample_ids_wgbs,'n_alleles_inactivated']==2

print(wilcox.test(hmrlengths[sample_ids_wgbs]~tosort$TP53))
print(wilcox.test(hmrlengths[sample_ids_wgbs]~tosort$MYC))

map = c('Mutation','SV','CN Gain','Mutation','Mutation','Mutation','SV','CN Loss','CN Loss')
gof = c(
    'activating_missense',
    'activating_sv',
    'CNA_amp')
lof = c(
    'inactivating_missense',
    'nonsense',
    'inactivating_germline',
    'inactivating_sv',
    'CNA_1',
    'CNA_2')
names(map) = c(gof,lof)
site = c('Bone','Lymph_node','Liver','Other')

height = 0.9

groups = c(10,20,30)
numgroups = length(groups)*2
groupnames = c(rev(paste('-',groups,'%',sep='')),paste('+',groups,'%',sep=''))

#What genes to use, and if DNA alterations are gain or loss of function
genes2use = rbind.data.frame(
    cbind.data.frame(gene='ERG',gof=TRUE),
    cbind.data.frame(gene='ETV1',gof=TRUE),
    cbind.data.frame(gene='ETV4',gof=TRUE),
    cbind.data.frame(gene='AR',gof=TRUE),
    cbind.data.frame(gene='KLK3',gof=TRUE),
    cbind.data.frame(gene='NKX3-1',gof=FALSE),
    cbind.data.frame(gene='FOLH1',gof=TRUE),
    cbind.data.frame(gene='SCHLAP1',gof=TRUE),
    #cbind.data.frame(gene='PCAT14',gof=TRUE),
    #cbind.data.frame(gene='PCGEM1',gof=TRUE), #PCAT9
    #cbind.data.frame(gene='CTBP1-AS',gof=TRUE), #PCAT10
    cbind.data.frame(gene='PIK3CA',gof=TRUE),
    cbind.data.frame(gene='MYC',gof=TRUE),
    #cbind.data.frame(gene='AXL',gof=TRUE),
    cbind.data.frame(gene='SPOP',gof=FALSE),
    cbind.data.frame(gene='CHD1',gof=FALSE),
    cbind.data.frame(gene='RB1',gof=FALSE),
    cbind.data.frame(gene='TP53',gof=FALSE),
    cbind.data.frame(gene='BRCA2',gof=FALSE),
    cbind.data.frame(gene='BRCA1',gof=FALSE)
    #cbind.data.frame(gene='ZFHX3',gof=FALSE),
    #cbind.data.frame(gene='GSTP1',gof=FALSE)
)

genes2use$gene = paste(genes2use$gene)


sumexpr = read.table(fn_hmr_summary_expr, sep='\t',stringsAsFactors=F,header=T,row.names=1)
sumexpr$anovafdr = p.adjust(sumexpr$anovap, method='fdr')
rowsfdr = sumexpr$anovafdr<=0.05
rowsfdr[is.na(rowsfdr)] = F
genesfdr = sumexpr[rowsfdr,'name']

numgenes = dim(genes2use)[1]

collvls = paste(tosort[order(tosort$ERG,tosort$ETV1,tosort$ETV4,tosort$ETV5,
                             tosort$CHD1,tosort$SPOP,
                             tosort$RB1,tosort$PTEN,tosort$MYC),
                       'sample_ids_wgbs'])

df = data.frame(
    x=sample_ids_wgbs,
    y=rep(-1,length(sample_ids_wgbs)),
    ht=height,
    fill=rep(NA,length(sample_ids_wgbs)),
    stringsAsFactors=F)
rownames(df) = sample_ids_wgbs
df[intersect(names( metastasis_locations )[metastasis_locations=='Bone'],sample_ids_wgbs),'fill'] = site[1]
df[intersect(names( metastasis_locations )[metastasis_locations=='Lymph_node'],sample_ids_wgbs),'fill'] = site[2]
df[intersect(names( metastasis_locations )[metastasis_locations=='Liver'],sample_ids_wgbs),'fill'] = site[3]
df[intersect(names( metastasis_locations )[metastasis_locations=='Other'],sample_ids_wgbs),'fill'] = site[4]

#Loop through genes and set up plotting data frames
for(i in 1:numgenes) {
    genei = paste(genes2use[i,'gene'])
    print(genei)
    gofi = genes2use[i,'gof']
    all = allele_effect(genei)$alleles
    
    df1 = data.frame(
        x=sample_ids_wgbs,
        y=rep((i*4)-3,length(sample_ids_wgbs)),
        ht=height,
        fill=rep(NA,length(sample_ids_wgbs)),
        stringsAsFactors=F)
    df2 = data.frame(
        x=sample_ids_wgbs,
        y=rep((i*4)-2,length(sample_ids_wgbs)),
        ht=height,
        fill=rep(NA,length(sample_ids_wgbs)),
        stringsAsFactors=F)
    df3 = data.frame(
        x=sample_ids_wgbs,
        y=rep((i*4)-1,length(sample_ids_wgbs)),
        ht=height,
        fill=rep(NA,length(sample_ids_wgbs)),
        stringsAsFactors=F)
    rownames(df1) = sample_ids_wgbs
    rownames(df2) = sample_ids_wgbs
    rownames(df3) = sample_ids_wgbs
    
    rowensembl = ensembl2sym$name==genei
    ensemblid = rownames(ensembl2sym)[rowensembl]
    stopifnot(sum(rowensembl)==1)
    chr = ensembl2sym[rowensembl,'chr']
    gene_start = ensembl2sym[rowensembl,'start']
    gene_end = ensembl2sym[rowensembl,'end']
    
    #Just plot the promoter
    if(ensembl2sym[ensemblid,'strand']=='+') {
        start = gene_start-1500
        end = gene_start+1500
        #end = gene_end
    } else {
        #start = gene_start
        start = gene_end-1500
        end = gene_end+1500
    }
    gr = GRanges(chr,IRanges(start,end))
    
    #get recurrent HMRs
    emr_gene = tracks[['HMRseggene']][tracks[['HMRseggene']]$gene==genei]
    emr_gene = emr_gene[countOverlaps(emr_gene,gr)>0 & !is.na(emr_gene$avg_cor)]
    #print(emr_gene)
    if(length(emr_gene)>0) {
        #Get the best eHMR
        #print(emr_gene)
        bestcor = min(emr_gene$avg_p)
        maxrow = emr_gene$avg_p==bestcor
        emr_best = emr_gene[maxrow]
        stopifnot(length(emr_best)==1)
        
        #Get the methylation in the eHMR
        emr_start = start(emr_best)
        emr_end = end(emr_best)
        methylmean = meanmethyl(ensemblid, dir_gene_output_mcrpc, emr_start,emr_end)[sample_ids_wgbs]
        methylmean_nt = meanmethyl(ensemblid, filepathctrl, emr_start,emr_end)
        names(methylmean_nt) = gsub('_R1_BISMARK_BT2_PE','',names(methylmean_nt),fixed=T)
        methylmean_nt = methylmean_nt[samples_ctrl]
        
        #Compare against the control
        if(ctrlgroup==1) {
            methyldiff = rep(NA,length(methylmean))
            for(loci in unique(locations)) {
                samplestosub = names(methylmean) %in% names(locations)[locations==loci]
                if(loci %in% locs_normal) {
                    tosub = mean(methylmean_nt[names(locs_normal)[locs_normal==loci]],na.rm=T)
                } else {
                    tosub = mean(methylmean_nt,na.rm=T)
                }
                methyldiff[samplestosub] = methylmean[samplestosub] - tosub
            }
        } else {
            methyldiff = methylmean-(mean(methylmean_nt,na.rm=T))
        }
        
        for(cutoff in groups) {
            rowspos = methyldiff>=cutoff
            rowspos[is.na(rowspos)] = FALSE
            df3[rowspos,'fill'] = paste('+',cutoff,'%',sep='')
            rowsneg = methyldiff<=(-cutoff)
            rowsneg[is.na(rowsneg)] = FALSE
            df3[rowsneg,'fill'] = paste('-',cutoff,'%',sep='')
        }
    } else {
        print(paste('No EMR for',genei))
    }
    
    #Add new AR enhancer amplification
    if(genei=='AR') { #Add all eHMR enhancers too
        chrx_1 = cnregion('chrX',66744900,66745900)[sample_ids_wgbs]
        df1[names(chrx_1)[chrx_1 >=GAIN_SEX],'fill'] = 'CNA_amp'
        chrx_2 = cnregion('chrX',66900400,66902700)[sample_ids_wgbs]
        df1[names(chrx_2)[chrx_2 >=GAIN_SEX],'fill'] = 'CNA_amp'
        chrx_3 = cnregion('chrX',66906200,66908000)[sample_ids_wgbs]
        df1[names(chrx_3)[chrx_3 >=GAIN_SEX],'fill'] = 'CNA_amp'
        chrx_4 = cnregion('chrX',67043300,67045400)[sample_ids_wgbs]
        df1[names(chrx_4)[chrx_4 >=GAIN_SEX],'fill'] = 'CNA_amp'
        chrx_5 = cnregion('chrX',67361700,67362400)[sample_ids_wgbs]
        df1[names(chrx_5)[chrx_5 >=GAIN_SEX],'fill'] = 'CNA_amp'
        chrx_6 = cnregion('chrX',67543100,67552000)[sample_ids_wgbs]
        df1[names(chrx_6)[chrx_6 >=GAIN_SEX],'fill'] = 'CNA_amp'
        chrx_7 = cnregion('chrX',67799500,67800400)[sample_ids_wgbs]
        df1[names(chrx_7)[chrx_7 >=GAIN_SEX],'fill'] = 'CNA_amp'
    }
    
    #Add DNA alterations
    if(gofi) {
        for(j in 1:length(gof)) {
            namegof = gof[j]
            samplesalt = sample_ids_wgbs[all[sample_ids_wgbs,namegof]]
            if(genei=='AR' && namegof=='CNA_amp') {
                df1[samplealt,'fill'] = namegof
            } else {
                for(samplealt in samplesalt) {
                    count = is.na(df1[samplealt,'fill']) + is.na(df2[samplealt,'fill'])
                    if(count==2) {
                        df1[samplealt,'fill'] = namegof
                    } else if(count==1) {
                        df2[samplealt,'fill'] = namegof
                    }
                }
            }
        }
    } else {
        for(j in 1:length(lof)) {
            namelof = lof[j]
            samplesalt = sample_ids_wgbs[all[sample_ids_wgbs,namelof]]
            for(samplealt in samplesalt) {
                if(namelof=='CNA_2' || (namelof=='inactivating_missense' && genei=='CDK12')) {
                    df1[samplealt,'fill'] = namelof
                    df2[samplealt,'fill'] = namelof
                } else {
                    count = is.na(df1[samplealt,'fill']) + is.na(df2[samplealt,'fill'])
                    if(count==2) {
                        df1[samplealt,'fill'] = namelof
                    } else if(count==1) {
                        df2[samplealt,'fill'] = namelof
                    } else {
                        print(paste(genei,samplealt,namelof))
                    }
                }
            }
        }
    }
    df = rbind.data.frame(df,df1,df2,df3)
}
df$x = factor( df$x, levels=collvls )
rowsna = is.na(df$fill)
df = df[!rowsna,]

rowsalt = df$fill %in% c(gof,lof)
df[rowsalt,'fill'] = map[df[rowsalt,'fill']]
df$fill = factor(df$fill,levels=c(site,' ',unique(map),'  ',groupnames))

reds = colorRampPalette(c('white','red'))(n=10)
blues = colorRampPalette(c('white','blue'))(n=10)
colpal = c(blues[10],blues[7],blues[5],reds[5],reds[7],reds[10])

xpos = -10
p = ggplot(df,aes(x,y,width=0.5,height=ht))+geom_tile(aes(fill=fill))+
    scale_fill_manual(drop=FALSE,name="",
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
                          "-30%"=colpal[1],
                          "-20%"=colpal[2],
                          "-10%"=colpal[3],
                          "+10%"=colpal[4],
                          "+20%"=colpal[5],
                          "+30%"=colpal[6]))+
    guides(fill=F)+
    theme_classic()+xlab('Samples')+
    theme(plot.margin=unit(c(0,3.2,0,0.2),"cm"))+
    theme(axis.text.x=element_text(angle=90, hjust=1,size=5))+
    scale_y_reverse(name='',minor_breaks=(0:numgenes)*4,
                    breaks=(c(1,(1:numgenes)*4)-2),labels=c('Biopsy Site',paste(genes2use$gene)))+
    coord_cartesian(ylim=c(84,2))+
    theme(panel.grid.minor=element_line(colour='black',size=0.5))+
    annotation_custom(grob=textGrob(expression(underline('  ETS  ')),
                                    hjust='center',rot=90,gp=gpar(cex=1)),
                      ymin=-6,ymax=-6,xmin=xpos,xmax=xpos)+
    annotation_custom(grob=textGrob(expression(underline('      AR      ')),
                                    hjust='center',rot=90,gp=gpar(cex=1)),
                      ymin=-20,ymax=-20,xmin=xpos,xmax=xpos)+
    annotation_custom(grob=textGrob(expression(underline('  lncRNA  ')),
                                    hjust='center',rot=90,gp=gpar(cex=1)),
                      ymin=-36,ymax=-36,xmin=xpos,xmax=xpos)+
    annotation_custom(grob=textGrob(expression(underline('Oncogene')),
                                    hjust='center',rot=90,gp=gpar(cex=1)),
                      ymin=-50,ymax=-50,xmin=xpos,xmax=xpos)+
    annotation_custom(grob=textGrob(expression(underline('      Tumor suppressor      ')),
                                    hjust='center',rot=90,gp=gpar(cex=1)),
                      ymin=-72,ymax=-72,xmin=xpos,xmax=xpos)+
    #annotation_custom(grob=textGrob(expression(underline('DNA repair')),
    #	hjust='center',rot=90,gp=gpar(cex=1)),
    #	ymin=-86,ymax=-86,xmin=xpos,xmax=xpos)+
    annotation_custom(grob=linesGrob(),xmin=101,xmax=103,ymin=-81,ymax=-79)+
    annotation_custom(grob=linesGrob(),xmin=101,xmax=103,ymin=-82,ymax=-81)+
    annotation_custom(grob=linesGrob(),xmin=101,xmax=103,ymin=-83,ymax=-83)+
    annotation_custom(grob=textGrob('DNA alteration 1',
                                    hjust=0,gp=gpar(cex=0.8)),
                      ymin=-79,ymax=-79,xmin=103.5,xmax=103.5)+
    annotation_custom(grob=textGrob('DNA alteration 2',
                                    hjust=0,gp=gpar(cex=0.8)),
                      ymin=-81,ymax=-81,xmin=103.5,xmax=103.5)+
    annotation_custom(grob=textGrob('Methylation',
                                    hjust=0,gp=gpar(cex=0.8)),
                      ymin=-83,ymax=-83,xmin=103.5,xmax=103.5)

#Add custom legends
legx = 105
legtx = 107
width = 5

ystart = 0
p = p +
    annotation_custom(grob=linesGrob(gp=gpar(col=wes_palette('Chevalier1')[1],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-5,ymax=ystart-7)+
    annotation_custom(grob=linesGrob(gp=gpar(col=wes_palette('Chevalier1')[2],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-8,ymax=ystart-10)+
    annotation_custom(grob=linesGrob(gp=gpar(col=wes_palette('Chevalier1')[4],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-11,ymax=ystart-13)+
    annotation_custom(grob=linesGrob(gp=gpar(col=wes_palette('Chevalier1')[3],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-14,ymax=ystart-16)+
    annotation_custom(grob=textGrob('Bone',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-6,ymax=ystart-6,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('LN',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-9,ymax=ystart-9,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('Liver',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-12,ymax=ystart-12,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('Other',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-15,ymax=ystart-15,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('Biopsy',hjust=0,gp=gpar(cex=1)),ymin=ystart,ymax=ystart,xmin=legx-2,xmax=legx-2)+
    annotation_custom(grob=textGrob('Site',hjust=0,gp=gpar(cex=1)),ymin=ystart-3,ymax=ystart-3,xmin=legx-2,xmax=legx-2)

ystart = -20
p = p +
    annotation_custom(grob=linesGrob(gp=gpar(col='purple',lwd=width)),xmin=legx,xmax=legx,ymin=ystart-5,ymax=ystart-7)+
    annotation_custom(grob=linesGrob(gp=gpar(col='orange',lwd=width)),xmin=legx,xmax=legx,ymin=ystart-8,ymax=ystart-10)+
    annotation_custom(grob=linesGrob(gp=gpar(col='green',lwd=width)),xmin=legx,xmax=legx,ymin=ystart-11,ymax=ystart-13)+
    annotation_custom(grob=linesGrob(gp=gpar(col='cyan4',lwd=width)),xmin=legx,xmax=legx,ymin=ystart-14,ymax=ystart-16)+
    annotation_custom(grob=textGrob('Mutation',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-6,ymax=ystart-6,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('SV',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-9,ymax=ystart-9,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('CN Gain',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-12,ymax=ystart-12,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('CN Loss',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-15,ymax=ystart-15,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('DNA',hjust=0,gp=gpar(cex=1)),ymin=ystart,ymax=ystart,xmin=legx-2,xmax=legx-2)+
    annotation_custom(grob=textGrob('Alteration',hjust=0,gp=gpar(cex=1)),ymin=ystart-3,ymax=ystart-3,xmin=legx-2,xmax=legx-2)

ystart = -40
p = p +
    annotation_custom(grob=linesGrob(gp=gpar(col=colpal[1],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-5,ymax=ystart-7)+
    annotation_custom(grob=linesGrob(gp=gpar(col=colpal[2],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-8,ymax=ystart-10)+
    annotation_custom(grob=linesGrob(gp=gpar(col=colpal[3],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-11,ymax=ystart-13)+
    annotation_custom(grob=linesGrob(gp=gpar(col=colpal[4],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-14,ymax=ystart-16)+
    annotation_custom(grob=linesGrob(gp=gpar(col=colpal[5],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-17,ymax=ystart-19)+
    annotation_custom(grob=linesGrob(gp=gpar(col=colpal[6],lwd=width)),xmin=legx,xmax=legx,ymin=ystart-20,ymax=ystart-22)+
    annotation_custom(grob=textGrob('-30%',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-6,ymax=ystart-6,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('-20%',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-9,ymax=ystart-9,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('-10%',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-12,ymax=ystart-12,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('+10%',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-15,ymax=ystart-15,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('+20%',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-18,ymax=ystart-18,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('+30%',hjust=0,gp=gpar(cex=0.8)),ymin=ystart-21,ymax=ystart-21,xmin=legtx,xmax=legtx)+
    annotation_custom(grob=textGrob('Promoter',hjust=0,gp=gpar(cex=1)),ymin=ystart,ymax=ystart,xmin=legx-2,xmax=legx-2)+
    annotation_custom(grob=textGrob('Methylation',hjust=0,gp=gpar(cex=1)),ymin=ystart-3,ymax=ystart-3,xmin=legx-2,xmax=legx-2)

g = ggplotGrob(p)
g$layout$clip[g$layout$name=="panel"] = "off"
g$layout$clip[g$layout$name=="legend"] = "off"
grid.draw(g)

#ggsave(fileout,g,width=10,height=6.66667)


