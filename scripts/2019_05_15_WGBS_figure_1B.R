# Plots Figure 1B

df_hmr_recurrent = read.delim(fn_hmr_recurrent, sep='\t', header=T, check.names=F)
samples_plot = intersect(sample_ids_wgbs,colnames(df_hmr_recurrent))

# restrict to the top 10% of variability in HMR by standard deviation
sds = apply(df_hmr_recurrent[,samples_plot],1,sd)
rows_top_decile = sds >= quantile( sds, 0.9, na.rm=T )
rows_top_decile[ is.na(rows_top_decile) ] = FALSE
df_hmr_recurrent = df_hmr_recurrent[rows_top_decile,]

mat_hmr_recurrent = as.matrix(df_hmr_recurrent[,samples_plot])
gr = makeGRangesFromDataFrame(df_hmr_recurrent)
cgi = rep('white',dim(mat_hmr_recurrent)[1])
rsidebar = t(matrix(cgi))

toplot = data.frame(samples_plot)
toplot$ERG = allele_effect('ERG')$alleles[samples_plot,'activating_sv']
toplot$ETV1 = allele_effect('ETV1')$alleles[samples_plot,'activating_sv']
toplot$ETV4 = allele_effect('ETV4')$alleles[samples_plot,'activating_sv']
toplot$ETV5 = allele_effect('ETV5')$alleles[samples_plot,'activating_sv']
toplot$CHD1 = allele_effect('CHD1')$alleles[samples_plot,'n_alleles_inactivated']==2
toplot$SPOP = allele_effect('SPOP')$alleles[samples_plot,'inactivating_missense']
toplot$RB1 = allele_effect('RB1')$alleles[samples_plot,'n_alleles_inactivated']==2
toplot$PTEN = allele_effect('PTEN')$alleles[samples_plot,'n_alleles_inactivated']==2
toplot$MYC = allele_effect('MYC')$alleles[samples_plot,'CNA_amp']
toplot$TP53 = allele_effect('TP53')$alleles[samples_plot,'n_alleles_inactivated']==2

CMP = samples_plot %in% samples_cmp
ETS = toplot$ERG|toplot$ETV1|toplot$ETV4|toplot$ETV5
CHD1_SPOP = toplot$CHD1 | toplot$SPOP
PTEN = toplot$PTEN
RB1 = toplot$RB1
TP53 = toplot$TP53
MYC = toplot$MYC
Mutation = samples_plot %in% c(samples_braf,samples_idh1,samples_tet2,samples_dnmt3b)
Sites = metastasis_locations[samples_plot]
print(fisher.test(CMP,ETS))
print(fisher.test(CMP,TP53))
print(fisher.test(CMP,MYC))
print(fisher.test(CMP,Mutation))
print(fisher.test(CMP,Sites))

colside = rep('white',dim(mat_hmr_recurrent)[2])
names(colside) = colnames(mat_hmr_recurrent)
stopifnot(length(colside)==length(samples_plot))

tSCNC = colside
tSCNC[samples_tscnc] = 'gray30'
stopifnot(length(tSCNC)==length(samples_plot))

ETS = colside
ETS[toplot$ERG] = 'gray30'
ETS[toplot$ETV1] = 'gray30'
ETS[toplot$ETV4] = 'gray30'
ETS[toplot$ETV5] = 'gray30'
stopifnot(length(ETS)==length(samples_plot))

MYC = colside
MYC[toplot$MYC] = 'gray30'
stopifnot(length(MYC)==length(samples_plot))

SPOP = colside
SPOP[toplot$SPOP] = 'gray30'
stopifnot(length(SPOP)==length(samples_plot))

RB1 = colside
RB1[toplot$RB1] = 'gray30'
stopifnot(length(RB1)==length(samples_plot))

PTEN = colside
PTEN[toplot$PTEN] = 'gray30'
stopifnot(length(PTEN)==length(samples_plot))

TP53 = colside
TP53[toplot$TP53] = 'gray30'
stopifnot(length(TP53)==length(samples_plot))

CMP = colside
CMP[samples_cmp] = 'mediumturquoise'#'lightsteelblue'
CMP[samples_braf] = 'orange'
CMP[samples_idh1] = 'green'
CMP[samples_tet2] = 'purple'
CMP[samples_dnmt3b] = 'deeppink'
#CMP['DTB-252-BL'] = 'black'
stopifnot(length(CMP)==length(samples_plot))

Site = colside
Site[intersect(names(metastasis_locations)[metastasis_locations=='Bone'],samples_plot)] = wes_palette('Chevalier1')[1] #darkgreen
Site[intersect(names(metastasis_locations)[metastasis_locations=='Lymph_node'],samples_plot)] = wes_palette('Chevalier1')[2] #gold
Site[intersect(names(metastasis_locations)[metastasis_locations=='Liver'],samples_plot)] = wes_palette('Chevalier1')[4] #brown
Site[intersect(names(metastasis_locations)[metastasis_locations=='Other'],samples_plot)] = wes_palette('Chevalier1')[3] #gray

purity = colside
purity[tumorpurity[colnames(mat),'Tumor.Purity.Histo']<50] = 'black'
stopifnot(length(purity)==length(samples_plot))

csidebar = cbind(Site,CMP,tSCNC,TP53,RB1,PTEN,MYC,ETS)
color = colorRampPalette(c('blue','white','red'))(n = 1000)

hclust2 = function(x, method="average") {
    hclust(x, method=method)
}
dist2 = function(x) {
    as.dist(1-cor(t(x), method="spearman"))
}

main_title = ""
heatmap.3(mat_hmr_recurrent, 
          hclustfun=hclust, 
          distfun=dist,
          na.rm=TRUE, 
          scale='none', 
          dendrogram="column", margins=c(4,4),
          Rowv=TRUE, Colv=TRUE, 
          ColSideColors=csidebar, 
          RowSideColors=rsidebar, 
          symbreaks=FALSE, 
          key=TRUE, 
          symkey=FALSE,
          density.info="none", trace="none",
          labCol=colnames(mat_hmr_recurrent), 
          labRow=rep("",(dim(mat_hmr_recurrent)[1])), 
          main=main_title, cexRow=1, cexCol=1, col=color,
          ColSideColorsSize=8, RowSideColorsSize=1)
