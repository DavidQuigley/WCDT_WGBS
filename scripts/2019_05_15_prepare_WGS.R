N_SAMPLES=101
dir_root = make_path(BASE_DIR, "secondary_data/WGS")
dir_results = paste( dir_root, "/results", sep='')
dir_meta = paste( dir_root, "/metadata", sep='')

fn_curated = paste( dir_meta,'/2018_05_07_curated_AR_SV.txt',sep='')
fn_sa = paste(dir_meta, '20180123_sample_attributes.txt', sep='/')
fn_clinical = paste(dir_meta, '2018_01_10_clinical_attributes.txt', sep='/')
fn_cyto = paste( dir_meta, '/HG38_cytoband_std_chrom.txt', sep='')
fn_gene_locs = paste(dir_meta, 'GRCh38Decoy_refseq_genelocs_from_refFlat.bed',sep='/')
fn_chrom_lengths = paste(dir_meta, 'HG38_chromosome_lengths.txt',sep='/')
fn_centromeres = paste(dir_meta, 'HG38_centromere_loci.txt',sep='/')
fn_curated_sv = paste(dir_meta, '/2018_03_15_curated_SV_table.txt', sep='')
fn_curated_missense = paste(dir_meta, '/2018_05_14_curated_missense_table.txt', sep='')
fn_curated_fs = paste(dir_meta, '/2018_05_14_curated_frameshift_table.txt', sep='')
fn_symbols_fig5 = paste(dir_meta, '/2018_05_09_symbols_figure_5.txt', sep='')

build_id='2018_04_15'
fn_cna_en = paste( dir_results, '/', build_id, '_matrix_AR_enhancer_CN.txt', sep='')
fn_sample_summary = paste(dir_results, '/', build_id, '_matrix_sample_summary.txt', sep='')
fn_binned_CN = paste(dir_results, '/', build_id, '_matrix_binned_weighted_CN_copycat.txt', sep='')
fn_cna = paste(dir_results, '/', build_id, '_matrix_CNA_symbol_copycat.txt', sep='')
fn_cna_int = paste(dir_results, '/', build_id, '_matrix_CN_integer_symbol_copycat.txt', sep='')
fn_germline = paste(dir_results, '/', build_id, '_matrix_germline_summary.txt', sep='')
fn_somatic = paste(dir_results, '/', build_id, '_list_somatic_strelka.txt', sep='')
fn_somatic_mutect = paste(dir_results, '/', build_id, '_list_somatic_mutect.txt', sep='')
fn_missense = paste(dir_results, '/', build_id, '_matrix_somatic_missense_strelka.txt', sep='')
fn_inactivating = paste(dir_results, '/', build_id, '_matrix_somatic_inactivating_strelka.txt', sep='')
fn_mutcount = paste(dir_results, '/', build_id, '_matrix_mutation_count_summary.txt', sep='')
fn_manta_sv = paste(dir_results, '/', build_id, '_list_manta_SV.txt', sep='')
fn_lumpy_sv = paste(dir_results, '/', build_id, '_list_lumpy_SV.txt', sep='')

fn_tpm = paste(dir_results, '/', build_id, '_matrix_rna_tpm.txt', sep='')
fn_tpm_full = paste(dir_results, '/', build_id, '_matrix_rna_tpm_full.txt', sep='')
fn_alignment_summary_t = paste(dir_results, '/', build_id, '_matrix_alignment_summary_tumor.txt', sep='')
fn_alignment_summary_n = paste(dir_results, '/', build_id, '_matrix_alignment_summary_normal.txt', sep='')
fn_genomewide_sliding_scan = paste(dir_results, '/genome_wide_SV_scan.txt', sep='')
fn_assignment = paste( dir_meta, '/sv_window_assignment.txt', sep='')

################################################################################
# Genome data
################################################################################
cyto_info = read.table(fn_cyto, stringsAsFactors=FALSE,header=TRUE, sep='\t')
chrom_lengths = read.table(fn_chrom_lengths, row.names=1, stringsAsFactors=FALSE)
centromeres = read.table( fn_centromeres, header=TRUE, stringsAsFactors = FALSE)
gene_locs = read.table( fn_gene_locs, stringsAsFactors = FALSE)
names(gene_locs) = c("chrom", "start", "end", "symbol", "score", "strand")
gene_locs$symbol = get.split.col(gene_locs$symbol, "~", first=TRUE)
seen = hsh_new()
chrom_valid = hsh_from_vectors(rownames(chrom_lengths), 1:dim(chrom_lengths)[1])
keep = rep(FALSE, dim(gene_locs)[1])
for(i in 1:dim(gene_locs)[1]){
    keep[i] = hsh_in( chrom_valid, gene_locs$chrom[i]) &
              !hsh_in( seen, gene_locs$symbol[i] )
    hsh_set( seen, gene_locs$symbol[i], 1 )
}
gene_locs = gene_locs[keep,]
rownames(gene_locs) = gene_locs$symbol
chrom_idx = get.split.col( gene_locs$chrom, "chr", last=TRUE)
chrom_idx[chrom_idx=="X"]=23
chrom_idx[chrom_idx=="Y"]=24
gene_locs$chrom_idx = as.numeric(chrom_idx)
rm(keep)
rm(seen)
rm(chrom_idx)


################################################################################
# Sample matrix file
################################################################################

matrix_samples = read.table(fn_sample_summary, header=TRUE, sep='\t',
                            stringsAsFactors=FALSE, check.names = FALSE, 
                            row.names = 1)

#canonical list and order set here
sample_ids = rownames(matrix_samples)  
n_samples = length(sample_ids)
sa = load.matrix( fn_sa )


################################################################################
# CNA: binned in 3mb windows across genome
################################################################################

cnbin = read.table(fn_binned_CN,
                   header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                   check.names = FALSE)

cols_cnbin = rep("blue", dim(cnbin)[1])
cols_cnbin[ cnbin$chrom %in% c("chr2","chr4","chr6","chr8","chr10","chr12",
                               "chr14","chr16","chr18","chr20","chr22",
                               "chrY")] = "coral"

matrix_CNA = read.table(fn_cna,
                        header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                        check.names = FALSE, row.names = 1)
matrix_CNA_int = read.table(fn_cna_int,
                            header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                            check.names = FALSE, row.names = 1)

cnamp = rowSums( cnbin[,3:dim(cnbin)[2]] >=3 ) / dim(cnbin)[2]
cndel = rowSums( cnbin[,3:dim(cnbin)[2]] <=1 ) / dim(cnbin)[2]
sexchrom = cnbin$chrom=="chrX" | cnbin$chrom=="chrY"
cnamp[sexchrom] = rowSums( cnbin[sexchrom,3:dim(cnbin)[2]] >=1.5 ) / dim(cnbin)[2]
cndel[sexchrom] = rowSums( cnbin[sexchrom,3:dim(cnbin)[2]] <=0.75 ) / dim(cnbin)[2]
cndel = -1 * cndel

# make sure genes in matrix_CNA, matirx_CNA_int match gene_locs 
m = match.idx( dimnames(matrix_CNA)[[1]], rownames(gene_locs))
matrix_CNA = matrix_CNA[m$idx.A,]
matrix_CNA_int = matrix_CNA_int[m$idx.A,]
gene_locs = gene_locs[m$idx.B,]

matrix_CNA_int_ploidy=matrix_CNA_int

# canonical symbol list and order set here
#-------------------------------------------------------------------------------
symbols = dimnames(matrix_CNA)[[1]]
n_symbols = length(symbols)


# recalculate percent altered 
is_sex = cnbin$chrom=="chrX" | cnbin$chrom=="chrY"
n_alt_notsex = colSums( cnbin[!is_sex,3:103] >= GAIN_NONSEX, na.rm=TRUE ) +
               colSums( cnbin[!is_sex,3:103] <= LOSS_SINGLE_NONSEX, na.rm=TRUE )
n_alt_sex = colSums( cnbin[is_sex,3:103] >= GAIN_SEX, na.rm=TRUE ) +
            colSums( cnbin[is_sex,3:103] <= LOSS_SEX, na.rm=TRUE )

percent_alt = (n_alt_notsex + n_alt_sex) / dim(cnbin)[1]
matrix_samples$percent_CNA_ref = 1-percent_alt


################################################################################
# Germline mutations that potentially are pathogenic
################################################################################

germ = read.table(fn_germline, header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                  check.names = FALSE, row.names = 1)
germ = germ[ match.idx( sample_ids, rownames(germ))$idx.B,]
germ = t(germ)
matrix_germline = matrix('.', nrow=n_symbols, ncol=n_samples, 
                         dimnames=list( y=symbols, x=sample_ids))

m_symbol = match.idx( rownames(matrix_germline), rownames(germ))
m_sample = match.idx( dimnames(matrix_germline)[[2]], dimnames(germ)[[2]] )
matrix_germline[ m_symbol$idx.A, m_sample$idx.B] = germ[m_symbol$idx.B, m_sample$idx.B]


###############################################
# Somatic mutation data
###############################################

list_somatic_data = read.table(fn_somatic,
                          header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                          check.names = FALSE)

# load inactivating somatic mutations; row is sample, col is symbol
som_inactive = read.table(fn_inactivating, header=TRUE, sep='\t', 
                          stringsAsFactors=FALSE, check.names = FALSE,
                          row.names = 1)
m = match.idx( sample_ids, rownames(som_inactive) )
som_inactive = som_inactive[m$idx.B,]
matrix_inactive = matrix('.', nrow=n_symbols, ncol=n_samples, 
                  dimnames=list( y=symbols, x=sample_ids))
som_inactive=t(som_inactive)
m = match.idx( rownames(matrix_inactive), rownames(som_inactive))
matrix_inactive[ m$idx.A,] = som_inactive[m$idx.B,]
rm(som_inactive)

# load missense somatic mutations; row is sample, col is symbol
som_missense = read.table(fn_missense,
                          header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names = FALSE,
                          row.names = 1)
matrix_missense = matrix('.', nrow=n_symbols, ncol=n_samples, 
                  dimnames=list( y=symbols, x=sample_ids))
m = match.idx( sample_ids, rownames(som_missense) )
som_missense = som_missense[m$idx.B,]
som_missense=t(som_missense)
m = match.idx( rownames(matrix_missense), rownames(som_missense))
matrix_missense[ m$idx.A,] = som_missense[m$idx.B,]
rm(som_missense)

matrix_mutcount  = read.table(fn_mutcount, header=TRUE, sep='\t', 
                              stringsAsFactors=FALSE, check.names = FALSE,
                            row.names = 1)
matrix_mutcount = cbind(matrix_mutcount, dummy=rep('.', dim(matrix_mutcount)[1]), stringsAsFactors=FALSE)
m = match.idx( sample_ids, rownames(matrix_mutcount))
matrix_mutcount = matrix_mutcount[m$idx.A,]

list_somatic_data_mutect = read.table(fn_somatic_mutect,
                               header=TRUE, sep='\t', stringsAsFactors=FALSE, 
                               check.names = FALSE)

# load inactivating somatic mutations; row is sample, col is symbol
som_inactive_mutect = read.table(fn_inactivating, header=TRUE, sep='\t', 
                          stringsAsFactors=FALSE, check.names = FALSE,
                          row.names = 1)
m = match.idx( sample_ids, rownames(som_inactive_mutect) )
som_inactive_mutect = som_inactive_mutect[m$idx.B,]
matrix_inactive_mutect = matrix('.', nrow=n_symbols, ncol=n_samples, 
                         dimnames=list( y=symbols, x=sample_ids))
som_inactive_mutect=t(som_inactive_mutect)
m = match.idx( rownames(matrix_inactive_mutect), rownames(som_inactive_mutect))
matrix_inactive_mutect[ m$idx.A,] = som_inactive_mutect[m$idx.B,]
rm(som_inactive_mutect)

# load missense somatic mutations; row is sample, col is symbol
som_missense_mutect = read.table(fn_missense,
                          header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names = FALSE,
                          row.names = 1)
matrix_missense_mutect = matrix('.', nrow=n_symbols, ncol=n_samples, 
                         dimnames=list( y=symbols, x=sample_ids))
m = match.idx( sample_ids, rownames(som_missense_mutect) )
som_missense_mutect = som_missense_mutect[m$idx.B,]
som_missense_mutect=t(som_missense_mutect)
m = match.idx( rownames(matrix_missense_mutect), rownames(som_missense_mutect))
matrix_missense_mutect[ m$idx.A,] = som_missense_mutect[m$idx.B,]
rm(som_missense_mutect)

matrix_mutcount  = read.table(fn_mutcount, header=TRUE, sep='\t', 
                              stringsAsFactors=FALSE, check.names = FALSE,
                              row.names = 1)
matrix_mutcount = cbind(matrix_mutcount, dummy=rep('.', dim(matrix_mutcount)[1]), stringsAsFactors=FALSE)
m = match.idx( sample_ids, rownames(matrix_mutcount))
matrix_mutcount = matrix_mutcount[m$idx.A,]


# compare mutect and strelka
mu = paste(list_somatic_data_mutect$sample_id,list_somatic_data_mutect$symbol, sep='|')
st = paste(list_somatic_data$sample_id,list_somatic_data$symbol, sep='|')



################################################################################
# Structural variants from manta (list_sv_m) 
################################################################################
list_sv_m = read.table(fn_manta_sv,
                  header=TRUE, sep='\t', 
                  stringsAsFactors=FALSE, check.names = FALSE)
list_sv_m = list_sv_m[list_sv_m$chrom_start != "chrM",]

# eliminate artifact on chromosome 22
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3996660/
list_sv_m= list_sv_m[!(list_sv_m$chrom_start=="chr22" & list_sv_m$svtype=="BND" & list_sv_m$pos_start>=28597176 & list_sv_m$pos_start<=29037622),]

chrom_idx = get.split.col( list_sv_m$chrom_start, "chr", last=TRUE)
chrom_idx[chrom_idx=="X"] = 23
chrom_idx[chrom_idx=="Y"] = 24
chrom_idx = as.numeric(chrom_idx)
list_sv_m = cbind( list_sv_m, chrom_idx, stringsAsFactors=FALSE)

# remove duplicate BND events
keep = rep(TRUE, dim(list_sv_m)[1])
idx_BND = which(list_sv_m$svtype=="BND")
seen = hsh_new()

as_written = paste(list_sv_m$chrom_start, list_sv_m$pos_start, list_sv_m$chrom_end, list_sv_m$pos_end)
reversed = paste(list_sv_m$chrom_end, list_sv_m$pos_end, list_sv_m$chrom_start, list_sv_m$pos_start)

for(i in 1:length(idx_BND)){
    if( !hsh_in( seen, as_written[idx_BND[i]] ) & 
        !hsh_in( seen, reversed[idx_BND[i]] ) ){
        hsh_set( seen, as_written[idx_BND[i]], 1)
        keep[idx_BND[i]] = TRUE
    }else{
        keep[idx_BND[i]] = FALSE
    }
}
list_sv_m = list_sv_m[keep,]


##############################################
# Summarize counts of SV events
##############################################
matrix_sv_m = data.frame( n_dels=matrix_samples$n_del_manta,
                          n_tandems=matrix_samples$n_dup_manta,
                          n_mh = matrix_samples$n_del_MHgt2_manta,
                          n_ins=matrix_samples$n_ins_manta,
                          n_inv=matrix_samples$n_inv_manta,
                          per=round(matrix_samples$n_del_MHgt2_manta/matrix_samples$n_del_manta,2), 
                          row.names=dimnames(matrix_samples)[[1]], 
                          stringsAsFactors=FALSE)
matrix_sv_m$n_bnd = as.numeric(table(list_sv_m$sample_id[list_sv_m$svtype=="BND"]))



################################################################################
# RNA
################################################################################

tpm = read.table(fn_tpm,
                 header=TRUE, sep='\t', stringsAsFactors=FALSE, check.names = FALSE,
                 row.names = 1)
m=match.idx( sample_ids, names(tpm) )
tpm = data.matrix(tpm)
matrix_tpm = matrix( NA, nrow=dim(tpm)[1], ncol=length(sample_ids) )
matrix_tpm[,m$idx.A] = tpm[,m$idx.B ]
matrix_tpm = data.frame(matrix_tpm)
names(matrix_tpm) = sample_ids
rownames(matrix_tpm) = dimnames(tpm)[[1]]





################################################################################
# Gene activation/ inactivation
################################################################################

matrix_somatic = matrix_inactive != '.' | matrix_missense != '.' 
sums_inactive = rowSums( matrix_inactive != '.' )
sums_missense = rowSums( matrix_missense != '.' )
sums_somatic = rowSums( matrix_missense != '.' | matrix_inactive != '.')

genes = rownames(matrix_CNA_int_ploidy)
n_genes = length(genes)
sexchrom = which( gene_locs$chrom == "chrX" | gene_locs$chrom=="chrY" )

curated_fs = read.table(fn_curated_fs,header=TRUE,sep='\t',stringsAsFactors=FALSE)
symbols_curated_fs = unique( curated_fs$symbol )
curated_sv = read.table(fn_curated_sv,header=TRUE,sep='\t',stringsAsFactors=FALSE)
curated_sv = curated_sv[curated_sv$reported=="yes",]
whitelist_SV = c( unique(curated_sv$threeprime), "AR" )
curated_missense = read.table(fn_curated_missense,header=TRUE,sep='\t',stringsAsFactors=FALSE)
symbols_curated_missense = unique( curated_missense$gene )
curated_missense = curated_missense[curated_missense$reported_pathogenic!="no",]


INACTIVE = matrix(0, nrow=length(genes), ncol=3)
ALLELES_INACTIVATED = matrix(0, nrow=length(genes), ncol=length(sample_ids))
matrix_SV = matrix('.', nrow=length(genes), ncol=length(sample_ids))
dimnames(ALLELES_INACTIVATED)[[1]] = genes
dimnames(ALLELES_INACTIVATED)[[2]] = sample_ids
dimnames(INACTIVE)[[1]] = genes
dimnames(INACTIVE)[[2]] = c("bi","mono", "none")
dimnames(matrix_SV)[[1]] = genes
dimnames(matrix_SV)[[2]] = sample_ids

# samples with tandem duplication hotspots identified in figure 2
matrix_SV["AR", assess_tandem_dup( chrom="chrX", locus=6692, 
                                   list_sv_m, chrom_lengths, sample_ids ) ] = 'tandem'
matrix_SV["FOXA1", assess_tandem_dup( chrom="chr14", locus=3754, 
                                      list_sv_m, chrom_lengths, sample_ids) ] = 'tandem'
matrix_SV["MYC", assess_tandem_dup( chrom="chr8", locus=12701, 
                                    list_sv_m, chrom_lengths, sample_ids ) ] = 'tandem'
matrix_SV["MYC", assess_tandem_dup( chrom="chr8", locus=12741, 
                                    list_sv_m, chrom_lengths, sample_ids ) ] = 'tandem'

# samples with curated SV
for(idx in 1:dim(curated_sv)[1]){
    sample_id = curated_sv$sample[ idx ]
    symbol = curated_sv$threeprime[ idx ]
    idx_in_samples = which( sample_ids==sample_id ) 
    consequence = curated_sv$consequence[ idx ]
    mechanism = curated_sv$mechanism[ idx ]
    matrix_SV[ symbol, sample_id ] = mechanism
}

matrix_CNA["CDK12", "DTB-183-BL"] = "LOH"

# Allele effects on expression

cors = rep(NA, dim(gene_locs)[1])
pvals = rep(NA, dim(gene_locs)[1])
n_0 = rep(NA, dim(gene_locs)[1])
n_1 = rep(NA, dim(gene_locs)[1])
n_2 = rep(NA, dim(gene_locs)[1])

s5 = read.table(fn_symbols_fig5,
                header=TRUE, sep='\t', stringsAsFactors = FALSE, row.names=1)
symbols = rownames(s5)
for(i in 1:length(symbols)){
    if( i%%100==0){
        print(i)
    }
    symbol=symbols[i]
    if( symbol != "AR enhancer"){
        ae=allele_effect( symbol, do_plot=FALSE, axis_override=NULL )
        INACTIVE[symbol, 1] = ae$n_bi
        INACTIVE[symbol, 2] = ae$n_mono
        INACTIVE[symbol, 3] = N_SAMPLES-ae$n_bi-ae$n_mono
        ALLELES_INACTIVATED[symbol,] = ae$alleles$n_alleles_inactivated
        idx = which(rownames(gene_locs)==symbols[i])
        cors[idx] = ae$cor_alleles_xp
        pvals[idx] = ae$pval_alleles_xp
        n_0[idx] = sum(ae$alleles$n_alleles_inactivated==2)
        n_1[idx] = sum(ae$alleles$n_alleles_inactivated==1)
        n_2[idx] = sum(ae$alleles$n_alleles_inactivated==0)
    }
}

matrix_samples$has_HRD = allele_effect("BRCA2")$alleles$n_alleles_inactivated==2
matrix_samples$has_ETS = allele_effect("ERG")$alleles$activating_sv |
                         allele_effect("ETV1")$alleles$activating_sv |
                         allele_effect("ETV4")$alleles$activating_sv |
                         allele_effect("ETV5")$alleles$activating_sv 


non_sex = cnbin$chrom!="chrX" & cnbin$chrom!="chrY"
percent_CNA = round( 
    (colSums(cnbin[non_sex,]>=GAIN_NONSEX|cnbin[non_sex,]<=LOSS_SINGLE_NONSEX, na.rm=TRUE) + 
    colSums(cnbin[!non_sex,]>=GAIN_SEX|cnbin[!non_sex,]<=GAIN_SEX,na.rm=TRUE) ) / 1042 , 2)[3:103]
matrix_samples$percent_CNA = percent_CNA

