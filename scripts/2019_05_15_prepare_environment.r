library(data.table)
library(forcats)
library(GenomicRanges)
library(ggpubr)
library(ggplot2)
library(gplots)
library(grid)
library(ggpubr)
library(plotrix)
library(plyr)
library(matrixStats)
library(reshape2)
library(rCGH) # Load hg38 chromosomal info
library(RColorBrewer)
library(Rtsne)
library(rtracklayer)
library(scales)
library(stringr)
library(Sushi)
library(wesanderson)

#BASE_DIR = '~/mnt/data1/projects/WCDT_WGBS_2019/WCDT_WGBS'
BASE_DIR = '/notebook/human_sequence_prostate_WGBS/reproduce//WCDT_WGBS'
make_path = function(BASE_DIR, s){
    paste(BASE_DIR, s, sep='/' )
}

#######################################################################################################
# Load all custom function calls
#######################################################################################################

source( make_path(BASE_DIR, 'scripts/2019_05_15_function_calls.R') )

#######################################################################################################
# File locations
#######################################################################################################

#dir_pmd = make_path( BASE_DIR, 'secondary_data/partially_methylated_domains/input/10k')
fn_pmd_mu_summary = make_path( BASE_DIR, 'secondary_data/partially_methylated_domains/output/pmd_mu_summary.txt')

dir_normal_wide_slice = make_path( BASE_DIR, 'secondary_data/methyl_slices_around_features/wide_slice/UPitt')
dir_mcrpc_wide_slice = make_path( BASE_DIR, 'secondary_data/methyl_slices_around_features/wide_slice/mCRPC')
dir_copycat = make_path( BASE_DIR, 'secondary_data/copycat' )
fn_tpm_wgbs = make_path(BASE_DIR, 'secondary_data/RNAseq/featurecounts_tpm_gene_gencodev28.tsv' )
fn_ensembl2sym = make_path(BASE_DIR, 'metadata/genes/ensembl2sym.txt' )
fn_genelists = make_path(BASE_DIR, 'metadata/genes/COSMIC_6_19_18.tsv' )
fn_annovar = make_path(BASE_DIR, 'secondary_data/somatic_mutations/annovar.hg38.txt' )
fn_annovar_full = make_path(BASE_DIR, 'secondary_data/somatic_mutations/annovar.hg38_full.txt' )
fns_chipseq = c(
    make_path(BASE_DIR, 'secondary_data/chipseq/crpc_are_sum_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chipseq/localized_are_sum_hg38.bed'  ),
    make_path(BASE_DIR, 'secondary_data/chipseq/pca100/AR_sum_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chipseq/pca100/H3K4me3_sum_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chipseq/pca100/H3K27ac_sum_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chipseq/pca100/H3K27me3_count_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chipseq/FOXA1_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chipseq/HOXB13_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chipseq/erg_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chiapet/AR_merged_hg38.bed' ),
    make_path(BASE_DIR, 'secondary_data/chiapet/ERG_merged_hg38.bed' )
)

fn_cpg = make_path(BASE_DIR, 'metadata/ucsc_cpg_island.txt' )
fn_DSSS_tscnc_adeno = make_path(BASE_DIR, 'secondary_data/DSS/DSS_tscnc_adeno.tsv' )
fn_DSS_localized_benign = make_path(BASE_DIR, 'secondary_data/DSS/DSS_localized_benign.tsv' )
fn_DSS_adeno_localized = make_path(BASE_DIR, 'secondary_data/DSS/DSS_metastatic_localized_SRA.tsv')
fn_DSS_adeno_benign = make_path(BASE_DIR, "secondary_data/DSS/DSS_mAdeno_benignprostate.tsv")

fn_hmrseg_gene = make_path(BASE_DIR, 'secondary_data/recurrent_HMRs/hmrseg_gene.txt' )
fn_hmrseg = make_path(BASE_DIR, 'secondary_data/recurrent_HMRs/hmrseg.txt' )
fn_hmr_recurrent = make_path(BASE_DIR, 'secondary_data/recurrent_HMRs/hmr_recurrent.txt' )

fn_recurrent_hmr_mcrpc = make_path(BASE_DIR, 'secondary_data/recurrent_HMRs/rhmr.tsv' )
fn_recurrent_hmr_nt = make_path(BASE_DIR, 'secondary_data/recurrent_HMRs/rhmr_nt.tsv' )
fn_recurrent_hmr_upitt = make_path(BASE_DIR, 'secondary_data/recurrent_HMRs/rhmr_upitt.tsv' )

dir_methylseqR = make_path(BASE_DIR,'secondary_data/MethylSeekR/')
fns_methylseqR = list.files(dir_methylseqR, pattern="UMRsLMRs.txt", full.names=F)

fn_treatment = make_path(BASE_DIR, 'metadata/clinical/treatment.txt' )
fn_array_27k = make_path(BASE_DIR, 'secondary_data/illumina_methylation_array_coverage/hm27.hg38.slim.tsv' )
fn_array_450k = make_path(BASE_DIR, 'secondary_data/illumina_methylation_array_coverage/hm450.hg38.slim.tsv' )
fn_array_epic = make_path(BASE_DIR, 'secondary_data/illumina_methylation_array_coverage/EPIC.hg38.slim.tsv' )

fn_wgs_sample_summary = make_path(BASE_DIR, 'metadata/WGS_sample_summary.txt' )
fn_tumor_purity = make_path(BASE_DIR, 'metadata/WGS_tumor_purity_AF_012318.csv' )
fn_manta_slim = make_path(BASE_DIR, 'metadata/WGS_vcf_manta_slim.txt')
fn_genelist_housekeeping = make_path(BASE_DIR,'metadata/genes/housekeeping_genes.txt')

fn_hmr_summary_expr = make_path(BASE_DIR, 'secondary_data/recurrent_HMRs/summary_expr.txt')

fn_pathways = make_path(BASE_DIR, "metadata/genes/h.all.v6.2.symbols.gmt")
dir_gene_output_mcrpc_narrow = make_path(BASE_DIR, 'secondary_data/gene_output/gene_output')
dir_gene_output_mcrpc = make_path(BASE_DIR, 'secondary_data/gene_output/gene_selected_output')
dir_gene_output_NT = make_path(BASE_DIR, 'secondary_data/gene_output/gene_output_NT')
dir_gene_output_upitt = make_path(BASE_DIR, 'secondary_data/gene_output/gene_output_upitt')

fn_expr_covar_sum = make_path(BASE_DIR, 'secondary_data/recurrent_HMRs/summary_expr.txt')

fn_pca_specific = make_path(BASE_DIR, 'metadata/genes/pca_specific.tsv')

fn_lift_38_to_19 = make_path(BASE_DIR, 'metadata/genome_builds/hg38ToHg19.over.chain')
fn_lift_19_to_38 = make_path(BASE_DIR, 'metadata/genome_builds/hg19ToHg38.over.chain')

fn_gencode_exons = make_path(BASE_DIR, 'metadata/genome_builds/gencode28_exons.txt')
dir_ar_enhancer_bw = make_path( BASE_DIR, 'secondary_data/chipseq/H3K27ac_chrX' )
fn_h3k27ac_avg = make_path(BASE_DIR, "secondary_data/chipseq/H3K27ac_cellline/h3k27ac_avg.bw" )
fn_chipseq_hoxb13 = make_path(BASE_DIR, "secondary_data/chipseq/HOXB13_hg38.bed")
fn_chipseq_foxa1 = make_path(BASE_DIR, "secondary_data/chipseq/FOXA1_hg38.bed")
fn_chipseq_erg_hg38 = make_path(BASE_DIR, "secondary_data/chipseq/erg_hg38.bed")
fn_chipseq_erg_vcap = make_path(BASE_DIR, "secondary_data/chipseq/erg_vcap_hg38.bed")
fn_crpc_are = make_path(BASE_DIR, 'secondary_data/chipseq/crpc_are_sum_hg38.bed')
fn_loc_are = make_path(BASE_DIR, 'secondary_data/chipseq/localized_are_sum_hg38.bed')
fn_pc100_are = make_path(BASE_DIR, 'secondary_data/chipseq/pca100/AR_sum_hg38.bed')
fn_pc100_h3k4me3 = make_path(BASE_DIR, 'secondary_data/chipseq/pca100/H3K4me3_sum_hg38.bed')
fn_pc100_h3k27ac = make_path(BASE_DIR, 'secondary_data/chipseq/pca100/H3K27ac_sum_hg38.bed')

fn_expressed_genes = make_path(BASE_DIR, 'metadata/genes/expressed_genes.txt')
fn_chia_ar_bind = make_path(BASE_DIR, 'secondary_data/chiapet/AR_merged_hg38.bed' )
fn_chiapet_ar = make_path( BASE_DIR, "secondary_data/chiapet/ar_chiapet_hg38.txt")
fn_chiapet_erg = make_path(BASE_DIR, 'secondary_data/chiapet/erg_chiapet_hg38.txt')
fn_chiapet_erg_bind = make_path(BASE_DIR, 'secondary_data/chiapet/ERG_merged_hg38.bed')

fn_hic_lncap = make_path(BASE_DIR, 'secondary_data/HiC/hic_lncap_hg38.txt' )
fn_hic_pc3 = make_path(BASE_DIR, 'secondary_data/HiC/hic_pc3_hg38.txt' )
fn_hic_prec = make_path(BASE_DIR, 'secondary_data/HiC/hic_prec_hg38.txt' )

fn_valley = make_path(BASE_DIR, 'secondary_data/valleys/valley_5kb.bed')
fn_centromeres = make_path(BASE_DIR, 'metadata/grch38_centromeres.txt')
fn_telomeres = make_path(BASE_DIR, 'metadata/grch_telomeres.txt')
fn_gap = make_path(BASE_DIR, 'metadata/grch38_gap.txt')
fn_surv = make_path(BASE_DIR, 'metadata/clinical/dt.101.full.csv')

fn_spline_foxa1 = make_path(BASE_DIR, 'secondary_data/tfbs_splines/dssmat.tfbs.foxa1.RData')
fn_spline_hoxb13 = make_path(BASE_DIR, 'secondary_data/tfbs_splines/dssmat.tfbs.foxa1.RData')
fn_spline_ar = make_path(BASE_DIR, 'secondary_data/tfbs_splines/dssmat.tfbs.foxa1.RData')
fn_spline_erg = make_path(BASE_DIR, 'secondary_data/tfbs_splines/dssmat.tfbs.foxa1.RData')
fn_spline_h3k27ac = make_path(BASE_DIR, 'secondary_data/tfbs_splines/dssmat.tfbs.h3k27ac.RData')

fn_figure6c_foxa1 = make_path(BASE_DIR, 'figures/figure6c_foxa1.pdf')
fn_figure6c_hoxb13 = make_path(BASE_DIR, 'figures/figure6c_hoxb13.pdf')
fn_figure6c_ar = make_path(BASE_DIR, 'figures/figure6c_ar.pdf')
fn_figure6c_erg = make_path(BASE_DIR, 'figures/figure6c_erg.pdf')
fn_figure6c_h3k27ac = make_path(BASE_DIR, 'figures/figure6c_h3k27ac.pdf')

fn_fig1a = make_path(BASE_DIR, 'figures/figure_1a.pdf')
fn_figure1b = make_path(BASE_DIR, 'figures/figure_1b.pdf')
fn_figure1c = make_path(BASE_DIR, 'figures/figure_1c.pdf')
fn_figure1d = make_path(BASE_DIR, 'figures/figure_1d.pdf')
fn_figure1e = make_path(BASE_DIR, 'figures/figure_1e.pdf')
fn_figure1f = make_path(BASE_DIR, 'figures/figure_1f.pdf')
fn_figure2a = make_path(BASE_DIR, 'figures/figure_2a.pdf')
fn_figure2b = make_path(BASE_DIR, 'figures/figure_2b.pdf')
fn_figure3a = make_path(BASE_DIR, 'figures/figure_3a.pdf')
fn_figure3b = make_path(BASE_DIR, 'figures/figure_3b.pdf')
fn_figure4a = make_path(BASE_DIR, 'figures/figure_4a.pdf')
fn_figure4b = make_path(BASE_DIR, 'figures/figure_4b.pdf')
fn_figure4c = make_path(BASE_DIR, 'figures/figure_4c.pdf')
fn_figure5a = make_path(BASE_DIR, 'figures/figure_5a.pdf')
fn_figure5b = make_path(BASE_DIR, 'figures/figure_5b.pdf')
fn_figure5c = make_path(BASE_DIR, 'figures/figure_5c.pdf')
fn_figure5d = make_path(BASE_DIR, 'figures/figure_5d.pdf')
fn_figure6a = make_path(BASE_DIR, 'figures/figure_6a.png')
fn_figure6a_colorbar = make_path(BASE_DIR, 'figures/figure_6a_colorbar.pdf')
fn_figure6b = make_path(BASE_DIR, 'figures/figure_6b.pdf')
fn_figure6d_h3k27ac
fn_figure_pmd_chr16 = make_path(BASE_DIR, 'figures/figure_pmd_chr16.pdf')
fn_figure_pmd_chrom_summary = make_path(BASE_DIR, 'figures/figure_pmd_chromosome_summary.pdf')
fn_figure_tsne_a= make_path(BASE_DIR, 'figures/supplementary_figure_tsne_a.pdf')

fn_figure_neg_ehmr_cor_allgenes = make_path(BASE_DIR, 'figures/figure_neg_ehmr_cor_allgenes.pdf')
fn_figure_neg_ehmr_cordensity_allgenes = make_path(BASE_DIR, 'figures/figure_neg_ehmr_cordensity_allgenes.pdf')
fn_figure_pos_ehmr_cor_allgenes = make_path(BASE_DIR, 'figures/figure_pos_ehmr_cor_allgenes.pdf')
fn_figure_pos_ehmr_cordensity_allgenes = make_path(BASE_DIR, 'figures/figure_pos_ehmr_cordensity_allgenes.pdf')
fn_figure_surv = make_path(BASE_DIR, 'figures/figure_km.pdf')

fn_table_ehmr = make_path(BASE_DIR, 'tables/table_eHMR.txt')

#######################################################################################################
# Hardcoded variables
#######################################################################################################

# Plot colors

colors_GOF <- c(
    'purple', #activating_missense
    'orange', #activating_sv
    'green') #CNA_amp
names(colors_GOF) = c(
    'activating_missense',
    'activating_sv',
    'CNA_amp')
colors_LOF <- c(
    'purple', #inactivating_missense
    'purple', #nonsense
    'purple', #inactivating_germline
    'orange',#inactivating_sv
    'cyan4', #CNA_1
    'cyan4') #CNA_2
names(colors_LOF) = c(
    'inactivating_missense',
    'nonsense',
    'inactivating_germline',
    'inactivating_sv',
    'CNA_1',
    'CNA_2')


# Copy gain/loss hardcoded limits
#------------------------------------------------------------------------------------------------------
GAIN_NONSEX = 3
LOSS_SINGLE_NONSEX = 1.65
LOSS_DOUBLE_NONSEX = 0.6
GAIN_SEX = 1.4
LOSS_SEX = 0.6

# hg38 chromosome dimensions
#------------------------------------------------------------------------------------------------------
chrsom = c(paste('chr', 1:22, sep=''))
chrsex = c('chrX', 'chrY')
chromosomes = c(paste('chr', 1:22, sep=''), 'chrX', 'chrY')
names(chromosomes) = c("NC_000001.11","NC_000002.12","NC_000003.12","NC_000004.12",
                       "NC_000005.10","NC_000006.12","NC_000007.14","NC_000008.11","NC_000009.12",
                       "NC_000010.11","NC_000011.10","NC_000012.12","NC_000013.11","NC_000014.9",
                       "NC_000015.10","NC_000016.10","NC_000017.11","NC_000018.10","NC_000019.10",
                       "NC_000020.11","NC_000021.9","NC_000022.11","NC_000023.11","NC_000024.10")
chrsizes = c(248956422,242193529,198295559,190214555,181538259,
             170805979,159345973,145138636,138394717,133797422,135086622,
             133275309,114364328,107043718,101991189,90338345,83257441,
             80373285,58617616,64444167,46709983,50818468,156040895,57227415)
names(chrsizes) = chromosomes

# Genomic range limits
#------------------------------------------------------------------------------------------------------
PROMOTER_RANGE = 1500
RANGE_UP_DOWN = 10000
HG38_SIZE = 3099734149

# physical bounds of genes
ensembl2sym = read.table(fn_ensembl2sym, header=TRUE, row.names=1, 
                         stringsAsFactors=FALSE, strip.white=TRUE)
gencode_protein = c('protein_coding','TEC')
gencode_pseudogene = c('transcribed_processed_pseudogene','transcribed_unitary_pseudogene',
                       'transcribed_unprocessed_pseudogene','translated_processed_pseudogene','unitary_pseudogene',
                       'unprocessed_pseudogene','TR_V_pseudogene','TR_J_pseudogene','IG_C_pseudogene',
                       'IG_J_pseudogene','IG_pseudogene','IG_V_pseudogene','polymorphic_pseudogene',
                       'processed_pseudogene','pseudogene')
gencode_immune = c('TR_V_gene','TR_J_gene','TR_D_gene','TR_C_gene',
                   'IG_C_gene','IG_D_gene','IG_J_gene','IG_V_gene')
gencode_noncoding = c('vaultRNA','3prime_overlapping_ncRNA','antisense',
                      'bidirectional_promoter_lncRNA','processed_transcript',
                      'non_coding','lincRNA','macro_lncRNA','misc_RNA','sense_intronic',
                      'sense_overlapping','miRNA','ribozyme','rRNA','scaRNA','scRNA',
                      'snoRNA','snRNA','sRNA','Mt_rRNA','Mt_tRNA')

# Patient treatment and PSA responses
#------------------------------------------------------------------------------------------------------
txresponse = read.table(fn_treatment, header=TRUE, sep='\t',
                        row.names=1, stringsAsFactors=FALSE, strip.white=TRUE)
rownames(txresponse) = convertClinSamplenames(rownames(txresponse))

# Locations of each metastatic tumor in the patient's body
#------------------------------------------------------------------------------------------------------
metastasis_locations = txresponse$Biopsy.site
names(metastasis_locations) = rownames(txresponse)
metastasis_locations[!(metastasis_locations %in% c('Bone','Liver','Lymph node'))] = 'Other'
metastasis_locations[metastasis_locations=='Lymph node'] = 'Lymph_node'
metastasis_locations['DTB-167-PRO'] = 'Bone'

# Sample identifiers
#------------------------------------------------------------------------------------------------------
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

sample_blacklist = c('DTB-004-BL', 'DTB-053-RP', 'DTB-265-BL')
samples_hypermut = c('DTB-083-BL','DTB-126-BL')

samples_sra_benign_prostate = c('SRR6156035','SRR6156036','SRR6156037')
samples_sra_localized_tumor = c('SRR6156038','SRR6156039','SRR6156040',
                                'SRR6156041','SRR6156042','SRR6156043','SRR6156044',
                                'SRR6156045','SRR6156046','SRR6156047','SRR6156048')
#samples_upitt_benign_prostate = toupper(samples_upitt_benign_prostate)
#samples_upitt_localized_tumor = toupper(samples_upitt_localized_tumor)

samples_normal = c('NT-DTB-003-BL','NT-DTB-021-BL','NT-DTB-035-BL','NT-DTB-037-BL',
                   'NT-DTB-074-BL','NT-DTB-077-Pro','NT-DTB-104-BL','NT-DTB-124-BL', 
                   'NT-DTB-188-BL','NT-DTB-194-Pro')
samples_normal_filecase = samples_normal
samples_normal = toupper(samples_normal)
locs_normal = metastasis_locations[gsub('NT-','',samples_normal)]
names(locs_normal) = paste('NT-',names(locs_normal),sep='')

samples_upitt_benign_prostate = c('Bis49AT',          'Bis159AT','Bis165AT','Bis171AT')
samples_upitt_localized_tumor <- c('Bis49T','Bis158T','Bis159T', 'Bis165T', 'Bis171T')
samples_upitt_benign_prostate <- toupper(samples_upitt_benign_prostate)
#samples_upitt_localized_tumor <- toupper(samples_upitt_localized_tumor)

samples_braf = "DTB-022-BL"
samples_idh1 = "DTB-112-BL"
samples_tet2 <- c("DTB-023-BL", "DTB-188-BL", "DTB-202-BL", "DTB-252-BL") #252 is the one out of CMP
samples_dnmt3b <- c("DTB-024-PRO", "DTB-094-BL")


samples_tscnc = c('DTB-003-BL','DTB-036-BL','DTB-040-BL','DTB-135-PRO','DTB-205-BL')

samples_cmp = c(
    "DTB-018-BL",
    "DTB-021-BL",
    "DTB-022-BL",  #BRAF
    "DTB-023-BL",  #TET2
    "DTB-024-PRO", #DNMT3B
    "DTB-037-BL",
    "DTB-042-BL",
    "DTB-074-BL",
    "DTB-089-BL",
    "DTB-094-BL",  #DNMT3B
    "DTB-102-PRO",
    "DTB-112-BL",  #IDH1
    "DTB-124-BL",
    "DTB-137-PRO",
    "DTB-141-BL",
    "DTB-151-BL",
    "DTB-188-BL",  #TET2
    "DTB-190-BL",
    "DTB-194-PRO",
    "DTB-202-BL",# TET2
    "DTB-204-BL",
    "DTB-260-BL")

# Gene lists
#------------------------------------------------------------------------------------------------------
cosmic = read.table( fn_genelists, sep='\t', header=TRUE, stringsAsFactors=FALSE, strip.white=TRUE)
cosmic[cosmic$Gene.Symbol=='C2orf44','Gene.Symbol'] = 'WDCP'
cosmic[cosmic$Gene.Symbol=='CASC5','Gene.Symbol'] = 'KNL1'
cosmic[cosmic$Gene.Symbol=='KIAA1598','Gene.Symbol'] = 'SHTN1'
cosmic[cosmic$Gene.Symbol=='LHFP','Gene.Symbol'] = 'LHFPL6'
cosmic[cosmic$Gene.Symbol=='MLLT4','Gene.Symbol'] = 'AFDN'
cosmic[cosmic$Gene.Symbol=='RUNDC2A','Gene.Symbol'] = 'SNX29'
cosmic[cosmic$Gene.Symbol=='WHSC1','Gene.Symbol'] = 'NSD2'
cosmic[cosmic$Gene.Symbol=='WHSC1L1','Gene.Symbol'] = 'NSD3'
cosmic[cosmic$Gene.Symbol=='ZNF198','Gene.Symbol'] = 'ZMYM2'
cosmic[cosmic$Gene.Symbol=='ZNF278','Gene.Symbol'] = 'PATZ1'
genelist_cosmic = cosmic$Gene.Symbol
genelist_oncogenes = cosmic[grepl('oncogene',cosmic$Role.in.Cancer,fixed=TRUE),'Gene.Symbol']
genelist_tumor_suppressors = cosmic[grepl('TSG',cosmic$Role.in.Cancer,fixed=TRUE),'Gene.Symbol']
genelist_housekeeping = read.delim(fn_genelist_housekeeping, sep='\t',
                                   header=FALSE,stringsAsFactors=FALSE,strip.white=TRUE)[,1]

pca_celldf = rbind.data.frame(
    cbind.data.frame(gene='ERG',gof=TRUE),
    cbind.data.frame(gene='ETV1',gof=TRUE),
    cbind.data.frame(gene='ETV4',gof=TRUE),
    cbind.data.frame(gene='ETV5',gof=TRUE),
    cbind.data.frame(gene='BRAF',gof=TRUE),
    cbind.data.frame(gene='HRAS',gof=TRUE),
    cbind.data.frame(gene='MYC',gof=TRUE),
    cbind.data.frame(gene='CHD1',gof=FALSE),
    cbind.data.frame(gene='SPOP',gof=FALSE),
    cbind.data.frame(gene='IDH1',gof=FALSE),
    cbind.data.frame(gene='AR',gof=TRUE),
    cbind.data.frame(gene='FOXA1',gof=TRUE),
    cbind.data.frame(gene='NCOR1',gof=FALSE),
    cbind.data.frame(gene='NCOR2',gof=FALSE),
    cbind.data.frame(gene='ASXL2',gof=FALSE),
    cbind.data.frame(gene='ZBTB16',gof=FALSE),
    cbind.data.frame(gene='TP53',gof=FALSE),
    cbind.data.frame(gene='PTEN',gof=FALSE),
    cbind.data.frame(gene='AKT1',gof=TRUE),
    cbind.data.frame(gene='PIK3CA',gof=TRUE),
    cbind.data.frame(gene='RB1',gof=FALSE),
    cbind.data.frame(gene='CDKN2A',gof=FALSE),
    cbind.data.frame(gene='CCND1',gof=TRUE),
    cbind.data.frame(gene='APC',gof=FALSE),
    cbind.data.frame(gene='CTNNB1',gof=TRUE),
    cbind.data.frame(gene='ZNRF3',gof=FALSE),
    cbind.data.frame(gene='KMT2C',gof=FALSE),
    cbind.data.frame(gene='KMT2D',gof=FALSE),
    cbind.data.frame(gene='KDM6A',gof=FALSE),
    cbind.data.frame(gene='HDAC4',gof=FALSE),
    cbind.data.frame(gene='NKX3-1',gof=FALSE),
    cbind.data.frame(gene='AXL',gof=TRUE),
    cbind.data.frame(gene='MED12',gof=TRUE),
    cbind.data.frame(gene='ZFHX3',gof=FALSE),
    cbind.data.frame(gene='GNAS',gof=FALSE),
    cbind.data.frame(gene='BRCA2',gof=FALSE),
    cbind.data.frame(gene='BRCA1',gof=FALSE),
    cbind.data.frame(gene='ATM',gof=FALSE),
    cbind.data.frame(gene='CDK12',gof=FALSE),
    cbind.data.frame(gene='ERCC2',gof=FALSE),
    cbind.data.frame(gene='PRKDC',gof=FALSE),
    cbind.data.frame(gene='MLH1',gof=FALSE),
    cbind.data.frame(gene='MSH2',gof=FALSE),
    cbind.data.frame(gene='MSH6',gof=FALSE),
    cbind.data.frame(gene='PIK3CA',gof=TRUE),
    cbind.data.frame(gene='MAP2K4',gof=TRUE),
    cbind.data.frame(gene='SCHLAP1',gof=TRUE))
pca_celldf$gene = paste(pca_celldf$gene)
genelist_pca_cell_genes = pca_celldf$gene

data.pathways = lapply(scan(fn_pathways, "character", sep = '\n'), function(line){
    line = strsplit(line, '\t')[[1]]
    line = line[-2]
    line = unique(line)
    return(line)
})
for(i in 1:length(data.pathways)) {
    names(data.pathways)[i] <- data.pathways[[i]][1]
    data.pathways[[i]] <- data.pathways[[i]][-1]
}
names(data.pathways) <- sub("HALLMARK_", "", names(data.pathways))
data.pathways[['ADIPOGENESIS']] <- NULL
data.pathways[['MYOGENESIS']] <- NULL
data.pathways[['PROTEIN_SECRETION']] <- NULL
data.pathways[['UV_RESPONSE_UP']] <- NULL
data.pathways[['UV_RESPONSE_DN']] <- NULL
data.pathways[['HEME_METABOLISM']] <- NULL
data.pathways[['COAGULATION']] <- NULL
data.pathways[['BILE_ACID_METABOLISM']] <- NULL
data.pathways[['PEROXISOME']] <- NULL
data.pathways[['ALLOGRAFT_REJECTION']] <- NULL
data.pathways[['SPERMATOGENESIS']] <- NULL
data.pathways[['PANCREAS_BETA_CELLS']] <- NULL

#######################################################################################################
# Load WGS variables
#######################################################################################################

source( make_path(BASE_DIR, 'scripts/2019_05_15_prepare_WGS.R') )

#######################################################################################################
# Load secondary (calculated) data
#######################################################################################################

# properties of the WGS sequences
#------------------------------------------------------------------------------------------------------------------
cellsummary = read.delim(fn_wgs_sample_summary, row.names=1, header=TRUE, sep='\t',
                         stringsAsFactors=FALSE, strip.white=TRUE)
tumorpurity = read.csv(fn_tumor_purity, row.names=2)

# RNAseq data load
#------------------------------------------------------------------------------------------------------------------
tpm = read.table(fn_tpm_wgbs, header=TRUE, row.names=1, stringsAsFactors=FALSE, 
                 strip.white=TRUE, check.names=F)
tpm = modify_RNAseq_names(tpm)

# Manta Structural Variation calls
#------------------------------------------------------------------------------------------------------------------
manta_counts = matrix_sv_m
names(manta_counts) = c("DEL", "DUP", "n_mh", "INS", "INV", "per", "BND")
manta_counts$SV = rowSums( manta_counts[, c("DEL", "DUP", "INS","INV", "BND")] )
manta_counts = manta_counts[,c(8,7,1,2,4,5)]


# Copy number data load
#------------------------------------------------------------------------------------------------------------------
copycat_bed = NULL
filenames_cc = list.files(dir_copycat)
for(i in 1:length(filenames_cc)) {
    filename_bed = paste(dir_copycat, filenames_cc[i], sep='/')
    data_i = read.table(filename_bed, header=FALSE, stringsAsFactors=FALSE, sep='\t')
    sample_i = gsub('_copycat.bed','',filenames_cc[i],fixed=TRUE)
    data_i$sample_id = sample_i
    copycat_bed = rbind.data.frame(copycat_bed, data_i)
}
colnames(copycat_bed) = c('chrom','start','end', 'call','cn','strand','sample_id')
copycat_bed$call = 'REF'
copycat_bed[copycat_bed$chr %in% chrsom & copycat_bed$cn >= GAIN_NONSEX,'call'] = 'GAIN'
copycat_bed[copycat_bed$chr %in% chrsex & copycat_bed$cn >= GAIN_SEX,'call'] = 'GAIN'
copycat_bed[copycat_bed$chr %in% chrsom & copycat_bed$cn <= LOSS_DOUBLE_NONSEX,'call'] = 'LOSS'
copycat_bed[copycat_bed$chr %in% chrsex & copycat_bed$cn <= LOSS_SEX,'call'] = 'LOSS'

copycat_gr <- makeGRangesFromDataFrame(copycat_bed[,-6],keep.extra.columns=TRUE)
copycat_gr <- copycat_gr[copycat_gr$call != 'REF']

# Annovar mutation annotations
#------------------------------------------------------------------------------------------------------------------
annovar = read.delim(fn_annovar, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=TRUE)
annovar = annotate_annovar(annovar)
annovar_full = read.delim( fn_annovar_full, 
                           header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=TRUE)
annovar_full = annotate_annovar(annovar_full)

# Chipseq track load
#------------------------------------------------------------------------------------------------------------------

colnamesbed = c('chrom','start','end','score')
namelist_chipseq = c(
    'CRPC_ARE',
    'LOC_ARE',
    'pca100_ARE',
    'pca100_H3K4me3',
    'pca100_H3K27ac',
    'pca100_H3K27me3',
    'chip_FOXA1',
    'chip_HOXB13',
    'chip_ERG',
    'chiapet_AR',
    'chiapet_ERG'
)
scorecut = c(0,0,0,0,0,0,0,0,0,0,0)
stopifnot(length(fns_chipseq)==length(namelist_chipseq))
stopifnot(length(scorecut)==length(namelist_chipseq))

tracks = list()
for(i in 1:length(namelist_chipseq)) {
    filename = fns_chipseq[i]
    trackname = namelist_chipseq[i]
    track = read.delim(filename, header=FALSE, sep='\t', stringsAsFactors=FALSE, strip.white=TRUE)
    colnames(track) = colnamesbed
    rows2keep = track$score >= scorecut[i]
    stopifnot(sum(rows2keep)>0)
    tracks[[trackname]] = makeGRangesFromDataFrame(track[rows2keep,])
}

# GENE REGIONS
#------------------------------------------------------------------------------------------------------------------
rowspos = ensembl2sym$strand=='+'
rowsneg = ensembl2sym$strand=='-'

promoterpos = cbind.data.frame(ensembl2sym$chr,
                               ensembl2sym$start-PROMOTER_RANGE,
                               ensembl2sym$start+PROMOTER_RANGE,
                               ensembl2sym$strand,ensembl2sym$name)[rowspos,]
promoterneg = cbind.data.frame(ensembl2sym$chr,
                               ensembl2sym$end-PROMOTER_RANGE,
                               ensembl2sym$end+PROMOTER_RANGE,
                               ensembl2sym$strand,ensembl2sym$name)[rowsneg,]
colnames(promoterpos) = c('chr','start','end','strand','gene')
colnames(promoterneg) = c('chr','start','end','strand','gene')
promoterdf = rbind.data.frame(promoterpos,promoterneg)
tracks[['promoter']] = makeGRangesFromDataFrame(promoterdf,keep.extra.columns=TRUE)

genebody = cbind.data.frame(ensembl2sym$chr,
                            ensembl2sym$start,
                            ensembl2sym$end,
                            ensembl2sym$strand,
                            ensembl2sym$name)
colnames(genebody) = c('chr','start','end','strand','gene')
tracks[['genebody']] = makeGRangesFromDataFrame(genebody,keep.extra.columns=TRUE)

upstreampos = cbind.data.frame(ensembl2sym$chr,
                               ensembl2sym$start-RANGE_UP_DOWN,
                               ensembl2sym$start,
                               ensembl2sym$strand,ensembl2sym$name)[rowspos,]
upstreamneg = cbind.data.frame(ensembl2sym$chr,
                               ensembl2sym$end,
                               ensembl2sym$end+RANGE_UP_DOWN,
                               ensembl2sym$strand,ensembl2sym$name)[rowsneg,]
colnames(upstreampos) = c('chr','start','end','strand','gene')
colnames(upstreamneg) = c('chr','start','end','strand','gene')
upstreamdf = rbind.data.frame(upstreampos,upstreamneg)
tracks[['upstream']] = makeGRangesFromDataFrame(upstreamdf,keep.extra.columns=TRUE)

downstreampos = cbind.data.frame(ensembl2sym$chr,ensembl2sym$end,ensembl2sym$end+RANGE_UP_DOWN,
                                 ensembl2sym$strand,ensembl2sym$name)[rowspos,]
downstreamneg = cbind.data.frame(ensembl2sym$chr,ensembl2sym$start-RANGE_UP_DOWN,ensembl2sym$start,
                                 ensembl2sym$strand,ensembl2sym$name)[rowsneg,]
colnames(downstreampos) = c('chr','start','end','strand','gene')
colnames(downstreamneg) = c('chr','start','end','strand','gene')
downstreamdf = rbind.data.frame(downstreampos,downstreamneg)
tracks[['downstream']] = makeGRangesFromDataFrame(downstreamdf,keep.extra.columns=TRUE)

### CpG islands/shores/shelves

cpg = read.delim(fn_cpg, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=TRUE)
cpg = cpg[,c('chrom','chromStart','chromEnd','perCpg')]
colnames(cpg) = colnamesbed
tracks[['cpg_island']] = makeGRangesFromDataFrame(cpg)
cpg$start = cpg$start-2000
cpg$end = cpg$end+2000
tracks[['cpg_shore']] = makeGRangesFromDataFrame(cpg)
cpg$start = cpg$start-2000
cpg$end = cpg$end+2000
tracks[['cpg_shelf']] = makeGRangesFromDataFrame(cpg)

### DSS
#DSS tSCNC
tscncdf = read.delim(file=fn_DSSS_tscnc_adeno,sep='\t',header=TRUE,stringsAsFactors=FALSE)
tscncdf$chr = chromosomes[tscncdf$chr]
dss_tscnc = makeGRangesFromDataFrame(tscncdf,keep.extra.columns=TRUE)

#DSS localized vs benign
locdf = read.delim(file=fn_DSS_localized_benign,sep='\t',header=TRUE,stringsAsFactors=FALSE)
locdf$chr = chromosomes[locdf$chr]
dss_loc = makeGRangesFromDataFrame(locdf,keep.extra.columns=TRUE)

#DSS adeno vs localized
metdf = read.delim(file=fn_DSS_adeno_localized, sep='\t',header=TRUE,stringsAsFactors=FALSE)
metdf$chr = chromosomes[metdf$chr]
dss_met = makeGRangesFromDataFrame(metdf,keep.extra.columns=TRUE)

#DSS adeno vs benign
metdf = read.delim(file=fn_DSS_adeno_benign, sep='\t',header=TRUE,stringsAsFactors=FALSE)
metdf$chr = chromosomes[metdf$chr]
dss_metbp = makeGRangesFromDataFrame(metdf,keep.extra.columns=TRUE)

### HMR
hmrseggene = read.delim(fn_hmrseg_gene,header=TRUE,stringsAsFactors=FALSE,sep='\t',strip.white=TRUE)
tracks[['HMRseggene']] = makeGRangesFromDataFrame(hmrseggene,keep.extra.columns=TRUE)

hmrseg = read.delim(fn_hmrseg,header=TRUE,stringsAsFactors=FALSE,sep='\t',strip.white=TRUE)
tracks[['HMRseg']] = makeGRangesFromDataFrame(hmrseg,keep.extra.columns=TRUE)

for(i in 1:length(tracks)) {
    tracksnotin = tracks[[i]][countOverlaps(tracks[[i]],tracks[['promoter']])==0 & countOverlaps(tracks[[i]],tracks[['genebody']])==0]
    pct_rhmr = sum(countOverlaps(tracksnotin,tracks[['HMRseg']])>0)/length(tracksnotin)
    pct_cgi = sum(countOverlaps(tracksnotin,tracks[['cpg_shelf']])>0)/length(tracksnotin)
}

# centromeres and telomeres
centromeresdf <- read.delim(fn_centromeres, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=TRUE)
centromeres <- makeGRangesFromDataFrame(centromeresdf)
telomeresdf <- read.delim(fn_telomeres, header=TRUE, sep='\t', stringsAsFactors=FALSE, strip.white=TRUE)
telomeres <- makeGRangesFromDataFrame(telomeresdf)

# Load MethylseqR HMRs
#------------------------------------------------------------------------------------------------------------------
all_hmrs = list()
for(curfile in fns_methylseqR) {
    curfile_data = extract_methylseqR_metadata(curfile)
    curfile = paste(dir_methylseqR,curfile,sep='')
    curtab = read.delim(file=curfile,sep='\t',header=TRUE,stringsAsFactors=F)
    #all_hmrs[[curfile_data$sname]][[curfile_data$chr]] = curtab
    all_hmrs[[ toupper(curfile_data$sname)]][[curfile_data$chr]] = curtab
}
for(i in 1:length(all_hmrs)){
    if(length(all_hmrs[[i]])!=24) {
        missing = setdiff(chromosomes,names(all_hmrs[[i]]))
        print(paste(names(all_hmrs)[i],'missing',missing))
    }
}
# Merge all chromosomes into 1 table for each sample
final_hmr_list = lapply(all_hmrs, function(x) {
    rbindlist(x)
})
hmr = list()
hmrlengths = c()
hmrwidths = c()
for(i in 1:length(final_hmr_list)) {
    sample = names(final_hmr_list)[i]
    hmr[[sample]] = makeGRangesFromDataFrame(as.data.frame(final_hmr_list[[i]]),keep.extra.columns=TRUE)
    hmrlengths[sample] = length(hmr[[sample]])
    hmrwidths[sample] = sum(width(hmr[[sample]]))
}
names(hmr) = toupper(names(hmr))
names(hmrlengths) = toupper(names(hmrlengths))
names(hmrwidths) = toupper(names(hmrwidths))
hmr = hmr[toupper(c(sample_ids_wgbs,samples_normal,samples_upitt_benign_prostate,samples_upitt_localized_tumor))]

# Load methylation array data
#------------------------------------------------------------------------------------------------------------------

array27kdf = read.delim(fn_array_27k,header=TRUE,sep='\t')
array27k = makeGRangesFromDataFrame(array27kdf)
array450kdf = read.delim(fn_array_450k,header=TRUE,sep='\t')
array450k = makeGRangesFromDataFrame(array450kdf)
epicdf = read.delim(fn_array_epic,header=TRUE,sep='\t')
epic = makeGRangesFromDataFrame(epicdf)


# Calculate granges for oncogenes and tumor suppressors
#------------------------------------------------------------------------------------------------------------------
rows_tsg = ensembl2sym$name %in% setdiff(genelist_tumor_suppressors, genelist_oncogenes)
df_tsg = cbind.data.frame(ensembl2sym[rows_tsg,'chr'],
                          ensembl2sym[rows_tsg,'start'] - PROMOTER_RANGE,
                          ensembl2sym[rows_tsg,'end'])
colnames(df_tsg) = c('chr','start','end')
gr_tsg = makeGRangesFromDataFrame(df_tsg)

rows_oncogenes = ensembl2sym$name %in% setdiff(genelist_oncogenes, genelist_tumor_suppressors)
df_oncogenes = cbind.data.frame(ensembl2sym[rows_oncogenes,'chr'],
                                ensembl2sym[rows_oncogenes,'start'] - PROMOTER_RANGE,
                                ensembl2sym[rows_oncogenes,'end'])
colnames(df_oncogenes) = c('chr','start','end')
gr_oncogenes = makeGRangesFromDataFrame(df_oncogenes)
