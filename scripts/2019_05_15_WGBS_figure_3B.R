#source('/Volumes/datasets_1/human_sequence_prostate_WGBS/reproduce/scripts/2019_05_15_prepare_environment.r')

# Plots figure 3B

p1=plot_gene_methylation_across_samples("TMEFF2") + theme(legend.position = "none") 
p2=plot_gene_methylation_across_samples("PCAT14") + theme(legend.position = "none") 
p3=plot_gene_methylation_across_samples("SCHLAP1") + theme(legend.position = "none") 

p4=plot_gene_methylation_across_samples("SLC45A3") + theme(legend.position = "none") 
p5=plot_gene_methylation_across_samples("SPON2") + theme(legend.position = "none") 
p6=plot_gene_methylation_across_samples("TDRD1") + theme(legend.position = "none") 

ggarrange(plotlist=list( p1, p2, p3, p4, p5, p6),
          ncol = 3, nrow=2)
ggsave(fn_figure3b,width=8,height=4,useDingbats=FALSE)