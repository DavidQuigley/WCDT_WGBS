#source('/Volumes/datasets_1/human_sequence_prostate_WGBS/reproduce/scripts/2019_05_15_prepare_environment.r')

# Plots figure 4B
p1=plot_gene_methylation_across_samples("KLK3") + theme(legend.position = "none") 
p2=plot_gene_methylation_across_samples("NKX3-1") + theme(legend.position = "none") 
p3=plot_gene_methylation_across_samples("FOLH1") + theme(legend.position = "none") 

ggarrange(plotlist=list( p1, p2, p3 ),
          ncol = 1, nrow=3)
ggsave(fn_figure4b,width=8,height=4,useDingbats=FALSE)