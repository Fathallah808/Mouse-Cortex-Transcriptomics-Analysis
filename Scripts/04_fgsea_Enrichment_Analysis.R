#======================================
# Gene Set Enrichment Analysis (fgsea)
#======================================

#--------------------------------------
# Libraries
#--------------------------------------
library(fgsea)
library(msigdbr)
library(tidyverse)
library(ggplot2)
library(patchwork)

set.seed(123)


#--------------------------------------
# Paths & Directories
#--------------------------------------
tables_dir = "Results/Tables/"
fgsea_dir  = "Results/FGSEA/"


#--------------------------------------
# Load Differential Expression Data
#--------------------------------------
deg_results = list(
  
  IB_vs_UB = read.csv(
    file.path(tables_dir, "edgeR_results_IB_vs_UB.csv"),
    row.names = 1,
    check.names = FALSE
  ),
  
  RC_vs_IB = read.csv(
    file.path(tables_dir, "edgeR_results_RC_vs_IB.csv"),
    row.names = 1,
    check.names = FALSE
  ),
  
  RC_vs_UB = read.csv(
    file.path(tables_dir, "edgeR_results_RC_vs_UB.csv"),
    row.names = 1,
    check.names = FALSE
  )

)

comparisons = c(
  IB_vs_UB  = "Injured Cortex vs. Uninjured Cortex",
  RC_vs_IB  = "Reactive Culture vs. Injured Cortex",
  RC_vs_UB  = "Reactive Culture vs. Uninjured Cortex"
)


#--------------------------------------
# Ranked Gene Lists
#--------------------------------------
CreateRankedList = function(deg_df){
  
  ranked_df = deg_df %>%
    mutate(
      ranking_score = sign(logFC) * (-log10(PValue))
      ) %>%
    arrange(PValue) %>%
    filter(
      !is.na(ranking_score), !is.na(ENTREZID)
      ) %>%
    distinct(ENTREZID , .keep_all = T)
  
  ranked_vector = ranked_df$ranking_score
  
  names(ranked_vector) = ranked_df$ENTREZID
  
  ranked_vector = sort(ranked_vector, decreasing = T)
  
}

ranked_lists = list()

for (comparison in names(deg_results)){
  
  ranked_lists[[comparison]] = CreateRankedList(
    deg_results[[comparison]]
  )

}


#--------------------------------------
# Load MSigDB Gene Sets
#--------------------------------------

# Hallmark
hallmark_gene_sets = msigdbr(
  species = "Mus musculus",
  db_species = "MM",
  collection = "MH"
)

hallmark_pathways = split(
  x = hallmark_gene_sets$ncbi_gene,
  f = hallmark_gene_sets$gs_name
)


#--------------------------------------
# fgsea Analysis Run
#--------------------------------------
fgseaRUN = function(ranked_lists, pathways){
  
  results = list()
  
  for (comparison in names(ranked_lists)){
    fgsea_results = fgseaMultilevel(
      pathways = pathways,
      stats = ranked_lists[[comparison]],
      minSize = 15,
      maxSize = 500
)
    
    results[[comparison]] = fgsea_results %>%
      filter(padj < 0.05) %>%
      mutate(
        group = case_when(
          NES > 0 ~ "Up",
          NES < 0 ~ "Down"
          )
        )%>%
      arrange(desc(NES))
  }
  
  return(results)
}

# Hallmark Results
hallmark_res = fgseaRUN(
  ranked_lists, 
  hallmark_pathways
)


#--------------------------------------
# Visualizations
#--------------------------------------

#-----Pathways-Barplot-----------------
BarPlot = function(
  fgsea_results, 
  plot_title,
  db
  ){
  
  ggplot(
    fgsea_results, 
    aes(x = reorder(
      pathway, NES
      ),
        y = NES,
        fill = group
  ))+
  geom_col()+
  labs(
    title = plot_title,
    subtitle = paste0(db, " Pathway Enrichment"),
    x = "", 
    y = "Normalized Enrichment Score (NES)"
  )+
  scale_fill_manual(values = c(
    "Up" = "#213F95",
    "Down" = "#69C36F"
  ))+
  coord_flip()+
  theme_minimal()+
  theme(
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(face = "italic"),
    axis.text.y = element_text(face = "bold", size = 10)
    )
}

for (comparison in names(hallmark_res)) {
  
  hallmark_bar_plot = BarPlot(
    hallmark_res[[comparison]], 
    comparisons[[comparison]],
    "Hallmark"
    )
  
  ggsave(
    filename = paste0("Gsea_Hallmark_Barplot_", comparison, ".pdf"),
    plot = hallmark_bar_plot,
    width = 12, 
    height = 8 , 
    path = fgsea_dir
    )
}


#-----Gsea-Table Plot------------------
GseaTable = function(
    pathways, 
    ranked_list, 
    fgsea_res) {
  
  plotGseaTable(
    pathways = pathways[fgsea_res$pathway],
    stats = ranked_list,
    fgseaRes = fgsea_res
    )
}

for (comparison in names(hallmark_res)){
  
  pdf(paste0(
    fgsea_dir, "Gsea_Hallmark_Table_", comparison, ".pdf"
    ),
    height = 8,
    width = 14
  )   
  
  print(
    GseaTable(
      hallmark_pathways, 
      ranked_lists[[comparison]],
      hallmark_res[[comparison]] 
    )
  )

  dev.off()
}


#-----Enrichment Plot------------------
EnrichmentPlot = function(
    target_pathways,
    pathways,
    ranked_list,
    plot_title,
    file_name
    ){
  
  plots = list()
  
  for (pathway in target_pathways) {
    
    title = ifelse(
      str_detect(pathway, "HALLMARK_"),
      str_replace_all(str_remove(pathway, "HALLMARK_"), "_", " "),
      str_replace_all(pathway, "_", " "))

    plots[[pathway]] = plotEnrichment(
      pathway = pathways[[pathway]], 
      stats =  ranked_list)+
      labs(title = title)+
      theme_minimal()+
      theme(
        plot.title = element_text(size = 15))
    }

  panel = wrap_plots(plots) +
    plot_annotation(
      title = plot_title,
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 18)
      )
    )
  
  ggsave(
    filename = file_name, 
    plot = panel,
    path = fgsea_dir,width = 16,
    height = 8
    )
}

selected_pathways = list(
  
  IB_vs_UB = c(
    "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
    "HALLMARK_IL6_JAK_STAT3_SIGNALING", 
    "HALLMARK_ALLOGRAFT_REJECTION"
  ),
  
  RC_vs_IB = c(
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_DNA_REPAIR",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
    "HALLMARK_INFLAMMATORY_RESPONSE"
  ),
  
  RC_vs_UB = c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_GLYCOLYSIS", 
    "HALLMARK_MTORC1_SIGNALING", 
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
  )

)

for (comparison in names(selected_pathways)){
  
  EnrichmentPlot(
    selected_pathways[[comparison]],
    hallmark_pathways,
    ranked_lists[[comparison]],
    paste0(comparisons[[comparison]], " Enriched Hallmark Pathways"),
    paste0(
      "Gsea_Hallmark_Enrichment_Pathways_", comparison, ".pdf"
    )
  )
}