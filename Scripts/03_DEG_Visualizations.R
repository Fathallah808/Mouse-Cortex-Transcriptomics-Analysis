#======================================
# Differential Expression Visualizations
#======================================

#--------------------------------------
# Libraries
#--------------------------------------
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

set.seed(123)


#--------------------------------------
# Paths & Directories
#--------------------------------------
input_dir   = "Data/"
tables_dir  = "Results/Tables/"
figures_dir = "Results/Figures/"


#--------------------------------------
# Load Data & MetaData
#--------------------------------------
IB_vs_UB = read.csv(
  file.path(tables_dir, "edgeR_results_IB_vs_UB.csv"),
  row.names = 1,
  check.names = FALSE
)

RC_vs_IB = read.csv(
  file.path(tables_dir, "edgeR_results_RC_vs_IB.csv"),
  row.names = 1,
  check.names = FALSE
)

RC_vs_UB = read.csv(
  file.path(tables_dir, "edgeR_results_RC_vs_UB.csv"),
  row.names = 1,
  check.names = FALSE
)

metadata = read.csv(
  file.path(input_dir, "MetaData.csv"),
  row.names = 1,
  check.names = FALSE
)

logCPM = read.csv(
  file.path(tables_dir, "logCPM.csv"),
  row.names = 1,
  check.names = FALSE
)

comparisons = list(
  IB_vs_UB  = IB_vs_UB,
  RC_vs_IB  = RC_vs_IB,
  RC_vs_UB  = RC_vs_UB
)

plot_titles = list(
  IB_vs_UB  = "Injured Cortex vs. Uninjured Cortex",
  RC_vs_IB  = "Reactive Culture vs. Injured Cortex",
  RC_vs_UB  = "Reactive Culture vs. Uninjured Cortex"
)


#--------------------------------------
# PCA Analysis
#--------------------------------------
pca = prcomp(t(logCPM))

pca_variance = pca$sdev^2
variance_percent = round(
pca_variance / sum(pca_variance) * 100, 1
)

pca_results = pca$x %>%
  as.data.frame() %>%
  mutate(
  samples = metadata$title,
  treatment = metadata$treatment
)

pca_plot = ggplot(
  pca_results, 
  aes(PC1, PC2, colour = treatment, shape = treatment))+
  geom_point(size = 4, alpha = 0.9)+
  labs(
  title = "PCA of Reactive Culture, Injured and Uninjured Cortex Samples",
  x = paste0("PC1: ",variance_percent[1],"%"),
  y = paste0("PC2: ",variance_percent[2],"%"))+
  geom_text_repel(aes(label = samples, colour = treatment),
  size = 3.3,
  show.legend = FALSE)+
  scale_color_manual(values = c(
  "ReactiveCulture" = "#B40426",
  "Injured" = "#1F73C2",
  "Uninjured" = "#3D434A"))+
  theme_minimal(base_size = 10)+
  theme(
    plot.title = element_text(size = 15)
)

ggsave(
  filename = "PCA_All_Samples.pdf",
  path = figures_dir, 
  plot = pca_plot,
  height = 7, 
  width = 9
)


#--------------------------------------
# LogCPM & Z-Score Heatmaps
#--------------------------------------
TopDEGs = function(expr_data){
  expr_data %>%
    filter(FDR < 0.05, abs(logFC) > 1, !is.na(SYMBOL)) %>%
    arrange(FDR, desc(abs(logFC))) %>%
    head(50)
}

CreateExprMatrix = function(top_degs){
  expr_matrix = as.matrix(logCPM[rownames(top_degs),])
  rownames(expr_matrix) = make.unique(top_degs$SYMBOL)
  return(expr_matrix)
}

ZscoreMatrix = function(expr_matrix){
  t(scale(t(expr_matrix)))
}

for (comparison in names(comparisons)){
  
top_degs = TopDEGs(comparisons[[comparison]])
  
expr_matrix = CreateExprMatrix(top_degs)
  
Annotations = data.frame(
    row.names = colnames(expr_matrix),
    treatment = metadata$treatment,
    tissue    = metadata$tissue 
)
  
AnnotationsColors = list(
    treatment = c(
    "ReactiveCulture" = "mediumturquoise",
    "Injured" = "steelblue1",
    "Uninjured" = "plum"),
    tissue = c("cortex" = "#3D434A")
)

logcpm_heatmap_col = rev(
    brewer.pal(11, "RdYlBu")
)

zscore_col = colorRamp2(
  c(-2, -1, 0, 1, 2),
  c("royalblue", "skyblue3", "white", "lightcoral", "red3")
)

HMAnnotations = HeatmapAnnotation(
    df = Annotations, 
    col = AnnotationsColors
)
  
zscore_matrix = ZscoreMatrix(expr_matrix)
  
pdf(
    file.path(
    figures_dir, 
    paste0("Heatmap_Top50_DEGs_", comparison, ".pdf")),
    width = 12, 
    height = 10
)
  
logcpm_heatmap = Heatmap(
    expr_matrix,
    name = "logCPM",
    col = logcpm_heatmap_col,
    cluster_rows = T,
    cluster_columns = T,
    top_annotation = HMAnnotations,
    column_title = paste0("Top 50-DEGs: ", plot_titles[[comparison]]),
    column_title_gp = gpar(
      fontsize = 15, 
      fontface = "bold"),
    column_names_gp = gpar(
      fontface = "italic"),
    column_names_rot = 45,
    row_names_gp = gpar(
      fontface = "italic"),
    heatmap_legend_param = list(
      legend_height = unit(5, "cm")
      )
)
  
draw(
  logcpm_heatmap,
    merge_legend = TRUE
)
  
dev.off()

pdf(
  file.path(
  figures_dir, 
  paste0("Heatmap_Zscore_" , comparison, ".pdf")),
  width = 12, 
  height = 10
)
 
zscore_heatmap = Heatmap(
    zscore_matrix,
    name = "Z-Score",
    col = zscore_col,
    cluster_rows = T,
    cluster_columns = T,
    top_annotation = HMAnnotations,
    column_names_gp = gpar(
    fontface = "italic"),
    column_title = paste0(
      "Top 50-DEGs (Z-Score): ", plot_titles[[comparison]]
      ),
    column_title_gp = gpar(
      fontsize = 15, 
      fontface = "bold"),
    column_names_rot = 45,
    row_names_gp = gpar(
    fontface = "italic"),
    heatmap_legend_param = list(
      legend_height = unit(5, "cm")
      )
)
  
draw(
  zscore_heatmap,
  merge_legend = TRUE
)
  
dev.off()
  
}


#--------------------------------------
# Volcano Plot
#--------------------------------------
SignificantGenes = function(main_data){
  main_data = main_data %>%
  mutate(expression = 
        case_when(
        logFC > 1 & FDR < 0.05 ~ "Upregulated",
        logFC < -1 & FDR < 0.05 ~ "Downregulated",
        TRUE ~ "Non-Significant")
    )
}

TopUpandDownGenes = function(sig_data){
  up_genes = sig_data %>%
    filter(expression == "Upregulated", !is.na(SYMBOL)) %>%
    arrange(FDR, desc(abs(logFC))) %>%
    slice(1:10)
  down_genes = sig_data %>%
    filter(expression == "Downregulated", !is.na(SYMBOL)) %>%
    arrange(FDR, desc(abs(logFC))) %>%
    slice(1:10)
  top_genes = bind_rows(up_genes, down_genes)
}

VolcanoPlot = function(data, top_genes, plot_title){
  ggplot(data,
  aes(x = logFC,  y = -log10(FDR), colour = expression)
  )+
  geom_point()+
  labs(
  title = plot_title,
  x = expression(log[2]~FC), 
  y = expression(-log[10]~FDR),
  color = NULL
  )+
  geom_text_repel(
  data = top_genes,aes(label = SYMBOL),
  colour = "black",
  show.legend = FALSE)+
  scale_color_manual(values = c(
      "Upregulated" = "#B40426",
      "Downregulated" = "#1F73C2",
      "Non-Significant" = "#3D434A")
  )+
  xlim(-15, 15)+
  geom_vline(xintercept = c(-1, 1), 
               colour = "black", linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05) , 
               colour = "black", linetype = "dashed")+
  theme_minimal()+
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.position = "top")
}

for (comparison in names(comparisons)){
  
  data = SignificantGenes(comparisons[[comparison]])
  top_genes = TopUpandDownGenes(data)
  
  plots = VolcanoPlot(data, 
  top_genes, 
  paste0("Differentailly Expressed Genes in ", plot_titles[[comparison]])
  )
  
  ggsave(filename = paste0( "Volcano_", comparison, ".pdf"),
         plot = plots,
         path = figures_dir,
         height = 8,
         width = 12
         )
}