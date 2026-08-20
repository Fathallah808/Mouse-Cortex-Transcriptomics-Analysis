#======================================
# Differential Expression Analysis (edgeR) 
#======================================

#--------------------------------------
# Libraries
#--------------------------------------
library(edgeR)
library(org.Mm.eg.db)
library(tidyverse)

set.seed(123)


#--------------------------------------
# Paths & Directories
#--------------------------------------
input_dir  = "Data/"
figures_dir = "Results/Figures/"
tables_dir  = "Results/Tables/"


#--------------------------------------
# Load Raw Counts & MetaData
#--------------------------------------
counts = read.csv(
  file.path(input_dir, "Counts.csv"),
  row.names = 1,
  check.names = FALSE
  )

metadata = read.csv(
  file.path(input_dir, "MetaData.csv"), 
  row.names = 1, 
  check.names = FALSE
  )

#----Validation------------------------
stopifnot("Sample Order Mismatch between MetaData & Counts" =
  all(rownames(metadata) == colnames(counts))
  )


#--------------------------------------
# edgeR Differential Analysis Run
#--------------------------------------

#----1: Experimental Design------------
TreatmentFactor = factor(metadata$treatment, 
                  levels = c("Uninjured", "Injured", "ReactiveCulture")
)

design = model.matrix(~0 + TreatmentFactor)


#----2: Create DGEList-----------------
dge = DGEList(
  counts = counts, 
  group = TreatmentFactor
)


#-------Library Size Barplot-----------
pdf(paste0(
  figures_dir, "Library_Sizes.pdf"), 
  height = 6, 
  width = 12
)

barplot(dge$samples$lib.size / 1e6, 
        names.arg = colnames(dge), 
        main = "Library Sizes",
        ylab = "Total Read Counts (Million)",
        xlab = "Samples",
        col = "grey70",
        width = 0.5, space = 0.5,
        cex.main = 1.2, 
        cex.axis = 1, 
        cex.lab = 1, 
        cex.names = 0.7
)

dev.off()


#----3: Low Genes Filtering------------
keep = filterByExpr(
  dge, 
  design
)

dge = dge[keep , , keep.lib.sizes = FALSE]


#----4: Normalization of Lib. Sizes----
dge = normLibSizes(
  dge, 
  method = "TMM"
)

logCPM = cpm(
  dge, 
  log = TRUE, 
  prior.count = 1, 
  normalized.lib.sizes = TRUE
) 


#------------MDS plot------------------
pdf(
  paste0(figures_dir, "MDS_All_Samples.pdf"),
  height = 7, 
  width = 12
)

plotMDS(dge, 
        labels = metadata$title,
        col = c("#3D434A" , "#1F73C2" ,"#B40426")[TreatmentFactor],
        main = "MDS Plot",
        xlab = paste0("MDS1"),
        ylab = paste0("MDS2"),
        cex.main = 1.5, 
        cex.axis = 1, 
        cex.lab = 1, 
        cex = 1
)

dev.off()


#----5: Estimate Dispersion------------
dge = estimateDisp(
  dge, 
  design
)


#------------BCV plot------------------
pdf(
  paste0(figures_dir, "BCV_Plot.pdf"), 
  height = 6,
  width = 8
)

plotBCV(dge,
        cex.main = 1.2, 
        cex.axis = 1, 
        cex.lab = 1
)

dev.off()


#----6: Quasi Likelihood Model Fitting-
fit = glmQLFit(
  dge, 
  design
)


#----------QLDisp Plot-----------------
pdf(
  paste0(figures_dir, "QLDisp_Plot.pdf"),
  height = 6,
  width = 8
)

plotQLDisp(fit,
           main = "Quasi-Likelihood Dispersion Plot",
           cex.main = 1.2, 
           cex.axis = 1, 
           cex.lab = 1
)

dev.off()


#----7: Differential Testing-----------
contrast_matrix = makeContrasts(
  RC_vs_IB = TreatmentFactorReactiveCulture - TreatmentFactorInjured,
  IB_vs_UB = TreatmentFactorInjured - TreatmentFactorUninjured,
  RC_vs_UB = TreatmentFactorReactiveCulture - TreatmentFactorUninjured,
  levels = design
)

qlf_results = list()

for (contrast in colnames(contrast_matrix)) {
  qlf_results[[contrast]] = glmQLFTest(fit, contrast =  contrast_matrix[, contrast])
}


#-------------MA Plot------------------
plot_titles = list(
  IB_vs_UB  = "Injured Cortex vs. Uninjured Cortex",
  RC_vs_IB  = "Reactive Culture vs. Injured Cortex",
  RC_vs_UB  = "Reactive Culture vs. Uninjured Cortex"
)

pdf(
  paste0(figures_dir, "MA_Plots.pdf"),
  height = 9,
  width = 20
)

par(
  mfrow = c(1,3)
)

for (name in names(qlf_results)) {
plotMD(qlf_results[[name]],
       main = plot_titles[[name]],
       cex = 0.4,
       cex.main = 2, 
       cex.axis = 1.5, 
       cex.lab = 1.5
       )
}

dev.off()


#--------------------------------------
# Extraction of Results & Gene Annotation
#--------------------------------------
resultsANN = function(qlf_result){
  results = topTags(qlf_result, n = Inf)$table
  results = results %>%
    mutate(
      SYMBOL = mapIds
           (org.Mm.eg.db,
             keys = rownames(results), 
             keytype = "ENSEMBL" , 
             column = "SYMBOL", 
             multiVals = "first"),
      GENENAME = mapIds
           (org.Mm.eg.db,
             keys = rownames(results), 
             keytype = "ENSEMBL" , 
             column = "GENENAME", 
             multiVals = "first"),
      ENTREZID = mapIds
           (org.Mm.eg.db,
             keys = rownames(results), 
             keytype = "ENSEMBL" , 
             column = "ENTREZID", 
             multiVals = "first")
      )
  }

results = list()

for (name in names(qlf_results)) {
  results[[name]] = resultsANN(qlf_results[[name]])
}


#--------------------------------------
# Export Results
#--------------------------------------
for (comparison in names(results)){
  write.csv(results[[comparison]],
    file.path(
      tables_dir, paste0("edgeR_results_" , comparison, ".csv"))
    )
}

write.csv(logCPM, 
          file.path(tables_dir, "LogCPM.csv")
) 