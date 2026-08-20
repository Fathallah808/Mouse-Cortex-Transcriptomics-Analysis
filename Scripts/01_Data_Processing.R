#======================================
# Data Import & Processing
#======================================

#--------------------------------------
# Libraries
#--------------------------------------
library(GEOquery)
library(tidyverse)


#--------------------------------------
# Paths & Directories
#--------------------------------------
output_dir = "Data/"

#--------------------------------------
# Import Raw Counts
#--------------------------------------
data = read.csv(
  "Data/GSE330004_raw_counts.csv", 
  row.names = 1
)

counts = data %>%
  select(1:9)

colnames(counts) = gsub(
  " ", "", colnames(counts)
)


#--------------------------------------
# Import MetaData from GEO
#--------------------------------------
gse = getGEO(
  "GSE330004", 
  GSEMatrix = TRUE
)

metadata = pData(
  phenoData(gse[[1]])
)


#--------------------------------------
# MetaData Prepare
#--------------------------------------
metadata = metadata %>%
  mutate(description = gsub(
    "Library name: ", "", description
    )
)

stopifnot(
  "Raw Sample Order Mismatch between Counts Matrix & MetaData" =
  colnames(counts) == metadata$description
)

metadata_sub = metadata %>%
  select(title, `treatment:ch1`, 
         `tissue:ch1`, characteristics_ch1.2) %>%
  rename(treatment = `treatment:ch1`,
         tissue = `tissue:ch1`,
         characteristics = "characteristics_ch1.2") %>%
  mutate(
    treatment = case_when(
    grepl("reactive", treatment) ~ "ReactiveCulture",
    grepl("cortical", treatment) ~ "Injured",
    grepl("uninjured cortex", treatment) ~ "Uninjured"
    ),
    title = 
    str_replace_all(
    title, c(
        "In vitro, Reactive culture, rep " = "Vitro_RC_",
        "In vivo, Injured brain, rep " = "Vivo_IB_",
        "In vivo, Uninjured brain, rep " = "Vivo_UB_")
    ),
    characteristics = gsub(
      "treatment: ", "", characteristics
      )
)


#--------------------------------------
# Matching MetaData with Count Matrix
#--------------------------------------
rownames(metadata_sub) = metadata_sub$title
colnames(counts) = metadata_sub$title


#----Validation------------------------
stopifnot(
  "Sample Names Mismatch between MetaData & Counts Matrix" =
    all(rownames(metadata_sub) == colnames(counts))
)


#--------------------------------------
# Export Data
#--------------------------------------
write.csv(
  counts, 
  file.path(output_dir, "Counts.csv")
) 

write.csv(
  metadata_sub, 
  file.path(output_dir, "MetaData.csv")
) 