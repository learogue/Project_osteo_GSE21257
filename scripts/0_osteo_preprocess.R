#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------------------------------------------------
# R script : Pre-process osterosarcoma microarray data from GSE169077
# Auteur  : Léa ROGUE
# Date    : 07-03-2025
# Description : This R script uses the package GEOquery to download microarray data from the GEO database (GSE21257) 
# and performs preprocessing, quality control, and normalization steps on the raw data. 
# The script includes:
# 	- Downloading raw Illumina microarray data (GSE21257) from GEO repository.
#	- Reading the raw data and summarizing it to generate expression matrices.
#	- Normalizing the data using the NEQC (normal-exponential background correction and quantile normalization).
#	- Annotating the data with gene symbols.
#	- Merging expression data with additional annotation such as CTA genes and immune cell signatures.
#	- Calculating Z-scores for gene expression values to facilitate comparisons between samples.
#	- Averaging expression data for immune cell signatures.
#	- Writing results to multiple output files for downstream analysis.
#The final outputs include:
#	- A summarized expression matrix with annotations (CTA genes, immune cell signatures).
#	- Z-scores for the expression values of genes, allowing comparison of expression patterns.
# ------------------------------------------------------------------------------------------------------------------------

# Load packages
library(GEOquery)
library(beadarray)
library(arrayQualityMetrics)
library(dplyr)

# Function to calculate Z scores
calculate_z_scores <- function(df_input, col) {
  # Exclude col
  data_values <- df_input[, -c(col)]

  # Calculate Z-scores
  z_scores_row <- t(scale(t(data_values)))

  # Add columns
  df_z_scores <- cbind(df_input[, c(col)], z_scores_row)

  # Return df
  return(df_z_scores)
}

# Data directory path
dir_files <- "/home/ubuntu/Project_osteo_GSE21257/data/GSE21257_raw"

# Create the dir
dir.create(dir_files)

# Dowload data
getGEOSuppFiles("GSE21257", baseDir = dir_files, makeDirectory = FALSE)

# Read raw data
data <- readIllumina(dir="/home/ubuntu/Project_osteo_GSE21257/data/GSE21257_raw/", useImages=FALSE, illuminaAnnotation="Humanv2")

# Sample groups
sample_group <- seq(1, 53)
sample_group <- rep(sample_group, each = 2)

# Change this function to have non log2 transform values to perform neqc normalization
greenChannel@transFun[[1]] <- function (BLData, array) 
{       
    x = getBeadData(BLData, array = array, what = "Grn")
    return(x)
}

# Summarize data
data_sum <- beadarray::summarize(data, useSampleFac=TRUE, sampleFac=sample_group, channelList=list(greenChannel))

# Add sample IDs
data_sum$sampleID <- c("GSM530667", "GSM530899", "GSM531283", "GSM531284", "GSM531285", "GSM531286",
                       "GSM531287", "GSM531288", "GSM531289", "GSM531290", "GSM531291", "GSM531292",
                       "GSM531293", "GSM531294", "GSM531295", "GSM531296", "GSM531297", "GSM531298",
                       "GSM531299", "GSM531300", "GSM531301", "GSM531302", "GSM531303", "GSM531304",
                       "GSM531305", "GSM531306", "GSM531307", "GSM531308", "GSM531309", "GSM531310",
                       "GSM531311", "GSM531312", "GSM531313", "GSM531314", "GSM531319", "GSM531320",
                       "GSM531321", "GSM531322", "GSM531323", "GSM531324", "GSM531325", "GSM531326",
                       "GSM531327", "GSM531328", "GSM531329", "GSM531330", "GSM531331", "GSM531332",
                       "GSM531333", "GSM531334", "GSM531335", "GSM531351", "GSM531352")


# Normalization with neqc which combine normal-exponential background correction with quantile normalisation as suggested in the limma package
data_sum_norm <- normaliseIllumina(data_sum, method="neqc", transform="none")

# Quality control like affymetrix
#arrayQualityMetrics(expressionset = data_sum, outdir = "../results/qc/qc_raw", force = TRUE, do.logtransform = TRUE)
#arrayQualityMetrics(expressionset = data_sum_norm, outdir = "../results/qc/qc_norm", force = TRUE, do.logtransform = TRUE)

# Annotate
anno <- addFeatureData(data_sum_norm, toAdd = "SYMBOL")
anno <- fData(anno)
df_anno <- subset(anno, !is.na(SYMBOL))
df_anno <- df_anno %>% dplyr::select(IlluminaID, SYMBOL)

# Merge expression matrix
expr_data <- as.data.frame(exprs(data_sum_norm))
colnames(expr_data) <- data_sum$sampleID
expr_data$IlluminaID <- rownames(expr_data)
df_expr <- merge(df_anno, expr_data, by = "IlluminaID")
df_expr <- df_expr[, -1]

# Mean for the same genes
df_expr_avg <- df_expr %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# Read CTA file
df_expr$CTA <-  NA
df_CTA <- read.table("../data/CTA_list_clean.txt", header = FALSE)

# Rename col
colnames(df_CTA) <- c("SYMBOL")

# Df with CTA and whole genes by adding CTA in col CTA for CTA genes
df_CTA_whole <- df_expr_avg %>%
  left_join(df_CTA, by = "SYMBOL") %>%
  mutate(CTA = ifelse(SYMBOL %in% df_CTA$SYMBOL, "CTA", NA))

# Reorganize columns
df_CTA_whole <- df_CTA_whole %>%
  dplyr::select(SYMBOL, CTA, everything())

# Read immune cells genes
df_immune_sign <- read.table("../data/immune_cells_mhc_genes.tsv", header = FALSE, sep = "\t")
colnames(df_immune_sign) <- c("Signature", "Gene")

# Merge df_CTA_whole and df_immune_sign to associate signatures
df_CTA_immune_sign_whole <- merge(df_CTA_whole, df_immune_sign, by.x = "SYMBOL", by.y = "Gene", all.x = TRUE)

# Reorganize the columns
df_CTA_immune_sign_whole <- df_CTA_immune_sign_whole %>%
  dplyr::select(SYMBOL, CTA, Signature, everything())
write.table(df_CTA_immune_sign_whole, file = "../results/expr_matrix_int_CTA_sign_imm.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Z-scores on rows on whole genes to comparate expression between genes with a robust manner
df_CTA_immune_whole_clean_z_scores <- calculate_z_scores(df_CTA_immune_sign_whole, c(1,2,3))
write.table(df_CTA_immune_whole_clean_z_scores, file = "../results/expr_matrix_CTA_sign_imm_z_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Average the expression between same immune cells types
# Take rows with immune cells signature from normalized data
df_avg_immune_sign <- df_CTA_immune_sign_whole %>%
  filter(Signature != "NA")

# Group by signature and calculate mean of expression values
df_avg_immune_sign_final <- df_avg_immune_sign %>%
  dplyr::select(-c(SYMBOL, CTA)) %>% 
  group_by(Signature) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) 

# Z-scores
df_avg_immune_sign_z_scores <- calculate_z_scores(df_avg_immune_sign_final, 1)
write.table(df_avg_immune_sign_z_scores, file = "../results/imm_sign_avg_z_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

