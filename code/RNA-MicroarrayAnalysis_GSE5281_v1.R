---
title: "RNA-analysis_GSE5281"
author: "DG"
date: "2026-04-16"
output:
  html_document: default
toc: true
toc_depth: 2
theme: flatly
---
#The code below demostrates analysis of microarray expression intensities from the GSE5281 dataset 
---
## 1.1 Install Packages 
install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "annotate", "hgu133plus2.db"))

install.packages("pheatmap")

install.packages(c("ggplot2", "dplyr", "ggrepel"))

install.packages("UpSetR")

install.packages("matrixStats")

## 1.2 Load Libraries 
library(GEOquery)

library(limma)

library(hgu133plus2.db)

library(pheatmap)

library(AnnotationDbi)

library(ggplot2)

library(dplyr)

library(ggrepel)  # Optional but recommended for clean labeling

library(UpSetR)

library(matrixStats)

## 2. Load Data
gse <- getGEO("GSE5281", GSEMatrix = TRUE)

gse <- gse[[1]] #if multiple platforms select the first

expr <- exprs(gse) #expression matrix

meta <- pData(gse) #metadata

dim(expr) #check dimensions expression matrix

head(meta) #check expression metadata


## 3. Data Preprocessing

### 3.1. Ensure Consistent Labeling
colnames(meta) #look at metadata columns

lapply(meta[, grep("characteristics", colnames(meta))], unique) #check all the characteristics_c1.* columns to find the one with the group information details

meta$clean <- gsub("bio-source name:\\s*", "", 
                   meta$characteristics_ch1.1, 
                   ignore.case = TRUE) #case insensitive removal of "bio-source name:"

meta$clean <- tolower(meta$clean) #turn every letter to lowercase

meta$clean <- gsub(" ", "_", meta$clean) #replace spaces with underscores

meta$clean <- gsub("_+", "_", meta$clean) #fix messy underscores

unique(meta$clean) #check that the conversion was applied correctly

### 3.2. Rename variable labels to represent experimental factors
meta$group <- ifelse(grepl("affected", meta$clean),
                     "AD", "Control") #disease variable classification

meta$region <- sub("^([a-z]+)_.*", "\\1", meta$clean) #brain region variable string extraction

meta$region <- toupper(meta$region) #make brain region abbreviations uppercase to comply with neuroscience labeling convention

unique(meta$region) 

table(meta$region, useNA ="ifany") #check for any unexpected values

table(meta$group, meta$region) #check correct conversion with group and region table

prop.table(table(meta$group, meta$region), margin = 2) #check for biological balance,QC

barplot(table(meta$group, meta$region),
        beside = TRUE,
        legend = TRUE) # check for biological balance,visualQC
        
### 3.3. Investigate raw data values - ranges show if the original data are log transformed
range(expr)

summary(as.vector(expr))

hist(expr, breaks = 100) #histogram

### 3.4. Perform log2 transformation
expr_log <- log2(expr + 1)

range(expr_log)

hist(expr_log, breaks = 100) #histogram

meta_clean <- data.frame(
  sample_id = colnames(expr),
  group = meta$group,
  region = meta$region,
  label = meta$clean
) 

all(meta_clean$sample_id == colnames(expr_log)) #check that the order matches the matrix

write.csv(meta_clean, "data/clean/meta_clean.csv", row.names = FALSE) #save metadata

write.csv(expr_log, "data/clean/expr_log2.csv") #save clean expression matrix

### 3.5. Perform Data Cleaning
expr_log_no_control <- expr_log[!grepl("^AFFX", rownames(expr_log)), ] #remove control (AFFX) probes

sum(grepl("^AFFX", rownames(expr_log_no_control)))   #check it worked, should be 0

write.csv(expr_log_no_control, "data/clean/expr_log2_noAFFXcontrol.csv") #save as .csv

expr_highconf <- expr_log_no_control[
  !grepl("_x_at$|_f_at$", rownames(expr_log_no_control)),
] #remove low specificity probes, keep only the high confidence ones

sum(grepl("_x_at$|_f_at$", rownames(expr_highconf)))  #check it worked, should be 0

write.csv(expr_highconf, "data/clean/expr_log2_noAFFX_noX_noF.csv") #save as .csv

### 3.6. Dataset Validation before limma analysis
dim(expr_highconf)

dim(meta_clean) #check dataset dimensions

all(colnames(expr_highconf) == meta_clean$sample_id) #check sample alignment

sum(is.na(expr_highconf)) #check for missing values

### 3.7. Check expression distribution
hist(expr_highconf, breaks = 100) #histogram

boxplot(expr_highconf, outline = FALSE, las = 2) #boxplot

plot(rowMeans(expr_highconf),
     matrixStats::rowVars(as.matrix(expr_highconf)),
     xlab = "Mean", ylab = "Variance") #mean variance

table(meta_clean$group, meta_clean$region) #check group balance

meta_clean$group_region <- factor(
  paste(meta_clean$group, meta_clean$region, sep = "_")
)

design <- model.matrix(~ 0 + group_region, data = meta_clean)

colnames(design) <- levels(meta_clean$group_region)

qr(design)$rank == ncol(design) 

table(meta_clean$group_region) #check proper matrix formation

### 4. Limma Analysis

### 4.1. Load Data
load("data/clean/expr_log2_noAFFX_noX_noF.csv") #load cleaned data

### 4.2. Align metadata
meta <- meta[match(colnames(expr_highconf), rownames(meta)), ]

### 4.3. Create group_region
meta$group_region <- factor(paste(meta$group, meta$region, sep = "_"))

design <- model.matrix(~ 0 + group_region, data = meta)

colnames(design) <- levels(meta$group_region)

### 4.4. Define contrasts
contrasts <- makeContrasts(
  EC  = AD_EC  - Control_EC,
  HIP = AD_HIP - Control_HIP,
  MTG = AD_MTG - Control_MTG,
  PC  = AD_PC  - Control_PC,
  SFG = AD_SFG - Control_SFG,
  VCX = AD_VCX - Control_VCX,
  levels = design
) #define contrasts

contrasts #check contrast matrix

### 4.5. Fit the model
fit <- lmFit(expr_highconf, design)

fit2 <- contrasts.fit(fit, contrasts)

fit2 <- eBayes(fit2) #fit limma model with contrasts

### 4.6. Extract results
res_EC <- topTable(fit2, coef = "EC", number = Inf) #extract results per brain region

results_list <- lapply(colnames(contrasts), function(region) {
  topTable(fit2, coef = region, number = Inf)
}) #loop for all regions

names(results_list) <- colnames(contrasts)

### 4.7. Save results
for (region in names(results_list)) {
  write.csv(
    results_list[[region]],
    paste0("results/limma_", region, ".csv"),
    row.names = TRUE
  )
} #save results for all regions

results_list

results_list[["EC"]] #results list for each element

### 4.8. Annotate results
annotate_results <- function(res_table) {
  
  probe_ids <- rownames(res_table)
  
  res_table$SYMBOL <- mapIds(
    hgu133plus2.db,
    keys = probe_ids,
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  res_table$ENTREZID <- mapIds(
    hgu133plus2.db,
    keys = probe_ids,
    column = "ENTREZID",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  res_table$GENENAME <- mapIds(
    hgu133plus2.db,
    keys = probe_ids,
    column = "GENENAME",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  return(res_table) 
} #annotate results

annotated_results <- lapply(results_list, annotate_results) #apply annotation to all regions

dir.create("results/annotated", showWarnings = FALSE)

for (region in names(annotated_results)) {
  
  file_name <- paste0("results/annotated/limma_", region, "_annotated.csv")
  
  write.csv(
    annotated_results[[region]],
    file = file_name,
    row.names = TRUE
  )
} #save the annotated results per region

all_probes <- rownames(expr_highconf)

### 4.9. Generate standalone mapping file
probe_map <- data.frame(
  PROBEID = all_probes,
  SYMBOL = mapIds(hgu133plus2.db, all_probes, "SYMBOL", "PROBEID", multiVals = "first"),
  ENTREZID = mapIds(hgu133plus2.db, all_probes, "ENTREZID", "PROBEID", multiVals = "first"),
  GENENAME = mapIds(hgu133plus2.db, all_probes, "GENENAME", "PROBEID", multiVals = "first")
) #create a standalone mapping file

write.csv(probe_map, "data/clean/probe_to_gene_mapping.csv", row.names = FALSE) #save the file

### 4.10. Filter significant genes
sig_EC <- annotated_results[["EC"]]

sig_EC <- sig_EC[sig_EC$adj.P.Val < 0.05, ] #filter significant genes

dir.create("results/DE/filtered_significant", recursive = TRUE, showWarnings = FALSE)

dir.create("results/plots/volcano", recursive = TRUE, showWarnings = FALSE)

dir.create("results/plots/cross_region", recursive = TRUE, showWarnings = FALSE)

dir.create("results/summaries/top_genes", recursive = TRUE, showWarnings = FALSE) #create directories to store analyses

padj_cutoff <- 0.05 #save significant genes per region

save_significant <- function(df, region_name) {
  
  sig <- df[df$adj.P.Val < padj_cutoff, ]
  
  # order by significance
  sig <- sig[order(sig$adj.P.Val), ]
  
  write.csv(
    sig,
    file = paste0("results/DE/filtered_significant/sig_", region_name, ".csv"),
    row.names = TRUE
  )
  
  return(sig)
} #function to filter significnt genes

sig_list <- lapply(names(annotated_results), function(r) {
  save_significant(annotated_results[[r]], r)
})

names(sig_list) <- names(annotated_results) #apply to all regions

### 5. Visualisationscolnames(rna_results)

### 5.1. Draw Volcano Plots
plot_volcano <- function(df, region_name) {
  
  df$negLogP <- -log10(df$adj.P.Val)
  
  df$label <- ifelse(df$adj.P.Val < 0.05 & abs(df$logFC) > 1,
                     df$SYMBOL,
                     NA)
  
  p <- ggplot(df, aes(x = logFC, y = negLogP)) +
    geom_point(alpha = 0.5) +
    geom_text(aes(label = label), na.rm = TRUE, size = 3) +
    ggtitle(paste("Volcano Plot -", region_name)) +
    xlab("log2 Fold Change") +
    ylab("-log10 adj P-value")
  
  ggsave(
    filename = paste0("results/plots/volcano/volcano_", region_name, ".png"),
    plot = p,
    width = 6,
    height = 5
  )
}

lapply(names(annotated_results), function(r) {
  plot_volcano(annotated_results[[r]], r)
})

### 5.2. Cross-Region Comparison Plots
top_genes <- lapply(sig_list, function(df) {
  head(df$SYMBOL, 20)
})

all_genes <- unique(unlist(top_genes))

binary_matrix <- sapply(top_genes, function(g) all_genes %in% g)

class(binary_matrix) #check current matrix form

binary_df <- as.data.frame(binary_matrix) #force it to be a datafram

dim(binary_df) #check dimensions

colSums(binary_df)

upset(as.data.frame(binary_matrix),
      order.by = "freq")
png("results/plots/cross_region/region_overlap.png", width = 800, height = 600)

upset(as.data.frame(binary_matrix), order.by = "freq")

dev.off() #save plot

common_genes <- Reduce(intersect, lapply(sig_list, function(x) x$SYMBOL))

### 5.3. Filtered DE Results per region
sig_list   

common_genes <- Reduce(intersect, lapply(sig_list, function(df) df$SYMBOL)) # extract logFC per region #extract logFC per region

logFC_mat <- sapply(names(sig_list), function(region) {
  
  df <- sig_list[[region]]
  
  df$logFC[match(common_genes, df$SYMBOL)]
  
})

rownames(logFC_mat) <- common_genes #build logFC matrix

logFC_mat <- logFC_mat[rowSums(is.na(logFC_mat)) < ncol(logFC_mat) - 1, ] #clean logFC matrix

pheatmap(
  logFC_mat,
  scale = "row",           # very important for interpretation
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  main = "DE Gene Patterns Across Brain Regions (AD vs Control)"
) #draw heatmap with pheatmap

rowVars(logFC_mat)

top_genes <- names(sort(rowVars(logFC_mat), decreasing = TRUE))[1:50]

logFC_top <- logFC_mat[top_genes, ] #filter to include only top variable genes
pheatmap(logFC_top, scale = "row", na_col = "grey")

pheatmap(logFC_top, scale = "row", na_col = "grey") #draw filtered heatmap
dir.create("results/plots/heatmap", recursive = TRUE, showWarnings = FALSE) #save heatmap

### 5.4. Define AD-relevant genes
alz_genes <- c(
  "APP",      # amyloid precursor protein
  "MAPT",     # tau
  "PSEN1",
  "PSEN2",
  "APOE",
  "BACE1",
  "CLU",
  "SORL1",
  "TREM2",
  "NEAT1"     # lncRNA involved in neurodegeneration
)

alz_genes_present <- intersect(rownames(logFC_mat), alz_genes) #match AD-relevant genes to matrix

notation_df <- data.frame(
  Alzheimer_gene = ifelse(rownames(logFC_mat) %in% alz_genes_present,
                          "YES", "NO")
)

rownames(annotation_df) <- rownames(logFC_mat) #create highlited AD version
pheatmap(
  logFC_mat,
  scale = "row",
  na_col = "grey",
  annotation_row = annotation_df,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "AD vs Control LogFC (Alzheimer Genes Highlighted)"
) #plot highlighted matrix version

alz_matrix <- logFC_mat[alz_genes_present, , drop = FALSE] #zoom into only AD-relevant genes

pheatmap(
  alz_matrix,
  scale = "row",
  na_col = "grey",
  main = "Alzheimer Disease Core Genes Across Brain Regions"
) #plot heatmap with AD-selected genes
















  
 