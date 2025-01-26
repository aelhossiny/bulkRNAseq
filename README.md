---
title: "Bulk RNASeq Analysis Pipeline"
output: github_document
---

# Table of content

* Introduction
* Analysis
* Notes
* Recommended Readings

# Introduction

This is an analysis pipeline for Bulk RNAseq data. 


```{r}
# setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(PCAtools)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
})

dir.create("outputs", showWarnings = F, recursive = T)

# loading data ------------------------------------------------------------

hallmarks_pathways <- msigdbr(species = "Homo sapiens", category = 'H', subcategory = NULL)
hallmarks_pathways <- split(x = hallmarks_pathways$ensembl_gene, f = hallmarks_pathways$gs_name)

ct_mtx <- read.table("data/human_samples_gene_count.txt", quote = "", header = T, fill = T, sep = "\t")
gene_symbols <- data.frame(gene_symbol = ct_mtx$gene_name) %>% `rownames<-`(ct_mtx$gene_id)

ct_mtx <- ct_mtx %>%
  select(-c("gene_name", "gene_chr", "gene_start", "gene_end", "gene_strand", "gene_length", "gene_biotype", "gene_description", "tf_family")) %>%
  column_to_rownames('gene_id')

metadata <- read.table("data/human_rnaseq_samples_info.txt", sep = '\t', header = T) %>%
  column_to_rownames('sample_id')
metadata <- metadata[colnames(ct_mtx),]
metadata$condition <- gsub("NF", "Normal_fat", metadata$condition)
metadata$condition <- gsub("WD", "WDLPS", metadata$condition)
metadata$condition <- gsub("DD", "DDLPS", metadata$condition)
metadata$recurrence <- gsub("NR", "No_Rec", metadata$recurrence)
metadata$recurrence <- gsub("R", "Rec", metadata$recurrence)

# DESeq object all samples ------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = ct_mtx,
                              colData = metadata,
                              rowData = gene_symbols,
                              design = ~ 1)
dds$condition <- factor(dds$condition, levels = c('Normal_fat', 'WDLPS', 'DDLPS'))

saveRDS(dds, "outputs/all_human_samples_dds.rds")
write.csv(merge(assay(vst(dds)),gene_symbols, by = 0),
          "outputs/human_normalized_count_mtx.csv",
          row.names = F)

# WDLPS vs fat ------------------------------------------------------------

wdlps_idx <- dds$condition %in% c('Normal_fat', 'WDLPS')
dds_wdlps <- dds[,wdlps_idx]
dds_wdlps$condition <- factor(dds_wdlps$condition, levels = c('Normal_fat', 'WDLPS'))
design(dds_wdlps) <- model.matrix(~ condition, data = colData(dds_wdlps))
# dds_wdlps <- estimateSizeFactors(dds_wdlps)
# idx <- rowSums(counts(dds_wdlps, normalized=TRUE) >= 10 ) >= 3
# dds_wdlps <- dds_wdlps[idx,]
dds_wdlps <- DESeq(dds_wdlps)

# > DGE analysis ----------------------------------------------------------

wdlps_v_fat <- lfcShrink(dds_wdlps, coef = 'conditionWDLPS', type = 'ashr') %>%
  data.frame() %>%
  merge(gene_symbols, by = 0) %>%
  filter(!is.na(padj)) %>%
  rename(gene_id = Row.names) %>%
  mutate(status = ifelse(log2FoldChange > 0.5 & padj < 0.05, 1,
                         ifelse(log2FoldChange < -0.5 & padj < 0.05, -1, 0)))

write.csv(wdlps_v_fat, "outputs/human_wdlps_v_fat_dge.csv", row.names = F, quote = F)

human_ann <- data.frame(gene = wdlps_v_fat$gene_id,
                        gene_symbol = wdlps_v_fat$gene_symbol,
                        pval_adj = wdlps_v_fat$padj)
human_wdlps_v_fat <- glimmaXY(x = wdlps_v_fat$log2FoldChange,
         xlab = 'log2FoldChange',
         y = -log(wdlps_v_fat$padj),
         ylab = 'significance',
         status = wdlps_v_fat$status,
         counts = vst(counts(dds_wdlps))[wdlps_v_fat$gene_id,],
         transform.counts = 'none',
         main = 'WDLPS vs Normal Fat',
         groups = colData(dds_wdlps)$condition,
         anno = human_ann,
         display.columns = c('gene_symbol','pval_adj'))

htmlwidgets::saveWidget(human_wdlps_v_fat,'docs/human_wdlps_v_fat_dge_interactive_app.html')

# > GSEA ------------------------------------------------------------------

wdlps_v_fat_rank <- wdlps_v_fat %>%
  pull(log2FoldChange) %>%
  `names<-`(wdlps_v_fat$gene_id) %>%
  sort()

wdlps_v_fat_gsea <- fgsea(pathways = hallmarks_pathways,
                          stats = wdlps_v_fat_rank)

saveRDS(wdlps_v_fat_gsea, file = "outputs/human_wdlps_v_fat_gsea.rds")
write.csv(wdlps_v_fat_gsea[,-8], file = "outputs/human_wdlps_v_fat_gsea.csv", quote = FALSE, row.names = F)

# DDLPS vs fat ------------------------------------------------------------

ddlps_idx <- dds$condition %in% c("Normal_fat","DDLPS")
dds_ddlps <- dds[,ddlps_idx]
dds_ddlps$condition <- factor(dds_ddlps$condition, levels = c("Normal_fat","DDLPS"))
design(dds_ddlps) <- model.matrix(~ condition, data = colData(dds_ddlps))
# dds_wdlps <- estimateSizeFactors(dds_wdlps)
# idx <- rowSums(counts(dds_wdlps, normalized=TRUE) >= 10 ) >= 3
# dds_wdlps <- dds_wdlps[idx,]
dds_ddlps <- DESeq(dds_ddlps)

# > DGE analysis ----------------------------------------------------------

ddlps_v_fat <- lfcShrink(dds_ddlps, coef = 'conditionDDLPS', type = 'ashr') %>%
  data.frame() %>%
  merge(gene_symbols, by = 0) %>%
  filter(!is.na(padj)) %>%
  rename(gene_id = Row.names) %>%
  mutate(status = ifelse(log2FoldChange > 0.5 & padj < 0.05, 1,
                         ifelse(log2FoldChange < -0.5 & padj < 0.05, -1, 0)))

write.csv(ddlps_v_fat, "outputs/human_ddlps_v_fat_dge.csv", row.names = F, quote = F)

human_ann <- data.frame(gene = ddlps_v_fat$gene_id,
                        gene_symbol = ddlps_v_fat$gene_symbol,
                        pval_adj = ddlps_v_fat$padj)
human_ddlps_v_fat <- glimmaXY(x = ddlps_v_fat$log2FoldChange,
                              xlab = 'log2FoldChange',
                              y = -log(ddlps_v_fat$padj),
                              ylab = 'significance',
                              status = ddlps_v_fat$status,
                              counts = vst(counts(dds_ddlps))[ddlps_v_fat$gene_id,],
                              transform.counts = 'none',
                              main = 'DDLPS vs Normal Fat',
                              groups = colData(dds_ddlps)$condition,
                              anno = human_ann,
                              display.columns = c('gene_symbol','pval_adj'))

htmlwidgets::saveWidget(human_ddlps_v_fat,'docs/human_ddlps_v_fat_dge_interactive_app.html')

# > GSEA ------------------------------------------------------------------

ddlps_v_fat_rank <- ddlps_v_fat %>%
  pull(log2FoldChange) %>%
  `names<-`(ddlps_v_fat$gene_id) %>%
  sort()

ddlps_v_fat_gsea <- fgsea(pathways = hallmarks_pathways,
                          stats = ddlps_v_fat_rank)

saveRDS(ddlps_v_fat_gsea, file = "outputs/human_ddlps_v_fat_gsea.rds")
write.csv(ddlps_v_fat_gsea[,-8], file = "outputs/human_ddlps_v_fat_gsea.csv", quote = FALSE, row.names = F)

# setup -------------------------------------------------------------------

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(PCAtools)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(mMCPcounter)
  library(MCPcounter)
  library(pheatmap)
  library(nichenetr)
  library(GSVA)
  library(cowplot)
  library("AnnotationDbi")
  library("org.Hs.eg.db")
})

dir.create("outputs", showWarnings = F, recursive = T)

# loading data ------------------------------------------------------------

human_dds <- readRDS("outputs/all_human_samples_dds.rds")
sample_info <- read.table("data/human_rnaseq_samples_info.txt", sep = '\t', quote = "", header = T) %>%
  column_to_rownames("sample_id")
sample_info <- sample_info[colnames(human_dds),]
sample_info$TILS <- factor(sample_info$TILS, levels = c('low', 'mod', 'high', 'NF'))

MCP_genes <- read.table(curl:::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),
                   sep = "\t", stringsAsFactors = FALSE, header = TRUE,
                   colClasses = "character", check.names = FALSE) %>%
  dplyr::select(c("Cell population", "ENSEMBL ID")) %>%
  split(f = .$`Cell population`)
MCP_genes <- lapply(MCP_genes, function(x) x$`ENSEMBL ID`)

## https://www.nature.com/articles/s41586-019-1906-8

Petitprez_et_al_genesets <- list('T_cells' = c('CD28', 'CD3D', 'CD3G', 'CD5', 'CD6', 'CHRM3-AS2', 'CTLA4', 'FLT3LG', 'ICOS', 'MAL', 'PBX4', 'SIRPG', 'THEMIS', 'TNFRSF25', 'TRAT1'),
                                 'cytotoxic lymphocytes' = c('CD8A', 'EOMES', 'FGFBP2', 'GNLY', 'KLRC3', 'KLRC4', 'KLRD1'),
                                 'B lineage' = c('BANK1','CD19', 'CD22', 'CD79A', 'CR2', 'FCRL2', 'IGKC', 'MS4A1', 'PAX5'),
                                 'natural killer cells' = c('CD160', 'KIR2DL1', 'KIR2DL3', 'KIR2DL4', 'KIR3DL1', 'KIR3DS1', 'NCR1', 'PTGDR', 'SH2D1B'),
                                 'monocytic lineage' = c('ADAP2', 'CSF1R', 'FPR3', 'KYNU', 'PLA2G7', 'RASSF4', 'TFEC'),
                                 'myeloid dendritic cells' = c('CD1A', 'CD1B', 'CD1E', 'CLEC10A', 'CLIC2', 'WFDC21P'),
                                 'neutrophils' = c('CA4', 'CEACAM3', 'CXCR1', 'CXCR2', 'CYP4F3', 'FCGR3B', 'HAL', 'KCNJ15', 'MEGF9', 'SLC25A37', 'STEAP4', 'TECPR2', 'TLE3', 'TNFRSF10C', 'VNN3'),
                                 'endothelial cells' = c('ACVRL1', 'APLN', 'BCL6B', 'BMP6', 'BMX', 'CDH5', 'CLEC14A', 'CXorf36', 'EDN1', 'ELTD1','EMCN', 'ESAM', 'ESM1', 'FAM124B', 'HECW2', 'HHIP', 'KDR', 'MMRN1', 'MMRN2', 'MYCT1', 'PALMD','PEAR1', 'PGF', 'PLXNA2', 'PTPRB', 'ROBO4', 'SDPR', 'SHANK3', 'SHE', 'TEK', 'TIE1', 'VEPH1', 'VWF'))

Petitprez_et_al_genesets <- lapply(Petitprez_et_al_genesets, function(x){
  genes <- mapIds(org.Hs.eg.db,
                  keys=x,
                  column="ENSEMBL",
                  keytype="SYMBOL",
                  multiVals="first")

  return(unlist(genes))
})


# GSVA Enrichment ---------------------------------------------------------

gaussian_gsva_enrich <- gsva(gsvaParam(log(cpm(counts(human_dds)) + 1), MCP_genes, kcdf = 'Gaussian'))
color_palette <- colorRampPalette(c('blue','white','red'))
pheatmap(gaussian_gsva_enrich[, order(sample_info$TILS)],
         color = color_palette(100),
         annotation_col = sample_info,
          gaps_col = c(9, 14, 16),
         cluster_rows = T, cluster_cols = F,
         main = 'GSVA Enrichemt (using MCP Signatures)')


gaussian_gsva_enrich <- gsva(gsvaParam(log(cpm(counts(human_dds)) + 1), Petitprez_et_al_genesets, kcdf = 'Gaussian'))
color_palette <- colorRampPalette(c('blue','white','red'))
pheatmap(gaussian_gsva_enrich[, order(sample_info$TILS)],
         color = color_palette(100),
         annotation_col = sample_info,
         gaps_col = c(9, 14, 16),
         cluster_rows = T, cluster_cols = F,
         main = 'GSVA Enrichemt (using Petitprez et al Signatures)')



```