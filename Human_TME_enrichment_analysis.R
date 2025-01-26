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


