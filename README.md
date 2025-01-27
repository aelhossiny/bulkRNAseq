Bulk RNASeq Analysis Pipeline
================

# Table of content

- Introduction
- Setup
- Exploratory Analysis
- Analysis
- Notes
- Recommended Readings

# Introduction

This is an analysis pipeline for Bulk RNAseq data. We will re-analyze
the mouse RNAseq data from [Donahue et
al](https://aacrjournals.org/cancerdiscovery/article/14/10/1964/748588/Oncogenic-KRAS-Dependent-Stromal-Interleukin-33)
it can be downloaded from
[here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269889)

This tutorial will not include essentials for RNAseq statistics. Your
can read more from the following resources:

- [DESeq2
  Vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [STAT 555 PSU course](https://online.stat.psu.edu/stat555/node/13/)
- [Design Matrix
  Guidlines](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)

# Setup

First we will load essential libraries and import the data. Note that we
will need to install the following packages (DESeq2, PCAtools,
EnhancedVolcano, fgsea, msigdbr, decoupleR)

The goal of this section is to construct a numeric matrix where the rows
represent gene names/ids and the columns represent samples and metadata
dataframe to create DESeq2 object.

``` r
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(PCAtools)
  library(fgsea)
  library(msigdbr)
  library(EnhancedVolcano)
  library(decoupleR)
})

## Create output directory
dir.create("outputs", showWarnings = F, recursive = T)

## loading count data
## NOTES: here we will store the gene symbols in a dataframe and keep the gene id as the rowname
## This is because in the count matrix that we have, some genes ids have the same gene symbols
## so we will keep the gene id as rownames to avoid duplicates and add the gene symbols in the end
## for visualization.
ct_mtx <- read.table("data/gene_expected_count.annot.txt", quote = "", header = T, fill = T, sep = "\t")
gene_symbols <- data.frame(gene_symbol = ct_mtx$external_gene_name) %>% `rownames<-`(ct_mtx$gene_id)
ct_mtx <- ct_mtx %>%
  select(-c("entrezgene_id", "external_gene_name", "description")) %>%
  column_to_rownames('gene_id')

## Remove the X in the begining of each sample name and make them match the sample names in metadata file
colnames(ct_mtx) <- gsub("^X", "", colnames(ct_mtx))
colnames(ct_mtx) <- gsub("\\.", "-", colnames(ct_mtx))


## loading metadata
metadata <- read.csv("data/metadata.csv", header = T) %>%
  column_to_rownames('sample')
metadata <- metadata[colnames(ct_mtx),]

# DESeq object all samples ------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = ct_mtx,
                              colData = metadata,
                              rowData = gene_symbols,
                              design = ~ 1)
```

# Exploratory Analysis

## PCA Analysis

Here we are visualizing the PCA plots for the samples. You can change
the shape and color based on the `color_by` and `shape_by` variables

``` r
vst <- assay(vst(dds))
rv <- rowVars(vst)
select <- order(rv, decreasing=TRUE)[seq_len(min(2000, length(dds)))]
pca <- pca(vst[select,], metadata = colData(dds))
color_by <- 'IL33.Status'
shape_by <- 'Cell.Line'
biplot(pca,
       showLoadings = F,
       colby = color_by,
       shape = shape_by,
       labSize = 2, pointSize = 5, drawConnectors = T,
       legendPosition = 'right',
       legendLabSize = 7, legendIconSize = 4, legendTitleSize = 10)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
