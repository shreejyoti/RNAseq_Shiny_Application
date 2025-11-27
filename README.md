# RNAseq Shiny Application  
### Differential Expression • GSEA • Interactive Dashboard

![R](https://img.shields.io/badge/R-4.2.0-blue)
![Shiny](https://img.shields.io/badge/Shiny-Application-orange)
![Status](https://img.shields.io/badge/Status-Completed-success)

## Project Overview 
This project is an end-to-end RNA-seq analysis workflow built in R and deployed as an interactive Shiny application
The goal of this project is to take a publicly available transcriptomic dataset, perform differential gene expression analysis (DGEA), run gene set enrichment analysis (GSEA), and provide an accessible dashboard for exploring results

## Dataset
The data comes from GEO accession GSE64810, a well-known transcriptomic dataset containing normalized FPKM counts and extensive clinical/biological metadata.
The dataset includes:
- Metadata (clinical variables) such as diagnosis, age of death, RIN, disease onset, Vonsattel grade, striatal score, cortical score, and more
- Normalized FPKM gene expression values across multiple brain tissue samples
- Differential expression results generated using DESeq2
- Pathway enrichment results generated using FGSEA with MSigDB Hallmark pathways
- You can access full dataset directly from GEO:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810

## Pre-Processing & Statistical Analysis
- All core analysis steps were performed in R before building the Shiny app, using a standard RNA-seq workflow. Metadata was cleaned using dplyr to correct formatting, unify clinical fields, and prepare consistent sample annotations. Normalized FPKM counts were reshaped and filtered using variance and zero-count thresholds to reduce noise and focus on biologically meaningful genes.
- Differential expression was carried out using DESeq2, which models RNA-seq counts with a negative binomial framework to estimate stable log2 fold changes and significance values. The ranked genes were then used as input for FGSEA, where Hallmark gene sets from msigdbr were evaluated across 1,000 permutations. 
- This preprocessing pipeline produced all cleaned metadata, filtered expression matrices, DE results, and pathway enrichment outputs—which the Shiny application then renders interactively for exploration.

## Shiny Application
1. Sample Tab:
This tab uses readr, dplyr, tidyr, and DT to read and clean metadata, generate structured summaries, and display interactive tables. I coded plotting using ggplot2, allowing the user to generate histograms, density plots, and violin plots by selecting any numeric variable and grouping column. The tab reacts instantly because all summaries and plots are built with reactive expressions in Shiny.

3. Counts Tab:
Normalized counts are processed using dplyr, tidyr, and base R functions to compute gene-wise variance and zero counts for filtering. Diagnostic visualizations are coded using ggplot2, while clustered heatmaps are generated through pheatmap. PCA is computed with prcomp, and I used both ggplot2 and ggbeeswarm to provide scatter and beeswarm PCA plots. These components use reactive filtering so users see updated heatmaps and PCA results based on their settings.

4. Differential Expression Tab:
I used DESeq2 results that were precomputed in R and then processed in the Shiny app using dplyr, ggplot2, and DT. The volcano plot, p-value histogram, and log2FC histogram were coded in ggplot2, with dynamic filtering based on user-selected padj thresholds. A jitter plot highlights top genes by plotting their normalized counts (log-transformed) across control and disease groups, providing a clear visualization of differential signal.

5. GSEA Tab:
GSEA was performed using fgsea and gene sets from msigdbr. I wrote code to rank genes, run enrichment with Hallmark pathways, flatten leading-edge lists, and save cleaned GSEA outputs. The Shiny tab displays these results using ggplot2 for barplots and NES–p-value scatter plots, and DT for sortable pathway tables. Users can filter pathways by NES direction or padj, and download the processed subset directly from the app.

You can find the Shiny app code inside Code folder named as updated_version_app.R

## Results

The final analysis revealed a set of clear transcriptional shifts between disease and control groups. Several genes—including PPP1R1B, HTR2A, and RGS4—showed strong differential expression patterns, supported by distinct log2 fold-change distributions and visually separated jitter plots. After filtering, PCA demonstrated clearer clustering of samples, indicating improved structure in the dataset.

Pathway enrichment further highlighted biologically meaningful signatures. Hallmark pathways such as Inflammatory Response, TNFα Signaling via NF-κB, Apoptosis, and Oxidative Phosphorylation appeared significantly enriched, reflecting expected stress and neurodegenerative processes in this dataset.

## How to Run the App

You can run this Shiny application locally by cloning this repository or downloading all files as a folder. Once the files are on your computer, open the project in RStudio and set your working directory to the main project folder.

Before running the app, install all required packages:

```r
install.packages(c(
  "shiny", "tidyverse", "DT", "ggplot2", "ggbeeswarm"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery", "fgsea", "msigdbr"
))
```

Load the libraries (optional but recommended to verify installation):

```r
library(shiny)
library(tidyverse)
library(DT)
library(ggplot2)
library(ggbeeswarm)
library(GEOquery)
library(fgsea)
library(msigdbr)
```

Finally, run the application using:
```r
shiny::runApp("Code/updated_version_app.R")
```

This will launch the full interactive Shiny interface, where you can upload metadata, counts, DESeq2 results, and fgsea outputs and explore the analyses.
