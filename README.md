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

## Analysis & Shiny Application
1. Sample Tab:
This tab uses readr, dplyr, tidyr, and DT to read and clean metadata, generate structured summaries, and display interactive tables. I coded plotting using ggplot2, allowing the user to generate histograms, density plots, and violin plots by selecting any numeric variable and grouping column. The tab reacts instantly because all summaries and plots are built with reactive expressions in Shiny.

3. Counts Tab:
Normalized counts are processed using dplyr, tidyr, and base R functions to compute gene-wise variance and zero counts for filtering. Diagnostic visualizations are coded using ggplot2, while clustered heatmaps are generated through pheatmap. PCA is computed with prcomp, and I used both ggplot2 and ggbeeswarm to provide scatter and beeswarm PCA plots. These components use reactive filtering so users see updated heatmaps and PCA results based on their settings.

4. Differential Expression Tab:
I used DESeq2 results that were precomputed in R and then processed in the Shiny app using dplyr, ggplot2, and DT. The volcano plot, p-value histogram, and log2FC histogram were coded in ggplot2, with dynamic filtering based on user-selected padj thresholds. A jitter plot highlights top genes by plotting their normalized counts (log-transformed) across control and disease groups, providing a clear visualization of differential signal.

5. GSEA Tab:
GSEA was performed using fgsea and gene sets from msigdbr. I wrote code to rank genes, run enrichment with Hallmark pathways, flatten leading-edge lists, and save cleaned GSEA outputs. The Shiny tab displays these results using ggplot2 for barplots and NES–p-value scatter plots, and DT for sortable pathway tables. Users can filter pathways by NES direction or padj, and download the processed subset directly from the app.
