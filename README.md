# Integrated-Differential-Gene-Expression-DEG-GO-Enrichment-Analysis-of-Microarray-Data-GDS2778
R Script for Integrated Differential Gene Expression DEG and Gene Ontology GO Enrichment Analysis of Microarray Data GDS2778
# 🧬 Comprehensive Microarray Analysis: DEG and GO Enrichment (GDS2778)

This repository contains a robust R script that executes a complete bioinformatics workflow for analyzing publicly available microarray data from the GEO database. The script uses the **GDS2778** dataset (Treatment vs. Control) as an example to demonstrate best practices for differential gene expression (DEG) and subsequent functional enrichment analysis.

The workflow leverages Bioconductor packages (`GEOquery`, `limma`, `clusterProfiler`) to provide statistical rigor and rich biological context.

---

## 🚀 Getting Started

### Prerequisites

You must have **R** and an R development environment (like RStudio) installed. The script automatically checks for and installs the necessary Bioconductor and CRAN packages, but you must have Bioconductor installed initially.

The core packages required are:

* `GEOquery`
* `limma`
* `clusterProfiler`
* `org.Hs.eg.db` (for Human annotation)
* `pheatmap`, `RColorBrewer`

### Installation and Setup

1.  **Clone the repository:**
    ```bash
    git clone [Your Repository URL]
    cd [Your Repository Folder]
    ```

2.  **Set Working Directory:**
    **Before running**, you must customize the `setwd()` path in the R script to a directory on your machine where you want all output plots and data tables to be saved.

    ```R
    # --- 1. SETUP ENVIRONMENT ---
    setwd("D:/DOWNLOADS") # <--- REPLACE THIS WITH YOUR DESIRED PATH
    ```

3.  **Run the script:**
    Execute the R script in R/RStudio:

    ```R
    source("R Script for Integrated Differential Gene Expression DEG and Gene Ontology GO Enrichment Analysis of Microarray Data GDS2778.R")
    ```

---

## ⚙️ Workflow and Output Files

The script performs a sequential four-part analysis, generating multiple output files:

### 1. Differential Gene Expression (DEG) Analysis

The `limma` package is used to fit a linear model and calculate moderated t-statistics, comparing the **Control** group to the **Treatment** group as defined by the dataset's phenotype data (`targets$agent`).

| Output File | Description |
| :--- | :--- |
| `GDS2778_DEG_results.csv` | **Full DEG Table:** Contains log-fold change (`logFC`), average expression (`AveExpr`), t-statistic, raw p-value, and Benjamini-Hochberg (BH) adjusted p-value (`adj.P.Val`) for all genes. |
| `GDS2778_Filtered_DEGs.csv` | **Significant DEGs:** A subset of the full DEG table, filtered by strict cutoffs. |

### 2. DEG Filtering Criteria

The script applies the following conventional thresholds to identify biologically relevant genes:

* **Statistical Significance:** Adjusted P-value (`adj.P.Val`) $\le 0.01$
* **Fold Change:** Absolute Log2 Fold Change (`|logFC|`) $\ge 1$ (equivalent to a 2-fold change).

### 3. Gene Ontology (GO) Enrichment Analysis

The `clusterProfiler` package is used to determine which biological functions or pathways are overrepresented among the significantly regulated genes (using a P-value cutoff of 0.05 for filtering the input gene list).

| Output File | Analysis Type | Description |
| :--- | :--- | :--- |
| `GDS2778_GO_Enrichment_MF_Barplot.png` | **Molecular Function (MF)** | Bar plot of the top 10 enriched MF terms. |
| `GDS2778_GO_Enrichment_MF_DAG.png` | **Molecular Function (MF)** | Directed Acyclic Graph (DAG) plot visualizing the hierarchical relationship between the enriched MF terms. |
| `GDS2778_GO_Enrichment_BP_Barplot.png` | **Biological Process (BP)** | Bar plot of the top 10 enriched BP terms. |

### 4. Clustering Heatmap

A heatmap is generated using the expression values of all genes that meet the statistical significance cutoff (`adj.P.Val <= 0.05`).

| Output File | Analysis Type | Description |
| :--- | :--- | :--- |
| `GDS2778_DEG_Heatmap.png` | **Clustering Heatmap** | Visualizes the relative expression levels of significant genes across all samples, clustered by gene and sample similarity. |

---

## 🔗 Dataset Information

* **GEO Accession:** **GDS2778**
* **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array (GPL570).
* **Context:** This dataset is a Genome Database (GDS), which is typically curated and provides a well-defined context, making it suitable for comparative analysis.
