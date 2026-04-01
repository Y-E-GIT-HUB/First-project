# README: m6A‑SNPs identified by integrating genomic data are associated with the occurrence and prognosis of HCC

This repository contains the complete R scripts and input data used to perform differential expression analyses for the manuscript above. The analyses include:

- **TCGA‑LIHC** (RNA‑seq)
- **GSE45114** (microarray)
- **GSE60502** and **GSE36376** (two additional microarray datasets, processed with the same workflow as GSE45114)

All scripts are fully self‑contained and, after adjusting file paths, can be run end‑to‑end to reproduce the results presented in the paper.

---

## Repository contents

```
The repository is organized into separate folders for each dataset. Each folder contains its own input files, analysis script, and output results.

.
├── TCGA-LIHC/
│   ├── input/
│   │   ├── TCGA-LIHC_MergeCOUNT.txt   # raw count matrix
│   │   └── id_symble_type.Rda         # gene ID mapping
│   ├── scripts/
│   │   └── TCGA_DEG_analysis.R        # DESeq2 analysis script
│   └── output/
│       └── TCGA-Differentially expressed genes.csv   # results (pre‑generated)
│
├── GSE45114/
│   ├── input/
│   │   ├── GPL5918.txt                # platform annotation
│   │   ├── GSE45114_series_matrix.txt # expression matrix
│   │   └── sample_Type.txt            # sample grouping
│   ├── scripts/
│   │   └── GSE45114_DEG_analysis.R    # limma analysis script
│   └── output/
│       └── GSE45114-Differentially expressed genes.csv
│
├── GSE60502/      # second GEO dataset (similar structure)
├── GSE36376/      # third GEO dataset (similar structure)
│
└── README.md      # this file
```

*Note:* The exact GEO accession numbers for the second and third datasets should be filled in according to your study.

---

## Requirements

- **R** version 4.2 or higher (tested with 4.3.0)

### Required R packages

The analyses use both CRAN and Bioconductor packages. The scripts automatically check for and install missing packages, but you may also pre‑install them using the commands below.

| Package        | Source       | Used in                           |
|----------------|--------------|-----------------------------------|
| DESeq2         | Bioconductor | TCGA analysis                     |
| edgeR          | Bioconductor | TCGA analysis (dependency)        |
| datTraits      | Bioconductor | TCGA analysis (dependency)        |
| ggplotify      | Bioconductor | TCGA analysis (graphic conversion)|
| patchwork      | CRAN         | TCGA analysis                     |
| cowplot        | CRAN         | TCGA analysis                     |
| stringr        | CRAN         | Both                              |
| tidyr          | CRAN         | Both                              |
| dplyr          | CRAN         | Both                              |
| data.table     | CRAN         | GEO analyses (fast file reading)  |
| pheatmap       | CRAN         | Both                              |
| RColorBrewer   | CRAN         | Both                              |
| ggplot2        | CRAN         | Both                              |
| ggrepel        | CRAN         | Both                              |
| limma          | Bioconductor | GEO analyses                      |

#### Installation commands


**CRAN packages:** 
```r
install.packages(c("patchwork", "cowplot", "stringr", "tidyr", "dplyr", "data.table",
                   "pheatmap", "RColorBrewer", "ggplot2", "ggrepel"))
```

**Bioconductor packages:**
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "datTraits", "ggplotify", "limma"))
```

**1. TCGA‑LIHC analysis**
Navigate to the TCGA‑LIHC folder and open scripts/TCGA_DEG_analysis.R.

Replace the file paths with the relative paths inside the repository:

read.table("C:/TCGA-LIHC_MergeCOUNT.txt") → read.table("input/TCGA-LIHC_MergeCOUNT.txt")

load("C:/id_symble_type.Rda") → load("input/id_symble_type.Rda")

write.csv(..., "C:/TCGA-Differentially expressed genes.csv") → write.csv(..., "output/TCGA-Differentially expressed genes.csv")

Run the script from R (set working directory to the TCGA‑LIHC folder first):

```r
source("scripts/TCGA_DEG_analysis.R")
```

The script will perform DESeq2 differential expression (tumor vs. normal), save the results as a CSV file, and generate a heatmap and volcano plot in the graphics window.

2. GSE45114 analysis (and the two additional GEO datasets)
The same workflow (using limma) is applied to each GEO dataset. For each dataset, follow these steps:

Create a folder for the dataset (e.g., GSE45114/, GSEXXXXX/, GSEYYYYY/) and place the corresponding input files inside the input/ subfolder.

Copy the script template (GSE45114_DEG_analysis.R) into the scripts/ subfolder and rename it appropriately.

Modify the file paths inside the script to use relative paths:

For GSE45114:

data.table::fread("C:/GSE45114/GPL5918.txt") → data.table::fread("input/GPL5918.txt")

read.table("C:/GSE45114/GSE45114_series_matrix.txt", ...) → read.table("input/GSE45114_series_matrix.txt", ...)

read.table("C:/GSE45114/sample_Type.txt", ...) → read.table("input/sample_Type.txt", ...)

write.csv(..., "C:/GSE45114/GSE45114-Differentially expressed genes.csv") → write.csv(..., "output/GSE45114-Differentially expressed genes.csv")

For the other GEO datasets, adjust the file names accordingly.

Run the script from R (set working directory to the dataset’s root folder):

```r
source("scripts/GSE45114_DEG_analysis.R")
```

The script will perform limma analysis, save the differentially expressed genes (adjusted p < 0.05 and |logFC| ≥ 1) as a CSV file, and display a volcano plot and heatmap.

Important notes for all GEO analyses
The script includes a log‑transformation check that automatically applies log2 transformation if the data are not already log‑normalized. This is essential for microarrays.

The volcano plot and heatmap are shown in the R graphics window but are not automatically saved. You may save them manually (e.g., with ggsave()) if you wish to keep them.

The three GEO datasets are processed identically; only the input files differ. The provided script for GSE45114 serves as a template for all GEO analyses.

Expected outputs
After successfully running all scripts, the output/ folders will contain:

TCGA-Differentially expressed genes.csv – DESeq2 results (log2FoldChange, pvalue, padj, change status)

GSE45114-Differentially expressed genes.csv – limma results (logFC, P.Value, adj.P.Val, etc.)

Similar CSV files for the other two GEO datasets.

These CSV files contain the exact results used in the manuscript. You can compare them with the pre‑generated files in the repository to verify reproducibility.

The scripts also produce figures (volcano plots and heatmaps) that match those shown in the paper.

Additional notes
Inline comments have been added to all scripts to explain each major step. Please refer to them for methodological details.

All analyses were performed on a standard desktop/laptop; no high‑performance computing cluster is required.

The input files are provided exactly as they were used in the original analysis. If you run the scripts on a different operating system, please ensure that file paths are correctly adjusted (use forward slashes / or file.path() for cross‑platform compatibility).
