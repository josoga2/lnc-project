# lncRNA Differential Expression and Correlation Analysis Pipeline

This repository contains a reusable DESeq2-based workflow for differential expression analysis of long non-coding RNAs (lncRNAs) and downstream correlation analysis. The pipeline was developed for hepatocellular carcinoma (HCC) transcriptomic datasets and supports subgroup-specific comparisons based on treatment response and disease stage.

The workflow combines differential expression analysis, quality control, lncRNA prioritization, and co-expression network discovery into a single framework.

---

## Repository Structure

```text
.
├── cleanup_lnc.R
├── deseq2_subset_pipeline.R
├── counts.csv
├── Metadata.xlsx
├── mart_export_all_lnc.txt
└── cleanup_pass/
    └── qc/
```

| File | Description |
|--------|-------------|
| `cleanup_lnc.R` | Main analysis script containing all biological comparisons. |
| `deseq2_subset_pipeline.R` | Core DESeq2 and correlation analysis functions. |
| `counts.csv` | Raw gene count matrix. |
| `Metadata.xlsx` | Sample metadata. |
| `mart_export_all_lnc.txt` | Ensembl-derived lncRNA annotation file. |

---

## Features

- Differential expression analysis using DESeq2
- Automatic sample subsetting based on metadata
- lncRNA-specific filtering and annotation
- Variance-stabilized expression transformation
- PCA visualization
- Volcano plots
- MA plots
- Sample distance heatmaps
- Top DEG heatmaps
- lncRNA-mRNA correlation network generation
- Flexible metadata-driven comparisons

---

## Installation

### CRAN Packages

```r
install.packages(c(
  "pheatmap",
  "ggplot2",
  "ggrepel",
  "matrixStats",
  "RColorBrewer",
  "readxl"
))
```

### Bioconductor Packages

```r
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "apeglm"
))
```

---

## Input Files

### Count Matrix

The count matrix should contain raw counts with genes in rows and samples in columns.

Example:

```csv
geneId,HCC001,HCC002,HCC003
ENSG000001,100,120,98
ENSG000002,5,7,3
```

### Metadata

The metadata file must contain at least the following columns:

| Column | Description |
|----------|-------------|
| `SampleID` | Sample identifier matching count matrix columns |
| `Condition` | PRE or POST treatment |
| `Group` | Responder or Non-responder |
| `BCLC_Stage` | Stage B or Stage C |

Example:

| SampleID | Condition | Group | BCLC_Stage |
|-----------|------------|---------|------------|
| HCC001 | PRE | Responder | Stage B |
| HCC002 | POST | Non-responder | Stage C |

---

## Running Differential Expression Analysis

Load the pipeline:

```r
source("deseq2_subset_pipeline.R")
```

Run a comparison:

```r
results <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Condition = "POST"
  ),
  contrast_col = "Group",
  contrast_levels = c(
    "Responder",
    "Non-responder"
  ),
  outdir = "results"
)
```

---

## Comparisons Performed in `cleanup_lnc.R`

### 1. Stage B vs Stage C Before Treatment

```r
Condition = PRE
Stage B vs Stage C
```

### 2. Stage B vs Stage C After Treatment

```r
Condition = POST
Stage B vs Stage C
```

### 3. Responder vs Non-responder Before Treatment

```r
Condition = PRE
Responder vs Non-responder
```

### 4. Responder vs Non-responder After Treatment

```r
Condition = POST
Responder vs Non-responder
```

### 5. Responder vs Non-responder Within Stage B

```r
Condition = POST
BCLC_Stage = Stage B
```

### 6. Responder vs Non-responder Within Stage C

```r
Condition = POST
BCLC_Stage = Stage C
```

### 7. Treatment Effect in Responders

```r
POST vs PRE
Group = Responder
```

---

## Output Objects

The function returns a list containing:

| Object | Description |
|----------|-------------|
| `UP` | Top upregulated genes |
| `DOWN` | Top downregulated genes |
| `ALL` | Ordered DEG results |
| `FULL_DEGS` | Complete DESeq2 output |
| `UP_ACT` | Significant upregulated genes |
| `DOWN_ACT` | Significant downregulated genes |
| `FOR_CORR` | Variance-stabilized expression matrix |
| `upLNC` | Upregulated lncRNAs |
| `downLNC` | Downregulated lncRNAs |

Example:

```r
results$UP
results$DOWN
results$FULL_DEGS
```

---

## Quality Control Outputs

For each analysis, the pipeline automatically generates quality control plots under:

```text
results/qc/
```

Generated visualizations include:

| Plot | Description |
|---------|-------------|
| PCA Plot | Sample clustering |
| Volcano Plot | Differential expression overview |
| MA Plot | Expression change visualization |
| Dispersion Plot | DESeq2 model diagnostics |
| Sample Distance Heatmap | Sample similarity assessment |
| Top 50 Upregulated Heatmap | Expression patterns of highly induced genes |
| Top 50 Downregulated Heatmap | Expression patterns of highly repressed genes |

---

## lncRNA Correlation Analysis

The pipeline includes a correlation module for identifying highly correlated lncRNA-mRNA pairs.

### Example

```r
corr <- makeCorrelation(
  direction = results$UP,
  for_corr = results$FOR_CORR,
  cutoff = 0.9,
  n_samp_min = 7
)
```

### Filtering Strategy

Genes are retained only if:

- Expression ≥ 10 in at least `n_samp_min` samples
- Variance > 0.75
- Absolute Pearson correlation exceeds the chosen cutoff

### Output

```r
head(corr)
```

| gene1 | gene2 | correlation |
|--------|--------|-------------|
| ENSG... | ENSG... | 0.95 |

Only lncRNA-protein coding gene pairs are retained.

---

## Example Workflow

```r
source("deseq2_subset_pipeline.R")

results <- run_deseq_subset(
  counts_file = "counts.csv",
  metadata_file = "Metadata.xlsx",
  subset_filters = list(
    Condition = "POST"
  ),
  contrast_col = "Group",
  contrast_levels = c(
    "Responder",
    "Non-responder"
  ),
  outdir = "results"
)

corr <- makeCorrelation(
  direction = results$UP,
  for_corr = results$FOR_CORR,
  cutoff = 0.9,
  n_samp_min = 7
)
```

---

## Biological Applications

This workflow can be used to identify:

- Treatment response-associated lncRNAs
- Stage-specific lncRNA signatures
- Treatment-induced transcriptional changes
- Co-expression networks linking lncRNAs to coding genes
- Candidate regulatory lncRNAs involved in chromatin regulation, transcriptional control, RNA processing, and signaling pathways

---

## Notes

- Sample IDs are automatically cleaned and standardized.
- Common metadata naming inconsistencies are automatically repaired.
- lncRNA annotation is performed using the supplied Ensembl annotation file.
- Variance-stabilized expression values are used for downstream correlation analysis.
- The pipeline is optimized for RNA-seq count data but can be adapted to similar count-based transcriptomic datasets.

---

## Citation

If this pipeline contributes to a publication, please cite:

- Love MI, Huber W, Anders S. DESeq2.
- Zhu A, Ibrahim JG, Love MI. apeglm.
- Ensembl Gene Annotation Database.
