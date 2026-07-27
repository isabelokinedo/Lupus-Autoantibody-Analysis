# Family-Level Spearman Correlation and Bipartite Network Analysis

This repository contains the R scripts used to perform family-level correlation analyses between autoantibody reactivity and inflammatory protein expression measured by the Olink® Target 96 Inflammation panel.

Pairwise Spearman correlations were calculated between individual autoantibodies and Olink proteins. Autoantibodies and proteins were then grouped into biologically relevant families to summarize correlation patterns and visualize interactions using bipartite network graphs.

---

## Repository Contents

| File | Description |
|------|-------------|
| `1. Autoantibody x Olink NPX Spearman.R` | Calculates pairwise Spearman correlation coefficients, p values, and Benjamini–Hochberg false discovery rate (FDR) adjusted p values between autoantibody scores and Olink NPX protein values. Generates correlation matrices and a master results table. |
| `2. Family summaries batched.R` | Aggregates individual antibody–protein correlations into family-level summaries using antibody and protein annotation tables. Calculates summary statistics for each antibody family–protein family interaction. |
| `3. Bipartite plots.R` | Generates family-level bipartite network visualizations from the summarized correlation data. Nodes represent antibody and protein families, while edges represent retained family-level correlations. |
| `4. Spearman rho matrix dotplots.R` | Generates dot-plot heatmaps of individual Spearman correlation coefficients using the correlation matrix generated in the first step. Provides an alternative visualization of individual antibody–protein correlations and is independent of the bipartite network analysis. |
| `Autoantibody_annotations.csv` | Annotation table assigning each autoantibody to a biological family and subfamily. |
| `Olink_protein_annotations.csv` | Annotation table assigning each Olink analyte to a biological protein family and subfamily. |
| `Olink NPX Data BF Plasma CSF.csv` | Example input dataset containing Olink NPX values used for correlation analysis. |

---

## Software Requirements

- R (version 4.4 or later recommended)
- Required packages include:

- tidyverse
- readr
- dplyr
- tidyr
- stringr
- purrr
- igraph
- tidygraph
- ggraph
- ggrepel

Additional packages may be required depending on the plotting options used.

---

## Analysis Workflow

Run the scripts in the following order:

### 1. Calculate pairwise correlations

```
Autoantibody x Olink NPX Spearman.R
```

Input:

- Olink NPX dataset
- R-ready autoantibody dataset

Output includes:

- Spearman rho matrix
- Raw p-value matrix
- FDR-adjusted p-value matrix
- Correlation master table

---

### 2. Generate family-level summaries

```
Family summaries batched.R
```

Input:

- Correlation master tables
- Autoantibody annotation table
- Olink protein annotation table

Output includes:

- Family summary tables
- Network family summary tables

---

### 3. Generate bipartite networks

```
Bipartite plots.R
```

Input:

- Network family summary tables

Output:

- Publication-quality bipartite network figures
- Retained edge tables
- Network summary tables

---

### Optional

```
Spearman rho matrix dotplots.R
```

Produces dot-plot visualizations of individual antibody–protein Spearman correlations.

---

## Input Data

This repository includes example annotation tables and an example Olink NPX dataset.

To reproduce the analyses with new data, users should provide:

- Olink NPX data
- R-ready autoantibody score data
- Autoantibody annotation table
- Olink protein annotation table

---

## Statistical Methods

Pairwise associations between autoantibody scores and Olink NPX protein values were evaluated using Spearman rank correlation coefficients (ρ).

Raw p values were adjusted for multiple comparisons using the Benjamini–Hochberg false discovery rate (FDR).

For family-level analyses, individual antibody–protein correlations were grouped according to biologically defined antibody and protein families. Summary statistics describing each family-family interaction were calculated and used to construct bipartite correlation networks.

---

## Citation

If you use this code, please cite the associated manuscript describing the methodology.
