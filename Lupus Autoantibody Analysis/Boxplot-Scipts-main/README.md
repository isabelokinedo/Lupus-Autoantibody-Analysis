# Lupus Autoantibody Boxplot Plot Analysis

This repository contains the R script used to analyze antibody-specific boxplots and statistical output tables frp, R-ready autoantibody data:
  CLE compartment analysis: CLE blister fluid (BF) vs plasma, displayed across IgA, IgG, and IgM
  CSF diagnostic-group analysis: NPSLE vs Other Inflammatory Neurological Disease (OIND), displayed across IgA and IgG

Both of the scripts included import separate isotype-specific CSV files, combine/validate the data, rin antibody-level statistical analysis, and generate one figure per antibody.

## Files
- CSF Boxplot.R AND BF Plasma Boxplot Script.R – Main analysis scripts
- Raw R Ready data/ – Input datasets
- Output/ – Output figures and stats

## Requirements
- R 4.4.2
- tidyverse
Install required packages: install.packages("tidyverse")

## Input data
Script requires R-ready CSV in long format. Raw data provided by UTSW. The following columns are required:
Sample_ID: unique sample identifier
Antibody: autoantibody or antigen name
Isotype: antibody isotype (IgA, IgG, or IgM)
Diagnosis: clinical group, such as NPSLE, OIND, Lupus, or Healthy
Sample_Type: biological compartment, such as CSF, BF, or Plasma
Score: numerical autoantibody reactivity score provided by UTSW

Each input CSV file should only have data for ONE isotype. Control antigens beginning with HIgControl or MIgControl are automatically excluded.

This script also:
- trims extra spaces
- standardize isotype labels
- standardize diagnosis and sample-type labels
- remove missing scores
- check for duplicate sample–antibody–isotype rows
- verifies that each selected file contains the expected isotype

# Analysis 1: CLE blister fluid versus plasma

## Purpose

The CLE script compares autoantibody reactivity between lupus blister fluid and lupus plasma while displaying IgA, IgG, and IgM in separate facets. It imports one R-ready CSV for each isotype.

## Required input files

When prompted, select the files in this order:

1. IgA R-ready CSV
2. IgG R-ready CSV
3. IgM R-ready CSV

The script combines the three files into one long-format dataset

## Study groups

The plotting and post hoc analysis are restricted to:

```text
Diagnosis = Lupus
Sample_Type = Blister or Plasma
```

The factor order is:

```text
Blister
Plasma
```

## Statistical analysis

### Two-way ANOVA

For each antibody, the script fits:

```r
aov(Score ~ Sample_Type * Isotype)
```

This evaluates:

- the main effect of sample type;
- the main effect of isotype;
- the sample type × isotype interaction.

The resulting ANOVA table includes:

- degrees of freedom;
- sum of squares;
- mean square;
- F statistic;
- p-value.

### Tukey HSD comparisons

Within each antibody and isotype, lupus blister and lupus plasma scores are analyzed using:

```r
aov(Score ~ Sample_Type)
TukeyHSD(model, "Sample_Type")
```

The relevant contrast is reported as:

```text
Lupus Plasma-Lupus Blister
```

The Tukey output contains:

| Column | Meaning |
|---|---|
| `diff` | Mean plasma score minus mean blister score |
| `lwr` | Lower 95% confidence limit |
| `upr` | Upper 95% confidence limit |
| `p adj` | Tukey-adjusted p-value |

Interpretation of `diff`:

- positive: higher mean score in plasma;
- negative: higher mean score in blister fluid.

Comparisons with `p adj < 0.05` are retained as significant.

## Figures

One figure is generated per antibody.

Each figure contains:

- IgA, IgG, and IgM facets;
- blister and plasma boxplots;
- individual sample points;
- white box interiors;
- red blister outlines and points;
- blue plasma outlines and points;
- independently scaled y-axes;
- significance brackets only in significant isotype panels.

Asterisks are defined as:

```text
*     p < 0.05
**    p < 0.01
***   p < 0.001
****  p < 0.0001
```

The default plot mode is:

```r
plot_mode <- "significant_only"
```

This generates figures only for antibodies with at least one significant blister-versus-plasma comparison.

To generate every antibody figure, change it to:

```r
plot_mode <- "all"
```

## CLE output files

Outputs are saved beside the selected IgA file in:

```text
Boxplot_Analysis_Outputs/
```

Statistical and quality-control files include:

```text
anova_results.csv
anova_results_labelled.csv
All_Tukey_Comparisons_LupusOnly.csv
Significant_Tukey_Comparisons_LupusOnly.csv
group_counts.csv
isotype_check.csv
duplicate_rows_check.csv
plot_manifest.csv
boxplot_R_sessionInfo_version3.txt
```

Figures are saved in:

```text
Boxplot_Analysis_Outputs/
└── Antibody_Plots/
    ├── AQP4_plot.pdf
    ├── AQP4_plot.png
    └── ...
```

PDF files preserve vector graphics. PNG files are exported at 600 dpi.

---

# Analysis 2: CSF NPSLE versus OIND

## Purpose

The CSF script compares CSF autoantibody scores between patients with NPSLE and patients with OIND. It imports separate IgA and IgG R-ready files and displays the two isotypes in adjacent facets.

## Required input files

When prompted, select the files in this order:

1. CSF IgA R-ready CSV
2. CSF IgG R-ready CSV

## Study groups

The script retains:

```text
Sample_Type = CSF
Diagnosis = OIND or NPSLE
```

The factor order is:

```text
OIND
NPSLE
```

These diagnostic groups are treated as independent groups.

## Statistical analysis

### Two-way ANOVA

For each antibody, the script fits:

```r
aov(Score ~ Diagnosis * Isotype)
```

This evaluates:

- the main effect of diagnosis;
- the main effect of isotype;
- the diagnosis × isotype interaction.

### Tukey HSD comparisons

Within each antibody and isotype, the script fits:

```r
aov(Score ~ Diagnosis)
TukeyHSD(model, "Diagnosis")
```

The resulting contrast is:

```text
NPSLE-OIND
```

The Tukey result contains:

| Column | Meaning |
|---|---|
| `diff` | Mean NPSLE score minus mean OIND score |
| `lwr` | Lower 95% confidence limit |
| `upr` | Upper 95% confidence limit |
| `p adj` | Tukey-adjusted p-value |
| `Direction` | Group with the higher mean score |

Interpretation of `diff`:

- positive: higher mean score in NPSLE;
- negative: higher mean score in OIND.

Comparisons with `p adj < 0.05` are retained as significant.

## Figures

Each antibody figure contains:

- IgA and IgG facets;
- OIND and NPSLE boxplots;
- individual sample points;
- white box interiors;
- blue OIND outlines and points;
- red NPSLE outlines and points;
- independently scaled y-axes;
- significance brackets only in significant isotype panels.

The default setting is:

```r
plot_mode <- "significant_only"
```

To generate all antibody figures:

```r
plot_mode <- "all"
```

## CSF output files

Outputs are saved beside the selected CSF IgA file in:

```text
CSF_Boxplot_Analysis_Outputs/
```

Statistical and quality-control files include:

```text
CSF_anova_results.csv
CSF_anova_results_labelled.csv
CSF_All_Tukey_Comparisons_NPSLE_vs_OIND.csv
CSF_Significant_Tukey_Comparisons_NPSLE_vs_OIND.csv
sample_counts.csv
isotype_check.csv
duplicate_rows_check.csv
CSF_plot_manifest.csv
CSF_boxplot_R_sessionInfo.txt
```

Figures are saved in:

```text
CSF_Boxplot_Analysis_Outputs/
└── Antibody_Plots/
    ├── AQP4_CSF_plot.pdf
    ├── AQP4_CSF_plot.png
    └── ...
```

PDF files preserve vector graphics. PNG files are exported at 600 dpi.

---

## Running the scripts

1. Open the relevant `.R` file in RStudio.
2. Run the entire script.
3. Select each requested isotype-specific CSV in the exact order shown in the file-selection prompts.
4. Wait for the console message confirming that the analysis is complete.
5. Open the automatically generated output folder beside the selected IgA file.

The console reports:

- number of antibodies analyzed;
- number of significant comparisons;
- number of figures generated;
- output-directory location.

---

## Figure interpretation

The boxplots display the median, interquartile range, and whiskers using the default `ggplot2::geom_boxplot()` definition. Individual observations are overlaid with minimal horizontal jitter.

Because facets use:

```r
scales = "free_y"
```

each isotype panel has its own y-axis range. This improves visibility when IgA, IgG, and IgM scores occur on different numerical scales. Absolute vertical positions should therefore not be compared visually across facets without reading the axis values.

---

## Important statistical considerations

### Multiple antibodies

The Tukey adjustment applies to the contrasts within each fitted model. Scripts do not perform an additional false-discovery-rate correction across the full set of antibodies.

Significant results should be described as Tukey-adjusted within-antibody comparisons rather than globally FDR-significant findings.

### CLE sample pairing

Blister-fluid and plasma specimens may come from the same patients. The current CLE script fits ordinary ANOVA models and does not include patient ID as a repeated-measures or blocking factor. Therefore, it treats blister and plasma observations as independent.

For a primary within-patient compartment analysis, a paired test or repeated-measures model may be more appropriate when reliable patient-level matching is available. The present script documents and reproduces the unpaired ANOVA/Tukey workflow.

### Parametric assumptions

ANOVA assumes approximately normally distributed residuals, homogeneous variance, and independent observations. With small groups, these assumptions can be difficult to assess. Diagnostic plots or nonparametric sensitivity analyses may be useful.

---

## Reproducibility

Each script saves the active R session information using:

```r
sessionInfo()
```

The saved file records:

- R version;
- operating system;
- attached packages;
- loaded namespaces;
- package versions.

This information should be retained with the statistical outputs when archiving or publishing the analysis.
