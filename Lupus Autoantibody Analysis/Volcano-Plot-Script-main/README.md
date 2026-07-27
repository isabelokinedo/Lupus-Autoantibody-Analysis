# Lupus Autoantibody Volcano Plot Analysis

This repository contains the R script used to analyze autoantibody specificity across isotype and compartment (CSF, blister fluid (BF), or plasma).
Script supports these comparisons:
  CSF: NPSLE CSF vs Other Inflammatory Neurological Disease (OIND) CSF
  CLE BF: CLE BF vs healthy BF
  CLE compartment: CLE BF vs CLE plasma

## Files
- UTSW Volcano Plot Script.R – Main analysis script
- R ready raw data – Input datasets
- BF vs healthy plots and stats/ – Output figures and statistics

## Requirements
- R 4.4.2
- tidyverse
- ggrepel
Install required packages: install.packages(c("tidyverse", "ggrepel"))

## Input data
Script requires R-ready CSV in long format. Raw data provided by UTSW. The following columns are required:
Sample_ID: unique sample identifier
Antibody: autoantibody or antigen name
Isotype: antibody isotype (IgA, IgG, or IgM)
Diagnosis: clinical group, such as NPSLE, OIND, Lupus, or Healthy
Sample_Type: biological compartment, such as CSF, BF, or Plasma
Score: numerical autoantibody reactivity score provided by UTSW

Each input CSV file should only have data for ONE isotype. Control antigens beginning with HIgControl or MIgControl are automatically excluded.

## Select the analysis

Set comparison variable to one of the following:

comparison <- "CSF"
comparison <- "CLE_blister"
comparison <- "CLE_compartment"

## Comparison specifics

comparison <- "CSF"
  - compares reference group (OIND CSF) to comparison group (NPSLE CSF)
  - median difference is NPSLE median - OIND median
  - positive values indicate higher autoantibody reactivity in NPSLE


comparison <- "CLE_blister"
  - compares reference group (Healthy BF) to comparison group (CLE BF)
  - median difference is CLE median - Healthy median
  - positive values indicate higher autoantibody reactivity in CLE BF


comparison <- "CLE_compartment"
  - compares reference group (CLE Plasma) to comparison group (CLE BF)
  - median difference is BF median - Plasma median
  - positive values indicate higher autoantibody reactivity in CLE BF

### Paired analysis
      - `CLE_compartment` analysis compares paired blister-fluid and plasma samples obtained from the same patients with cutaneous lupus erythematosus (CLE)
      - Unlike the `CSF` and `CLE_blister` comparisons, which use unpaired Mann–Whitney U tests because they compare independent groups, the `CLE_compartment` comparison uses a paired Wilcoxon signed-rank test (`paired = TRUE`) to account for the matched study design
      - Only patients with both blister-fluid and plasma samples are included in the paired analysis. Samples without a corresponding paired specimen are excluded automatically
      - For paired analysis, the volcano plot x-axis represents the median within-patient difference in antibody reactivity (blister fluid − plasma), while p-values are calculated using the paired Wilcoxon signed-rank test. Benjamini–Hochberg false discovery rate correction is then applied to account for multiple testing.

## Running the script
- When prompted by file.choose(), select the R-ready CSV file containing the data for the desired isotype comparison
- Output files are saved automatically in the same folder as the selected input CSV file

## Statistical Analysis
For each autoantibody, the script calculates:
- Number of non-missing samples in each group
- Median reactivity score in the reference group
- Median reactivity score in the comparison group
- Difference between group medians
- Nonparametric comparison appropriate for selected analysis:
    - CSF and CLE_blister use unpaired Mann-Whitney U test, as mentioned above
    - CLE_compartment uses paired Wilcoxon signed-rank test for matched BF and plasma samples from the same patients
- Benjamini–Hochberg false discovery rate
- Negative log10-transformed p-value
- Negative log10-transformed FDR

Statistical testing is performed using the appropriate Wilcoxon rank-based test for the selected analysis:

**CSF and CLE_blister (independent groups):**

```r
wilcox.test(
  comparison_group_scores,
  reference_group_scores,
  paired = FALSE,
  exact = FALSE
)
```

**CLE_compartment (paired samples):**

```r
wilcox.test(
  blister_scores,
  plasma_scores,
  paired = TRUE,
  exact = FALSE
)
```

Multiple-testing correction is performed using:

p.adjust(
  p_value,
  method = "BH"
)

Significance categories are defined as:

FDR < 0.05
Nominal p < 0.05
Not significant

Autoantibodies with nominal p-values below 0.05 are labeled on the volcano plot

Before the analysis, this script:
- Confirms that all required columns are present
- Confirms that the selected file contains only one antibody isotype
- Reports the number of observations in each group
- Reports the number of unique samples
- Reports the number of unique autoantibodies
- Checks for duplicate sample–antibody combinations
- Confirms that both required comparison groups are present

## Volcano Plot Interpretation

X-axis shows the difference between the group median autoantibody scores.

Positive values indicate higher reactivity in the comparison group.
Negative values indicate higher reactivity in the reference group.
The vertical dashed line at zero indicates no difference between group medians.

The y-axis shows: −log10(p-value)
The horizontal dashed line represents: p = 0.05

Filled red points indicate nominally significant or FDR-significant results. Open black points indicate results that are not statistically significant.

Because both nominally significant and FDR-significant points are displayed in red, the corresponding statistics CSV should be consulted to distinguish between these categories.

## Output Files

The script generates three files for each analysis:

<comparison>_<isotype>_statistics.csv
<comparison>_<isotype>_volcano.pdf
<comparison>_<isotype>_volcano.png

The statistics CSV contains the complete results for every autoantibody.
The PDF is saved as a vector-format figure suitable for editing and publication.
The PNG is saved at 600 dpi for high-resolution raster output.

## Extra Considerations
- Missing values: Missing scores are excluded from group sample counts and median calculations; autoantibodies with insufficient observations may return missing p-values
- Multiple testing: nominal p-values below 0.05 should be interpreted cautiously when many autoantibodies are tested. FDR-adjusted values should be used to evaluate significance after correction for multiple comparisons.
- Small sample sizes: results from comparisons with small numbers of samples should be described as exploratory and interpreted based on effect direction and magnitude in addition to p-values.

Citation

When using this script in a publication, report that autoantibody reactivity scores were compared using Wilcoxon rank-based tests and that p-values were adjusted using the Benjamini–Hochberg false discovery rate procedure.
