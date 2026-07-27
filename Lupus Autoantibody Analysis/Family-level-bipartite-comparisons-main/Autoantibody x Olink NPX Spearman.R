############################################################
# BF IgM autoantibody × Olink NPX Spearman correlations
############################################################

library(tidyverse)


############################################################
# 1. Analysis settings
############################################################

sample_type_to_use <- "Blister"
isotype_to_use <- "IgM"
output_prefix <- "BF_IgM_Olink_Spearman"


############################################################
# 2. Select input files
############################################################

message("Select the cleaned Olink blister-fluid CSV file")
olink_file <- file.choose()

message("Select the UTSW autoantibody CSV file")
antibody_file <- file.choose()


############################################################
# 3. Read files
############################################################

olink <- read.csv(
  olink_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

antibody_raw <- read.csv(
  antibody_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove any unnamed columns created by Excel
antibody_raw <- antibody_raw |>
  select(-matches("^$"))

############################################################
# 4. Check required columns
############################################################

required_olink_columns <- c(
  "Sample_ID",
  "Clinical_Label",
  "OlinkID"
)

required_antibody_columns <- c(
  "Sample_ID",
  "Sample_Type",
  "Isotype",
  "Antibody",
  "Score"
)

missing_olink_columns <- setdiff(
  required_olink_columns,
  names(olink)
)

missing_antibody_columns <- setdiff(
  required_antibody_columns,
  names(antibody_raw)
)

if (length(missing_olink_columns) > 0) {
  stop(
    "Missing required Olink columns: ",
    paste(missing_olink_columns, collapse = ", ")
  )
}

if (length(missing_antibody_columns) > 0) {
  stop(
    "Missing required antibody columns: ",
    paste(missing_antibody_columns, collapse = ", ")
  )
}


############################################################
# 5. Clean identifiers and text fields
############################################################

olink <- olink |>
  mutate(
    Sample_ID = trimws(Sample_ID),
    Clinical_Label = trimws(Clinical_Label),
    OlinkID = trimws(OlinkID)
  )

antibody_raw <- antibody_raw |>
  mutate(
    Sample_ID = trimws(Sample_ID),
    Sample_Type = trimws(Sample_Type),
    Isotype = trimws(Isotype),
    Antibody = trimws(Antibody),
    Score = as.numeric(Score)
  )


############################################################
# 6. Keep the selected compartment and isotype
############################################################

antibody_filtered <- antibody_raw |>
  filter(
    Sample_Type == sample_type_to_use,
    Isotype == isotype_to_use
  )

if (nrow(antibody_filtered) == 0) {
  stop(
    "No antibody rows were found for Sample_Type = ",
    sample_type_to_use,
    " and Isotype = ",
    isotype_to_use
  )
}


############################################################
# 7. Check for duplicate sample-antibody measurements
############################################################

duplicate_measurements <- antibody_filtered |>
  count(Sample_ID, Antibody) |>
  filter(n > 1)

if (nrow(duplicate_measurements) > 0) {
  stop(
    "Duplicate Sample_ID–Antibody measurements were found. ",
    "Review the antibody input file before continuing."
  )
}


############################################################
# 8. Convert antibody data from long to wide
############################################################

antibody_wide <- antibody_filtered |>
  select(
    Sample_ID,
    Antibody,
    Score
  ) |>
  pivot_wider(
    names_from = Antibody,
    values_from = Score
  )


############################################################
# 9. Identify Olink protein columns
############################################################

olink_metadata_columns <- c(
  "Sample_ID",
  "Clinical_Label",
  "OlinkID",
  "Plate ID",
  "QC Warning"
)

protein_cols <- setdiff(
  names(olink),
  olink_metadata_columns
)

if (length(protein_cols) == 0) {
  stop("No Olink protein columns were detected.")
}

olink <- olink |>
  mutate(
    across(
      all_of(protein_cols),
      as.numeric
    )
  )


############################################################
# 10. Identify antibody columns
############################################################

antibody_cols <- setdiff(
  names(antibody_wide),
  "Sample_ID"
)

if (length(antibody_cols) == 0) {
  stop("No antibody columns were detected.")
}

antibody_wide <- antibody_wide |>
  mutate(
    across(
      all_of(antibody_cols),
      as.numeric
    )
  )


############################################################
# 11. Check matched sample IDs
############################################################

matched_ids <- intersect(
  olink$Sample_ID,
  antibody_wide$Sample_ID
)

message("Matched samples: ", length(matched_ids))
print(matched_ids)

message("\nOlink samples without matching antibody data:")
print(
  setdiff(
    olink$Sample_ID,
    antibody_wide$Sample_ID
  )
)

message("\nAntibody samples without matching Olink data:")
print(
  setdiff(
    antibody_wide$Sample_ID,
    olink$Sample_ID
  )
)

if (length(matched_ids) < 3) {
  stop(
    "Fewer than three paired samples were found. ",
    "Spearman correlations cannot be calculated reliably."
  )
}


############################################################
# 12. Prefix measurement columns before merging
############################################################

olink_prefixed <- olink |>
  rename_with(
    ~ paste0("Olink__", .x),
    all_of(protein_cols)
  )

antibody_prefixed <- antibody_wide |>
  rename_with(
    ~ paste0("Ab__", .x),
    all_of(antibody_cols)
  )


############################################################
# 13. Merge paired samples
############################################################

merged <- olink_prefixed |>
  inner_join(
    antibody_prefixed,
    by = "Sample_ID"
  )

message(
  "\nNumber of paired samples after merging: ",
  nrow(merged)
)

protein_cols_prefixed <- paste0(
  "Olink__",
  protein_cols
)

antibody_cols_prefixed <- paste0(
  "Ab__",
  antibody_cols
)

stopifnot(
  all(protein_cols_prefixed %in% names(merged)),
  all(antibody_cols_prefixed %in% names(merged))
)


############################################################
# 14. Initialize result matrices
############################################################

rho_matrix <- matrix(
  NA_real_,
  nrow = length(antibody_cols),
  ncol = length(protein_cols),
  dimnames = list(
    antibody_cols,
    protein_cols
  )
)

p_matrix <- matrix(
  NA_real_,
  nrow = length(antibody_cols),
  ncol = length(protein_cols),
  dimnames = list(
    antibody_cols,
    protein_cols
  )
)

n_matrix <- matrix(
  NA_integer_,
  nrow = length(antibody_cols),
  ncol = length(protein_cols),
  dimnames = list(
    antibody_cols,
    protein_cols
  )
)


############################################################
# 15. Calculate Spearman correlations
############################################################

for (i in seq_along(antibody_cols)) {
  
  for (j in seq_along(protein_cols)) {
    
    x <- merged[[antibody_cols_prefixed[i]]]
    y <- merged[[protein_cols_prefixed[j]]]
    
    complete_rows <- complete.cases(x, y)
    
    x_complete <- x[complete_rows]
    y_complete <- y[complete_rows]
    
    n_complete <- length(x_complete)
    
    n_matrix[i, j] <- n_complete
    
    if (
      n_complete >= 3 &&
      length(unique(x_complete)) > 1 &&
      length(unique(y_complete)) > 1
    ) {
      
      test_result <- suppressWarnings(
        cor.test(
          x_complete,
          y_complete,
          method = "spearman",
          exact = FALSE
        )
      )
      
      rho_matrix[i, j] <-
        unname(test_result$estimate)
      
      p_matrix[i, j] <-
        test_result$p.value
    }
  }
}


############################################################
# 16. Benjamini-Hochberg FDR correction
############################################################

fdr_matrix <- matrix(
  p.adjust(
    as.vector(p_matrix),
    method = "BH"
  ),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix),
  dimnames = dimnames(p_matrix)
)


############################################################
# 17. Create long-format master results table
############################################################

results_long <- expand_grid(
  Protein = colnames(rho_matrix),
  Antibody = rownames(rho_matrix)
) |>
  mutate(
    Spearman_rho = as.vector(rho_matrix),
    Raw_p = as.vector(p_matrix),
    FDR = as.vector(fdr_matrix),
    N = as.vector(n_matrix)
  ) |>
  select(
    Antibody,
    Protein,
    Spearman_rho,
    Raw_p,
    FDR,
    N
  )


############################################################
# 18. Save matched sample list
############################################################

matched_sample_table <- merged |>
  select(
    Sample_ID,
    Clinical_Label,
    OlinkID
  )


############################################################
# 19. Save outputs to the working directory
############################################################

write.csv(
  rho_matrix,
  paste0(output_prefix, "_rho_matrix.csv")
)

write.csv(
  p_matrix,
  paste0(output_prefix, "_raw_p_matrix.csv")
)

write.csv(
  fdr_matrix,
  paste0(output_prefix, "_FDR_matrix.csv")
)

write.csv(
  n_matrix,
  paste0(output_prefix, "_n_matrix.csv")
)

write.csv(
  results_long,
  paste0(output_prefix, "_master_table.csv"),
  row.names = FALSE
)

write.csv(
  matched_sample_table,
  paste0(output_prefix, "_matched_samples.csv"),
  row.names = FALSE
)


############################################################
# 20. Final checks
############################################################

message("\nAnalysis complete.")
message("Matched samples: ", nrow(merged))
message("Antibodies analyzed: ", length(antibody_cols))
message("Olink proteins analyzed: ", length(protein_cols))
message("Results saved to: ", getwd())

message("\nSpearman rho range:")
print(
  range(
    rho_matrix,
    na.rm = TRUE
  )
)

message("\nPairwise sample-size distribution:")
print(
  table(
    n_matrix,
    useNA = "ifany"
  )
)