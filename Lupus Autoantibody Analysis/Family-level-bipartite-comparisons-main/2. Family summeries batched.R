############################################################
# BATCH FAMILY SUMMARY GENERATOR
#
# To generate standard and network-ready family summaries for
# every Spearman master table in a selected folder.
#
# Required inputs:
#   1. Folder containing all *_master_table.csv files
#   2. Autoantibody annotation CSV
#   3. Olink protein annotation CSV
############################################################


# ==========================================================
# 1. INSTALL AND LOAD PACKAGES
# ==========================================================

required_packages <- c(
  "tidyverse"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)


# ==========================================================
# 2. SELECT INPUT FILES
# ==========================================================

cat(
  "\nSelect any Spearman master-table file located in the\n",
  "folder containing all of your master tables.\n\n"
)

example_master_file <- file.choose()

master_table_directory <- dirname(
  example_master_file
)


cat(
  "\nSelect the corrected autoantibody annotation CSV.\n\n"
)

antibody_annotation_file <- file.choose()


cat(
  "\nSelect the Olink protein annotation CSV.\n\n"
)

protein_annotation_file <- file.choose()


# ==========================================================
# 3. CREATE OUTPUT DIRECTORY
# ==========================================================

output_directory <- file.path(
  master_table_directory,
  "Family_summary_results"
)

if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE
  )
}

cat(
  "\nResults will be saved in:\n",
  output_directory,
  "\n\n"
)


# ==========================================================
# 4. FIND ALL MASTER TABLES
# ==========================================================

master_files <- list.files(
  path = master_table_directory,
  pattern = "_master_table\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Exclude previously generated summary files, if present
master_files <- master_files[
  !str_detect(
    basename(master_files),
    "family_summary|network_family_summary"
  )
]

if (length(master_files) == 0) {
  stop(
    "No files ending in '_master_table.csv' were found."
  )
}

cat(
  "Master tables found:",
  length(master_files),
  "\n\n"
)

print(
  basename(master_files)
)


# ==========================================================
# 5. IMPORT ANNOTATION FILES
# ==========================================================

antibody_annotations <- read.csv(
  antibody_annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

protein_annotations <- read.csv(
  protein_annotation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ==========================================================
# 6. STANDARDIZE ANNOTATION COLUMN NAMES
# ==========================================================

find_first_column <- function(
    data,
    possible_names,
    table_name
) {
  
  matched_name <- possible_names[
    possible_names %in% names(data)
  ]
  
  if (length(matched_name) == 0) {
    stop(
      paste0(
        "Could not identify the required column in ",
        table_name,
        ". Expected one of: ",
        paste(possible_names, collapse = ", ")
      )
    )
  }
  
  matched_name[[1]]
}


antibody_name_column <- find_first_column(
  antibody_annotations,
  c(
    "Antibody",
    "Autoantibody",
    "Antigen",
    "Target"
  ),
  "autoantibody annotation file"
)

antibody_family_column <- find_first_column(
  antibody_annotations,
  c(
    "Family",
    "AntibodyFamily",
    "Antibody_family",
    "AutoantibodyFamily"
  ),
  "autoantibody annotation file"
)

antibody_subfamily_column <- find_first_column(
  antibody_annotations,
  c(
    "Subfamily",
    "AntibodySubfamily",
    "Antibody_subfamily",
    "AutoantibodySubfamily"
  ),
  "autoantibody annotation file"
)


protein_name_column <- find_first_column(
  protein_annotations,
  c(
    "Protein",
    "Olink",
    "Cytokine",
    "Analyte",
    "Assay"
  ),
  "Olink protein annotation file"
)

protein_family_column <- find_first_column(
  protein_annotations,
  c(
    "Family",
    "ProteinFamily",
    "Protein_family",
    "CytokineFamily"
  ),
  "Olink protein annotation file"
)

protein_subfamily_column <- find_first_column(
  protein_annotations,
  c(
    "Subfamily",
    "ProteinSubfamily",
    "Protein_subfamily",
    "CytokineSubfamily"
  ),
  "Olink protein annotation file"
)


antibody_annotations_clean <- antibody_annotations %>%
  transmute(
    Antibody = trimws(
      as.character(
        .data[[antibody_name_column]]
      )
    ),
    
    AntibodyFamily = trimws(
      as.character(
        .data[[antibody_family_column]]
      )
    ),
    
    AntibodySubfamily = trimws(
      as.character(
        .data[[antibody_subfamily_column]]
      )
    )
  ) %>%
  filter(
    !is.na(Antibody),
    Antibody != ""
  ) %>%
  distinct(
    Antibody,
    .keep_all = TRUE
  )


protein_annotations_clean <- protein_annotations %>%
  transmute(
    Protein = trimws(
      as.character(
        .data[[protein_name_column]]
      )
    ),
    
    ProteinFamily = trimws(
      as.character(
        .data[[protein_family_column]]
      )
    ),
    
    ProteinSubfamily = trimws(
      as.character(
        .data[[protein_subfamily_column]]
      )
    )
  ) %>%
  filter(
    !is.na(Protein),
    Protein != ""
  ) %>%
  distinct(
    Protein,
    .keep_all = TRUE
  )


# ==========================================================
# 7. CHECK ANNOTATIONS FOR DUPLICATE TARGET NAMES
# ==========================================================

duplicate_antibodies <- antibody_annotations_clean %>%
  count(Antibody) %>%
  filter(n > 1)

duplicate_proteins <- protein_annotations_clean %>%
  count(Protein) %>%
  filter(n > 1)

if (nrow(duplicate_antibodies) > 0) {
  warning(
    "Duplicate antibody names were found in the annotation file."
  )
}

if (nrow(duplicate_proteins) > 0) {
  warning(
    "Duplicate protein names were found in the annotation file."
  )
}


# ==========================================================
# 8. FUNCTION TO PROCESS ONE MASTER TABLE
# ==========================================================

process_master_table <- function(
    master_file,
    antibody_annotations_clean,
    protein_annotations_clean,
    output_directory
) {
  
  cat(
    "\nProcessing:",
    basename(master_file),
    "\n"
  )
  
  
  # --------------------------------------------------------
  # 8A. Import master table
  # --------------------------------------------------------
  
  master <- read.csv(
    master_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  
  # --------------------------------------------------------
  # 8B. Identify required master-table columns
  # --------------------------------------------------------
  
  antibody_column <- find_first_column(
    master,
    c(
      "Antibody",
      "Autoantibody",
      "Antigen"
    ),
    basename(master_file)
  )
  
  protein_column <- find_first_column(
    master,
    c(
      "Protein",
      "Olink",
      "Cytokine",
      "Analyte"
    ),
    basename(master_file)
  )
  
  rho_column <- find_first_column(
    master,
    c(
      "rho",
      "Rho",
      "Spearman_rho",
      "Spearman rho",
      "SpearmanRho",
      "r"
    ),
    basename(master_file)
  )
  
  p_column <- find_first_column(
    master,
    c(
      "p",
      "P",
      "p_value",
      "p.value",
      "P_value",
      "P.value",
      "Raw_p"
    ),
    basename(master_file)
  )
  
  fdr_column <- find_first_column(
    master,
    c(
      "FDR",
      "fdr",
      "FDR_value",
      "Adjusted_p",
      "adjusted_p",
      "p_adjusted",
      "p_adj"
    ),
    basename(master_file)
  )
  
  n_column <- find_first_column(
    master,
    c(
      "N",
      "n",
      "Pairwise_N",
      "pairwise_n",
      "Sample_N"
    ),
    basename(master_file)
  )
  
  
  # --------------------------------------------------------
  # 8C. Clean master table
  # --------------------------------------------------------
  
  master_clean <- master %>%
    transmute(
      Antibody = trimws(
        as.character(
          .data[[antibody_column]]
        )
      ),
      
      Protein = trimws(
        as.character(
          .data[[protein_column]]
        )
      ),
      
      rho = suppressWarnings(
        as.numeric(
          .data[[rho_column]]
        )
      ),
      
      p = suppressWarnings(
        as.numeric(
          .data[[p_column]]
        )
      ),
      
      FDR = suppressWarnings(
        as.numeric(
          .data[[fdr_column]]
        )
      ),
      
      N = suppressWarnings(
        as.numeric(
          .data[[n_column]]
        )
      )
    ) %>%
    filter(
      !is.na(Antibody),
      Antibody != "",
      !is.na(Protein),
      Protein != "",
      is.finite(rho)
    )
  
  
  # --------------------------------------------------------
  # 8D. Join annotations
  # --------------------------------------------------------
  
  annotated <- master_clean %>%
    left_join(
      antibody_annotations_clean,
      by = "Antibody"
    ) %>%
    left_join(
      protein_annotations_clean,
      by = "Protein"
    )
  
  
  # --------------------------------------------------------
  # 8E. Save unmatched annotation names
  # --------------------------------------------------------
  
  unmatched_antibodies <- annotated %>%
    filter(
      is.na(AntibodyFamily)
    ) %>%
    distinct(Antibody) %>%
    arrange(Antibody)
  
  unmatched_proteins <- annotated %>%
    filter(
      is.na(ProteinFamily)
    ) %>%
    distinct(Protein) %>%
    arrange(Protein)
  
  
  base_name <- basename(master_file) %>%
    str_remove(
      regex(
        "\\.csv$",
        ignore_case = TRUE
      )
    )
  
  
  if (nrow(unmatched_antibodies) > 0) {
    
    write.csv(
      unmatched_antibodies,
      file.path(
        output_directory,
        paste0(
          base_name,
          "_UNMATCHED_ANTIBODIES.csv"
        )
      ),
      row.names = FALSE
    )
    
    warning(
      paste0(
        basename(master_file),
        ": ",
        nrow(unmatched_antibodies),
        " antibody names did not match the annotation file."
      )
    )
  }
  
  
  if (nrow(unmatched_proteins) > 0) {
    
    write.csv(
      unmatched_proteins,
      file.path(
        output_directory,
        paste0(
          base_name,
          "_UNMATCHED_PROTEINS.csv"
        )
      ),
      row.names = FALSE
    )
    
    warning(
      paste0(
        basename(master_file),
        ": ",
        nrow(unmatched_proteins),
        " protein names did not match the annotation file."
      )
    )
  }
  
  
  # Use only fully annotated rows for family summaries
  annotated_complete <- annotated %>%
    filter(
      !is.na(AntibodyFamily),
      AntibodyFamily != "",
      !is.na(ProteinFamily),
      ProteinFamily != ""
    )
  
  
  if (nrow(annotated_complete) == 0) {
    warning(
      paste0(
        "No fully annotated rows remained for ",
        basename(master_file)
      )
    )
    
    return(
      tibble(
        File = basename(master_file),
        Total_rows = nrow(master_clean),
        Annotated_rows = 0,
        Family_pairs = 0,
        Status = "No annotated rows"
      )
    )
  }
  
  
  # --------------------------------------------------------
  # 8F. Generate standard family summary
  # --------------------------------------------------------
  
  family_summary <- annotated_complete %>%
    group_by(
      AntibodyFamily,
      ProteinFamily
    ) %>%
    summarise(
      Number_of_pairs = n(),
      
      Number_of_antibodies = n_distinct(
        Antibody
      ),
      
      Number_of_proteins = n_distinct(
        Protein
      ),
      
      Mean_rho = mean(
        rho,
        na.rm = TRUE
      ),
      
      Median_rho = median(
        rho,
        na.rm = TRUE
      ),
      
      Mean_abs_rho = mean(
        abs(rho),
        na.rm = TRUE
      ),
      
      Median_abs_rho = median(
        abs(rho),
        na.rm = TRUE
      ),
      
      Max_positive = if_else(
        any(rho > 0, na.rm = TRUE),
        max(
          rho[rho > 0],
          na.rm = TRUE
        ),
        NA_real_
      ),
      
      Max_negative = if_else(
        any(rho < 0, na.rm = TRUE),
        min(
          rho[rho < 0],
          na.rm = TRUE
        ),
        NA_real_
      ),
      
      Significant_raw_count = sum(
        p < 0.05,
        na.rm = TRUE
      ),
      
      Significant_FDR_count = sum(
        FDR < 0.05,
        na.rm = TRUE
      ),
      
      Percent_significant_raw =
        100 *
        Significant_raw_count /
        Number_of_pairs,
      
      Percent_significant_FDR =
        100 *
        Significant_FDR_count /
        Number_of_pairs,
      
      Median_pairwise_N = median(
        N,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    arrange(
      desc(Median_abs_rho)
    )
  
  
  standard_output_name <- paste0(
    base_name,
    "_family_summary.csv"
  )
  
  write.csv(
    family_summary,
    file.path(
      output_directory,
      standard_output_name
    ),
    row.names = FALSE
  )
  
  
  # --------------------------------------------------------
  # 8G. Generate network-ready family summary
  # --------------------------------------------------------
  
  network_family_summary <- annotated_complete %>%
    group_by(
      AntibodyFamily,
      ProteinFamily
    ) %>%
    summarise(
      Number_of_pairs = n(),
      
      Number_of_antibodies = n_distinct(
        Antibody
      ),
      
      Number_of_proteins = n_distinct(
        Protein
      ),
      
      Mean_rho = mean(
        rho,
        na.rm = TRUE
      ),
      
      Median_rho = median(
        rho,
        na.rm = TRUE
      ),
      
      Mean_abs_rho = mean(
        abs(rho),
        na.rm = TRUE
      ),
      
      Median_abs_rho = median(
        abs(rho),
        na.rm = TRUE
      ),
      
      Max_positive = if_else(
        any(rho > 0, na.rm = TRUE),
        max(
          rho[rho > 0],
          na.rm = TRUE
        ),
        NA_real_
      ),
      
      Max_negative = if_else(
        any(rho < 0, na.rm = TRUE),
        min(
          rho[rho < 0],
          na.rm = TRUE
        ),
        NA_real_
      ),
      
      Significant_raw_count = sum(
        p < 0.05,
        na.rm = TRUE
      ),
      
      Significant_FDR_count = sum(
        FDR < 0.05,
        na.rm = TRUE
      ),
      
      Number_positive = sum(
        rho > 0,
        na.rm = TRUE
      ),
      
      Number_negative = sum(
        rho < 0,
        na.rm = TRUE
      ),
      
      Number_zero = sum(
        rho == 0,
        na.rm = TRUE
      ),
      
      # Strong pair defined as |rho| >= 0.50
      Strong_pairs = sum(
        abs(rho) >= 0.50,
        na.rm = TRUE
      ),
      
      Median_pairwise_N = median(
        N,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) %>%
    mutate(
      Percent_significant_raw =
        100 *
        Significant_raw_count /
        Number_of_pairs,
      
      Percent_significant_FDR =
        100 *
        Significant_FDR_count /
        Number_of_pairs,
      
      Proportion_positive =
        Number_positive /
        Number_of_pairs,
      
      Proportion_negative =
        Number_negative /
        Number_of_pairs,
      
      Direction_consistency = pmax(
        Proportion_positive,
        Proportion_negative
      ),
      
      Dominant_direction = case_when(
        Number_positive > Number_negative ~
          "Positive",
        
        Number_negative > Number_positive ~
          "Negative",
        
        TRUE ~
          "Mixed"
      ),
      
      Proportion_strong =
        Strong_pairs /
        Number_of_pairs,
      
      Ranking_score =
        Median_abs_rho *
        Direction_consistency *
        Proportion_strong
    ) %>%
    arrange(
      desc(Ranking_score),
      desc(Median_abs_rho)
    )
  
  
  network_output_name <- paste0(
    base_name,
    "_network_family_summary.csv"
  )
  
  write.csv(
    network_family_summary,
    file.path(
      output_directory,
      network_output_name
    ),
    row.names = FALSE
  )
  
  
  # --------------------------------------------------------
  # 8H. Return processing log
  # --------------------------------------------------------
  
  tibble(
    File = basename(master_file),
    
    Total_rows =
      nrow(master_clean),
    
    Annotated_rows =
      nrow(annotated_complete),
    
    Unmatched_antibodies =
      nrow(unmatched_antibodies),
    
    Unmatched_proteins =
      nrow(unmatched_proteins),
    
    Family_pairs =
      nrow(family_summary),
    
    Status =
      "Completed"
  )
}


# ==========================================================
# 9. PROCESS EVERY MASTER TABLE
# ==========================================================

processing_log <- map_dfr(
  master_files,
  process_master_table,
  antibody_annotations_clean =
    antibody_annotations_clean,
  protein_annotations_clean =
    protein_annotations_clean,
  output_directory =
    output_directory
)


# ==========================================================
# 10. SAVE PROCESSING LOG
# ==========================================================

write.csv(
  processing_log,
  file.path(
    output_directory,
    "family_summary_processing_log.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 11. PRINT FINAL SUMMARY
# ==========================================================

cat(
  "\n============================================\n",
  "BATCH FAMILY SUMMARY GENERATION COMPLETE\n",
  "============================================\n\n"
)

print(
  processing_log
)

cat(
  "\nFiles were saved in:\n",
  output_directory,
  "\n\n"
)
