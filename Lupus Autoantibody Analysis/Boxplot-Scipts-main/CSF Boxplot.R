# ============================================================
# CSF autoantibody boxplot analysis: NPSLE versus OIND
#
# Designed for separate R-ready IgA and IgG CSV files.
#
# Expected columns:
#   Sample_ID
#   Antibody
#   Isotype
#   Diagnosis
#   Sample_Type
#   Score
#
# Analysis:
#   1. Combine IgA and IgG CSF data.
#   2. For each antibody, fit a two-way ANOVA:
#
#          Score ~ Diagnosis * Isotype
#
#   3. Within each antibody and isotype, compare NPSLE with OIND
#      using Tukey HSD after a one-way model:
#
#          Score ~ Diagnosis
#
#   4. Create one figure per antibody, with IgA and IgG shown
#      side by side.
# ============================================================


# ------------------------------------------------------------
# 1. INSTALL AND LOAD TIDYVERSE
# ------------------------------------------------------------

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(tidyverse)


# ------------------------------------------------------------
# 2. SELECT THE TWO R-READY INPUT FILES
# ------------------------------------------------------------

message("Select the CSF IgA R-ready CSV file.")
iga_file <- file.choose()

message("Select the CSF IgG R-ready CSV file.")
igg_file <- file.choose()

input_files <- c(
  IgA = iga_file,
  IgG = igg_file
)


# Outputs are saved beside the IgA input file.
output_dir <- file.path(
  dirname(iga_file),
  "CSF_Boxplot_Analysis_Outputs"
)

plot_dir <- file.path(
  output_dir,
  "Antibody_Plots"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. IMPORT AND COMBINE DATA
# ------------------------------------------------------------

dat <- purrr::imap_dfr(
  input_files,
  function(file_path, expected_isotype) {

    readr::read_csv(
      file_path,
      show_col_types = FALSE
    ) |>
      mutate(
        Source_File = basename(file_path),
        Expected_Isotype = expected_isotype
      )
  }
)


# ------------------------------------------------------------
# 4. VALIDATE REQUIRED COLUMNS
# ------------------------------------------------------------

required_columns <- c(
  "Sample_ID",
  "Antibody",
  "Isotype",
  "Diagnosis",
  "Sample_Type",
  "Score"
)

missing_columns <- setdiff(
  required_columns,
  names(dat)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------
# 5. CLEAN AND STANDARDIZE DATA
# ------------------------------------------------------------

dat <- dat |>
  mutate(
    Sample_ID = str_squish(as.character(Sample_ID)),
    Antibody = str_squish(as.character(Antibody)),
    Isotype = str_squish(as.character(Isotype)),
    Diagnosis = str_squish(as.character(Diagnosis)),
    Sample_Type = str_squish(as.character(Sample_Type)),
    Score = as.numeric(Score)
  ) |>
  mutate(
    Isotype = case_when(
      str_to_upper(Isotype) == "IGA" ~ "IgA",
      str_to_upper(Isotype) == "IGG" ~ "IgG",
      TRUE ~ Isotype
    ),
    Diagnosis = case_when(
      str_to_upper(Diagnosis) == "OIND" ~ "OIND",
      str_to_upper(Diagnosis) == "NPSLE" ~ "NPSLE",
      TRUE ~ Diagnosis
    ),
    Sample_Type = case_when(
      str_to_upper(Sample_Type) == "CSF" ~ "CSF",
      TRUE ~ Sample_Type
    )
  ) |>
  filter(
    !is.na(Score),
    !is.na(Antibody),
    Antibody != "",
    Diagnosis %in% c("OIND", "NPSLE"),
    Sample_Type == "CSF"
  ) |>
  filter(
    !str_detect(
      Antibody,
      regex(
        "^(HIgControl|MIgControl|Control)",
        ignore_case = TRUE
      )
    )
  ) |>
  mutate(
    Isotype = factor(
      Isotype,
      levels = c("IgA", "IgG")
    ),
    Diagnosis = factor(
      Diagnosis,
      levels = c("OIND", "NPSLE")
    )
  )


# ------------------------------------------------------------
# 6. QUALITY-CONTROL OUTPUTS
# ------------------------------------------------------------

isotype_check <- dat |>
  distinct(
    Source_File,
    Expected_Isotype,
    Isotype
  )

write_csv(
  isotype_check,
  file.path(
    output_dir,
    "isotype_check.csv"
  )
)

if (
  any(
    as.character(isotype_check$Isotype) !=
      isotype_check$Expected_Isotype
  )
) {
  warning(
    "At least one selected file may contain the wrong isotype. ",
    "Review isotype_check.csv."
  )
}


sample_counts <- dat |>
  distinct(
    Sample_ID,
    Diagnosis,
    Isotype
  ) |>
  count(
    Diagnosis,
    Isotype,
    name = "n_samples"
  )

write_csv(
  sample_counts,
  file.path(
    output_dir,
    "sample_counts.csv"
  )
)


duplicate_check <- dat |>
  count(
    Sample_ID,
    Antibody,
    Isotype,
    name = "n"
  ) |>
  filter(n > 1)

write_csv(
  duplicate_check,
  file.path(
    output_dir,
    "duplicate_rows_check.csv"
  )
)

if (nrow(duplicate_check) > 0) {
  warning(
    "Duplicate Sample_ID–Antibody–Isotype rows were found. ",
    "Review duplicate_rows_check.csv."
  )
}


# ------------------------------------------------------------
# 7. TWO-WAY ANOVA FOR EACH ANTIBODY
# ------------------------------------------------------------

# Model terms:
#   Diagnosis
#   Isotype
#   Diagnosis:Isotype
#   Residuals

run_two_way_anova <- function(x) {

  model_data <- x |>
    drop_na(
      Score,
      Diagnosis,
      Isotype
    ) |>
    droplevels()

  if (
    n_distinct(model_data$Diagnosis) < 2 ||
    n_distinct(model_data$Isotype) < 2
  ) {
    return(
      tibble(
        Term = NA_character_,
        Df = NA_real_,
        `Sum Sq` = NA_real_,
        `Mean Sq` = NA_real_,
        `F value` = NA_real_,
        `Pr(>F)` = NA_real_
      )
    )
  }

  model <- aov(
    Score ~ Diagnosis * Isotype,
    data = model_data
  )

  as.data.frame(
    summary(model)[[1]]
  ) |>
    rownames_to_column("Term") |>
    as_tibble()
}


anova_results_labelled <- dat |>
  group_by(Antibody) |>
  group_modify(
    ~ run_two_way_anova(.x)
  ) |>
  ungroup() |>
  relocate(
    Antibody,
    Term
  )

write_csv(
  anova_results_labelled,
  file.path(
    output_dir,
    "CSF_anova_results_labelled.csv"
  )
)


anova_results <- anova_results_labelled |>
  select(
    Df,
    `Sum Sq`,
    `Mean Sq`,
    `F value`,
    `Pr(>F)`,
    Antibody
  )

write_csv(
  anova_results,
  file.path(
    output_dir,
    "CSF_anova_results.csv"
  )
)


# ------------------------------------------------------------
# 8. NPSLE-VERSUS-OIND TUKEY HSD TESTS
# ------------------------------------------------------------

# Testing is performed separately for each antibody and isotype.
# With two diagnostic groups, there is one Tukey contrast:
#
#   NPSLE - OIND
#
# Positive difference means higher mean score in NPSLE.
# Negative difference means higher mean score in OIND.

run_tukey <- function(x) {

  model_data <- x |>
    drop_na(
      Score,
      Diagnosis
    ) |>
    droplevels()

  if (n_distinct(model_data$Diagnosis) < 2) {
    return(
      tibble(
        Comparison = NA_character_,
        diff = NA_real_,
        lwr = NA_real_,
        upr = NA_real_,
        `p adj` = NA_real_
      )
    )
  }

  model <- aov(
    Score ~ Diagnosis,
    data = model_data
  )

  TukeyHSD(
    model,
    "Diagnosis"
  )$Diagnosis |>
    as.data.frame() |>
    rownames_to_column("Comparison") |>
    as_tibble() |>
    select(
      Comparison,
      diff,
      lwr,
      upr,
      `p adj`
    )
}


all_tukey_results <- dat |>
  group_by(
    Antibody,
    Isotype
  ) |>
  group_modify(
    ~ run_tukey(.x)
  ) |>
  ungroup() |>
  select(
    Comparison,
    diff,
    lwr,
    upr,
    `p adj`,
    Antibody,
    Isotype
  ) |>
  arrange(
    Isotype,
    `p adj`,
    Antibody
  )

write_csv(
  all_tukey_results,
  file.path(
    output_dir,
    "CSF_All_Tukey_Comparisons_NPSLE_vs_OIND.csv"
  )
)


significant_tukey <- all_tukey_results |>
  filter(
    Comparison == "NPSLE-OIND",
    !is.na(`p adj`),
    `p adj` < 0.05
  ) |>
  mutate(
    Direction = case_when(
      diff > 0 ~ "Higher in NPSLE",
      diff < 0 ~ "Higher in OIND",
      TRUE ~ "No difference"
    )
  )

write_csv(
  significant_tukey,
  file.path(
    output_dir,
    "CSF_Significant_Tukey_Comparisons_NPSLE_vs_OIND.csv"
  )
)


# ------------------------------------------------------------
# 9. PREPARE SIGNIFICANCE BRACKETS
# ------------------------------------------------------------

annotation_data <- significant_tukey |>
  mutate(
    significance = case_when(
      `p adj` < 0.0001 ~ "****",
      `p adj` < 0.001 ~ "***",
      `p adj` < 0.01 ~ "**",
      `p adj` < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) |>
  left_join(
    dat |>
      group_by(
        Antibody,
        Isotype
      ) |>
      summarise(
        y_min = min(Score, na.rm = TRUE),
        y_max = max(Score, na.rm = TRUE),
        y_range = y_max - y_min,
        .groups = "drop"
      ),
    by = c("Antibody", "Isotype")
  ) |>
  mutate(
    bracket_height = case_when(
      y_range > 0 ~ y_max + 0.08 * y_range,
      TRUE ~ y_max + 0.5
    ),
    tick_height = case_when(
      y_range > 0 ~ 0.025 * y_range,
      TRUE ~ 0.15
    ),
    label_height = bracket_height +
      case_when(
        y_range > 0 ~ 0.035 * y_range,
        TRUE ~ 0.2
      )
  )


# ------------------------------------------------------------
# 10. FIGURE STYLE
# ------------------------------------------------------------

# OIND is blue and NPSLE is red.
diagnosis_colors <- c(
  "OIND" = "#4575B4",
  "NPSLE" = "#D73027"
)


safe_file_name <- function(x) {

  x |>
    str_replace_all(
      "[^A-Za-z0-9_-]+",
      "_"
    ) |>
    str_replace_all(
      "_+",
      "_"
    ) |>
    str_remove("^_") |>
    str_remove("_$")
}


make_antibody_plot <- function(
    antibody_name,
    save_plot = TRUE
) {

  plot_data <- dat |>
    filter(
      Antibody == antibody_name
    ) |>
    droplevels()

  if (nrow(plot_data) == 0) {
    warning(
      "No CSF data were found for: ",
      antibody_name
    )
    return(NULL)
  }

  current_annotations <- annotation_data |>
    filter(
      Antibody == antibody_name
    )

  p <- ggplot(
    plot_data,
    aes(
      x = Diagnosis,
      y = Score,
      color = Diagnosis
    )
  ) +
    geom_boxplot(
      width = 0.58,
      fill = "white",
      linewidth = 1.05,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.055,
      height = 0,
      size = 2.6,
      alpha = 0.80
    ) +
    facet_wrap(
      ~ Isotype,
      scales = "free_y",
      nrow = 1,
      drop = FALSE
    ) +
    scale_color_manual(
      values = diagnosis_colors,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "Anti–",
        antibody_name,
        " Antibody Scores in CSF"
      ),
      x = "Diagnosis",
      y = "Antibody Score"
    ) +
    theme_bw(
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        size = 21,
        face = "plain",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      axis.title.x = element_text(
        size = 17,
        margin = margin(t = 8)
      ),
      axis.title.y = element_text(
        size = 17,
        margin = margin(r = 8)
      ),
      axis.text.x = element_text(
        size = 14,
        color = "grey20"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "grey20"
      ),
      strip.background = element_rect(
        fill = "grey85",
        color = "grey20",
        linewidth = 0.7
      ),
      strip.text = element_text(
        size = 14,
        face = "bold",
        color = "grey10"
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "grey20",
        fill = NA,
        linewidth = 0.7
      ),
      legend.position = "none",
      plot.margin = margin(
        t = 10,
        r = 12,
        b = 10,
        l = 12
      )
    )

  # Adds significance brackets only to significant isotype panels.
  if (nrow(current_annotations) > 0) {

    p <- p +

      geom_segment(
        data = current_annotations,
        aes(
          x = 1,
          xend = 2,
          y = bracket_height,
          yend = bracket_height
        ),
        inherit.aes = FALSE,
        linewidth = 0.7,
        color = "black"
      ) +

      geom_segment(
        data = current_annotations,
        aes(
          x = 1,
          xend = 1,
          y = bracket_height,
          yend = bracket_height - tick_height
        ),
        inherit.aes = FALSE,
        linewidth = 0.7,
        color = "black"
      ) +

      geom_segment(
        data = current_annotations,
        aes(
          x = 2,
          xend = 2,
          y = bracket_height,
          yend = bracket_height - tick_height
        ),
        inherit.aes = FALSE,
        linewidth = 0.7,
        color = "black"
      ) +

      geom_text(
        data = current_annotations,
        aes(
          x = 1.5,
          y = label_height,
          label = significance
        ),
        inherit.aes = FALSE,
        size = 5,
        fontface = "bold",
        color = "black"
      )
  }

  if (save_plot) {

    file_stub <- safe_file_name(
      antibody_name
    )

    ggsave(
      filename = file.path(
        plot_dir,
        paste0(
          file_stub,
          "_CSF_plot.pdf"
        )
      ),
      plot = p,
      width = 8.5,
      height = 8,
      units = "in",
      device = cairo_pdf
    )

    ggsave(
      filename = file.path(
        plot_dir,
        paste0(
          file_stub,
          "_CSF_plot.png"
        )
      ),
      plot = p,
      width = 8.5,
      height = 8,
      units = "in",
      dpi = 600,
      bg = "white"
    )
  }

  p
}


# ------------------------------------------------------------
# 11. CHOOSE WHICH FIGURES TO GENERATE
# ------------------------------------------------------------

# Options:
#
#   "significant_only"
#       Plot antibodies with at least one significant NPSLE-versus-
#       OIND comparison in IgA or IgG.
#
#   "all"
#       Plot every antibody.

plot_mode <- "significant_only"


if (plot_mode == "significant_only") {

  antibodies_to_plot <- significant_tukey |>
    distinct(Antibody) |>
    arrange(Antibody) |>
    pull(Antibody)

} else if (plot_mode == "all") {

  antibodies_to_plot <- dat |>
    distinct(Antibody) |>
    arrange(Antibody) |>
    pull(Antibody)

} else {

  stop(
    'plot_mode must be "significant_only" or "all".'
  )
}


# ------------------------------------------------------------
# 12. GENERATE FIGURES
# ------------------------------------------------------------

message(
  "Generating ",
  length(antibodies_to_plot),
  " CSF figures..."
)

walk(
  antibodies_to_plot,
  ~ make_antibody_plot(
    antibody_name = .x,
    save_plot = TRUE
  )
)


# ------------------------------------------------------------
# 13. SAVE PLOT MANIFEST
# ------------------------------------------------------------

plot_manifest <- tibble(
  Antibody = antibodies_to_plot
) |>
  mutate(
    File_Stub = map_chr(
      Antibody,
      safe_file_name
    ),
    PDF_File = paste0(
      File_Stub,
      "_CSF_plot.pdf"
    ),
    PNG_File = paste0(
      File_Stub,
      "_CSF_plot.png"
    )
  )

write_csv(
  plot_manifest,
  file.path(
    output_dir,
    "CSF_plot_manifest.csv"
  )
)


# ------------------------------------------------------------
# 14. SAVE SESSION INFORMATION
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "CSF_boxplot_R_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 15. PRINT SUMMARY
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "CSF analysis complete\n",
  "============================================================\n",
  "Output directory:\n",
  output_dir,
  "\n\n",
  "Antibodies analyzed: ",
  n_distinct(dat$Antibody),
  "\n",
  "Significant NPSLE-versus-OIND comparisons: ",
  nrow(significant_tukey),
  "\n",
  "Figures generated: ",
  length(antibodies_to_plot),
  "\n",
  sep = ""
)

print(
  significant_tukey,
  n = Inf
)
