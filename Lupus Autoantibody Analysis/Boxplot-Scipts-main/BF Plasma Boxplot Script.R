# ============================================================
# Lupus autoantibody boxplot analysis

# This script:
#   1. Imports separate IgA, IgG, and IgM R-ready CSV files.
#   2. Combines and validates the data.
#   3. Runs a two-way ANOVA for each antibody:
#
#        Score ~ Sample_Type * Isotype
#
#   4. Runs lupus-only Tukey HSD testing separately for each
#      antibody and isotype:
#
#        Score ~ Sample_Type
#
#   5. Saves significant plasma-versus-blister comparisons.
#   6. Creates one figure per antibody with IgA, IgG, and IgM
#      displayed side by side.
# ============================================================


# ------------------------------------------------------------
# 1. INSTALL AND LOAD REQUIRED PACKAGE
# ------------------------------------------------------------

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(tidyverse)


# ------------------------------------------------------------
# 2. SELECT INPUT FILES
# ------------------------------------------------------------

message("Select the IgA R-ready CSV file.")
iga_file <- file.choose()

message("Select the IgG R-ready CSV file.")
igg_file <- file.choose()

message("Select the IgM R-ready CSV file.")
igm_file <- file.choose()

input_files <- c(
  IgA = iga_file,
  IgG = igg_file,
  IgM = igm_file
)

output_dir <- file.path(
  dirname(iga_file),
  "Boxplot_Analysis_Outputs"
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
      str_to_upper(Isotype) == "IGM" ~ "IgM",
      TRUE ~ Isotype
    ),
    Diagnosis = case_when(
      str_to_lower(Diagnosis) == "lupus" ~ "Lupus",
      str_to_lower(Diagnosis) == "healthy" ~ "Healthy",
      TRUE ~ Diagnosis
    ),
    Sample_Type = case_when(
      str_to_lower(Sample_Type) %in%
        c("blister", "blisterfluid", "blister fluid") ~ "Blister",
      str_to_lower(Sample_Type) == "plasma" ~ "Plasma",
      TRUE ~ Sample_Type
    )
  ) |>
  filter(
    !is.na(Score),
    !is.na(Antibody),
    Antibody != ""
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
      levels = c("IgA", "IgG", "IgM")
    ),
    Sample_Type = factor(
      Sample_Type,
      levels = c("Blister", "Plasma")
    )
  )


# ------------------------------------------------------------
# 6. BASIC QUALITY-CONTROL OUTPUTS
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
    "One or more selected files may contain the wrong isotype. ",
    "Review isotype_check.csv."
  )
}

group_counts <- dat |>
  count(
    Diagnosis,
    Sample_Type,
    Isotype,
    name = "n_observations"
  )

write_csv(
  group_counts,
  file.path(
    output_dir,
    "group_counts.csv"
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

# The structure of the supplied ANOVA table is reproduced by:
#
#   aov(Score ~ Sample_Type * Isotype)
#
# This model creates rows for:
#   Sample_Type
#   Isotype
#   Sample_Type:Isotype
#   Residuals

run_two_way_anova <- function(x) {

  model_data <- x |>
    drop_na(
      Score,
      Sample_Type,
      Isotype
    ) |>
    droplevels()

  if (
    n_distinct(model_data$Sample_Type) < 2 ||
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
    Score ~ Sample_Type * Isotype,
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
    "anova_results_labelled.csv"
  )
)


# Save a second version matching the supplied original layout.
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
    "anova_results.csv"
  )
)


# ------------------------------------------------------------
# 8. LUPUS-ONLY TUKEY HSD TESTS
# ------------------------------------------------------------

# Each antibody and isotype is analyzed separately.
# Only lupus blister and lupus plasma observations are included.

lupus_data <- dat |>
  filter(
    Diagnosis == "Lupus",
    Sample_Type %in% c("Blister", "Plasma")
  ) |>
  droplevels()


run_tukey <- function(x) {

  model_data <- x |>
    drop_na(
      Score,
      Sample_Type
    ) |>
    droplevels()

  if (n_distinct(model_data$Sample_Type) < 2) {
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
    Score ~ Sample_Type,
    data = model_data
  )

  TukeyHSD(
    model,
    "Sample_Type"
  )$Sample_Type |>
    as.data.frame() |>
    rownames_to_column("Raw_Comparison") |>
    as_tibble() |>
    mutate(
      Comparison = case_when(
        Raw_Comparison == "Plasma-Blister" ~
          "Lupus Plasma-Lupus Blister",
        TRUE ~ Raw_Comparison
      )
    ) |>
    select(
      Comparison,
      diff,
      lwr,
      upr,
      `p adj`
    )
}


all_tukey_results <- lupus_data |>
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
    "All_Tukey_Comparisons_LupusOnly.csv"
  )
)


significant_tukey <- all_tukey_results |>
  filter(
    Comparison == "Lupus Plasma-Lupus Blister",
    !is.na(`p adj`),
    `p adj` < 0.05
  )

write_csv(
  significant_tukey,
  file.path(
    output_dir,
    "Significant_Tukey_Comparisons_LupusOnly.csv"
  )
)


# ------------------------------------------------------------
# 9. PREPARE SIGNIFICANCE-BRACKET DATA
# ------------------------------------------------------------

# Brackets are drawn manually with geom_segment(), avoiding the
# need for ggpubr or ggsignif.

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
    lupus_data |>
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

sample_colors <- c(
  "Blister" = "#D73027",
  "Plasma" = "#4575B4"
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

  plot_data <- lupus_data |>
    filter(
      Antibody == antibody_name
    ) |>
    droplevels()

  if (nrow(plot_data) == 0) {
    warning(
      "No lupus data were found for: ",
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
      x = Sample_Type,
      y = Score,
      color = Sample_Type
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
      values = sample_colors,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "Anti–",
        antibody_name,
        " Antibody Scores in Lupus Patients"
      ),
      x = "Sample Type",
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

  # Add significance brackets only where p < 0.05.
  if (nrow(current_annotations) > 0) {

    p <- p +

      # Horizontal bracket line
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

      # Left vertical bracket tick
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

      # Right vertical bracket tick
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

      # Asterisk label
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
          "_plot.pdf"
        )
      ),
      plot = p,
      width = 12,
      height = 8,
      units = "in",
      device = cairo_pdf
    )

    ggsave(
      filename = file.path(
        plot_dir,
        paste0(
          file_stub,
          "_plot.png"
        )
      ),
      plot = p,
      width = 12,
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
#   plot_mode <- "significant_only"
#       Generates figures only for antibodies with at least one
#       significant lupus blister-versus-plasma comparison.
#
#   plot_mode <- "all"
#       Generates a figure for every antibody.

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
  " figures..."
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
      "_plot.pdf"
    ),
    PNG_File = paste0(
      File_Stub,
      "_plot.png"
    )
  )

write_csv(
  plot_manifest,
  file.path(
    output_dir,
    "plot_manifest.csv"
  )
)


# ------------------------------------------------------------
# 14. SAVE SESSION INFORMATION
# ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "boxplot_R_sessionInfo_version3.txt"
  )
)


# ------------------------------------------------------------
# 15. PRINT SUMMARY
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "Analysis complete\n",
  "============================================================\n",
  "Output directory:\n",
  output_dir,
  "\n\n",
  "Antibodies analyzed: ",
  n_distinct(dat$Antibody),
  "\n",
  "Significant Tukey comparisons: ",
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
