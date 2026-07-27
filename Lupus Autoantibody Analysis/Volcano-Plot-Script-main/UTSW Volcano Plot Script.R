# ============================================================
# Reusable autoantibody volcano-plot script
#
# Available comparisons:
#   "CSF"             = NPSLE CSF versus OIND CSF
#   "CLE_blister"     = Lupus blister versus healthy blister
#   "CLE_compartment" = Lupus blister versus lupus plasma
# ============================================================

library(tidyverse)
library(ggrepel)

# ------------------------------------------------------------
# 1. CHOOSE THE ANALYSIS
# ------------------------------------------------------------

comparison <- "CLE_blister"

# Other options:
# comparison <- "CSF"
# comparison <- "CLE_compartment"


# ------------------------------------------------------------
# 2. MANUALLY SELECT R-READY CSV FILE
# ------------------------------------------------------------

file_path <- file.choose()
output_dir <- dirname(file_path)

dat <- read_csv(
  file_path,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 3. CLEAN SLASH VALIDATE DATA
# ------------------------------------------------------------

required_columns <- c(
  "Sample_ID",
  "Antibody",
  "Isotype",
  "Diagnosis",
  "Sample_Type",
  "Score"
)

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
  stop(
    paste(
      "The following required columns are missing:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

dat <- dat |>
  mutate(
    across(
      c(Sample_ID, Antibody, Isotype, Diagnosis, Sample_Type),
      str_squish
    ),
    Score = as.numeric(Score)
  ) |>
  filter(
    !str_detect(Antibody, "^HIgControl"),
    !str_detect(Antibody, "^MIgControl")
  )

# Confirm that the selected file contains only one isotype
isotypes_present <- unique(dat$Isotype)

if (length(isotypes_present) != 1) {
  stop(
    paste(
      "The selected file contains more than one isotype:",
      paste(isotypes_present, collapse = ", ")
    )
  )
}

isotype_name <- isotypes_present[[1]]


# ------------------------------------------------------------
# 4. DEFINE SELECTED COMPARISON
# ------------------------------------------------------------

if (comparison == "CSF") {
  
  dat_analysis <- dat |>
    filter(
      Sample_Type == "CSF",
      Diagnosis %in% c("OIND", "NPSLE")
    ) |>
    mutate(
      Group = factor(
        Diagnosis,
        levels = c("OIND", "NPSLE")
      )
    )
  
  reference_group <- "OIND"
  comparison_group <- "NPSLE"
  
  left_label <- "Higher in OIND"
  right_label <- "Higher in NPSLE"
  
  plot_title <- paste0(
    "CSF ",
    isotype_name,
    " autoantibody reactivity in NPSLE versus OIND"
  )
  
  x_axis_label <- "Median antibody score difference (NPSLE − OIND)"
  
  output_prefix <- paste0(
    "CSF_",
    isotype_name,
    "_NPSLE_vs_OIND"
  )
  
} else if (comparison == "CLE_blister") {
  
  dat_analysis <- dat |>
    filter(
      Sample_Type == "Blister",
      Diagnosis %in% c("Healthy", "Lupus")
    ) |>
    mutate(
      Group = factor(
        Diagnosis,
        levels = c("Healthy", "Lupus")
      )
    )
  
  reference_group <- "Healthy"
  comparison_group <- "Lupus"
  
  left_label <- "Higher in healthy blister"
  right_label <- "Higher in CLE blister"
  
  plot_title <- paste0(
    "Blister-fluid ",
    isotype_name,
    " autoantibody reactivity in CLE versus healthy controls"
  )
  
  x_axis_label <- "Median antibody score difference (CLE − healthy)"
  
  output_prefix <- paste0(
    "Blister_",
    isotype_name,
    "_CLE_vs_Healthy"
  )
  
} else if (comparison == "CLE_compartment") {
  
  dat_analysis <- dat |>
    filter(
      Diagnosis == "Lupus",
      Sample_Type %in% c("Plasma", "Blister")
    ) |>
    mutate(
      Group = factor(
        Sample_Type,
        levels = c("Plasma", "Blister")
      )
    )
  
  reference_group <- "Plasma"
  comparison_group <- "Blister"
  
  left_label <- "Higher in plasma"
  right_label <- "Higher in blister fluid"
  
  plot_title <- paste0(
    isotype_name,
    " autoantibody reactivity in CLE blister fluid versus plasma"
  )
  
  x_axis_label <- "Median antibody score difference (blister − plasma)"
  
  output_prefix <- paste0(
    "CLE_",
    isotype_name,
    "_Blister_vs_Plasma"
  )
  
} else {
  
  stop(
    paste0(
      "Invalid comparison. Choose one of: ",
      "'CSF', 'CLE_blister', or 'CLE_compartment'."
    )
  )
}


# ------------------------------------------------------------
# 5. QUALITY-CONTROL CHECKS
# ------------------------------------------------------------

cat("\nSelected comparison:", comparison, "\n")
cat("Isotype:", isotype_name, "\n\n")

print(table(dat_analysis$Group))

cat(
  "\nUnique samples:",
  n_distinct(dat_analysis$Sample_ID),
  "\n"
)

cat(
  "Unique antibodies:",
  n_distinct(dat_analysis$Antibody),
  "\n"
)

duplicates <- dat_analysis |>
  count(Sample_ID, Antibody) |>
  filter(n > 1)

if (nrow(duplicates) > 0) {
  warning("Duplicate sample-antibody combinations were found.")
  print(duplicates)
}

if (n_distinct(dat_analysis$Group) != 2) {
  stop(
    "Both required comparison groups were not found in the selected file."
  )
}


# ------------------------------------------------------------
# 6. MANN-WHITNEY TESTS FOR EACH ANTIBODY
# ------------------------------------------------------------

volcano_results <- dat_analysis |>
  group_by(Antibody) |>
  summarise(
    n_reference = sum(
      Group == reference_group & !is.na(Score)
    ),
    
    n_comparison = sum(
      Group == comparison_group & !is.na(Score)
    ),
    
    median_reference = median(
      Score[Group == reference_group],
      na.rm = TRUE
    ),
    
    median_comparison = median(
      Score[Group == comparison_group],
      na.rm = TRUE
    ),
    
    median_difference =
      median_comparison - median_reference,
    
    p_value = tryCatch(
      wilcox.test(
        Score[Group == comparison_group],
        Score[Group == reference_group],
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    ),
    
    .groups = "drop"
  ) |>
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    ),
    
    minus_log10_p = -log10(p_value),
    
    minus_log10_FDR = -log10(FDR),
    
    Direction = case_when(
      median_difference > 0 ~
        paste("Higher in", comparison_group),
      
      median_difference < 0 ~
        paste("Higher in", reference_group),
      
      TRUE ~ "No difference"
    ),
    
    Significance = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      p_value < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "Not significant"
    ),
    
    Label = if_else(
      p_value < 0.05,
      Antibody,
      NA_character_
    )
  ) |>
  mutate(
    Significance = factor(
      Significance,
      levels = c(
        "FDR < 0.05",
        "Nominal p < 0.05",
        "Not significant"
      )
    )
  ) |>
  arrange(p_value)


# ------------------------------------------------------------
# 7. VIEW SLASH EXPORT THE STATS
# ------------------------------------------------------------

print(
  volcano_results |>
    select(
      Antibody,
      median_reference,
      median_comparison,
      median_difference,
      Direction,
      p_value,
      FDR
    ),
  n = 20
)

write_csv(
  volcano_results,
  file.path(
    output_dir,
    paste0(
      output_prefix,
      "_statistics.csv"
    )
  )
)


# ------------------------------------------------------------
# 8. CREATE ACTUAL VOLCANO PLOT
# ------------------------------------------------------------

# Set vertical positions for directional headings and antibody labels
y_data_max <- max(volcano_results$minus_log10_p, na.rm = TRUE)

y_plot_max <- y_data_max * 1.22
y_heading <- y_data_max * 1.16
y_label_max <- y_data_max * 1.08

volcano_plot <- ggplot(
  volcano_results,
  aes(
    x = median_difference,
    y = minus_log10_p
  )
) +
  
  geom_point(
    aes(
      shape = Significance,
      color = Significance
    ),
    size = 2.8,
    alpha = 0.8
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  annotate(
    "text",
    x = -Inf,
    y = y_heading,
    label = left_label,
    hjust = -0.05,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  annotate(
    "text",
    x = Inf,
    y = y_heading,
    label = right_label,
    hjust = 1.05,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  
  geom_text_repel(
    aes(label = Label),
    color = "black",
    segment.color = "red",
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    ylim = c(0, y_label_max),
    show.legend = FALSE
  ) +
  
  scale_shape_manual(
    values = c(
      "FDR < 0.05" = 16,
      "Nominal p < 0.05" = 16,
      "Not significant" = 1
    ),
    drop = TRUE
  ) +
  
  scale_color_manual(
    values = c(
      "FDR < 0.05" = "red",
      "Nominal p < 0.05" = "red",
      "Not significant" = "black"
    ),
    drop = TRUE
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  scale_y_continuous(
    limits = c(0, y_plot_max),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  labs(
    title = plot_title,
    x = x_axis_label,
    y = expression(
      -log[10](italic(p))
    ),
    color = NULL,
    shape = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )

print(volcano_plot)


# ------------------------------------------------------------
# 9. SAVE THE PLOT
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    paste0(
      output_prefix,
      "_volcano.pdf"
    )
  ),
  plot = volcano_plot,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(
    output_dir,
    paste0(
      output_prefix,
      "_volcano.png"
    )
  ),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 600
)

cat(
  "\nFiles saved to:\n",
  output_dir,
  "\n"
)