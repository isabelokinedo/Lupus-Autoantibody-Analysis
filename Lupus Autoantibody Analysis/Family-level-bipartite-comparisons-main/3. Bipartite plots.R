############################################################
# COMPARATIVE THREE-PANEL BIPARTITE NETWORK
#
# Panels:
#   CSF | Blister fluid | Plasma
#
# Input:
#   Three network-family-summary CSV files for the same
#   antibody isotype.
#
# Example:
#   CSF_IgG_..._network_family_summary.csv
#   BF_IgG_..._network_family_summary.csv
#   Plasma_IgG_..._network_family_summary.csv
#
# Edge width = absolute median Spearman rho
# Edge color = direction of median Spearman rho
# Edge alpha = direction consistency
############################################################


# ==========================================================
# 1. INSTALL AND LOAD PACKAGES
# ==========================================================

required_packages <- c(
  "tidyverse",
  "igraph",
  "tidygraph",
  "ggraph",
  "patchwork"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(patchwork)


# ==========================================================
# 2. SELECT THE THREE INPUT FILES
# ==========================================================

cat(
  "\nSelect the CSF network-family-summary file.\n"
)

csf_file <- file.choose()


cat(
  "\nSelect the blister-fluid network-family-summary file.\n"
)

bf_file <- file.choose()


cat(
  "\nSelect the plasma network-family-summary file.\n"
)

plasma_file <- file.choose()


# ==========================================================
# 3. USER SETTINGS
# ==========================================================
#
# Adjust these thresholds if the networks are too sparse
# or too crowded.
# ==========================================================

minimum_pairs <- 3

minimum_median_abs_rho <- 0.40

minimum_direction_consistency <- 0.65

minimum_proportion_strong <- 0.20

# Maximum number of retained edges per tissue
maximum_edges_per_tissue <- 15

# Enter the isotype represented by the three selected files
isotype_label <- "IgG"

# Set TRUE to require at least one FDR-significant
# individual antibody-protein correlation within the
# family pair.
require_fdr_support <- FALSE

minimum_fdr_significant_pairs <- 1


# ==========================================================
# 4. CREATE OUTPUT DIRECTORY
# ==========================================================

output_directory <- file.path(
  dirname(csf_file),
  paste0(
    "Comparative_",
    isotype_label,
    "_bipartite_network"
  )
)

if (!dir.exists(output_directory)) {
  dir.create(
    output_directory,
    recursive = TRUE
  )
}

cat(
  "\nOutput directory:\n",
  output_directory,
  "\n\n"
)


# ==========================================================
# 5. FUNCTION TO IMPORT AND STANDARDIZE ONE SUMMARY
# ==========================================================

read_network_summary <- function(
    file,
    tissue_label
) {
  
  data <- read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  
  required_columns <- c(
    "AntibodyFamily",
    "ProteinFamily",
    "Number_of_pairs",
    "Median_rho",
    "Median_abs_rho",
    "Direction_consistency",
    "Proportion_strong"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  
  
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "\nThe following columns are missing from ",
        basename(file),
        ":\n",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
  
  
  # FDR column is optional unless FDR filtering is requested
  if (
    require_fdr_support &&
    !"Significant_FDR_count" %in% names(data)
  ) {
    stop(
      paste0(
        basename(file),
        " does not contain Significant_FDR_count, ",
        "but require_fdr_support is TRUE."
      )
    )
  }
  
  
  if (!"Significant_FDR_count" %in% names(data)) {
    data$Significant_FDR_count <- NA_real_
  }
  
  
  data %>%
    transmute(
      Tissue = tissue_label,
      
      AntibodyFamily = trimws(
        as.character(AntibodyFamily)
      ),
      
      ProteinFamily = trimws(
        as.character(ProteinFamily)
      ),
      
      Number_of_pairs = suppressWarnings(
        as.numeric(Number_of_pairs)
      ),
      
      Median_rho = suppressWarnings(
        as.numeric(Median_rho)
      ),
      
      Median_abs_rho = suppressWarnings(
        as.numeric(Median_abs_rho)
      ),
      
      Direction_consistency = suppressWarnings(
        as.numeric(Direction_consistency)
      ),
      
      Proportion_strong = suppressWarnings(
        as.numeric(Proportion_strong)
      ),
      
      Significant_FDR_count = suppressWarnings(
        as.numeric(Significant_FDR_count)
      )
    ) %>%
    filter(
      !is.na(AntibodyFamily),
      AntibodyFamily != "",
      !is.na(ProteinFamily),
      ProteinFamily != "",
      is.finite(Median_rho),
      is.finite(Median_abs_rho)
    ) %>%
    mutate(
      Correlation_direction = case_when(
        Median_rho > 0 ~ "Positive",
        Median_rho < 0 ~ "Negative",
        TRUE ~ "Zero"
      ),
      
      # Used only to rank retained edges after the
      # interpretable thresholds are applied
      Ranking_score =
        Median_abs_rho *
        Direction_consistency *
        Proportion_strong
    )
}


# ==========================================================
# 6. IMPORT ALL THREE TISSUES
# ==========================================================

csf_summary <- read_network_summary(
  csf_file,
  "CSF"
)

bf_summary <- read_network_summary(
  bf_file,
  "Blister fluid"
)

plasma_summary <- read_network_summary(
  plasma_file,
  "Plasma"
)


all_summaries <- bind_rows(
  csf_summary,
  bf_summary,
  plasma_summary
)


# ==========================================================
# 7. FILTER FAMILY-LEVEL ASSOCIATIONS
# ==========================================================

filtered_summaries <- all_summaries %>%
  filter(
    Number_of_pairs >= minimum_pairs,
    Median_abs_rho >= minimum_median_abs_rho,
    Direction_consistency >=
      minimum_direction_consistency,
    Proportion_strong >=
      minimum_proportion_strong
  )


if (require_fdr_support) {
  
  filtered_summaries <- filtered_summaries %>%
    filter(
      Significant_FDR_count >=
        minimum_fdr_significant_pairs
    )
}


# Keep the highest-ranked edges separately within each tissue
retained_edges <- filtered_summaries %>%
  group_by(Tissue) %>%
  arrange(
    desc(Ranking_score),
    desc(Median_abs_rho),
    desc(Number_of_pairs),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = maximum_edges_per_tissue
  ) %>%
  ungroup()


if (nrow(retained_edges) == 0) {
  stop(
    paste0(
      "No family relationships passed the current filters. ",
      "Try lowering minimum_median_abs_rho, ",
      "minimum_direction_consistency, or ",
      "minimum_proportion_strong."
    )
  )
}


cat(
  "\nRetained edges by tissue:\n"
)

print(
  retained_edges %>%
    count(Tissue)
)


# ==========================================================
# 8. SAVE RETAINED EDGE TABLES
# ==========================================================

write.csv(
  retained_edges,
  file.path(
    output_directory,
    paste0(
      isotype_label,
      "_all_retained_family_edges.csv"
    )
  ),
  row.names = FALSE
)


retained_edges %>%
  group_split(Tissue) %>%
  walk(
    function(x) {
      
      tissue_name <- unique(x$Tissue) %>%
        str_replace_all(
          "[^A-Za-z0-9]+",
          "_"
        )
      
      write.csv(
        x,
        file.path(
          output_directory,
          paste0(
            isotype_label,
            "_",
            tissue_name,
            "_retained_edges.csv"
          )
        ),
        row.names = FALSE
      )
    }
  )


# ==========================================================
# 9. DETERMINE THE COMMON NODE SET
# ==========================================================
#
# A union of nodes is used so that a family seen in any
# tissue receives the same position in every panel.
# ==========================================================

all_antibody_families <- retained_edges %>%
  distinct(AntibodyFamily) %>%
  pull(AntibodyFamily) %>%
  sort()


all_protein_families <- retained_edges %>%
  distinct(ProteinFamily) %>%
  pull(ProteinFamily) %>%
  sort()


# ==========================================================
# 10. ORDER NODES USING THEIR OVERALL EVIDENCE
# ==========================================================
#
# Each family receives a score across all tissues.
# This gives one fixed vertical ordering for every panel.
# ==========================================================

antibody_order <- retained_edges %>%
  group_by(AntibodyFamily) %>%
  summarise(
    Total_evidence = sum(
      Ranking_score,
      na.rm = TRUE
    ),
    
    Strongest_edge = max(
      Median_abs_rho,
      na.rm = TRUE
    ),
    
    Tissue_count = n_distinct(Tissue),
    
    .groups = "drop"
  ) %>%
  arrange(
    desc(Tissue_count),
    desc(Total_evidence),
    desc(Strongest_edge),
    AntibodyFamily
  ) %>%
  pull(AntibodyFamily)


protein_order <- retained_edges %>%
  group_by(ProteinFamily) %>%
  summarise(
    Total_evidence = sum(
      Ranking_score,
      na.rm = TRUE
    ),
    
    Strongest_edge = max(
      Median_abs_rho,
      na.rm = TRUE
    ),
    
    Tissue_count = n_distinct(Tissue),
    
    .groups = "drop"
  ) %>%
  arrange(
    desc(Tissue_count),
    desc(Total_evidence),
    desc(Strongest_edge),
    ProteinFamily
  ) %>%
  pull(ProteinFamily)


# ==========================================================
# 11. CREATE FIXED MANUAL COORDINATES
# ==========================================================
#
# Antibody families are placed on the left.
# Olink protein families are placed on the right.
# ==========================================================

antibody_coordinates <- tibble(
  name = paste0(
    "Antibody__",
    antibody_order
  ),
  
  label = antibody_order,
  
  node_type = "Autoantibody family",
  
  x = -1,
  
  y = seq(
    from = 1,
    to = -1,
    length.out = length(
      antibody_order
    )
  )
)


protein_coordinates <- tibble(
  name = paste0(
    "Protein__",
    protein_order
  ),
  
  label = protein_order,
  
  node_type = "Olink protein family",
  
  x = 1,
  
  y = seq(
    from = 1,
    to = -1,
    length.out = length(
      protein_order
    )
  )
)


common_node_coordinates <- bind_rows(
  antibody_coordinates,
  protein_coordinates
)


# ==========================================================
# 12. CREATE ONE TISSUE-SPECIFIC BIPARTITE PLOT
# ==========================================================

create_tissue_plot <- function(
    tissue_name,
    show_left_labels = TRUE,
    show_right_labels = TRUE,
    show_legend = FALSE
) {
  
  tissue_edges <- retained_edges %>%
    filter(
      Tissue == tissue_name
    ) %>%
    transmute(
      from = paste0(
        "Antibody__",
        AntibodyFamily
      ),
      
      to = paste0(
        "Protein__",
        ProteinFamily
      ),
      
      median_rho = Median_rho,
      
      abs_median_rho =
        Median_abs_rho,
      
      direction_consistency =
        Direction_consistency,
      
      proportion_strong =
        Proportion_strong,
      
      pair_count =
        Number_of_pairs,
      
      correlation_direction =
        Correlation_direction,
      
      ranking_score =
        Ranking_score
    )
  
  
  # Keep all common nodes, including isolated nodes.
  # This preserves fixed positions across panels.
  tissue_nodes <- common_node_coordinates %>%
    mutate(
      node_present =
        name %in%
        c(
          tissue_edges$from,
          tissue_edges$to
        )
    )
  
  
  tissue_graph <- graph_from_data_frame(
    d = tissue_edges,
    directed = FALSE,
    vertices = tissue_nodes
  )
  
  
  V(tissue_graph)$node_degree <- as.numeric(
    igraph::degree(
      tissue_graph
    )
  )
  
  
  graph_tbl <- as_tbl_graph(
    tissue_graph
  )
  
  
  vertex_order <- match(
    V(tissue_graph)$name,
    common_node_coordinates$name
  )
  
  
  x_coordinates <-
    common_node_coordinates$x[
      vertex_order
    ]
  
  y_coordinates <-
    common_node_coordinates$y[
      vertex_order
    ]
  
  
  label_data <- common_node_coordinates %>%
    mutate(
      node_present =
        name %in%
        c(
          tissue_edges$from,
          tissue_edges$to
        ),
      
      show_label = case_when(
        node_type ==
          "Autoantibody family" ~
          show_left_labels,
        
        node_type ==
          "Olink protein family" ~
          show_right_labels,
        
        TRUE ~ FALSE
      ),
      
      label_x = case_when(
        node_type ==
          "Autoantibody family" ~
          x - 0.08,
        
        node_type ==
          "Olink protein family" ~
          x + 0.08
      ),
      
      label_hjust = case_when(
        node_type ==
          "Autoantibody family" ~ 1,
        
        node_type ==
          "Olink protein family" ~ 0
      )
    ) %>%
    filter(
      show_label
    )
  
  
  ggraph(
    graph_tbl,
    layout = "manual",
    x = x_coordinates,
    y = y_coordinates
  ) +
    
    geom_edge_link(
      aes(
        edge_width =
          abs_median_rho,
        
        edge_color =
          correlation_direction,
        
        edge_alpha =
          direction_consistency
      ),
      lineend = "round"
    ) +
    
    geom_node_point(
      aes(
        size = node_degree,
        fill = node_type,
        alpha = node_present
      ),
      shape = 21,
      color = "black",
      stroke = 0.3
    ) +
    
    geom_text(
      data = label_data,
      aes(
        x = label_x,
        y = y,
        label = label,
        hjust = label_hjust
      ),
      inherit.aes = FALSE,
      size = 3.1,
      fontface = "bold"
    ) +
    
    scale_fill_manual(
      values = c(
        "Autoantibody family" =
          "#D95F02",
        
        "Olink protein family" =
          "#1B9E77"
      ),
      name = "Node"
    ) +
    
    scale_edge_color_manual(
      values = c(
        "Positive" = "#B2182B",
        "Negative" = "#2166AC"
      ),
      breaks = c(
        "Negative",
        "Positive"
      ),
      limits = c(
        "Negative",
        "Positive"
      ),
      drop = FALSE,
      name = "Direction"
    ) +
    
    # Fixed limits are crucial for comparing tissues
    scale_edge_width_continuous(
      limits = c(
        minimum_median_abs_rho,
        1
      ),
      range = c(
        0.7,
        4.5
      ),
      name = "Median |ρ|"
    ) +
    
    scale_edge_alpha_continuous(
      limits = c(
        minimum_direction_consistency,
        1
      ),
      range = c(
        0.30,
        0.95
      ),
      name = "Direction consistency"
    ) +
    
    scale_size_continuous(
      range = c(
        3,
        8
      ),
      guide = "none"
    ) +
    
    scale_alpha_manual(
      values = c(
        "TRUE" = 1,
        "FALSE" = 0.35
      ),
      guide = "none"
    ) +
    
    guides(
      fill = guide_legend(
        order = 1,
        title.position = "top",
        override.aes = list(
          size = 6,
          alpha = 1
        )
      ),
      
      edge_color = guide_legend(
        order = 2,
        title.position = "top",
        override.aes = list(
          edge_width = 3,
          edge_alpha = 1
        )
      ),
      
      edge_width = guide_legend(
        order = 3,
        title.position = "top",
        nrow = 1
      ),
      
      edge_alpha = guide_legend(
        order = 4,
        title.position = "top",
        nrow = 1,
        override.aes = list(
          edge_width = 3
        )
      )
    ) +
    
    labs(
      title = tissue_name
    ) +
    
    theme_graph(
      base_family = "sans"
    ) +
    
    theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      ),
      
      legend.position = if (show_legend) "bottom" else "none",
      legend.justification = "center",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.direction = "horizontal",
      
      plot.margin = margin(
        t = 15,
        r = ifelse(
          show_right_labels,
          150,
          15
        ),
        b = 15,
        l = ifelse(
          show_left_labels,
          150,
          15
        )
      )
    ) +
    
    coord_cartesian(
      xlim = c(
        -1.15,
        1.15
      ),
      ylim = c(
        -1.08,
        1.08
      ),
      clip = "off"
    )
}



# ==========================================================
# 13. GENERATE THE THREE PANELS, THEN COMBINE
# ==========================================================
#
# The first panel shows antibody-family labels.
# The final panel shows Olink-family labels.
# ==========================================================

csf_plot <- create_tissue_plot(
  tissue_name = "CSF",
  show_left_labels = TRUE,
  show_right_labels = FALSE,
  show_legend = FALSE
)

bf_plot <- create_tissue_plot(
  tissue_name = "Blister fluid",
  show_left_labels = FALSE,
  show_right_labels = FALSE,
  show_legend = TRUE
)

plasma_plot <- create_tissue_plot(
  tissue_name = "Plasma",
  show_left_labels = FALSE,
  show_right_labels = TRUE,
  show_legend = FALSE
)


combined_plot <- (
  csf_plot |
    bf_plot |
    plasma_plot
) +
  
  plot_layout(
    widths = c(
      1.15,
      1,
      1.15
    )
  ) +
  
  plot_annotation(
    title = paste0(
      isotype_label,
      " autoantibody-family and Olink-family associations across tissues"
    ),
    
    subtitle = paste0(
      "Edges retained when: number of pairs ≥ ",
      minimum_pairs,
      "; median |Spearman ρ| ≥ ",
      minimum_median_abs_rho,
      "; direction consistency ≥ ",
      minimum_direction_consistency,
      "; proportion strong ≥ ",
      minimum_proportion_strong,
      ". Maximum ",
      maximum_edges_per_tissue,
      " edges per tissue."
    ),
    
    caption = paste0(
      "Edge width represents median absolute Spearman ρ; ",
      "red indicates positive and blue indicates negative median correlations.\n",
      "Edge transparency represents direction consistency. ",
      "Faded nodes preserve identical family positions across tissues but have no retained edge."
    ),
    
    theme = theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      
      plot.caption = element_text(
        size = 9,
        hjust = 0.5,
        margin = margin(t = 12)
      ),
      
      plot.margin = margin(
        t = 15,
        r = 20,
        b = 35,
        l = 20
      )
    )
  )


# ==========================================================
# 15. SAVE COMPARATIVE FIGURE
# ==========================================================

ggsave(
  filename = file.path(
    output_directory,
    paste0(
      isotype_label,
      "_comparative_bipartite_network.pdf"
    )
  ),
  plot = combined_plot,
  width = 20,
  height = 12,
  units = "in",
  device = grDevices::pdf
)


ggsave(
  filename = file.path(
    output_directory,
    paste0(
      isotype_label,
      "_comparative_bipartite_network.png"
    )
  ),
  plot = combined_plot,
  width = 20,
  height = 12,
  units = "in",
  dpi = 600
)


# ==========================================================
# 16. SAVE COMMON NODE ORDER
# ==========================================================

write.csv(
  common_node_coordinates,
  file.path(
    output_directory,
    paste0(
      isotype_label,
      "_common_node_positions.csv"
    )
  ),
  row.names = FALSE
)


# ==========================================================
# 17. CREATE CROSS-TISSUE MODULE SUMMARY
# ==========================================================
#
# This identifies whether each family pair is retained in:
#   one tissue only,
#   two tissues,
#   or all three tissues.
# ==========================================================

cross_tissue_summary <- retained_edges %>%
  group_by(
    AntibodyFamily,
    ProteinFamily
  ) %>%
  summarise(
    Tissues_retained = paste(
      sort(
        unique(Tissue)
      ),
      collapse = "; "
    ),
    
    Number_of_tissues =
      n_distinct(Tissue),
    
    Mean_median_rho =
      mean(
        Median_rho,
        na.rm = TRUE
      ),
    
    Maximum_abs_median_rho =
      max(
        Median_abs_rho,
        na.rm = TRUE
      ),
    
    Mean_direction_consistency =
      mean(
        Direction_consistency,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Module_class = case_when(
      Number_of_tissues == 3 ~
        "Shared across all tissues",
      
      Number_of_tissues == 2 ~
        "Shared across two tissues",
      
      Number_of_tissues == 1 ~
        "Tissue-specific",
      
      TRUE ~
        "Other"
    )
  ) %>%
  arrange(
    desc(Number_of_tissues),
    desc(Maximum_abs_median_rho)
  )


write.csv(
  cross_tissue_summary,
  file.path(
    output_directory,
    paste0(
      isotype_label,
      "_cross_tissue_module_summary.csv"
    )
  ),
  row.names = FALSE
)


cat(
  "\n============================================\n",
  "COMPARATIVE BIPARTITE NETWORK COMPLETE\n",
  "============================================\n\n"
)

cat(
  "Results saved in:\n",
  output_directory,
  "\n"
)
