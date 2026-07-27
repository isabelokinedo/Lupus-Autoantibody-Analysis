############################################################
# Correlation dot plot from Spearman rho matrix
############################################################

library(tidyverse)


############################################################
# 1. Select files
############################################################

message("Select the Spearman rho matrix CSV file")
rho_file <- file.choose()

message("Select Autoantibody_annotations.csv")
ab_annot_file <- file.choose()

message("Select Olink_protein_annotations.csv")
protein_annot_file <- file.choose()


############################################################
# 2. Read files
############################################################

cor_mat <- read.csv(
  rho_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ab_annot <- read.csv(
  ab_annot_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

protein_annot <- read.csv(
  protein_annot_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


############################################################
# 3. Fix first column of correlation matrix
############################################################

# write.csv() saves matrix row names as the first column
# with a blank header or a name such as X
names(cor_mat)[1] <- "Antibody"


############################################################
# 4. Clean names
############################################################

cor_mat <- cor_mat |>
  mutate(
    Antibody = trimws(Antibody)
  )

ab_annot <- ab_annot |>
  mutate(
    Antibody = trimws(Antibody)
  )

protein_annot <- protein_annot |>
  mutate(
    Protein = trimws(Protein)
  )


############################################################
# 5. Convert matrix to long format and add annotations
############################################################

long <- cor_mat |>
  pivot_longer(
    cols = -Antibody,
    names_to = "Protein",
    values_to = "rho"
  ) |>
  left_join(
    ab_annot,
    by = "Antibody"
  ) |>
  rename(
    AntibodyFamily = Family,
    AntibodySubfamily = Subfamily,
    AntibodyOrder = PlotOrder
  ) |>
  left_join(
    protein_annot,
    by = "Protein"
  ) |>
  rename(
    ProteinFamily = Family,
    ProteinSubfamily = Subfamily,
    ProteinOrder = PlotOrder
  )


############################################################
# 6. Check unmatched annotations
############################################################

unmatched_antibodies <- long |>
  filter(is.na(AntibodyFamily)) |>
  distinct(Antibody)

unmatched_proteins <- long |>
  filter(is.na(ProteinFamily)) |>
  distinct(Protein)

message("Antibodies missing annotations:")
print(unmatched_antibodies)

message("Proteins missing annotations:")
print(unmatched_proteins)


############################################################
# 7. Preserve annotation order
############################################################

long <- long |>
  mutate(
    Antibody = factor(
      Antibody,
      levels = ab_annot$Antibody[
        order(ab_annot$PlotOrder)
      ]
    ),
    Protein = factor(
      Protein,
      levels = protein_annot$Protein[
        order(protein_annot$PlotOrder)
      ]
    )
  )


############################################################
# 8. Filter correlations for plotting
############################################################

rho_cutoff <- 0.50

plot_data <- long |>
  filter(
    !is.na(rho),
    !is.na(AntibodyFamily),
    !is.na(ProteinFamily),
    abs(rho) >= rho_cutoff
  )


############################################################
# 9. Generate dot-plot heatmap
############################################################

p <- ggplot(
  plot_data,
  aes(
    x = Protein,
    y = Antibody,
    size = abs(rho),
    color = rho
  )
) +
  geom_point(alpha = 0.9) +
  
  scale_size_continuous(
    name = "|Spearman rho|",
    range = c(1.5, 7),
    limits = c(rho_cutoff, 1)
  ) +
  
  scale_color_gradient2(
    name = "Spearman rho",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  
  facet_grid(
    rows = vars(AntibodyFamily),
    cols = vars(ProteinFamily),
    scales = "free",
    space = "free",
    switch = "y",
    labeller = labeller(
      AntibodyFamily = label_wrap_gen(width = 22)
    )
  ) +
  labs(
    x = "Olink protein",
    y = "Autoantibody target"
  ) +
  theme_bw(base_size = 10) +
  theme(
    
    ##########################################################
    # Individual protein and antibody labels
    ##########################################################
    
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 7
    ),
    
    axis.text.y = element_text(
      size = 7
    ),
    
    ##########################################################
    # Protein-family labels across the top
    ##########################################################
    
    strip.text.x = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      size = 8,
      face = "bold",
      margin = margin(
        t = 6,
        r = 4,
        b = 6,
        l = 4
      )
    ),
    
    ##########################################################
    # Antibody-family labels along the left
    ##########################################################
    
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 8,
      face = "bold",
      margin = margin(
        t = 4,
        r = 8,
        b = 4,
        l = 8
      )
    ),
    
    ##########################################################
    # Facet-strip formatting
    ##########################################################
    
    strip.background = element_rect(
      fill = "grey95",
      color = "grey50",
      linewidth = 0.4
    ),
    
    strip.placement = "outside",
    
    # Prevent long facet labels from being clipped
    strip.clip = "off",
    
    ##########################################################
    # Grid and layout
    ##########################################################
    
    panel.grid.major = element_line(
      linewidth = 0.2
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.spacing.x = unit(
      0.08,
      "in"
    ),
    
    panel.spacing.y = unit(
      0.08,
      "in"
    ),
    
    legend.position = "right",
    
    # More room for top and left facet labels
    plot.margin = margin(
      t = 80,
      r = 20,
      b = 20,
      l = 120,
      unit = "pt"
    )
  )

############################################################
# 10. Save outputs
############################################################

ggsave(
  "CSF_IgA_Spearman_dotplot.pdf",
  p,
  width = 36,
  height = 24,
  units = "in",
  limitsize = FALSE
)

ggsave(
  "CSF_IgA_Spearman_dotplot.png",
  p,
  width = 36,
  height = 24,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)
