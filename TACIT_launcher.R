devtools::install("../TACIT-bkotzen")
library(TACIT.bkotzen)
library(Seurat)
library(class)
library(segmented)
library(readr)
library(dplyr)
library(tidyr)
library(future)
library(ggplot2)

plan(multisession)  # or multicore on Linux/Mac
options(future.globals.maxSize = 2*1024^3)
suffix = "Mean"



######################
# CELLxFEATURE
# Read CELLxFEATURE and filter based on size and NaNs
CELLxFEATURE_raw=read_csv("S1_Protein_Quantification_250116.csv")
id_col <- "Object ID"
feature_cols <- colnames(CELLxFEATURE_raw)[
  startsWith(colnames(CELLxFEATURE_raw), "ROI:")
]

size_cutoff <- 30

CELLxFEATURE_raw <- CELLxFEATURE_raw %>%
  mutate(
    `Removed, size` = `Area µm^2` < size_cutoff,
    `Removed, nan`  = if_any(all_of(feature_cols), ~ is.na(.)),
    Removed = case_when(
      `Removed, nan`  ~ "Contains NaN",
      `Removed, size` ~ "Size filter",
      TRUE            ~ 'Not removed'
    )
  )

stopifnot(sum(is.na((CELLxFEATURE_raw %>% filter('Removed' == 'Size filter'))))==0)
CELLxFEATURE <- CELLxFEATURE_raw %>% filter(Removed == "Not removed")
stopifnot(sum(is.na(CELLxFEATURE[feature_cols]))==0)

cols <- grep(paste0(suffix, "$"), names(CELLxFEATURE), value = TRUE)
  
# Subset table with Object ID
sub_tbl <- CELLxFEATURE[, c(id_col, cols)]
  
# Clean column names: extract protein name from full text
# Assumes column format: "ROI: 2.00 µm per pixel: ProteinName: Aggregate"
new_colnames <- sapply(cols, function(x) {
  parts <- strsplit(x, ":")[[1]]
  # protein name is the second-to-last part
  protein <- trimws(parts[length(parts) - 1])
  return(protein)
})
  
colnames(sub_tbl) <- c(id_col, new_colnames)
  
# Store table
CELLxFEATURE_df <- sub_tbl




######################
# TYPExMARKER
# Create TYPExFEATURE matrix
# Exclude the Object ID column
features <- setdiff(colnames(CELLxFEATURE_df), id_col)

# Create identity matrix
TYPExMARKER <- diag(length(features))
colnames(TYPExMARKER) <- features
rownames(TYPExMARKER) <- paste0(features, "+")




######################
# RUN TACIT
# Set parameters
r <- 10   # (resolution) Depends on the data but aims to create microcluster cell communities with sizes averaging between 0.1–0.5% of cells per microcluster. 
# The higher the resolution, the greater the number of microclusters.
p <- 10  # (dimension) Number of dimensions used for microclusters.

# For now, run in series (switch to parallel for Orion)
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3)

TACIT = TACIT(
  CELLxFEATURE_df %>% select(-all_of(id_col)),
  TYPExMARKER,
  r=r,p=p,
  suffix=suffix)
final_threshold <- TACIT$threshold
TACIT <- TACIT$df




######################
# UMAP PLOT POSITIVES
library(RColorBrewer)
library(ggrastr)
library(forcats)

# Post-processing: plots
Signature <- as.data.frame(TYPExMARKER)

final_threshold_df <- tibble(
  Signature = colnames(Signature),
  Threshold = final_threshold$value) %>%
  arrange(desc(Threshold)) %>%
  mutate(Signature = fct_reorder(Signature, Threshold))

markers <- final_threshold_df$Signature

p <- ggplot(
  final_threshold_df,
  aes(x = Threshold, y = Signature, fill = Threshold)
) +
  geom_col() +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,      # darker = smaller, lighter = larger
    guide = "none"
  ) +
  labs(
    x = "Final Threshold",
    y = NULL,
    title = "Final Thresholds by Signature"
  ) +
  theme_minimal(base_size = 12)
# print(p)
ggsave('final_threshold.png', plot = p,
       width = 10, height = 6, units = c("in"), dpi = 300)




# Get feature values from CELLxFEATURE
X <- CELLxFEATURE_df %>%
  select(all_of(markers)) %>%
  as.matrix()
# Get threshold values from final_threshold_df and reshape to match X
T <- matrix(
  final_threshold_df$Threshold,
  nrow = nrow(X),
  ncol = ncol(X),
  byrow = TRUE
)
# Calculate positivity
CELL_POSITIVITY <- as.data.frame((X > T) * 1L)
colnames(CELL_POSITIVITY) <- markers
# Add in metadata
meta_cols <- setdiff(colnames(CELLxFEATURE), feature_cols)
stopifnot(
  dim(T)[1] == nrow(CELLxFEATURE[,meta_cols])
)
CELL_POSITIVITY <- bind_cols(
  CELLxFEATURE %>% select(all_of(meta_cols)),
  CELL_POSITIVITY
)

# Extract UMAP info and create plot_df
UMAP_df <- CELLxFEATURE %>%
  dplyr::select(all_of(id_col)) %>%
  dplyr::bind_cols(
    TACIT[c('UMAP1', 'UMAP2')]
  )


plot_df <- CELL_POSITIVITY %>%
  dplyr::select(all_of(c(id_col, as.character(markers)))) %>%
  dplyr::left_join(
    UMAP_df,
    by = id_col
  ) %>%
  tidyr::pivot_longer(
    cols = all_of(as.character(markers)),
    names_to = "Feature",
    values_to = "Positive"
  )


# Plot with downsampling
feature_colors <- setNames(
  colorRampPalette(brewer.pal(11, "Spectral"))(length(markers)),
  sort(as.character(markers))
)

plot_df <- plot_df %>%
  mutate(
    PlotColor = ifelse(
      Positive == 1,
      feature_colors[Feature],  # positive: use feature color
      "grey80"                  # negative: light gray
    )
  )

ds <- 1/50  # keep 10% of cells
p <- ggplot(
  plot_df %>%
    group_by(Feature) %>%
    slice_sample(prop = ds) %>%
    ungroup() %>%
    arrange(Positive),  # negatives first, positives last
  aes(x = UMAP1, y = UMAP2)
) +
  geom_point_rast(aes(color = PlotColor), size = 0.01, alpha = 0.75) +
  facet_wrap(~ Feature, ncol = 6) +
  scale_color_identity() +   # uses exact colors in PlotColor
  theme_void() +
  theme(
    strip.text = element_text(size = 10),
    legend.position = "none"
  )
# print(p)
ggsave('positives_all_UMAP.png', plot = p, background='white',
       width = 14, height = 10, units = c("in"), dpi = 300)

# Now some single-feature plots
ds <- 1/30
# Loop over features and create plots
single_feature_plots <- lapply(c("NeuN", "DAPI", "ApoE", "Tau4R"), function(feat) {
  ggplot(
    plot_df %>%
      filter(Feature == feat) %>%
      slice_sample(prop = ds) %>%
      arrange(Positive),   # negatives first, positives last
    aes(x = UMAP1, y = UMAP2)
  ) +
    geom_point_rast(aes(color = PlotColor), size = 0.01, alpha = 0.75) +
    scale_color_identity() +
    theme_void() +
    theme(
      strip.text = element_text(size = 12),
      legend.position = "none"
    ) +
    ggtitle(feat)
  
})
# Print plots individually
ggsave(paste0('NeuN_positivity_', suffix, '_UMAP.png'), plot = single_feature_plots[[1]],
       width = 7, height = 4, units = c("in"), dpi = 300)
ggsave(paste0('DAPI_positivity_', suffix, '_UMAP.png'), plot = single_feature_plots[[2]],
       width = 7, height = 4, units = c("in"), dpi = 300)
ggsave(paste0('ApoE_positivity_', suffix, '_UMAP.png'), plot = single_feature_plots[[3]],
       width = 7, height = 4, units = c("in"), dpi = 300)
ggsave(paste0('Tau4R_positivity_', suffix, '_UMAP.png'), plot = single_feature_plots[[4]],
       width = 7, height = 4, units = c("in"), dpi = 300)




