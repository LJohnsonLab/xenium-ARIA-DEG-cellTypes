# =============================================================================
# Characterisation of geometric outliers in astrocyte cluster 0
# =============================================================================
# Cluster 0 of the reclustered/cleaned astrocyte object contains a population
# of cells that sit spatially apart from the main blob in UMAP space (a
# trailing "tail"). These cells likely represent transcriptional bleed-in from
# adjacent non-astrocyte cell types rather than a true astrocyte substate.
#
# Approach:
#   1. Compute Mahalanobis distance from the centroid of cluster 0 in UMAP
#      space — this respects the shape and orientation of the blob.
#   2. Label cells beyond the 97.5th-percentile distance as "bleed" and the
#      rest as "core".
#   3. Run a Wilcoxon differential expression test (bleed vs core) to identify
#      genes that are enriched in the outlier tail, which can inform whether
#      to exclude these cells from downstream analysis.
#
# Input:  seurat_objects/20260318-astro_cleaned2_0.3.qs2
# Output: bleed_vs_core_markers_cluster0.csv
# =============================================================================

library(Seurat)
library(tidyverse)

obj <- qs2::qs_read("seurat_objects/20260318-astro_cleaned2_0.3.qs2")

# UMAP coordinates provide the geometric definition of the main blob and its outer tail.
umap_tbl <-
  Embeddings(obj, "umap") |>
  as.data.frame() |>
  rownames_to_column("cell") |>
  as_tibble()

# the cluster labels identify which cells belong to cluster 0 before applying
# the boundary-based classification.
cluster_tbl <-
  tibble(
    cell = colnames(obj),
    cluster = as.character(Idents(obj)))

c0_umap_tbl <-
  umap_tbl |>
  inner_join(cluster_tbl, by = "cell") |>
  filter(cluster == "0")

# Why: Mahalanobis distance is more appropriate than raw Euclidean distance because it
# respects the shape and orientation of the main cluster-0 blob.
c0_coords <-
  c0_umap_tbl |>
  select(umap_1, umap_2) |>
  as.matrix()

center <- colMeans(c0_coords)
cov_mat <- cov(c0_coords)

c0_umap_tbl <-
  c0_umap_tbl |>
  mutate(
    maha_dist = mahalanobis(c0_coords, center = center, cov = cov_mat)
  )

# Why: the cutoff defines the outer boundary of the main blob. Cells beyond that
# boundary are treated as geometric outliers, which here are the bleeding cells.
cutoff <- quantile(c0_umap_tbl$maha_dist, 0.975, na.rm = TRUE)

c0_umap_tbl <-
  c0_umap_tbl |>
  mutate(
    bleed_group = if_else(maha_dist > cutoff, "bleed", "core")
  )

# Why: only the categorical label is needed downstream, so we store just that metadata.
obj$geom_bleed <- NA_character_
obj$geom_bleed[match(c0_umap_tbl$cell, colnames(obj))] <- c0_umap_tbl$bleed_group

# Why: this plot is the visual sanity check that the chosen cutoff is capturing the
# boundary cells outside the main cluster-0 blob.
bleed_cells <-
  c0_umap_tbl |>
  filter(bleed_group == "bleed") |>
  pull(cell)

DimPlot(
  obj,
  group.by = "seurat_clusters",
  cells.highlight = bleed_cells,
  cols = "lightgray",
  cols.highlight = "red"
)

# Why: gene characterization should be done only within cluster 0, otherwise the result
# would mostly reflect differences between major clusters rather than the bleed/core contrast.
sub <-
  subset(obj, cells = c0_umap_tbl$cell)

Idents(sub) <- sub$geom_bleed

table(Idents(sub))

# Why: this directly tests which genes are enriched in bleeding cells compared with the
# core of cluster 0. The default Wilcoxon test is usually a reasonable first pass.
markers_bleed_vs_core <-
  FindMarkers(
    sub,
    ident.1 = "bleed",
    ident.2 = "core",
    assay = DefaultAssay(sub),
    logfc.threshold = 0.1,
    min.pct = 0.05,
    test.use = "wilcox"
  ) |>
  rownames_to_column("gene") |>
  as_tibble() |>
  arrange(desc(avg_log2FC))

markers_core_vs_bleed <-
  markers_bleed_vs_core |>
  arrange(avg_log2FC)

markers_bleed_vs_core |> slice_head(n = 30)
markers_core_vs_bleed |> slice_head(n = 30)

# Optional export
write_csv(markers_bleed_vs_core, "bleed_vs_core_markers_cluster0.csv")