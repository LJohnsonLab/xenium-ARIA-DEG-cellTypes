# seurat_objects/ (DEG & cell-type analysis)

Seurat objects consumed by this repository. Ddownload them from the Dropbox link in the top-level
`README.md` into this folder. 

## Runnable set (loaded by current code)

| File | Size | Produced by | Consumed by |
|------|------|-------------|-------------|
| `20260217_Part2_fullobj.qs2` | ~9.0 GB | Upstream full object, all cell types, post-RCTD/SPLIT (Part2 pipeline, corrected 20260217 version) | all `deg_analyses/*.qmd`, `20260302-cell_types.qmd`, `20260305-astrocytes_reclustered.qmd` |
| `20260225-mg_kai.qs2` | ~5.6 GB | Microglia subset with Kai's annotation (canonical) | `20260303-cellType_proportions_microglia.qmd`, `deg_analyses/20260303-DEG_MG_brainRegions.qmd`, `deg_analyses/20260304-DEG_allCellTypes_brainRegions.qmd` |
| `20260305-astro_0.3.qs2` | ~5.6 GB | `20260305-astrocytes_reclustered.qmd` (astrocyte subset of the full object, reclustered at resolution 0.3; SCT + PCA + UMAP) | `20260305-astrocytes_cluster_characterization.qmd`, `20260305-astrocytes_reclustered.qmd`, `20260309-astrocytes_reclustered_cleaned.qmd` |
| `20260318-astro_cleaned2_0.3.qs2` | ~5.5 GB | `20260309-astrocytes_reclustered_cleaned.qmd` (low-quality clusters/genes removed, reclustered at 0.3) | `20260309-astrocytes_reclustered_cleaned.qmd`, `20260310-cellType_proportions_astrocytes.qmd`, `deg_analyses/20260303-DEG_astrocytes_brainRegions.qmd` |

## Upstream provenance (NOT required to run; kept for lineage)

These precede or split the runnable objects above and are loaded by no current
notebook. They are retained for provenance only.

| File | Size | Notes |
|------|------|-------|
| `20260126_mergedobj_QC.rds` | ~8.7 GB | Merged, QC-filtered object (both slides). Referenced only by the collaborator provenance script, via an absolute path to that machine. |
| `20260204_RCTDobj_bothslides.rds` | ~0.4 GB | RCTD cell-type deconvolution output. |
| `20260209_aria_RCTDSPLIT.qs2` | ~9.1 GB | Post-RCTD + SPLIT-decontaminated object; precursor to `20260217_Part2_fullobj.qs2`. |
| `20260217_Part2_fullobj_slide1.qs2` | ~4.5 GB | Per-slide split of the full object (slide 1). |
| `20260217_Part2_fullobj_slide2.qs2` | ~4.5 GB | Per-slide split of the full object (slide 2). |
| `20260225-mg_kai_slide1.qs2` | ~2.8 GB | Per-slide microglia (slide 1). |
| `20260225-mg_kai_slide2.qs2` | ~2.8 GB | Per-slide microglia (slide 2). |
