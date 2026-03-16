# Results

Saved outputs from differential expression analyses. All comparisons are **Aducanumab vs IgG** using pseudobulk DESeq2, stratified by brain region unless otherwise noted.

------------------------------------------------------------------------

## DEG by Cell Type and Brain Region

Named lists (by brain region) of `FindMarkers` + DESeq2 result dataframes. Each list element is a dataframe with genes as rownames (RDS) or as a `gene` column (JSON).

| File | Cell type | Source script |
|-----------------|-----------------------|--------------------------------|
| `20260223-DEG_astrocytes_byRegion.rds/.json` | Astrocytes | `deg_analyses/20260303-DEG_astrocytes_brainRegions.qmd` |
| `20260223-DEG_oligodendrocytes_byRegion.rds/.json` | Oligodendrocytes | `deg_analyses/20260303-DEG_oligo_brainRegions.qmd` |
| `20260223-DEG_OPCs_byRegion.rds` | OPCs | `deg_analyses/20260303-DEG_OPC_brainRegions.qmd` |
| `20260223-DEG_endothelial_byRegion.rds` | Endothelial | `deg_analyses/20260303-DEG_ENDOTHELIAL_brainRegions.qmd` |
| `20260303-DEG_MG_byRegion.rds/.json` | Microglia | `deg_analyses/20260303-DEG_MG_brainRegions.qmd` |
| `20260309-DEG_BAM_byRegion.rds` | Border-associated macrophages | `deg_analyses/20260309-DEG_BAM_brainRegions.qmd` |
| `20260309-DEG_TCell_byRegion.rds` | T cells | `deg_analyses/20260309-DEG_TCell_brainRegions.qmd` |
| `20260310-DEG_CP_byRegion.rds` | Choroid plexus | `deg_analyses/20260310-DEG_CP_brainRegions.qmd` |
| `20260310-DEG_Ependymal_byRegion.rds` | Ependymal | `deg_analyses/20260310-DEG_Ependymal_brainRegions.qmd` |
| `20260310-DEG_Fibroblast_byRegion.rds` | Fibroblasts | `deg_analyses/20260310-DEG_Fibroblast_brainRegions.qmd` |
| `20260310-DEG_Pericyte_byRegion.rds` | Pericytes | `deg_analyses/20260310-DEG_Pericyte_brainRegions.qmd` |
| `20260310-DEG_VLMC_byRegion.rds` | VLMCs | `deg_analyses/20260310-DEG_VLMC_brainRegions.qmd` |
| `20260310-DEG_VSMC_byRegion.rds` | VSMCs | `deg_analyses/20260310-DEG_VSMC_brainRegions.qmd` |
| `20260310-DEG_GABAergicNeuron{1-4}_byRegion.rds` | GABAergic neuron subtypes 1–4 | `deg_analyses/20260310-DEG_GABAergicNeuron{1-4}_brainRegions.qmd` |
| `20260310-DEG_GlutamatergicNeuron{1-5}_byRegion.rds` | Glutamatergic neuron subtypes 1–5 | `deg_analyses/20260310-DEG_GlutamatergicNeuron{1-5}_brainRegions.qmd` |

------------------------------------------------------------------------

## Multi-Cell-Type Summary

| File | Contents |
|---------------------------|---------------------------------------------|
| `20260304-DEG_allCellTypes.rds` | Named list by cell type, each element is a dataframe of all DESeq2 results across regions |
| `20260304-DEG_allCellTypes.json` | Same, JSON format with `gene` column |

Source: `deg_analyses/20260304-DEG_allCellTypes_brainRegions.qmd`

------------------------------------------------------------------------

## Astrocyte Subcluster Analyses

Results from reclustering the astrocyte subset (SCT, resolution 0.3).

| File | Contents | Source script |
|-----------------|----------------------|---------------------------------|
| `20260305-astrocyte_all_markers_0.3.csv` | `FindAllMarkers` output for all astrocyte clusters (initial reclustering) | `20260305-astrocytes_reclustered.qmd` |
| `20260305-DEG_astrocytes_clusters0-1.rds` | Adu vs IgG DESeq2 results for astrocyte clusters 0 and 1 | `20260305-astrocytes_reclustered.qmd` |
| `20260305-pairwise_DE_astro_clusters.rds` | Pairwise DESeq2 between astrocyte clusters 0, 1, and 7 | `20260305-astrocytes_cluster_characterization.qmd` |
| `20260309-astrocyte_cleaned_all_markers_0.3.csv` | `FindAllMarkers` after removing marker genes from panel (cluster 6 retained) | `20260309-astrocytes_reclustered_cleaned.qmd` |
| `20260312-astrocyte_cleaned_all_markers_0.3.csv` | `FindAllMarkers` after removing marker genes and cluster 6 (canonical cleaned object) | `20260309-astrocytes_reclustered_cleaned.qmd` |

------------------------------------------------------------------------

## Reference

| File | Contents |
|---------------------------|---------------------------------------------|
| `Xenium_Final_Gene_List_Updated_20260220_corrected.json` | Curated 480-gene Xenium panel gene list |

------------------------------------------------------------------------

## Utility Scripts

| File | Purpose |
|-----------------------------|-------------------------------------------|
| `20260316-export_DEG_byRegion_toJSON.R` | Converts astrocyte, MG, and oligodendrocyte byRegion RDS files to JSON |