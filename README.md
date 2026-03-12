## ARIA: Xenium Spatial Transcriptomics Analysis (DEGs & Cell-type proportions analysis)

https://ljohnsonlab.github.io/xenium-ARIA-DEG-cellTypes/

### Overview

This project contains a reproducible Quarto website documenting the analysis of Xenium spatial transcriptomics data from the ARIA study, investigating the differential brain effects of Aducanumab (anti-Aβ antibody) treatment in a mouse model of Amyloid-Related Imaging Abnormalities (ARIA) expressing human APOE4.

The 480-gene targeted Xenium panel was applied to 6 mouse brains, enabling cell-type identification and spatially resolved gene expression across multiple brain regions. Analyses include:

- **Cell-type characterisation** — UMAP clustering, spatial mapping, and cell-type proportion comparisons across treatment groups
- **Pseudobulk differential expression** — DESeq2-based Adu vs IgG comparisons for 20+ cell types across brain regions, including major glia (microglia, astrocytes, oligodendrocytes, OPC), vascular cells (endothelial, pericyte, VSMC, VLMC), immune subsets (BAM, T-cells), neurons (GABAergic, glutamatergic), and other cell types (fibroblast, ependymal, choroid plexus)
- **Microglia subpopulation analysis** — Reclustering and characterization of microglia subclusters
- **Astrocyte subpopulation analysis** — Reclustering and characterization of astrocyte subclusters


### Study Design

-   Platform: 10x Genomics Xenium In Situ spatial transcriptomics
-   Treatment Groups:
    -   Aducanumab (Adu): Monoclonal antibody treatment
    -   IgG control: Control antibody
-   Samples: 6 biological replicates across 2 slides
    -   Slide 1 (0021991): KK4_465, KK4_492, KK4_504
    -   Slide 2 (0021998): KK4_496, KK4_502, KK4_464

### Raw data

**Seurat objects** in this [Dropbox](https://www.dropbox.com/scl/fo/vdnbvohbtzyjgbxeai0jv/AFgDiHe_zcHO3BOD07V7xsM?rlkey=dtnvry6f2pbrtzod3w5kafsnb&dl=0). Download the data into the seurat_objects folder.