## 1. Load packages ----

library(Seurat)
library(data.table)
library(dplyr)


## 2. Load expression matrix ----

# Data downloaded from the GEO:
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229169/suppl/

GSE229169_counts <- ReadMtx(
  mtx = "./GSE229169/GSE229169_mouse_filtered_matrix.mtx.gz",
  features = "./GSE229169/features.tsv.gz",
  cells = "./GSE229169/barcodes.tsv.gz",
  feature.column = 1
)


## 3. Load metadata ----

GSE229169_meta <- data.table::fread(
  "./GSE229169/mouse_filtered_metadata.tsv",
  header = TRUE,
  data.table = FALSE
)

rownames(GSE229169_meta) <- GSE229169_meta[[1]]

GSE229169_meta <- GSE229169_meta[
  ,
  -1,
  drop = FALSE
]


## 4. Match metadata to expression matrix ----

GSE229169_meta <- GSE229169_meta[
  colnames(GSE229169_counts),
  ,
  drop = FALSE
]

identical(
  rownames(GSE229169_meta),
  colnames(GSE229169_counts)
)


## 5. Create Seurat object ----

GSE229169_scedata <- CreateSeuratObject(
  counts = GSE229169_counts,
  meta.data = GSE229169_meta
)


## 6. Normalize data and identify variable features ----

GSE229169_scedata <- NormalizeData(
  GSE229169_scedata,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

GSE229169_scedata <- FindVariableFeatures(
  GSE229169_scedata,
  selection.method = "vst",
  nfeatures = 2000
)


## 7. Scale data ----

all_genes <- rownames(GSE229169_scedata)

GSE229169_scedata <- ScaleData(
  GSE229169_scedata,
  features = all_genes
)


## 8. Principal component analysis ----

GSE229169_scedata <- RunPCA(
  GSE229169_scedata,
  features = VariableFeatures(
    object = GSE229169_scedata
  ),
  npcs = 50
)


## 9. Construct neighbor graph and perform clustering ----

GSE229169_scedata <- FindNeighbors(
  GSE229169_scedata,
  dims = 1:20
)

GSE229169_scedata <- FindClusters(
  GSE229169_scedata,
  resolution = 0.7
)


## 10. UMAP dimensionality reduction ----

GSE229169_scedata <- RunUMAP(
  GSE229169_scedata,
  dims = 1:20
)


## 11. Standardize cell-type labels ----

GSE229169_scedata$subclass <- dplyr::recode(
  as.character(GSE229169_scedata$subclass),
  "ASC" = "Astrocyte",
  "MGC" = "Microglia",
  "ODC" = "Oligodendrocytes",
  "Endo" = "Endothelial cells"
)

GSE229169_Nus1_pos <- subset(
  GSE229169_scedata,
  subset = Nus1 > 0
)


## 12. Visualize ----

VlnPlot(
  GSE229169_Nus1_pos,
  features = "Nus1",
  group.by = "subclass",
  pt.size = 0.1,
  log = FALSE
) +
  xlab("") +
  theme(
    axis.title.y = element_text(
      colour = "black",
      size = 18
    ),
    axis.text = element_text(
      colour = "black",
      size = 16
    )
  )