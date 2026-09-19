## 1. Load packages ----

library(DoubletFinder)
library(Seurat)
library(data.table)
library(dplyr)


## 2. Load individual samples ----

counts_GSM9205080 <- Read10X(
  data.dir = "GSE306642_RAW/GSM9205080_FW1023M_matrix"
)

seu_GSM9205080 <- CreateSeuratObject(
  counts = counts_GSM9205080,
  project = "GSM9205080",
  min.cells = 5,
  min.features = 200
)

counts_GSM9205081 <- Read10X(
  data.dir = "GSE306642_RAW/GSM9205081_MW210116M_matrix"
)

seu_GSM9205081 <- CreateSeuratObject(
  counts = counts_GSM9205081,
  project = "GSM9205081",
  min.cells = 5,
  min.features = 200
)


## 3. Merge samples ----

GSE306642_scedata <- merge(
  x = seu_GSM9205080,
  y = seu_GSM9205081
)


## 4. Normalize data and identify variable features ----

GSE306642_scedata <- NormalizeData(
  GSE306642_scedata,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

GSE306642_scedata <- FindVariableFeatures(
  GSE306642_scedata,
  selection.method = "vst",
  nfeatures = 2000
)


## 5. Scale data ----

all_genes <- rownames(GSE306642_scedata)

GSE306642_scedata <- ScaleData(
  GSE306642_scedata,
  features = all_genes
)


## 6. Principal component analysis ----

GSE306642_scedata <- RunPCA(
  GSE306642_scedata,
  features = VariableFeatures(
    object = GSE306642_scedata
  )
)


## 7. Construct neighbor graph and perform clustering ----

GSE306642_scedata <- FindNeighbors(
  GSE306642_scedata,
  dims = 1:20
)

GSE306642_scedata <- FindClusters(
  GSE306642_scedata,
  resolution = 0.5
)


## 8. UMAP dimensionality reduction ----

GSE306642_scedata <- RunUMAP(
  GSE306642_scedata,
  dims = 1:20
)


## 9. Doublet detection ----

GSE306642_scedata <- JoinLayers(
  GSE306642_scedata
)

sweep_res <- paramSweep(
  GSE306642_scedata,
  PCs = 1:20,
  sct = FALSE
)

sweep_stats <- summarizeSweep(
  sweep_res,
  GT = FALSE
)

bcmvn <- find.pK(
  sweep_stats
)

optimal_pK <- as.numeric(
  as.vector(
    bcmvn$pK[
      which.max(bcmvn$BCmetric)
    ]
  )
)

optimal_pK

cluster_annotations <- GSE306642_scedata$seurat_clusters

homotypic_prop <- modelHomotypic(
  cluster_annotations
)

nExp_poi <- round(
  0.075 * nrow(GSE306642_scedata@meta.data)
)

nExp_poi_adj <- round(
  nExp_poi * (1 - homotypic_prop)
)

GSE306642_scedata <- doubletFinder(
  GSE306642_scedata,
  PCs = 1:20,
  pN = 0.25,
  pK = optimal_pK,
  nExp = nExp_poi,
  reuse.pANN = FALSE,
  sct = FALSE
)

GSE306642_scedata <- doubletFinder(
  GSE306642_scedata,
  PCs = 1:20,
  pN = 0.25,
  pK = optimal_pK,
  nExp = nExp_poi_adj,
  reuse.pANN = paste0(
    "pANN_0.25_",
    optimal_pK,
    "_",
    nExp_poi
  ),
  sct = FALSE
)

table(
  GSE306642_scedata$DF.classifications_0.25_0.005_978
)

GSE306642_scedata <- GSE306642_scedata[
  ,
  GSE306642_scedata$DF.classifications_0.25_0.005_978 == "Singlet"
]


## 10. Recreate Seurat object after doublet removal ----

GSE306642_scedata <- CreateSeuratObject(
  counts = GSE306642_scedata@assays$RNA$counts,
  meta.data = GSE306642_scedata@meta.data
)


## 11. Reprocess singlet cells ----

GSE306642_scedata <- NormalizeData(
  GSE306642_scedata,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

GSE306642_scedata <- FindVariableFeatures(
  GSE306642_scedata,
  selection.method = "vst",
  nfeatures = 2000
)

all_genes <- rownames(GSE306642_scedata)

GSE306642_scedata <- ScaleData(
  GSE306642_scedata,
  features = all_genes
)

GSE306642_scedata <- RunPCA(
  GSE306642_scedata,
  features = VariableFeatures(
    object = GSE306642_scedata
  )
)

GSE306642_scedata <- FindNeighbors(
  GSE306642_scedata,
  dims = 1:20
)

GSE306642_scedata <- FindClusters(
  GSE306642_scedata,
  resolution = 0.5
)

GSE306642_scedata <- RunUMAP(
  GSE306642_scedata,
  dims = 1:20
)


## 12. Remove low-quality cell clusters ----

GSE306642_scedata <- subset(
  GSE306642_scedata,
  idents = c(
    "12", "17", "13",
    "20", "21", "22"
  ),
  invert = TRUE
)


## 13. Annotate cell types ----

cluster_to_celltype <- c(
  "0" = "Microglia",
  "1" = "Pericyte",
  "2" = "Oligodendrocyte",
  "3" = "Endothelial cells",
  "4" = "Astrocyte",
  "5" = "Endothelial cells",
  "6" = "Endothelial cells",
  "7" = "Ependymal cells",
  "8" = "Oligodendrocyte",
  "9" = "Macrophage",
  "10" = "Oligodendrocyte",
  "11" = "Oligodendrocyte",
  "14" = "Macrophage",
  "15" = "Pericyte",
  "16" = "OPC",
  "18" = "Oligodendrocyte",
  "19" = "Endothelial cells"
)

GSE306642_scedata <- RenameIdents(
  GSE306642_scedata,
  cluster_to_celltype
)

GSE306642_scedata$celltype <- GSE306642_scedata@active.ident

GSE306642_Nus1_pos <- subset(
  GSE306642_scedata,
  subset = Nus1 > 0
)


## 14. Visualize ----

VlnPlot(
  GSE306642_Nus1_pos,
  features = "Nus1",
  group.by = "celltype",
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