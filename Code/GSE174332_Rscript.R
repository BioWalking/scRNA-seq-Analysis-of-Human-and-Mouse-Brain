## Data downloaded from Brain Cell Atlas (https://www.braincellatlas.org/dataSet)

library(reticulate)
library(anndata)
library(Seurat)
library(Matrix)

use_python("D:/Python/python.exe", required = TRUE)


## 1. Load AnnData object ----

adata <- read_h5ad(
  "./human_brain_motor_cortex_Pineda_2021_10x/processedData/annot.h5ad"
)


## 2. Extract expression matrix ----

# AnnData stores the expression matrix as [Cell x Gene],
# whereas Seurat requires [Gene x Cell].

counts_matrix <- t(adata$X)

colnames(counts_matrix) <- adata$obs_names
rownames(counts_matrix) <- adata$var_names


## 3. Extract metadata ----

metadata <- as.data.frame(
  adata$obs
)


## 4. Create Seurat object ----

pbmc <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = metadata
)


## 5. Retain selected samples ----

pbmc <- subset(
  pbmc,
  subset = Batch %in% c(
    "301MCX", "302MCX", "303MCX", "304MCX",
    "306MCX", "307MCX", "308MCX", "309MCX",
    "311MCX", "317MCX", "318MCX", "319MCX",
    "322MCX", "323MCX", "324MCX", "325MCX",
    "328MCX"
  )
)


## 6. Normalize data and identify variable features ----

pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)


## 7. Scale data ----

all_genes <- rownames(pbmc)

pbmc <- ScaleData(
  pbmc,
  features = all_genes
)


## 8. Principal component analysis ----

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(
    object = pbmc
  ),
  npcs = 50
)


## 9. Construct neighbor graph and perform clustering ----

pbmc <- FindNeighbors(
  pbmc,
  dims = 1:25
)

pbmc <- FindClusters(
  pbmc,
  resolution = 0.5
)


## 10. UMAP dimensionality reduction ----

pbmc <- RunUMAP(
  pbmc,
  dims = 1:25
)


GSE174332_Nus1_pos <- subset(
  pbmc,
  subset = NUS1 > 0
)


## 11. Visualize ----

VlnPlot(
  GSE174332_Nus1_pos,
  features = "NUS1",
  group.by = "CellType",
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