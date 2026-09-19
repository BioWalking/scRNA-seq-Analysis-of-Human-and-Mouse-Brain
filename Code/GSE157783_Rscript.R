## 1. Load packages ----

library(Seurat)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)


## 2. Load raw data ----

# Data downloaded from the GEO:
# https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157783/suppl/

meta.data <- data.table::fread(
  "GSE157783/IPDCO_hg_midbrain_cell.tsv",
  header = TRUE,
  data.table = FALSE
)

GSE157783_data <- Read10X(
  data.dir = "./GSE157783/"
)


## 3. Convert Ensembl IDs to gene symbols ----

gene_id <- bitr(
  rownames(GSE157783_data),
  fromType = "ENSEMBL",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

gene_symbol <- gene_id$SYMBOL[
  match(rownames(GSE157783_data), gene_id$ENSEMBL)
]

gene_symbol[is.na(gene_symbol)] <- rownames(GSE157783_data)[
  is.na(gene_symbol)
]

rownames(GSE157783_data) <- make.unique(gene_symbol)


## 4. Create Seurat object ----

GSE157783_scedata <- CreateSeuratObject(
  counts = GSE157783_data
)


## 5. Add metadata ----

rownames(meta.data) <- meta.data[[1]]

meta.data <- meta.data[
  ,
  -1,
  drop = FALSE
]

meta.data <- meta.data[
  colnames(GSE157783_scedata),
  ,
  drop = FALSE
]

identical(
  rownames(meta.data),
  colnames(GSE157783_scedata)
)

GSE157783_scedata <- AddMetaData(
  GSE157783_scedata,
  metadata = meta.data
)

GSE157783_scedata$Celltype <- ifelse(
  GSE157783_scedata$cell_ontology %in% c(
    "CADPS2+ neurons",
    "DaNs",
    "Excitatory",
    "GABA",
    "Inhibitory"
  ),
  "Neuron",
  GSE157783_scedata$cell_ontology
)


## 6. Retain control samples ----

GSE157783_scedata <- GSE157783_scedata[
  ,
  GSE157783_scedata$patient %in% c(
    "C1", "C2", "C3",
    "C4", "C5", "C6"
  )
]


## 7. Normalize data and identify variable features ----

GSE157783_scedata <- NormalizeData(
  GSE157783_scedata,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

GSE157783_scedata <- FindVariableFeatures(
  GSE157783_scedata,
  selection.method = "vst",
  nfeatures = 2000
)


## 8. Scale data ----

all.genes <- rownames(GSE157783_scedata)

GSE157783_scedata <- ScaleData(
  GSE157783_scedata,
  features = all.genes
)


## 9. Principal component analysis ----

GSE157783_scedata <- RunPCA(
  GSE157783_scedata,
  features = VariableFeatures(
    object = GSE157783_scedata
  ),
  npcs = 50
)


## 10. Construct neighbor graph and perform clustering ----

GSE157783_scedata <- FindNeighbors(
  GSE157783_scedata,
  dims = 1:10
)

GSE157783_scedata <- FindClusters(
  GSE157783_scedata,
  resolution = 0.5
)


## 11. UMAP dimensionality reduction ----

GSE157783_scedata <- RunUMAP(
  GSE157783_scedata,
  dims = 1:10
)

GSE157783_NUS1_pos <- subset(
  GSE157783_scedata,
  subset = NUS1 > 0
)


## 12. Visualize ----

VlnPlot(
  GSE157783_NUS1_pos,
  features = "NUS1",
  group.by = "Celltype",
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