# sc-rna
practice script


library(Matrix)

# 1. Load features (genes)
features <- read.delim(
  "C:/Users/christinel/Desktop/scRNA-seq/scRNAseq_analysis_vignette-master/data/DS1/features.tsv.gz",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 2. Load barcodes (cell IDs)
barcodes <- read.delim(
  "C:/Users/christinel/Desktop/scRNA-seq/scRNAseq_analysis_vignette-master/data/DS1/barcodes.tsv.gz",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)[,1]

# 3. Load matrix
counts <- readMM(
  "C:/Users/christinel/Desktop/scRNA-seq/scRNAseq_analysis_vignette-master/data/DS1/matrix.mtx.gz"
)

# 4. Assign row/column names
rownames(counts) <- make.unique(features[,2])   # gene names
colnames(counts) <- barcodes                    # cell barcodes


Use Seurat's Read10X

Instead of manually reading 3 files, do:

library(Seurat)

counts <- Read10X("C:/Users/christinel/Desktop/scRNA-seq/scRNAseq_analysis_vignette-master/data/DS1/")

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)

| QC metric        | Interpretation                                              | Recommended threshold                    |
| ---------------- | ----------------------------------------------------------- | ---------------------------------------- |
| **percent.mt**   | Some high-mt cells at low counts → likely stressed          | **< 5%**                                 |
| **nCount_RNA**   | Broad normal distribution                                   | Automatically filtered via gene cutoff   |
| **nFeature_RNA** | Strong correlation with nCount_RNA (0.93) → dataset is good | **500–5000 genes**                       |
| **Doublets**     | Likely around 5,500–6,000+ genes                            | Remove via upper cutoff or DoubletFinder |


seurat <- NormalizeData(seurat)

top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
plot1 + plot2

seurat <- ScaleData(seurat)

seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))
seurat <- SCTransform(seurat, variable.features.n = 3000)
seurat <- SCTransform(seurat,
                      vars.to.regress = c("nFeature_RNA", "percent.mt"),
                      variable.features.n = 3000)

seurat <- RunPCA(seurat, npcs = 50)

ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))
PCHeatmap(seurat, dims = 1:20, cells = 500, balanced = TRUE, ncol = 4)

seurat <- RunTSNE(seurat, dims = 1:20)
seurat <- RunUMAP(seurat, dims = 1:20)

plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2


plot1 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "tsne")
plot2 <- FeaturePlot(seurat, c("MKI67","NES","DCX","FOXG1","DLX2","EMX1","OTX2","LHX9","TFAP2A"),
                     ncol=3, reduction = "umap")
plot1 / plot2

cluster the cells
seurat <- FindNeighbors(seurat, dims = 1:20)

seurat <- FindClusters(seurat, resolution = 1)

The easiest to visualize expression of marker genes of interest across cell clusters is probably by a heatmap.
ct_markers <- c("MKI67","NES","DCX","FOXG1", # G2M, NPC, neuron, telencephalon
                "DLX2","DLX5","ISL1","SIX3","NKX2.1","SOX6","NR2F2", # ventral telencephalon related
                "EMX1","PAX6","GLI3","EOMES","NEUROD6", # dorsal telencephalon related
                "RSPO3","OTX2","LHX9","TFAP2A","RELN","HOXB2","HOXB5") # non-telencephalon related
DoHeatmap(seurat, features = ct_markers) + NoLegend()

cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
library(dplyr)
cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



top2 <- cl_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)


# Find markers
cl_markers <- FindAllMarkers(
  seurat, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = log(1.2)
)

library(dplyr)

# Get top 2 markers per cluster (using avg_log2FC)
top2 <- cl_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 2)

library(presto)
cl_markers_presto <- wilcoxauc(seurat)
cl_markers_presto %>%
    filter(logFC > log(1.2) & pct_in > 20 & padj < 0.05) %>%
    group_by(group) %>%
    arrange(desc(logFC), .by_group=T) %>%
    top_n(n = 2, wt = logFC) %>%
    print(n = 40, width = Inf)

    library(presto)

cl_markers_presto <- wilcoxauc(seurat)

top2_presto <- cl_markers_presto %>%
  dplyr::filter(logFC > log(1.2),
                pct_in > 20,
                padj < 0.05) %>%
  dplyr::group_by(group) %>%
  dplyr::arrange(dplyr::desc(logFC), .by_group = TRUE) %>%
  dplyr::slice_head(n = 2)


library(dplyr)

top10_presto <- cl_markers_presto %>%
  filter(logFC > log(1.2),        # optional: keep only up-regulated genes
         pct_in > 20,
         padj <  0.05) %>%
  group_by(group) %>%
  slice_max(order_by = logFC, n = 10)

DoHeatmap(seurat, features = top10_presto$feature) +
+     NoLegend()
Error in DoHeatmap(seurat, features = top10_presto$feature) : 
  could not find function "DoHeatmap"

  load seurat

  
DoHeatmap(seurat, features = top10_presto$feature) + NoLegend()
plot1 <- FeaturePlot(seurat, c("NEUROD2","NEUROD6"), ncol = 1)
plot2 <- VlnPlot(seurat, features = c("NEUROD2","NEUROD6"), pt.size = 0)
plot1 + plot2 + plot_layout(widths = c(1, 2))

library(patchwork)

plot1 <- FeaturePlot(seurat, c("NEUROD2","NEUROD6"), ncol = 1)
plot2 <- VlnPlot(seurat, features = c("NEUROD2","NEUROD6"), pt.size = 0)

(p <- plot1 + plot2 + plot_layout(widths = c(1, 2)))
