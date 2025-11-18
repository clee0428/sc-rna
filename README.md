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


