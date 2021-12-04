library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(corrplot)

blastoid.data <- Read10X(data.dir = "F:/scRNA_Seq_Data_Analysis/Blastoids/data/")
blastoid <- CreateSeuratObject(counts = blastoid.data, project = "iblastoid", min.cells = 3, min.features = 200)
blastoid
blastoid[["percent.mt"]] <- PercentageFeatureSet(blastoid, pattern = "^MT-")
head(blastoid@meta.data, 5)
VlnPlot(blastoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(blastoid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(blastoid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

blastoid <- subset(blastoid, subset = nFeature_RNA > 1300 & percent.mt < 15)
blastoid
blastoid <- NormalizeData(blastoid, normalization.method = "LogNormalize", scale.factor = 10000)
blastoid
blastoid <- FindVariableFeatures(blastoid, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(blastoid), 10)
top10
###Run PCA for dimensionality reduction ##
blastoid <- ScaleData(blastoid)
blastoid <- RunPCA(blastoid, features = VariableFeatures(object = blastoid))
print(blastoid[["pca"]], dims = 1:5, nfeatures = 5)

### Visualize PCA ###
VizDimLoadings(blastoid, dims = 1:2, reduction = "pca")
### Visualize PCA ###
VizDimLoadings(blastoid, dims = 1:2, reduction = "pca")
DimPlot(blastoid, reduction = "pca")
DimHeatmap(blastoid, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(blastoid, dims = 1:5, cells = 500, balanced = TRUE)
blastoid <- JackStraw(blastoid, num.replicate = 100)
blastoid <- ScoreJackStraw(blastoid, dims = 1:20)
JackStrawPlot(blastoid, dims = 1:15)
JackStrawPlot(blastoid, dims = 1:10)
JackStrawPlot(blastoid, dims = 1:15)
ElbowPlot(blastoid)
### Find clusters ###
blastoid <- FindNeighbors(blastoid, dims = 1:10)
blastoid <- FindClusters(blastoid, resolution = 0.5)
head(Idents(blastoid), 5)
### UMAP ###
blastoid <- RunUMAP(blastoid, dims = 1:10)
DimPlot(blastoid, reduction = "umap")
# Identify gene markers###
all_markers <-FindAllMarkers(blastoid,
                             min.pct =  0.25,
                             min.diff.pct = 0.25)
## Loading the annotation file
annotations<- read.csv("annotation.csv")
head(annotations)
## Merge gene annotations to marker results
all_markers<- left_join(all_markers,
                        annotations[, c(1:2, 4,5)],
                        by =c("gene" = "gene_name"))
all_markers<-all_markers[ , c(6:8, 1:5, 9:10)]
view(all_markers)
## Create a data frame of clusters 0-11 sorted by average log FC
gen_marker_table<- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n=100)
view(gen_marker_table)
is.data.frame(gen_marker_table)
## Marker Visualization
VlnPlot(blastoid, features = c("POU5F1", "TMEM54"))
FeaturePlot(blastoid, features = c("POU5F1", "NANOG"))
new.cluster.ids<- c("Primitive Endoderm 1",
                    "Trophectoderm 1",
                    "Epiblast 1",
                    "Trophectoderm 2",
                    " Trophoectoderm 3",
                    "Primitive Endoderm 2",
                    "Unknown Lineage 1",
                    "UnKnown Lineage 2",
                    "Epiblast 2",
                    "Unknown Lineage 3",
                    "Unknown Lineage 4",
                    "Intermediate State")
names(new.cluster.ids) <- levels(blastoid)
blastoid <- RenameIdents(blastoid, new.cluster.ids)
DimPlot(blastoid, reduction = "umap", label = TRUE, pt.size = 1.5)

## Feature plot for Canonical markers of Epiblast, Primitive Endoderm & Trophectoderm
FeaturePlot(blastoid, features = c("SOX17", "GATA6")) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlGn")))
FeaturePlot(blastoid, features = c("POU5F1", "NANOG")) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlGn")))
FeaturePlot(blastoid, features = c("CDX2", "GATA2")) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlGn")))

table(Idents(blastoid))
prop.table(table(Idents(blastoid)))
WhichCells(blastoid, idents = "Epiblast 1")

###Calculating the average gene expression within a cluster
cluster.averages <- AverageExpression(blastoid)

head(cluster.averages[["RNA"]][, 1:5])

### correlation plot within clusters

data<- read.csv("ClusterAverages.csv")
View(data)

data<-data[,-1]

View(data)

log.data<-log2(data +1)
View(log.data)

cor.data<-cor(log.data,method = c("spearman"))

cor.data
corrplot(cor.data,tl.col = "black")
dev.off()
