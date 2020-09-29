# Seurat v3.0.0 pipeline

library(Seurat)
library(dplyr)
library(magrittr)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(xlsx)

# import UMAP
# IMPORTANT: Have to set path before loading reticulate package
Sys.setenv(PATH= paste("C:/Users/Stephanie/Anaconda3/envs/r-reticulate/Library/bin",Sys.getenv()["PATH"],sep=";"))
Sys.setenv(RETICULATE_PYTHON = "C:/Users/Stephanie/Anaconda3/envs/r-reticulate/python.exe")
library(reticulate)
use_condaenv("r-reticulate")
py_config()
import("umap")
#py_install("umap-learn") # Only run if haven't installed umap python package

# PREPROCESSING
# using filtered for example
KOFIL.data <- Read10X("KO/KO/raw_gene_bc_matrices/mm10/")
WTFIL.data <- Read10X("WT/WT/raw_gene_bc_matrices/mm10/")

# Create KO and WT Seurat objects
KOFIL <- CreateSeuratObject(counts = KOFIL.data, project = 'KOFIL', min.cells = 3, min.features = 100)
WTFIL <- CreateSeuratObject(counts = WTFIL.data, project = 'WTFIL', min.cells = 3, min.features = 100)

# Add any meta.data needed
KOFIL@meta.data$group <- "KO"
WTFIL@meta.data$group <- "WT"

# If more than 2 objects, then need to merge
# Add IDs to cell names to distinguish groups
merge <- merge(x = KOFIL, y = WTFIL, add.cell.ids = c('KO', 'WT'))

# IMPORTANT
# Changing between meta.data for identities
Idents(object = merge) <- 'group'
# Check active identity
levels(merge)

# QC
# get percent of mitochondrial genes reads
merge[["percent.mt"]] <- PercentageFeatureSet(object = merge, pattern = "^mt-")

# Look at nFeature, nCount, percent.mt
VlnPlot(object = merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        pt.size = .5)

# Another way to look at nFeature, nCount, percent.mt
plot1 <- FeatureScatter(object = merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Filter cells based on nFeature
# filtering cells that have nCount < 50000 and percent.mt < 20%
merge <- subset(x = merge, subset = nCount_RNA < 50000 & percent.mt < 15)

# Normalize merged data
merge <- NormalizeData(object = merge, normalization.method = "LogNormalize", 
                             scale.factor = 10000)
# Can use code below if want default parameters (same thing as above)
# merge <- NormalizeData(object = merge)

# find subset of features that exhibit high cell-to-cell variation in the dataset
# vst --> loess.span, clip.max, nfeatures
# change top variable genes you want to keep w/ nfeatures
merge <- FindVariableFeatures(object = merge, selection.method = 'vst', loess.span = .3,
                              nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = merge), 10)

plot1 <- VariableFeaturePlot(object = merge)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))
# ignore warning for now

# mean.var.plot --> mean.function, dispersion.function, num.bin, mean/dispersion.cutoff
# merge <- FindVariableFeatures(object = merge, selection.method = 'mean.var.plot', ,
#                               mean.cutoff = c(0.0125,Inf), dispersion.cutoff = c(0.5,Inf))

# Scale the data
# if leave features paramter blank, then it will use variable genes
# if don't use all genes, then DoHeatmap function might not work
all.genes <- rownames(merge)
merge <- ScaleData(object = merge, features = all.genes, 
                   vars.to.regress = c("nCount", "percent.mt"))

# Dimensional Reduction: PCA
# Can control the number of pcs that you want to calculate w/ npcs
merge <- RunPCA(object = merge, features = VariableFeatures(object = merge), 
                nfeatures.print = 5)

# Choosing how many PCs for clustering
# visualizing pcs
ElbowPlot(object = merge, ndims = 50)

# JackStraw if you want
# merge <- JackStraw(object = merge, num.replicate = 100)
# merge <- ScoreJackStraw(object = merge, dims = 1:20)
# JackStrawPlot(object = merge, dims = 1:15)

# Calculating --> until get over 90% variance
stdev <- merge[['pca']]@stdev
var <- (merge[['pca']]@stdev)^2
sum(var[1:31])/sum(var) 

# CLUSTERING
merge <- FindNeighbors(object = merge, dims = 1:31, k.param = 20)
merge <- FindClusters(object = merge, resolution = 0.5)
head(x = Idents(object = merge), 5)

# VISUALIZATION: TSNE/UMAP
# TSNE example
merge <- RunTSNE(merge, dims = 1:31)
DimPlot(object = merge, reduction = "tsne")

# UMAP example
merge <- RunUMAP(merge, dims = 1:31)
DimPlot(object = merge, reduction = 'umap')

# LABELING CLUSTERS
# Visualizing expression 
FeaturePlot(merge, features = c("Krt18", "Cdh11", "Ptprc", "Hba-a2", "Cdh5"))
DotPlot(merge, features = c("Krt18", "Cdh11", "Ptprc", "Hba-a2", "Cdh5"))

# differential expression
# markers for cluster 1
cluster1.markers <- FindMarkers(object = merge, ident.1 = 1, min.pct = 0.25)
head(x = cluster1.markers, n = 5)

# markers for all clusters
merge.markers <- FindAllMarkers(object = merge, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25)
merge.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Putting new labels on clusters
new.cluster.ids <- c("A", "B", "C", "D", "E", "F","G", "H", "I", "J", "K", "L", "M", "N",
                     "O", "P", "Q")
names(x = new.cluster.ids) <- levels(x = merge)
merge <- RenameIdents(object = merge, new.cluster.ids)
DimPlot(object = merge, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
merge$labeled_clusters <- Idents(object = merge) # save identity to meta.data

# DON'T MOVE ONE TO LABELING UNTIL YOU'VE FILTERED ENOUGH
# otherwise --> repeat from filtering step

# VISUALIZATIONS
# DoHeatmap
top10 <- merge.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = merge, features = top10$gene, size = 2) + NoLegend()

# pheatmap
# make sure when fetching data that there is no repeated genes
cells <- FetchData(object = merge, vars = c(unique(top10$gene), 'labeled_clusters'), 
                   slot = 'scale.data')
cells <- cells[order(cells$labeled_clusters),]

anno <- data.frame(cluster = cells$labeled_clusters)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:147], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells[,1:147]), show_colnames = F,
         cluster_rows = F, cluster_col = F,
         annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#8A2BE2','#000000','#FFFF00'))(250))
