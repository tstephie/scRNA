# Seurat v4 pipeline

library(tidyverse)
library(Seurat)
library(scImpute)
library(harmony)
library(pheatmap)
library(RColorBrewer)
library(scales)

# Run through the analysis first and then if you see that you need to correct for batch effects
# Example of analysis between two groups

# PREPROCESSING
## using filtered as example
KOFIL.data <- Read10X("KO/KO/filtered_gene_bc_matrices/mm10/")
WTFIL.data <- Read10X("WT/WT/filtered_gene_bc_matrices/mm10/")

# Put label on cell names if you have more than one group
colnames(KOFIL.data) <- paste("KO_", colnames(KOFIL.data), sep = '')
colnames(WTFIL.data) <- paste("WT_", colnames(WTFIL.data), sep = '')

'----------------------------------------------------------------------------------------------------------------'

# # Dropout correction
# # THIS NEEDS A LOT OF MEMORY TO RUN!
# # Will take some time to run
# combined.data <- cbind(KOFIL.data, WTFIL.data) # merge data
# saveRDS(combined.data, file = 'raw_count_combined_KPC.rds')
# count_path <- "./raw_count_combined_KPC.rds"
# infile="rds"
# outfile="rds"
# drop_thre=0.5 # threshold to determine dropout values
# ncores=1 # if working on Windows, set to 1
# out_dir="./"
# combined_scimpute <- scimpute(count_path, infile, outfile, out_dir, 
#                                    labeled = FALSE,  drop_thre, ncores = ncores, 
#                                    type = "count", Kcluster = 5)
# scimpute_count <- readRDS("scimpute_count.rds")
# # Create seurat object with dropout corrected data
# merge <- CreateSeuratObject(counts = scimpute_count, project = 'MERGED',  min.cells = 3, min.features = 100)

'----------------------------------------------------------------------------------------------------------------'

## Create KO and WT Seurat objects
KOFIL <- CreateSeuratObject(counts = KOFIL.data, project = 'KOFIL', min.cells = 3, min.features = 100)
WTFIL <- CreateSeuratObject(counts = WTFIL.data, project = 'WTFIL', min.cells = 3, min.features = 100)

## Add any meta.data needed
### OPTION 1 (old):
KOFIL@meta.data$group <- "KO"
WTFIL@meta.data$group <- "WT"

## If more than 2 objects, then need to merge
## Add IDs to cell names to distinguish groups
merge <- merge(x = KOFIL, y = WTFIL, add.cell.ids = c('KO', 'WT'))

## OPTION 2:
### merge objects together and rename cells
merge <- merge(x = KOFIL, y = WIFIL, add.cell.ids = c('KO', 'WT'))

### Add meta.data needed
### get number of cells for each sample
n_cells <- c(length(grep('^KO', colnames(merge))), length(grep('^WT', colnames(merge))))

metadata <- data.frame(run = rep(c('1234', '5678'), times = n_cells),
                       state = rep(c('KO', 'WT'), times = n_cells), 
                       row.names = colnames(merge), stringsAsFactors = F)

### add metadata to merged object
merge <- AddMetaData(merge, metadata = metadata)

### OPTION 3: Add groups to meta.data by using current ident
merge$group <- Idents(object = merge)

'----------------------------------------------------------------------------------------------------------------'

# Cell cycle genes (if want to do cell cycle correction)
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

'----------------------------------------------------------------------------------------------------------------'

# IMPORTANT
# Changing between meta.data for identities
Idents(object = merge) <- 'group'
# Check active identity
levels(merge)

'----------------------------------------------------------------------------------------------------------------'

# QC
## get percent of mitochondrial genes reads
merge[["percent.mt"]] <- PercentageFeatureSet(object = merge, pattern = "^mt-")

## Look at nFeature, nCount, percent.mt
VlnPlot(object = merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        pt.size = .5)

## Another way to look at nFeature, nCount, percent.mt
plot1 <- FeatureScatter(object = merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

## Filter cells based on nFeature
### filtering cells that have nCount < 50000 and percent.mt < 20%
merge <- subset(x = merge, subset = nCount_RNA < 50000 & percent.mt < 25)

'----------------------------------------------------------------------------------------------------------------'

# NORMALIZATION
merge <- NormalizeData(object = merge, normalization.method = "LogNormalize", 
                       scale.factor = 10000)
## Can use code below if want default parameters (same thing as above)
## merge <- NormalizeData(object = merge)

'----------------------------------------------------------------------------------------------------------------'

# VARIABLE FEATURE
## find subset of features that exhibit high cell-to-cell variation in the dataset
## vst --> loess.span, clip.max, nfeaturesmerge 
merge <- FindVariableFeatures(object = merge, selection.method = 'vst', loess.span = .3, nfeatures = 2000)
# change top variable genes you want to keep w/ nfeatures

top10 <- head(x = VariableFeatures(object = merge), 10)

plot1 <- VariableFeaturePlot(object = merge)
plot2 <- LabelPoints(plot = plot1, points = top10)
CombinePlots(plots = list(plot1, plot2))
# ignore warning for now

## mean.var.plot --> mean.function, dispersion.function, num.bin, mean/dispersion.cutoff
merge <- FindVariableFeatures(object = merge, selection.method = 'mean.var.plot', ,
                              mean.cutoff = c(0.0125,Inf), dispersion.cutoff = c(0.5,Inf))

'----------------------------------------------------------------------------------------------------------------'

# SCALING
## Scale the data
## if leave features paramter blank, then it will use variable genes
## if don't use all genes, then DoHeatmap function might not work
## OPTION 1 (all genes):
all.genes <- rownames(merge)
merge <- ScaleData(object = merge, features = all.genes,
                   vars.to.regress = c("nCount", "percent.mt"))
## OPTION 2 (only variable genes):
merge <- ScaleData(object = merge,
                   vars.to.regress = c("nCount", "percent.mt"))
## OPTION 3 (all genes w/ cell cycling correction): 
all.genes <- rownames(merge)
# use if want all cell cycle signals regressed out
merge <- CellCycleScoring(object = merge, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
merge <- ScaleData(object = merge, features = all.genes, 
                   vars.to.regress = c("nCount", "percent.mt", "S.Score", "G2M.Score")) 
### use if want some cell cycle signals regressed out (cycling vs. non-cycling will save)
merge$CC.Difference <- merge$S.Score - merge$G2M.Score
merge <- ScaleData(object = merge, features = all.genes, 
                   vars.to.regress = c("nCount", "percent.mt", "CC.Difference"))

'----------------------------------------------------------------------------------------------------------------'

# DIMENSIONAL REDUCTION: PCA
## Can control the number of pcs that you want to calculate w/ npcs
merge <- RunPCA(object = merge, features = VariableFeatures(object = merge), 
                nfeatures.print = 5)

# Choosing how many PCs for clustering
# visualizing pcs
ElbowPlot(object = merge, ndims = 50)

# JackStraw if you want
merge <- JackStraw(object = merge, num.replicate = 100)
merge <- ScoreJackStraw(object = merge, dims = 1:20)
JackStrawPlot(object = merge, dims = 1:15)

# Calculating --> until get over 90% variance
stdev <- merge[['pca']]@stdev
var <- (merge[['pca']]@stdev)^2
sum(var[1:31])/sum(var) 

'----------------------------------------------------------------------------------------------------------------'
# BATCH CORRECTION
## if know that you have batch effect from either previous info or from running through clustering and visualization
## using harmony package
merge <- RunHarmony(merge, group.by.vars = c('run'), dims.use = 1:31, verbose = F)

'----------------------------------------------------------------------------------------------------------------'

# CLUSTERING
merge <- FindNeighbors(object = merge, dims = 1:31, k.param = 20)
merge <- FindClusters(object = merge, resolution = 0.5)
head(x = Idents(object = merge), 5)

## if batch corrected:
merge <- FindNeighbors(object = merge, dims = 1:31, k.param = 20, reduction = 'harmony')
merge <- FindClusters(object = merge, resolution = 0.5)
head(x = Idents(object = merge), 5)

'----------------------------------------------------------------------------------------------------------------'

# VISUALIZATION: TSNE/UMAP
## TSNE example
merge <- RunTSNE(merge, dims = 1:31)
DimPlot(object = merge, reduction = "tsne")

## UMAP example
merge <- RunUMAP(merge, dims = 1:31)
DimPlot(object = merge, reduction = 'umap')

### if batch corrected:
merge <- RunUMAP(merge, dims = 1:31, reduction = 'harmony')
DimPlot(object = merge, reduction = 'umap')

## Visualizing expression 
FeaturePlot(merge, features = c("Krt18", "Cdh11", "Ptprc", "Hba-a2", "Cdh5"))
DotPlot(merge, features = c("Krt18", "Cdh11", "Ptprc", "Hba-a2", "Cdh5"))

'----------------------------------------------------------------------------------------------------------------'

# DIFFERENTIAL EXPRESSION
## markers for cluster 1
cluster1.markers <- FindMarkers(object = merge, ident.1 = 1, min.pct = 0.25)
head(x = cluster1.markers, n = 5)

## markers for all clusters
merge.markers <- FindAllMarkers(object = merge, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25)
merge.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

'----------------------------------------------------------------------------------------------------------------'

# LABELING CLUSTERS
# Putting new labels on clusters
new.cluster.ids <- c("A", "B", "C", "D", "E", "F","G", "H", "I", "J", "K", "L", "M", "N",
                     "O")
names(x = new.cluster.ids) <- levels(x = merge)
merge <- RenameIdents(object = merge, new.cluster.ids)
DimPlot(object = merge, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
merge$labeled_clusters <- Idents(object = merge) # save identity to meta.data

# DON'T MOVE ONE TO LABELING UNTIL YOU'VE FILTERED ENOUGH
# otherwise --> repeat from filtering step


'----------------------------------------------------------------------------------------------------------------'

# VISUALIZATIONS
# Get top n genes for each cluster
top10 <- merge.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

# DoHeatmap
DoHeatmap(object = merge, features = top10$gene, size = 2) + NoLegend()

# pheatmap
## make sure when fetching data that there is no repeated genes
cells <- FetchData(object = merge, vars = c(unique(top10$gene), 'labeled_clusters'), 
                   slot = 'scale.data')
cells <- cells[order(cells$labeled_clusters),]

## creating annotation bars
anno <- data.frame(cluster = cells$labeled_clusters)
rownames(anno) <- rownames(cells)

## rescaling
cells <- apply(cells[,1:147], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells[,1:147]), show_colnames = F,
         cluster_rows = F, cluster_col = F,
         annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#8A2BE2','#000000','#FFFF00'))(250))

'--------------'
# EXTRA 

## Subset specific cell type: myeloid
Idents(merge) <- 'labeled_clusters'
myeloid <- subset(merge, subset = labeled_clusters == 'Myeloid')

## Creating metadata based on current idents
old_ids <- levels(merge)
new_ids <- c('A', 'B')

merge$test <- plyr::mapvalues(merge$group, from = old_ids, to = new_ids)

Idents(merge) <- 'test'
levels(merge)
table(merge$test, merge$group)


