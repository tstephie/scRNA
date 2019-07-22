# Key Visualization Code
require(Seurat)
require(dplyr)
require(magrittr)
require(data.table)
require(pheatmap)
require(scales)
require(RColorBrewer)

# Heatmap
# pheatmap w/ 1 annotation bar
## need data.frame w/ dim = (cells,genes/variables)
cells <- FetchData(object = merge_2, vars = c(unique(rownames(group_DE)), 'group'), 
                   slot = 'scale.data')
cells <- cells[order(cells$group),]

anno <- data.frame(cluster = cells$group)
rownames(anno) <- rownames(cells)

cells <- apply(cells[,1:63], MARGIN = 2, FUN = function(x) (squish(x, range = c(-2.5,2.5))))

pheatmap(mat = t(cells[,1:63]), show_colnames = F,
         cluster_rows = F, cluster_col = F,
         annotation_col = anno, fontsize = 7,
         color = colorRampPalette(colors = c('#8A2BE2','#000000','#FFFF00'))(250))

# Bar

# Volcano 

# Box

# Scatter

# Histogram

