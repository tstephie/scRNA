# Interactome template

library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)

source(file.choose()) # import functions

# THIS PART MUST BE DONE BEFORE USING THE FUNCTIONS

# Separate objects
Idents(Tissue_Merge) <- 'Patient_ID'

# object_0: norm_19_227
# object_1: pda_181429
object_0 <- subset(Tissue_Merge, subset = Patient_ID == '19_227')
object_1 <- subset(Tissue_Merge, subset = Patient_ID == '181429')

# make sure that it has 2 columns (ligand-receptor)
LR_Pairs_database <- as.data.frame(fread("Literature Supported Receptor Ligand Pairs Ramilowski.csv", header = T, stringsAsFactors = F))
LR_Pairs_database <- LR_Pairs_database[,1:2]

# list of wanted genes
# no CCL12, CCR12, CCL6, CCL9
# 181429: no CCL24, CCL25, CCL1, CCR8, CCR3
ligand_genes <- as.character(LR_Pairs_database$Ligand.ApprovedSymbol)
receptor_genes <- as.character(LR_Pairs_database$Receptor.ApprovedSymbol)

# Create merged Seurat object
meta_label <- 'Simple_Relabelling'
# Switch ids of objects to meta_label
Idents(object_0) <- meta_label
Idents(object_1) <- meta_label
merge <- merge(object_0, object_1)

# Find common cell type populations and rename if need
object0_cell_types <- as.character(levels(object_0))
object1_cell_types <- as.character(levels(object_1))
cell_types <- object0_cell_types[which(object0_cell_types %in% object1_cell_types)]
cell_types <- cell_types[c(1, 2, 3, 5, 6, 8)] # where you choose cell types want to look at

# Get info for only common cell types
Idents(merge) <- meta_label
merge <- subset(merge, idents = cell_types)

# Get gene list in object 
merge_genes <- rownames(merge)

# Create group Seurat objects
Idents(merge) <- 'State'
normal <- subset(merge, idents = 'DuoAd')
pda <- subset(merge, idents = 'PDA')
Idents(normal) <- meta_label
Idents(pda) <- meta_label

'-------------------------------------------------------------------------------'
not_found_ligands <- check_genes(ligand_genes, LR_Pairs_database, merge_genes) 
not_found_receptors <- check_genes(receptor_genes, LR_Pairs_database, merge_genes)

# separate genes into ligands and receptors and get rid of genes not found
not_in_lig <- unique(not_found_ligands$object)
wanted_ligands <- unique(ligand_genes[!ligand_genes %in% not_in_lig])

not_in_rec <- unique(not_found_receptors$object)
wanted_receptors <- unique(receptor_genes[!receptor_genes %in% not_in_rec])

genes <- unique(c(wanted_ligands, wanted_receptors))

# get need info from object into table
# change second vars to what meta.data refers to labels
normal_data <- FetchData(object = normal, vars = c(genes, meta_label))
pda_data <- FetchData(object = pda, vars = c(genes, meta_label))

# Make all possible LR pairs
LRs <- make_LR_pairs(wanted_ligands, wanted_receptors, LR_Pairs_database)

# Calculate average expression for all cell types for each group
avg_0 <- AverageExpression(object = normal, features = genes, verbose = F)
avg_1 <- AverageExpression(object = pda, features = genes, verbose = F)

# Create table w/ average expressions for every ligand/receptor and source/target cell combination
LR_table <- create_LR_table(wanted_ligands, wanted_receptors, cell_types, LRs, avg_0$RNA, avg_1$RNA)
LR_table[,1:4] <- data.frame(apply(LR_table[,1:4], 2, as.character), stringsAsFactors = FALSE)

# Choose threshold --> here I'm using median
summary(c(LR_table$avg_lig_0, LR_table$avg_lig_1, LR_table$avg_rec_0, LR_table$avg_rec_1))
threshold <- median(c(LR_table$avg_lig_0, LR_table$avg_lig_1, LR_table$avg_rec_0, LR_table$avg_rec_1))

# Filter LR pairs based on average expresions
LR_table <- avg_LR_filt(LR_table, threshold)

# Calculate p-value between groups (can add alpha = )
LR_table <- LR_diff(LR_table, normal_data, pda_data, genes, meta_label)

# Filter LR pairs where expression is constant between groups
LR_table <- LR_table[which(LR_table$lig_diff != 0 | LR_table$rec_diff != 0),]

# Sort by ligand
LR_table <- arrange(LR_table, desc(lig_diff != 0 & rec_diff != 0), lig_diff_p_adj)

# This goes into cytoscape
write.csv(LR_table, file = 'human_pda_61354830_61356750_LR_052719_lig_sort.csv', row.names = F)


# list of genes
genes_list <- c(LR_table$ligands, LR_table$receptors)
genes_list <- unique(genes_list)
write.csv(genes_list, file = 'human_pda_61354830_61356750_LR_052719_genes.csv')

# get supplemental data from Cytoscape (need to have at least shared name column)
cytoscape_nodes <- fread('human_pda_61354830_61356750_LR_052719_order.csv')
node_supplement <- generate_supplement(LR_table, cytoscape_nodes)
write.csv(node_supplement, file = 'human_pda_61354830_61356750_LR_names_052719.csv', row.names = F)
