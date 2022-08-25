library(monocle3)
library(ggplot2)
library(Seurat)

## input data
pbmc <- readRDS("/Volumes/ZZ/ZZ_scRNA/pbmc.robj")
DimPlot(pbmc, reduction = "tsne", label = T, pt.size = 3)

exp_matrix <- as.matrix(pbmc@assays$RNA@counts)

fd <- data.frame(gene_short_name = row.names(exp_matrix), 
                 row.names = row.names(exp_matrix))
pd <- data.frame(pbmc@meta.data)
x <- pd$seurat_clusters
levels(x) <- list(ChP_like = "0", Neuron = "1", Photoreceptor = "2", Unknown = "3",
                  RGC = "4", Unknown = "5", ChP = "6", Unknown = "7", Neuro_epithelial = "8")


pd$for_plot <- x
cds <- new_cell_data_set(exp_matrix, cell_metadata = pd, gene_metadata = fd)

##Determine the num_dim for monocle3
for(x in c(50:80)){
  cds <- preprocess_cds(cds, num_dim = x, norm_method = "size_only")
  cds <- align_cds(cds)
  cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  DIR_file <- paste0("Desktop/HBO/ZZ_scRNA/", x, ".png")
  png(filename = DIR_file)
  print(plot_cells(cds, genes = c("SOX2", "NES", "PAX6", "MKI67", "STMN2", "FOXG1", "CLIC6", "NEUROD4", "AUTS2"),
                   alpha = 0.8, cell_size = 1, 
                   reduction_method = "UMAP", scale = F) + 
          scale_color_continuous(low = "grey", high = "blue") )
  dev.off()
}

cds <- preprocess_cds(cds, num_dim = 51, norm_method = "size_only")
plot_pc_variance_explained(cds)
cds <- align_cds(cds)
#cds <- reduce_dimension(cds, reduction_method = "tSNE", preprocess_method = "PCA")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")

plot_cells(cds, genes = c("FOXG1", "AUTS2", "CLIC6", "NEUROD4", "POSTN", "MKI67"),
           trajectory_graph_segment_size = .3, 
           alpha = 1, cell_size = 1, 
           reduction_method = "UMAP") + 
  scale_color_continuous(low = "grey", high = "blue")

cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, 
           genes = c("SOX2", "NES", "PAX6", "MKI67", "STMN2","CLIC6", "NEUROD4", "AUTS2"),
           alpha = 0.8, cell_size = 1, 
           reduction_method = "UMAP", scale = F) + 
  scale_color_continuous(low = "grey", high = "blue")

#Set neuro-epithelial cells as root
cds <- learn_graph(cds)
cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1, 
           reduction_method = "UMAP", label_cell_groups = F, 
           trajectory_graph_segment_size = 0.5,
           labels_per_group = F, show_trajectory_graph = T, label_branch_points = F, 
           label_principal_points = F, label_leaves = F, label_roots = F 
           )


feature <- c("STMN2", "FOXG1", "GAP43", "NES", "MKI67", "TOP2A", "NEUROD1", "RCVRN", "GNB1",  "TTR", "KRT18", "CLIC6", "POSTN", "COL3A1", "COL1A1")

plot_cells(cds, genes = c("STMN2", "FOXG1", "GAP43"), 
           alpha = 1, cell_size = 1, 
           reduction_method = "UMAP", trajectory_graph_segment_size = 0.5,
           labels_per_group = F, show_trajectory_graph = T, label_branch_points = F, 
           label_principal_points = F, label_leaves = F, label_roots = F) + 
  scale_color_continuous(low = "grey", high = "blue")

plot_cells(cds, genes = c("NES", "MKI67", "TOP2A"), 
           alpha = 1, cell_size = 1, 
           reduction_method = "UMAP", trajectory_graph_segment_size = 0.5,
           labels_per_group = F, show_trajectory_graph = T, label_branch_points = F, 
           label_principal_points = F, label_leaves = F, label_roots = F) + 
  scale_color_continuous(low = "grey", high = "blue")

plot_cells(cds, genes = c("NEUROD1", "RCVRN", "GNB1"), 
           alpha = 1, cell_size = 1, 
           reduction_method = "UMAP", trajectory_graph_segment_size = 0.5,
           labels_per_group = F, show_trajectory_graph = T, label_branch_points = F, 
           label_principal_points = F, label_leaves = F, label_roots = F) + 
  scale_color_continuous(low = "grey", high = "blue")

plot_cells(cds, genes = c("KRT18", "CLIC6", "TTR"), 
           alpha = 1, cell_size = 1, 
           reduction_method = "UMAP", trajectory_graph_segment_size = 0.5,
           labels_per_group = F, show_trajectory_graph = T, label_branch_points = F, 
           label_principal_points = F, label_leaves = F, label_roots = F) + 
  scale_color_continuous(low = "grey", high = "blue")

plot_cells(cds, genes = c("POSTN", "COL3A1", "COL1A1"), 
           alpha = 1, cell_size = 1, 
           reduction_method = "UMAP", trajectory_graph_segment_size = 0.5,
           labels_per_group = F, show_trajectory_graph = T, label_branch_points = F, 
           label_principal_points = F, label_leaves = F, label_roots = F) + 
  scale_color_continuous(low = "grey", high = "blue")



plot_cells(cds, color_cells_by = "for_plot", 
           cell_size = 1, reduction_method = "UMAP", 
           label_cell_groups = F, 
           trajectory_graph_segment_size = 0.5,
           labels_per_group = F, show_trajectory_graph = T, label_branch_points = F, 
           label_principal_points = F, label_leaves = F, label_roots = F) + 
  facet_wrap(vars(for_plot), nrow = 3)

plot_cells(cds, color_cells_by = "genotype", 
           cell_size = 1, reduction_method = "UMAP", 
           label_cell_groups = F, 
           trajectory_graph_segment_size = 0.5,
           labels_per_group = F, show_trajectory_graph = T, label_branch_points = F, 
           label_principal_points = F, label_leaves = F, label_roots = F)


gene_name <- c("NES", "MAP2", "AUTS2", "CLIC6", "NEUROD1")
cds_filter <- cds[gene_name,]
plot_genes_in_pseudotime(cds_filter)





