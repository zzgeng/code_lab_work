library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(cowplot)
library(ggplot2)

##Combine our data with published data using seurat

BO.flu<-read.table("/Volumes/ZZ/scRNA_Seq/G268_reads/G268_all_reads.txt")
pbmc <- CreateSeuratObject(counts = BO.flu, project = "AUTS2_scRNA", min.cells = 40, min.features = 500)
pbmc@meta.data <- pbmc@meta.data[, !colnames(pbmc@meta.data) %in% c("genotype", "CellType")]
pbmc@meta.data[, "protocol"] <- "ours"

##published data, sample down to 2000 cells
sc_all <- read.csv("Desktop/HBO/scRNA_sampleDown//GSM2295946_dge.3m.csv")
rownames(sc_all) <- sc_all$X
sc_all <- sc_all[,-1]
set.seed(101)
sc_2000 <- sc_all[sample(1:15402, 2000)]
rm(sc_all)

#Import published data and label the meta.data as "published"
pbmc_2000 <- CreateSeuratObject(counts = sc_2000, project = "sample_2000", min.cells = 3, min.features = 200)
pbmc_2000@meta.data[, "protocol"] <- "published"

##Combine our data with data from publication
pbmc.combine <- merge(pbmc, pbmc_2000, add.cell.ids = c("ours", "published"), project = "protocol")

species.list <- SplitObject(pbmc.combine, split.by = "protocol")

species.list <- lapply(species.list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(species.list)

species.list <- lapply(X = species.list, FUN = function(x){
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

species.anchors <- FindIntegrationAnchors(object.list = species.list, anchor.features = features, reduction = "rpca")
species.integrated <- IntegrateData(anchorset = species.anchors)

DefaultAssay(species.integrated) <- "integrated"
species.integrated <- ScaleData(species.integrated, verbose = FALSE)
species.integrated <- RunPCA(species.integrated, npcs = 30, verbose = FALSE)
species.integrated <- JackStraw(species.integrated, num.replicate = 100)
species.integrated <- ScoreJackStraw(species.integrated, dim = 1:20)

JackStrawPlot(species.integrated, dims = 1:20)
ElbowPlot(species.integrated)
species.integrated <- RunUMAP(species.integrated, reduction = "pca",dims = 1:30)
species.integrated <- FindNeighbors(species.integrated, reduction = "pca", dims = 1:30)
species.integrated <- FindClusters(species.integrated, resolution = 0.6)

DimPlot(species.integrated, reduction = "umap", group.by = "protocol")
DimPlot(species.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)

##Label the cells based on the genotype
genotyping <- c()
for(x in 1:3111){
  if(grepl("homo", colnames(species.integrated)[x])){
    genotyping[x] <-  "homo"
  } 
  if(grepl("het", colnames(species.integrated)[x])){
    genotyping[x] <-  "het"}
  if(grepl("Part", colnames(species.integrated)[x])){
    genotyping[x] <-  "wt"
  }
}

species.integrated$genotype <- genotyping
species.integrated$CellType <- Idents(species.integrated)
Idents(species.integrated) <- "genotype"
DimPlot(species.integrated, reduction = "umap", group.by = "genotype")

DimPlot(species.integrated, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T)
p2 <- DimPlot(species.integrated, reduction = "umap", group.by = "RNA_snn_res.0.8")
plot_grid(p1, p2)

marker <- FindAllMarkers(species.integrated, group.by = "seurat_clusters", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- marker %>% 
  group_by(cluster) %>% 
  top_n(n = 3, wt = avg_log2FC)

DoHeatmap(species.integrated, features = c(top10$gene), draw.lines = T) 

FeaturePlot(species.integrated, features = c("CLIC6", "TTR", "KRT18", "MAP2", "STMN2", "NEFL", "NEUROD1", "NEUROD4", "RCVRN"), min.cutoff = c("q1", "q10"))

FeaturePlot(species.integrated, features = c("CLIC6", "KRT18", "STMN2", "NEFL", "NEUROD1", "NEUROD4", "MKI67", "TOP2A"), min.cutoff = c("q1", "q10"))


FeaturePlot(species.integrated, features = c("SLC1A3", "SOX3", "ID4", "APOE"), min.cutoff = c("q1", "q10"))

#Neural proliferative procursors
FeaturePlot(species.integrated, features = c("MKI67", "TOP2A", "CENPF", "UBE2C", "NUSAP1", "NES"), min.cutoff = c("q1", "q10"))

#neuroepithelial cells
FeaturePlot(species.integrated, features = c("CHGB", "PCP4", "CLDN11", "ENPP2", "HHIP", "PLAT"), min.cutoff = c("q1", "q10"))

FeaturePlot(species.integrated, features = c("POSTN", "FOXG1"), min.cutoff = c("q1", "q10"))

#Astroglia
FeaturePlot(species.integrated, features = c("GFAP", "AGT", "AQP4", "GJA1"), min.cutoff = c("q1", "q10"))



Charles_marker <- c("AUTS2", 
                    "MAP2", "SOX4", "TUBB3", "MLLT11", "DCX", "NEUROD6", "GAP43","STMN2",
                    "C1orf61", "FABP7", "MARCKSL1","TUBB2B","TUBA1A", "NFIB","POU3F3", "FOXG1", "FOXG1", "POU3F2", "PTPRZ1",
                    "TOP2A", "CENPF", "SGOL2",  "CCNB2", "UBE2C", "NUSAP1", "ASPM", "CDCA8", "GFAP", "VIM",
                    "NES", "SOX2", "PAX6","MKI67", "SOX9", "NOTCH1",  "SIX3","MSI1",
                    "RCVRN", "NEUROD1", "GNB3", "PDE6H", "PCBP4", "AIPL1", "TULP1", 
                    "POSTN", "COL3A1", "COL1A1", "COL1A2", "BGN", "IGFBP4", "DCN", "SFRP2", "PLAT", 
                    "TTR",  "IGFBP7", "CXCL14", "CA2", "KRT18", "CLIC6","SLC4A5","FOLR1")

DoHeatmap(species.integrated, features = Charles_marker, draw.lines = T, combine = F) 


##plot the proportion of different genotypes in each cluster
freq_table <- prop.table(x = table(species.integrated@active.ident, 
                                   species.integrated@meta.data[, "genotype"]),
                         margin = 2)
barplot(freq_table)

marker_results <- cluster_identity(marker)
marker_results <- Marker_plot(test, top_n = 2)
Heatmap(marker_results, cluster_rows = F, cluster_columns = F, column_split = paste0("C", 0:(ncol(marker_results)-1)), row_split = rep(paste0("C", 0:(ncol(marker_results)-1)), each = 2))

genotype_percentage <- species.integrated@meta.data
ggplot(genotype_percentage, aes(x = seurat_clusters, fill = genotype)) + geom_bar(position = "fill")
ggplot(genotype_percentage, aes(x = seurat_clusters, fill = genotype)) + geom_bar()




