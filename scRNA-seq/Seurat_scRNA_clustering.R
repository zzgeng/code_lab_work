library(Seurat)
library(tidyr)
library(ggplot2)

setwd("/Volumes/ZZ/scRNA_Seq/G268_reads/")

## input data
BO.flu<-read.table("G268_all_reads.txt")
row.names(BO.flu)<-gsub("-",".",row.names(BO.flu))   
BO.flu.het<-BO.flu[,grep("het",colnames(BO.flu))]
BO.flu.homo<-BO.flu[,grep("homo",colnames(BO.flu))]

##seurat analysis and version check
print(packageVersion('Seurat'))
{
##Create seurat object, with min.cells = 40, min.features = 500
pbmc <- CreateSeuratObject(counts = BO.flu, project = "AUTS2_scRNA", min.cells = 40, min.features = 500)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

##Normalizing teh data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

##Find the highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

##Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

##perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

##Determine the dimensionality of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 1000)
pbmc <- ScoreJackStraw(pbmc, dim = 1:20)

JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

##Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:7)
pbmc <- FindClusters(pbmc, resolution = 0.3)

##try non-linear regression
####Using tsne
pbmc <- RunTSNE(pbmc, dims = 1:7)
DimPlot(pbmc, reduction = "tsne", label = T)

####using ump
pbmc <- RunUMAP(pbmc, dims = 1:5)
DimPlot(pbmc, reduction = "umap",  label = T)

##save the data
saveRDS(pbmc, "/Volumes/ZZ/ZZ_scRNA/pbmc.robj")
}

##plot some of the selected genes and observe the distribution among clusters

#for cluster0/6
FeaturePlot(pbmc, features = c("CLIC6", "HTR2C", "MSX1", "TTR", "KRT18",  "CLDN5", "IGFBP7"), reduction = "tsne")

##for neural markers
FeaturePlot(pbmc, features = c("MAP2", "TUBB3", "FOXG1", "STMN2", "FABP7", "TUBA1A"), reduction = "tsne")

FeaturePlot(pbmc, features = c("NEUROD1", "NEUROD4", "HES1", "VSX2"), reduction = "tsne")

##violin plot
VlnPlot(pbmc, features = c( "AUTS2", "TTR", "CLIC6", "MKI67", "NEUROD1", "SOX9", "CSF1R"), log = T)


##find markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top10 <- pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

##run GO on each cluster
library(clusterProfiler)
library(org.Hs.eg.db)

cluster_marker <- pbmc.markers %>% 
  group_by(cluster) 

cluster_marker <- as.data.frame(cluster_marker)

for(x in 0:7){
  filename <- paste0("GO-cluster", x)
  go <- enrichGO(cluster_marker[cluster_marker$cluster == x, 7], ont = "BP", OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
  assign(filename, go)
}


merge_result(list(cluster_0 = `GO-cluster0`, cluster_1 = `GO-cluster1`, 
                  cluster_2 = `GO-cluster2`, cluster_3 = `GO-cluster3`, 
                  cluster_4 = `GO-cluster4`, cluster_5 = `GO-cluster5`, 
                  cluster_6 = `GO-cluster6`, cluster_7 = `GO-cluster7`)) %>%
  dotplot(.,showCategory = 5)

DoHeatmap(pbmc, features = c("AUTS2",top10$gene)) + NoLegend()

Charles_marker <- c("AUTS2", 
                    "MAP2", "SOX4", "TUBB3", "MLLT11", "DCX", "NEUROD6", "GAP43","STMN2",
                    "FABP7","TUBB2B","TUBA1A", "NFIB","FOXG1", 
                    "MKI67","NES", "SOX2", "PAX6","TOP2A",
                    "NEUROD1", "NEUROD4","OPN1SW", "RCVRN", "GNB3",
                    "TTR","CLIC6","SLC4A5","FOLR1", "IGFBP7", 
                    "POSTN", "COL3A1", "COL1A1", "COL1A2", "PLAT")
DoHeatmap(pbmc, features = Charles_marker) + NoLegend()

##Label genotypes
genotyping <- c()
for(x in 1:1111){
  ifelse(grepl("homo", colnames(pbmc)[x]), 
          genotyping[x] <-  "homo", 
          genotyping[x] <-  "het")
}

pbmc$genotype <- genotyping
pbmc$CellType <- Idents(pbmc)
Idents(pbmc) <- "genotype"
DimPlot(pbmc, reduction = "tsne")


save.image("/Volumes/ZZ/ZZ_scRNA/scRNA.RData")