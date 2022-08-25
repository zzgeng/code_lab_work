library(tidyr)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(ComplexHeatmap)
library(ggpubr)
library(annotables)

##import data from DepMap
exp_matrix <- read.csv("Desktop/DepMap/Data/CCLE_expression.csv")
rownames(exp_matrix) <- exp_matrix$X
crisp_score <- read.csv("Desktop/DepMap/Data/CRISPR_gene_effect.csv")
rownames(crisp_score) <- crisp_score$DepMap_ID
sample_info <- read.csv("Desktop/DepMap/Data/sample_info.csv")

source("Desktop/DepMap/fun_DEPMAP.R")
##Extract a selected gene's crispr score, and gene expression
BMI1 <- CRISPR_Exp(exp_matrix = exp_matrix, 
                   crisp_score = crisp_score, 
                   sample_info = sample_info, 
                   select_gene = "BMI1")

##identify the cell line dependent on RHOA
png("Desktop/DepMap/figures_BMI1//scatter_CRIPSR_exp.png", width = 1000, height = 600)
plot_scatter_YFG_crispr(BMI1, threhold_crispr = -0.5, threhold_exp = 5, col_by = "sample_collection_site")
dev.off()

png("Desktop/DepMap/figures_BMI1//cellline_info_BMI1.png", width = 1200, height = 600)
collection_site <- plot_hist_YFG_crispr(BMI1, group_by = "sample_collection_site")
pri_disease <- plot_hist_YFG_crispr(BMI1, group_by = "primary_disease")
ggarrange(collection_site, pri_disease, 
          labels = c("A", "B"), ncol = 2, nrow = 1)
dev.off()


###Download gene abundancy matrix from NCBI/GEO
GEO_exp <- read.delim("Desktop/DepMap/Data/GSE176520_readcount_genes_all.txt")

##convert EMSEMBL to gene symbol
GEO_exp <- merge(GEO_exp, grch38, by.x = "gene_id", by.y = "ensgene")
GEO_exp <- aggregate(GEO_exp[, 2:10], by = list(GEO_exp$gene_name), FUN = sum)
rownames(GEO_exp) <- GEO_exp$Group.1
colnames(GEO_exp)
GEO_exp <- GEO_exp[, c(2:10)]
colnames(GEO_exp)

meta_Data <- data.frame(condition = c(rep("control", 3), rep("treated_24hr", 3), rep("treated_48hr", 3)), cellType = rep("Rh28"))

dds <- DESeqDataSetFromMatrix(countData = GEO_exp, 
                              colData = meta_Data, 
                              design = ~condition)
dds <- DESeq(dds)
res1 <- results(dds, contrast = c("condition", "control", "treated_48hr"))
res2 <- results(dds, contrast = c("condition", "control", "treated_24hr"))

##Volcano plot
png("Desktop/DepMap/figures_BMI1/volcano_res1.png")
plot_volcano(res1, padj_cutoff = 0.1)
dev.off()

png("Desktop/DepMap/figures_BMI1/volcano_res2.png")
plot_volcano(res2, padj_cutoff = 0.1)
dev.off()


##perform GO analysis
up_list_1 <- row.names(res1[(res1$padj < 0.01 & res1$log2FoldChange > 0) %in% TRUE, ])
down_list_1 <- row.names(res1[(res1$padj < 0.01 & res1$log2FoldChange < 0) %in% TRUE, ])
up_GO_1 <- enrichGO(up_list_1, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
down_GO_1 <- enrichGO(down_list_1, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

png("Desktop/DepMap/figures_BMI1//GO_upregualted.png", width = 600, height = 800)
dotplot(up_GO_1, showCategory = 15)
dev.off()

png("Desktop/DepMap/figures_RHOA/GO_downregualted.png", width = 600, height = 800)
dotplot(down_GO_1, showCategory = 15)
dev.off()

up_list_2 <- row.names(res2[(res2$padj < 0.01 & res2$log2FoldChange > 0) %in% TRUE, ])
down_list_2 <- row.names(res2[(res2$padj < 0.01 & res2$log2FoldChange < 0) %in% TRUE, ])
up_GO_2 <- enrichGO(up_list_2, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
down_GO_2 <- enrichGO(down_list_2, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

png("Desktop/DepMap/figures_BMI1//GO_upregualted_2.png", width = 600, height = 800)
dotplot(up_GO_2, showCategory = 15)
dev.off()

png("Desktop/DepMap/figures_RHOA/GO_downregualted_2.png", width = 600, height = 800)
dotplot(down_GO_2, showCategory = 15)
dev.off()

##Heatmap for DEG
GEO_exp_normalized <- counts(dds, normalized = TRUE)
GEO_exp_DEG <- GEO_exp_normalized[rownames(GEO_exp_normalized) %in% c(up_list_1, down_list_1, up_list_2, down_list_2),]
hmp <- Heatmap(t(scale(t(GEO_exp_DEG))), cluster_columns = F, show_row_names = FALSE, km = 4)
draw(hmp);RowOrder <- row_order(hmp)

for(x in 1:4){
  Gene_List <- rownames(GEO_exp_DEG[RowOrder[[x]],])
  temp <- enrichGO(Gene_List, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
  file_name <- paste0("Cluster_", x)
  assign(file_name, temp)
}

png("Desktop/DepMap/figures_RHOA/GO_each_cluster.png", width = 600, height = 800)
merge_result(list(cluster_1 = Cluster_1, cluster_2 = Cluster_2, cluster_3 = Cluster_3, cluster_4 = Cluster_4)) %>% 
  dotplot(., showCategory = 5)
dev.off()

##test the expression changes of BMI1 target genes
HOX <- GEO_exp_normalized[grep("^HOX",rownames(GEO_exp_normalized)), ]
HOX <- HOX[rowMeans(HOX) > 0,]

png("Desktop/DepMap/figures_BMI1/HOX_heatmap.png", width = 800, height = 800)
Heatmap(t(scale(t(HOX))), cluster_columns = F)
dev.off()

##Celligner
##no information for cellline 1321N1 from GSE111571. The the parternal line of 1321N1, U118MG
Celligner <- read.csv("Desktop/DepMap/Data/Celligner_info.csv")

png("Desktop/DepMap/figures_BMI1//celligner_hist.png", width = 1000, height = 1000)
plot_hist_Cellligner(input = Celligner, cell_line = "RH28",group_by = "disease")
dev.off()

