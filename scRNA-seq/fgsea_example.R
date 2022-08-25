library(fgsea)
library(DESeq2)

##gmt files can be downloaded from 
##http://www.gsea-msigdb.org/gsea/downloads.jsp
hm_fgsea <- gmtPathways("Desktop/HBO/GSEA dataset/h.all.v7.5.1.symbols.gmt")

##DEG_res is the DEG file from DESeq2 
rank <- DEG_res[!is.na(DEG_res$log2FoldChange), 2]
names(rank) <- rownames(DEG_res[!is.na(DEG_res$log2FoldChange),])
fgseaRes <- fgsea(pathways = hm_fgsea, 
                  stats = rank, 
                  minSize = 15, 
                  maxSize = 500)

##plot fgsea table 
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(hm_fgsea[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)
