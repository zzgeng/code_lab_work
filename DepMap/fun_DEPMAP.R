##Extract a selected gene's crispr score, and gene expression
CRISPR_Exp <- function(exp_matrix = exp_matrix, 
                       crisp_score = crisp_score, 
                       sample_info = sample_info, 
                       select_gene){
  gene_name <- paste0("^", select_gene)
  YFG_exp <- exp_matrix[, c(1, grep(gene_name, colnames(exp_matrix)))]
  YFG_crispr <- crisp_score[, c(1, grep(gene_name, colnames(crisp_score)))]
  YFG <- merge(YFG_exp, YFG_crispr, by = 0)
  YFG <- YFG[, c(1,3,5)]
  colnames(YFG) <- c("DepMap_ID", "Expression", "CRISPR_score")
  YFG <- merge(YFG, sample_info, by = "DepMap_ID")
  return(YFG)
}

##plot scatter
plot_scatter_YFG_crispr <- function(input, threhold_crispr = -0.5, threhold_exp = 0, col_by = NULL){
  ggplot(input, aes_string(x = "CRISPR_score", y = "Expression", color = col_by)) + 
    geom_point(position = "identity") + 
    geom_hline(yintercept = threhold_exp, col = "red", linetype = "dashed") + 
    geom_vline(xintercept = c(threhold_crispr, 0), col = c("red", "black"), linetype = "dashed")
}

#plot bar-plot to show the statics of cancers (cancer type, sample collection site, etc) that are affected by KO YFG
plot_hist_YFG_crispr <- function(input, threhold_crispr = -0.5, threhold_exp = 0, group_by = "primary_disease"){
  input_filtered <- input[input$CRISPR_score < -0.5 & input$Expression > threhold_exp, ]
  input_filtered$count <- 1
  input_filtered <- aggregate(input_filtered$count, by = list(input_filtered[, group_by]), FUN = sum)
  colnames(input_filtered) <- c("A", "count")
  ggplot(input_filtered, aes(x = reorder(A, count), y = count)) + 
    geom_bar(stat = "identity") +
    coord_flip() + 
    xlab(label = NULL) + 
    ggtitle(label = group_by)
}

##Volcano plot for DEG
plot_volcano <- function(input, FC_cutoff = 1, padj_cutoff = 0.01){
  xlim <- abs(input$log2FoldChange[order(abs(input$log2FoldChange), decreasing = T)])[1]
  with(input, plot(x = log2FoldChange, y = -log10(pvalue), 
                   pch = 20, col = "grey", xlim = c(-xlim, xlim)))
  with(subset(input, padj < padj_cutoff), points(x = log2FoldChange, y = -log10(pvalue), pch = 20, col = "blue"))
  with(subset(input, abs(log2FoldChange) > FC_cutoff & padj < padj_cutoff), points(x = log2FoldChange, y = -log10(pvalue), pch = 20, col = "red"))
  abline(v = c(-FC_cutoff, FC_cutoff), col = "red", lty = 2)
  abline(h = -log10(padj_cutoff), col = "red", lty = 2)
}

##Cellligner hist
plot_hist_Cellligner <- function(input, cell_line, group_by = "disease"){
  cluster_N <- input[grep(cell_line, input$sampleID_CCLE_Name), "cluster"]
  input_cluster <- input[input$cluster == cluster_N,]
  input_cluster$count <- 1
  input_cluster <- aggregate(input_cluster$count, by = list(input_cluster[, group_by]), FUN = sum)
  colnames(input_cluster) <- c("A", "count")
  ggplot(input_cluster, aes(x = reorder(A, count), y = count)) + 
    geom_bar(stat = "identity") +
    coord_flip() + 
    xlab(label = group_by) + ggtitle(label = paste0("Cluster ", cluster_N))
}

##Cellligner scatterplot
plot_hist_Cellligner <- function(input, cell_line, group_by = "disease"){
  cluster_N <- input[grep(cell_line, input$sampleID_CCLE_Name), "cluster"]
  input_cluster <- input[input$cluster == cluster_N,]
  input_cluster$count <- 1
  input_cluster <- aggregate(input_cluster$count, by = list(input_cluster[, group_by]), FUN = sum)
  colnames(input_cluster) <- c("A", "count")
  ggplot(input_cluster, aes(x = reorder(A, count), y = count)) + 
    geom_bar(stat = "identity") +
    coord_flip() + 
    ylab(label = paste0("cluster_", cluster_N))
}


plot_crispr_Cellligner <- function(input, cluster_N, gene, crisp_score_input = crisp_score, threhold_crispr = -0.5){
  input_cluster <- input[input$cluster == cluster_N,]
  Gene_name <- paste0("^", gene)
  crisp_score_selected <- crisp_score[, c(1,grep(Gene_name, colnames(crisp_score_input)))]
  crisp_score_selected <- crisp_score_selected[rownames(crisp_score_selected) %in% input_cluster$sampleID, ]
  colnames(crisp_score_selected) <- c("SAMPLE_ID", "select_gene")
  crisp_score_selected$cluster <- cluster_N
  
  xlabel <- paste0(nrow(crisp_score_selected[crisp_score_selected$select_gene < threhold_crispr,]),  " / ", nrow(crisp_score_selected))
  ggplot(crisp_score_selected, aes(x = cluster, y = select_gene)) + 
    geom_violin() + geom_jitter() + 
    geom_hline(yintercept = threhold_crispr, col = "red", linetype = "dashed") + 
    xlab(label = xlabel)
}







