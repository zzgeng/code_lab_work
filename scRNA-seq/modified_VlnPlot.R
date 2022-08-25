library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)

x <- NormalizeData(pbmc)
x_exp <- as.matrix(x@assays$RNA@data)
pd <- x@meta.data$seurat_clusters

feature <- c("AUTS2", "FOXG1", "MAP2", "STMN2", "MKI67", "TOP2A", "CENPF", "RCVRN", "NEUROD1", "NEUROD4", "NES", "PAX6", "SOX2", "DCX", "TTR", "CLIC6", "KRT18")


plot_vio <- data.frame(auts2 = t(x_exp[rownames(x_exp) %in% feature,]), pd = pd)
plot_test <- melt(plot_vio, id.vars = "pd")

for(x in 0:8){
  temp <- plot_test[plot_test$pd == x, ]
  temp <- temp %>% dplyr::group_by(variable) %>%
    mutate(med = median(value))
  temp <- tidyr::separate(temp, variable, into = c("a", "variable"))
  temp$variable <- factor(temp$variable, levels = feature)
  temp_plot <- ggplot(temp, aes(x = variable, y = value, fill = med)) + 
    geom_violin(scale = "width", show.legend = FALSE) + 
    scale_fill_gradientn(limits = c(0, 8), colors = c("lemonchiffon", "coral2")) +
    theme(axis.title.y=element_blank()) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  
  file_name <- paste0("Plot_C_",x)
  assign(value = temp_plot, x = file_name)
}

rm(temp)
rm(temp_plot)

ggarrange(Plot_C_0, Plot_C_1, Plot_C_2, Plot_C_3, Plot_C_4, Plot_C_5, Plot_C_6, Plot_C_7, Plot_C_8, ncol = 1, nrow = 9)

temp_plot <- ggplot(temp, aes(x = variable, y = value, fill = med)) + 
  geom_violin(scale = "width") + 
  scale_fill_gradientn(limits = c(0, 8), colors = c("lemonchiffon", "coral2")) +
  theme(axis.title.y=element_blank()) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
temp_plot



FeaturePlot(pbmc, features = c("SOX2", "DCX", "PCNA", "CR", "CB"))
