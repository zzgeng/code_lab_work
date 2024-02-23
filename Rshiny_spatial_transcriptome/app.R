#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(Seurat)

#load("/System/Volumes/Data/PSU/Visium_female_Rshiny/gene_exp_SCT.RData")
load("/Volumes/ZZ/aflah/spatial rna/R/Rshiny/gene_exp_SCT.RData")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Visium"),

    # Sidebar with a slider input for number of bins 
    textInput(inputId = "Gene",label = "Genes") ,

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("Exp_hist"),
           plotOutput("Exp_hist_2"),
           plotOutput("Spatial")
        )
    )


# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$Exp_hist <- renderPlot({
      gene_list <- input$Gene
      gene_list <- tolower(gene_list)
      gene_list <- tools::toTitleCase(gene_list)
      x <- c(WT_F@assays[["SCT"]]@scale.data[rownames(WT_F@assays[["SCT"]]@scale.data) == gene_list, ], KO_F@assays[["SCT"]]@scale.data[rownames(KO_F@assays[["SCT"]]@scale.data) == gene_list, ])
      Temp_data <- data.frame(Gene = x,
                                genotype = x_meta$orig.ident,
                                cluster = x_meta$manual)
        
        ggplot(Temp_data, aes(x = genotype, y = Gene, fill = genotype)) + 
          geom_violin() + geom_boxplot() +
          facet_grid(~cluster, scales = "free_x", space = "free_x")
          
    })
    output$Exp_hist_2 <- renderPlot({
      gene_list <- input$Gene
      gene_list <- tolower(gene_list)
      gene_list <- tools::toTitleCase(gene_list)
      x <- c(WT_F@assays[["SCT"]]@scale.data[rownames(WT_F@assays[["SCT"]]@scale.data) == gene_list, ], KO_F@assays[["SCT"]]@scale.data[rownames(KO_F@assays[["SCT"]]@scale.data) == gene_list, ])
      Temp_data <- data.frame(Gene = x,
                              genotype = x_meta$orig.ident,
                              cluster = x_meta$seurat_clusters)
      
      ggplot(Temp_data, aes(x = genotype, y = Gene, fill = genotype)) + 
        geom_violin() + geom_boxplot() +
        facet_grid(~cluster, scales = "free_x", space = "free_x")
      
    })
    
    output$Spatial <- renderPlot({
      gene_list <- input$Gene
      gene_list <- tolower(gene_list)
      gene_list <- tools::toTitleCase(gene_list)
      Visium_max = max(max(as.matrix(GetAssayData(WT_F, slot = "scale.data"))[rownames(as.matrix(GetAssayData(WT_F, slot = "scale.data"))) == gene_list,]), max(as.matrix(GetAssayData(KO_F, slot = "scale.data"))[rownames(as.matrix(GetAssayData(KO_F, slot = "scale.data"))) == gene_list,]))
      Visium_min = min(min(as.matrix(GetAssayData(WT_F, slot = "scale.data"))[rownames(as.matrix(GetAssayData(WT_F, slot = "scale.data"))) == gene_list,]), min(as.matrix(GetAssayData(KO_F, slot = "scale.data"))[rownames(as.matrix(GetAssayData(KO_F, slot = "scale.data"))) == gene_list,]))
      color_range <- c(Visium_min, Visium_max)
      color_midpoint = (min(color_range) + max(color_range))/2
      
      p1 <- SpatialFeaturePlot(WT_F, features = gene_list, alpha = 1, slot = "scale.data") + 
        ggtitle("WT_F") + 
        scale_fill_gradient2(limits = color_range, low = "blue4", high = "red2", mid = "khaki", midpoint = color_midpoint) +
        theme(legend.position = "right")
      
      p3 <- SpatialFeaturePlot(KO_F, features = gene_list, alpha = 1, slot = "scale.data") + 
        ggtitle("KO_F") + 
        scale_fill_gradient2(limits = color_range, low = "blue4", high = "red2", mid = "khaki", midpoint = color_midpoint) +
        theme(legend.position = "right")
      
      p3+p1
      
    })
    
     
}

# Run the application 
shinyApp(ui = ui, server = server)
