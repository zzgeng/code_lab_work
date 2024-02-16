#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#

#options(shiny.reactlog = 100000)
#Sys.setenv(VROOM_CONNECTION_SIZE = 5000000)

library(shiny)
library(tidyverse)
library(ggplot2)
library(patchwork)

## read files
exp_matrix_P0 <- readr::read_tsv("./data/P0.tsv.gz")
exp_matrix_E14 <- readr::read_tsv("./data/E14.tsv.gz")

# Genes
bait <- c("Rnf2")
PRC1_component <- c("Pcgf1", "Pcgf2", "Pcgf3", "Bmi1", "Pcgf5", "Pcgf6")
canonical_component <- c("Cbx2", "Cbx4", "Cbx6", "Cbx7", "Cbx8", "Phc1", "Phc2")
non_canonical <- c("Rybp", "Dcaf7", "Auts2", "Yaf2", "Pogz", "Kdm2b", "Bcor", "Bcorl1", "Fbrsl1", "Fbrs", "Max", "Mga", "L3mbtl2", "Wdr5", "Hdac1", "Hdac2", "E2f6")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("scRNA-seq mouse brain"),
    
    sidebarLayout(
      sidebarPanel(
        # Select box with terms
        selectInput("gene_set", "Select Promoter:", 
                    choices = c("PRC1 core component", "Canonical", "Non-canonical", "All"),
                    selected = "All")  # Optional: set default selected term
      ,
      
    textInput(inputId = "Gene",label = "Promoter"),
      ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
    output$distPlot <- renderPlot({
      promoter_name <- str_to_title(input$Gene)
      bait <- c(bait, promoter_name)
      exp_matrix_E14 <- exp_matrix_E14[exp_matrix_E14[,promoter_name]>0, ]
      exp_matrix_P0 <- exp_matrix_P0[exp_matrix_P0[,promoter_name]>0, ]
      
        if(input$gene_set == "All"){
          
          list_gene <- c(bait, PRC1_component, canonical_component, non_canonical)
          P0_df <- exp_matrix_P0 %>% 
            select(list_gene) %>% 
            reshape2::melt() %>% 
            mutate(label = case_when(variable %in% bait ~ "bait", 
                                        variable %in% PRC1_component ~ "PRC1_component", 
                                        variable %in% canonical_component ~ "canonical_component",
                                        variable %in% non_canonical ~ "non_canonical"))
          E14_df <- exp_matrix_E14 %>% 
            select(list_gene) %>% 
            reshape2::melt() %>% 
            mutate(label = case_when(variable %in% bait ~ "bait", 
                                        variable %in% PRC1_component ~ "PRC1_component", 
                                        variable %in% canonical_component ~ "canonical_component",
                                        variable %in% non_canonical ~ "non_canonical"))
          
          P0_plt <- ggplot(P0_df, aes(x = factor(variable, levels = c(bait, PRC1_component, canonical_component, non_canonical)), 
                                               y = value, color = label)) + 
            geom_jitter(alpha = 0.3, width = 0.3) + 
            geom_boxplot(alpha = 0.5, width = 0.3, outlier.shape = NA) + 
            xlab("") + 
            ylab("") + 
            theme(
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(color = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            ) + 
            ggtitle("P0 mouse brain")
          
          E14_plt <- ggplot(E14_df, aes(x = factor(variable, levels = c(bait, PRC1_component, canonical_component, non_canonical)), 
                                               y = value, color = label)) + 
            geom_jitter(alpha = 0.3, width = 0.3) + 
            geom_boxplot(alpha = 0.5, width = 0.3, outlier.shape = NA) + 
            xlab("") + 
            ylab("") + 
            theme(
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(color = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            ) + 
            ggtitle("E14 mouse brain")
          
          P0_plt / E14_plt
          
        } else if (input$gene_set == "PRC1 core component"){
          
          list_gene <- c(bait, PRC1_component)
          P0_df <- exp_matrix_P0 %>% 
            select(list_gene) %>% 
            reshape2::melt() %>% 
            mutate(label = case_when(variable %in% bait ~ "bait", 
                                     variable %in% PRC1_component ~ "PRC1_component", 
                                     #variable %in% canonical_component ~ "canonical_component",
                                     #variable %in% non_canonical ~ "non_canonical"
                                     ))
          E14_df <- exp_matrix_E14 %>% 
            select(list_gene) %>% 
            reshape2::melt() %>% 
            mutate(label = case_when(variable %in% bait ~ "bait", 
                                     variable %in% PRC1_component ~ "PRC1_component", 
                                     #variable %in% canonical_component ~ "canonical_component",
                                     #variable %in% non_canonical ~ "non_canonical"
                                     ))
          
          P0_plt <- ggplot(P0_df, aes(x = factor(variable, levels = c(bait, PRC1_component)), 
                                      y = value, color = label)) + 
            geom_jitter(alpha = 0.3, width = 0.3) + 
            geom_boxplot(alpha = 0.5, width = 0.3, outlier.shape = NA) + 
            xlab("") + 
            ylab("") + 
            theme(
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(color = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            ) + 
            ggtitle("P0 mouse brain")
          
          E14_plt <- ggplot(E14_df, aes(x = factor(variable, levels = c(bait, PRC1_component)), 
                                        y = value, color = label)) + 
            geom_jitter(alpha = 0.3, width = 0.3) + 
            geom_boxplot(alpha = 0.5, width = 0.3, outlier.shape = NA) + 
            xlab("") + 
            ylab("") + 
            theme(
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(color = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            ) + 
            ggtitle("E14 mouse brain")
          
          P0_plt / E14_plt
          
        } else if (input$gene_set == "Canonical"){
          
          list_gene <- c(bait, canonical_component)
          P0_df <- exp_matrix_P0 %>% 
            select(list_gene) %>% 
            reshape2::melt() %>% 
            mutate(label = case_when(variable %in% bait ~ "bait", 
                                     variable %in% PRC1_component ~ "PRC1_component", 
                                     variable %in% canonical_component ~ "canonical_component",
                                     variable %in% non_canonical ~ "non_canonical"
                                     ))
          E14_df <- exp_matrix_E14 %>% 
            select(list_gene) %>% 
            reshape2::melt() %>% 
            mutate(label = case_when(variable %in% bait ~ "bait", 
                                     variable %in% PRC1_component ~ "PRC1_component", 
                                     variable %in% canonical_component ~ "canonical_component",
                                     variable %in% non_canonical ~ "non_canonical"
                                     ))
          
          P0_plt <- ggplot(P0_df, aes(x = factor(variable, levels = c(bait, canonical_component)), 
                                      y = value, color = label)) + 
            geom_jitter(alpha = 0.3, width = 0.3) + 
            geom_boxplot(alpha = 0.5, width = 0.3, outlier.shape = NA) + 
            xlab("") + 
            ylab("") + 
            theme(
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(color = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            ) + 
            ggtitle("P0 mouse brain")
          
          E14_plt <- ggplot(E14_df, aes(x = factor(variable, levels = c(bait, canonical_component)), 
                                        y = value, color = label)) + 
            geom_jitter(alpha = 0.3, width = 0.3) + 
            geom_boxplot(alpha = 0.5, width = 0.3, outlier.shape = NA) + 
            xlab("") + 
            ylab("") + 
            theme(
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(color = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            ) + 
            ggtitle("E14 mouse brain")
          
          P0_plt / E14_plt
          
        } else if (input$gene_set == "Non-canonical"){
          
          list_gene <- c(bait, non_canonical)
          P0_df <- exp_matrix_P0 %>% 
            select(list_gene) %>% 
            reshape2::melt() %>% 
            mutate(label = case_when(variable %in% bait ~ "bait", 
                                     variable %in% PRC1_component ~ "PRC1_component", 
                                     variable %in% canonical_component ~ "canonical_component",
                                     variable %in% non_canonical ~ "non_canonical"))
          E14_df <- exp_matrix_E14 %>% 
            select(list_gene) %>% 
            reshape2::melt() %>% 
            mutate(label = case_when(variable %in% bait ~ "bait", 
                                     variable %in% PRC1_component ~ "PRC1_component", 
                                     variable %in% canonical_component ~ "canonical_component",
                                     variable %in% non_canonical ~ "non_canonical"))
          
          P0_plt <- ggplot(P0_df, aes(x = factor(variable, levels = c(bait, non_canonical)), 
                                      y = value, color = label)) + 
            geom_jitter(alpha = 0.3, width = 0.3) + 
            geom_boxplot(alpha = 0.5, width = 0.3, outlier.shape = NA) + 
            xlab("") + 
            ylab("") + 
            theme(
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(color = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            ) + 
            ggtitle("P0 mouse brain")
          
          E14_plt <- ggplot(E14_df, aes(x = factor(variable, levels = c(bait, non_canonical)), 
                                        y = value, color = label)) + 
            geom_jitter(alpha = 0.3, width = 0.3) + 
            geom_boxplot(alpha = 0.5, width = 0.3, outlier.shape = NA) + 
            xlab("") + 
            ylab("") + 
            theme(
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(color = "black"), 
              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)
            ) + 
            ggtitle("E14 mouse brain")
          
          P0_plt / E14_plt
          
        } 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
