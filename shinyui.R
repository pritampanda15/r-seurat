library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)
library(hdf5r)
library(ggsci)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(glmGamPoi)
library(cowplot)
library(pheatmap)
library(scran)
library(SingleR)
library(presto)
library(ggrepel)
library(dplyr)

# Define UI
ui <- fluidPage(
  titlePanel("Seurat Processing App"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload h5 file")
    ),
    mainPanel(
      plotOutput("vlnplot"),
      plotOutput("dimplot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  output$vlnplot <- renderPlot({
    req(input$file)
    
    # Read uploaded h5 file
    h5_file <- input$file$datapath
    adj_matrix <- Read10X_h5(h5_file, use.names = TRUE)
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = adj_matrix, project = 'TISCH2', min.cells = 3, min.features = 200)
    
    # Generate vlnplot
    plot1 <- VlnPlot(object = seurat_obj, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt', 'percent.rb'), ncol = 4)
    plot1
  })
  
  output$dimplot <- renderPlot({
    req(input$file)
    
    # Read uploaded h5 file
    h5_file <- input$file$datapath
    adj_matrix <- Read10X_h5(h5_file, use.names = TRUE)
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = adj_matrix, project = 'TISCH2', min.cells = 3, min.features = 200)
    
    # Run PCA
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    
    # Generate dimplot
    plot2 <- DimPlot(object = seurat_obj, label = TRUE)
    plot2
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
