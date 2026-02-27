# Server code for shinyapp, currently using template output file.
# Last edited: 27/02/26. To
# With 'server.R' and 'ui.R' files in 'ShinyApp/' subdirectory of wd, run command >runApp("ShinyApp") in console.


library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(rsconnect)
library(dplyr)
library(Seurat)
library(patchwork)

# Issues with consistently having this command work, change to correctly denote location of output.rds if needed
# This has to be loaded in both ui and server to correctly work
loadedpbmc <- LoadSeuratRds("./output.rds")

server <- function(input, output) {
  
  # Original pipeline graph 1
  # Selects specific clusters to show
  # Allows user to input specific cells to highlight
  output$originalOne<- renderPlotly(
    DimPlot(loadedpbmc, cells.highlight=input$cellsID ,cells=WhichCells(loadedpbmc, idents=input$SelectCluster), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  )
  
  # Original pipeline graph 2
  # Selects specific clusters to show
  output$originalTwo <- renderPlot({ 
    DimPlot(loadedpbmc, cells=WhichCells(loadedpbmc, idents=input$SelectCluster), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  })
  
}