# Server code for shinyapp, currently using template output file.
# Last edited: 26/02/26. To
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
loadedpbmc <- LoadSeuratRds("./output.rds")

server <- function(input, output) {
  output$page1<- renderPlotly(
    DimPlot(loadedpbmc, cells=WhichCells(loadedpbmc, idents=input$SelectCluster), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  )
  
  output$page2 <- renderPlot({ 
    DimPlot(loadedpbmc, cells=WhichCells(loadedpbmc, idents=input$SelectCluster), reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  })
  
}