# UI code for shinyapp, currently using template output file.
# Last edited: 26/02/26. Tom
# With 'server.R' and 'ui.R' files in 'ShinyApp/' subdirectory of wd, run command >runApp("ShinyApp") in console.


library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(rsconnect)
library(dplyr)
library(Seurat)
library(patchwork)

# Currently not in output.rds, could be beneficial to add this as a column to avoid coding it here.
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
ui <- page_fluid(
  theme = bs_theme(version = 5, bootswatch = "lux"),
  sidebarLayout(
    sidebarPanel(
      
      checkboxGroupInput("SelectCluster", 
                         "Select Clusters", 
                         new.cluster.ids,
      ), 
    ),
    
    mainPanel(
      navset_card_underline( 
        nav_panel("Page 1", plotlyOutput(outputId = "page1")), 
        nav_panel("Page 2", plotOutput("page2"))
      ),
    ),
  ),
)
