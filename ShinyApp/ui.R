# UI code for shinyapp, currently using template output file.
# Last edited: 27/02/26. Tom
# With 'server.R' and 'ui.R' files in 'ShinyApp/' subdirectory of wd, run command >runApp("ShinyApp") in console.

library(shiny)
library(bslib)
library(ggplot2)
library(plotly)
library(rsconnect)
library(dplyr)
library(Seurat)
library(patchwork)

# This has to be loaded in both ui and server to correctly work
loadedpbmc <- LoadSeuratRds("./output.rds")

# Currently not in output.rds, could be beneficial to add this as a column to avoid coding it here.
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")

# Changed to taglist, gives a nice page layout, may need to be changed again in the future
ui <- tagList(
  navbarPage(
    
    # Website title
    "Group B",
    
    # Original pipeline tab
    tabPanel("Original Pipeline",
             sidebarPanel(
               checkboxGroupInput("SelectCluster", 
                                  "Select Clusters", 
                                  new.cluster.ids
                                  ),
               selectizeInput(
                 "cellsID",
                 "Cell ID",
                 choices = loadedpbmc$nCount_RNA,
                 multiple = TRUE,
                 
                 ),

             ),
             mainPanel(
               tabsetPanel(
                 
                 # Original pipeline graph one goes here
                 tabPanel("Graph 1",
                          
                          plotlyOutput("originalOne"),
                          
                 ),
                 
                 # Original pipeline graph two goes here
                 tabPanel("Graph 2",
                          plotOutput("originalTwo")
                 ),
                 
                 # Original pipeline graph three goes here
                 tabPanel("Graph 3", "Graph 3 will go here")
               )
             )
    ),
    
    # Alternative pipeline page
    tabPanel("Alternative Pipeline",
             mainPanel(
               tabsetPanel(
                 
                 # All code to do with graph one of the alternative pipeline goes here
                 tabPanel("Graph 1",
                          
                          # Trialing showing checkbox menu inside tabPanel, so it only shows when a specific graph is selected
                          sidebarPanel(
                            checkboxGroupInput("SelectCluster", 
                                               "Select Clusters", 
                                               new.cluster.ids), 
                          ),
                          
                          
                 ),
                 
                 # All code to do with graph 2 of the alternative pipeline goes here
                 tabPanel("Graph 2",
                          
                 ),
                 
                 # All code to do with graph 3 of the alternative pipeline goes here
                 tabPanel("Graph 3",
                          
                 )
               ),
             ),
    )
  )
)