
library(shiny)
library(bslib)

# Define UI for application that draws the plots
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "minty"),
  
  navbarPage("RNA Seq Visualization",
             tabPanel("Mapping",
                      sidebarLayout(
                        sidebarPanel(
                          fileInput("file", "Upload a .bam or .sam file", buttonLabel = "Upload..."),
                          downloadButton("downloadData", "Download"),
                          ),
                        mainPanel(
                          
                          )
                        )
                      ),
             tabPanel("Count Matrix",
                      sidebarLayout(
                        sidebarPanel(
                          fileInput("file", "Upload a .txt file", buttonLabel = "Upload..."),
                          downloadButton("downloadData", "Download"),
                        ),
                        mainPanel(
                     
                          )
                        )
                      ),
             navbarMenu("Differential Gene Expression",
                        tabPanel("Heatmap",
                                 sidebarLayout(
                                   sidebarPanel(
                                     fileInput("file", "Upload normalized counts", buttonLabel = "Upload..."),
                                     downloadButton("downloadData", "Download"),
                                     ),
                                   mainPanel(
                                     
                                   )
                                   )
                                 ),
                        tabPanel("PCA Analysis",
                                 sidebarLayout(
                                   sidebarPanel(
                                     fileInput("file", "Upload normalized counts", buttonLabel = "Upload..."),
                                     downloadButton("downloadData", "Download"),
                                   ),
                                   mainPanel(
                                     
                                   )
                                 )
                        )
                        )
             )
  )

# Define server logic required to display the plots
server <- function(input, output) {
  

}

# Run the application 
shinyApp(ui = ui, server = server)
