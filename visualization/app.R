# Load necessary libraries
library(shiny)
library(shinydashboard)
library(limma)
library(edgeR)
library(pheatmap)
library(ggrepel)
library(ggplot2)
library(plotly)

# Creating the front end of the application
ui <- dashboardPage(
  # Setting the layout for the title, file upload and normalization button
  dashboardHeader(title = "RNAgenie"),
  dashboardSidebar(
    fileInput("count_table", "Upload count table"),
    actionButton("normalization", "Normalize Data")
    ),
  
  # Setting the layout for the plots that will get displayed
  dashboardBody(
    tabsetPanel(
      tabPanel("Counts Table", dataTableOutput("data_table")
               ),
      tabPanel("Plots",
               fluidRow(
                 box(plotOutput("raw_plot", height = 400), width = 6),
                 box(plotOutput("normalized_plot", height = 400), width = 6)
                 ),
               fluidRow(
                 box(plotlyOutput("pca_plot", height = 400), width = 6),
                 box(plotlyOutput("box_plot", height = 400), width = 6)
                 )
               ),
      )
    )
  )

# Initializing server to create plots
server <- function(input, output) {
  # Loading count table from user input file and processing it to use it
  # in the next steps
  data <- reactive({
    req(input$count_table)
    inFile <- input$count_table
    data <- read.table(inFile$datapath, header = TRUE, row.names = 1)

    # Log transforming the data to avoid the insertion of negative numbers
    data <- data + 1
    data <- log2(data)
    return(data)
  })
  
  output$data_table <- renderDataTable({
    req(data())
    data()
    })

  # Creating a density plot for raw count table
  output$raw_plot <- renderPlot({
    counts <- data()
    plot(density(colSums(counts)), main = "Raw Count Data", xlab = "Total Counts", ylab = "Density")
  })

  # Normalizing count table and creating the plot
  normalized_count <- eventReactive(input$normalization, {
    counts <- data()
    y <- DGEList(counts)
    y <- calcNormFactors(y)
    norm_factors <- y$samples$norm.factors
    norm_factors[norm_factors == Inf | is.na(norm_factors)] <- 1
    norm_counts <- counts / norm_factors
    v <- voom(norm_counts, design = NULL, plot = FALSE)
    n <- data.frame(v$E)
    colnames(n) <- make.names(colnames(n))
    row.names(n) <- row.names(norm_counts)
    list(norm_counts = norm_counts, E = n)
  })

  # Creating a density plot for normalized count table
  output$normalized_plot <- renderPlot({
    norm_counts <- normalized_count()$norm_counts
    v <- voom(norm_counts, design = NULL, plot = FALSE)
    norm_data <- data.frame(v$E)
    rownames(norm_data) <- rownames(norm_counts)
    norm_data
    plot(density(colSums(norm_data)), main = "Normalized Count Data", xlab = "Total Counts", ylab = "Density")
  })

  # Creating PCA plot using plotly for an interactive view
  output$pca_plot <- renderPlotly({
    norms <- normalized_count()$n
    pca <- prcomp(t(norms))
    variance <- round(summary(pca)$percentage[2,1]*100, 1)
    p <- plot_ly(data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=factor(colnames(norms)))) %>%
      add_markers(x = ~PC1, y = ~PC2, color = ~sample,
                  marker = list(size = 8, opacity = 0.8)) %>%
      layout(title = paste("PCA Plot: Variance between PC1 and PC2 = ", variance, "%", sep=""),
             xaxis = list(title = paste("PC1 (", summary(pca)$percentage[1, 1], "%)", sep="")),
             yaxis = list(title = paste("PC2 (", summary(pca)$percentage[2, 1], "%)", sep="")))
    p
  })

  # Creating boxplot using plotly for an interactive view
  output$box_plot <- renderPlotly({
    norm_counts <- normalized_count()$norm_counts
    boxplot_df <- data.frame(Sample = colnames(norm_counts), Counts = as.numeric(as.vector(t(norm_counts))))
    b <- plot_ly(boxplot_df, x = ~Sample, y = ~log10(Counts), type = "box") %>%
      layout(title = "Box plot", xaxis = list(title = "Sample"), yaxis = list(title = "log10 Counts"))
    b
  })
}

# Running application
shinyApp(ui, server)
