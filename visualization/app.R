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
  dashboardHeader(title = "RNA Genie"),
  dashboardSidebar(
    fileInput("count_table", "Upload count data"),
    actionButton("normalize_button", "Normalize Data")
  ),
  dashboardBody(
    fluidRow(
      box(plotOutput("raw_counts_plot", height = 400), width = 6),
      box(plotOutput("normalized_counts_plot", height = 400), width = 6)
    ),
    fluidRow(
      box(plotlyOutput("pca_plot"), height = 400), width = 6,
      box(plotlyOutput("box_plot"), height = 400), width = 6
    )  
  )
)

# Initializing server to create plots
server <- function(input, output) {
  
  # Load count data
  data <- reactive({
    req(input$count_table)
    
    inFile <- input$count_table
    data <- read.table(inFile$datapath, header = TRUE, row.names = 1)
    
    data <- data + 1
    data <- log2(data)
    
    return(data)
  })
  
 # Generate raw count data plot
  output$raw_counts_plot <- renderPlot({
    counts <- data()
    plot(density(colSums(counts)), main = "Raw Count Data", xlab = "Total Counts", ylab = "Density")
  })

  # Normalize count data and generate plot
  normalized_data <- eventReactive(input$normalize_button, {
    counts <- data()
    y <- DGEList(counts)
    y <- calcNormFactors(y)
    norm_factors <- y$samples$norm.factors
    norm_factors[norm_factors == Inf | is.na(norm_factors)] <- 1
    norm_counts <- counts / norm_factors
    v <- voom(norm_counts, design = NULL, plot = FALSE)
    E <- data.frame(v$E)
    colnames(E) <- make.names(colnames(E))
    row.names(E) <- row.names(norm_counts)
    list(norm_counts = norm_counts, E = E)
  })

  # Generate normalized count data plot
  output$normalized_counts_plot <- renderPlot({
    norm_counts <- normalized_data()$norm_counts
    plot(density(colSums(norm_counts)), main = "Normalized Count Data", xlab = "Total Counts", ylab = "Density")
  })

  # Generate PCA plot
  output$pca_plot <- renderPlotly({
    norms <- normalized_data()$E
    pca <- prcomp(t(norms))
    variance <- round(summary(pca)$importance[2,1]*100, 1)
    fig <- plot_ly(data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=factor(colnames(norms)))) %>%
      add_markers(x = ~PC1, y = ~PC2, color = ~sample,
                  marker = list(size = 8, opacity = 0.8)) %>%
      layout(title = paste("PCA Plot: Variance explained by PC1 and PC2 = ", variance, "%", sep=""),
             xaxis = list(title = paste("PC1 (", summary(pca)$importance[1, 1], "%)", sep="")),
             yaxis = list(title = paste("PC2 (", summary(pca)$importance[2, 1], "%)", sep="")))
    fig
  })
  
  # Generate boxplot
  output$box_plot <- renderPlotly({
    counts <- data()
    y <- DGEList(counts)
    y <- calcNormFactors(y)
    norm_factors <- y$samples$norm.factors
    norm_factors[norm_factors == Inf | is.na(norm_factors)] <- 1
    norm_counts <- counts / norm_factors
    data_df <- data.frame(Sample = colnames(norm_counts), Counts = as.numeric(as.vector(t(norm_counts))))
    fig <- plot_ly(data_df, x = ~Sample, y = ~log10(Counts), type = "box") %>%
      layout(title = "Box plot", xaxis = list(title = "Sample"), yaxis = list(title = "log10 Counts"))
    fig
  })
}

# Running application
shinyApp(ui, server)
