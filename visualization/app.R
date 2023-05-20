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
  dashboardHeader(title = "RNAgenie"),
  dashboardSidebar(
    fileInput("count_table", "Upload count table"),
    actionButton("normalization", "Normalize Data")
  ),
  dashboardBody(
    fluidRow(
      box(plotOutput("raw_plot", height = 400), width = 6),
      box(plotOutput("normalized_plot", height = 400), width = 6)
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
  
 # Create raw count data plot using a density plot
  output$raw_plot <- renderPlot({
    counts <- data()
    plot(density(colSums(counts)), main = "Raw Count Data", xlab = "Total Counts", ylab = "Density")
  })

  # Normalize count data and create plot
  normalized_count <- eventReactive(input$normalization, {
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

  # Create normalized count data plot using a density plot
  output$normalized_plot <- renderPlot({
    norm_counts <- normalized_count()$norm_counts
    plot(density(colSums(norm_counts)), main = "Normalized Count Data", xlab = "Total Counts", ylab = "Density")
  })

  # Create PCA plot using plotly for an interactive view
  output$pca_plot <- renderPlotly({
    norms <- normalized_count()$E
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
  
  # Create boxplot using plotly for an interactive view
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
