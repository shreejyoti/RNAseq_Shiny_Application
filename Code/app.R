# Loading libraries to perform required set of task in R shiny
# Library DT allows for the creation of interactive and customizable data tables in R applications, especially in Shiny apps
# Library ggbeeswarm is required for creating beeswarm plots using ggplot2
# Library fgsea is used for Gene Set Enrichment Analysis (GSEA) on large datasets
# Library(msigdbr) prvides access to the Molecular Signatures Database (MSigDB), a resource of gene sets that is required in enrichment analyses like GSEA
library(readr)
library(GEOquery)
library(tibble)
library(dplyr)
library(shiny)
library(tidyr)
library(ggplot2)
library(DT)
library(ggbeeswarm)
library(tibble)
library(fgsea)
library(msigdbr)



# Loading required data file for analysis 
metadata <- read_csv("metadata.csv")
norm_counts <- read_csv("norm_counts.csv")
DESeq2_results <- read_csv("DESeq2_results.csv")
fgsea_results<- read_csv("fgsea_results.csv")



# Increase the maximum upload size to 70 MB
options(shiny.maxRequestSize = 70*1024^2)



#Comments for user layout  
#fluidPage(): Creates a web page layout for the app with a flexible UI design
#titlePanel(): Sets the title of the app displayed at the top of the page, so i have set "Final Project App"
#tabsetPanel(): Creates a tabbed interface to organize content into multiple tabs, in this case we have multiple tabs such as Sample, Counts etc
#tabPanel("Sample", ...): This defines a single tab labeled "Sample" to contain input and output elements related to metadata processing
#fluidRow(): Arranges UI components horizontally within a single row
#fileInput(): This function allows to upload input file and accept = c(".csv",".tsv") is used to restrict file input only to csv and tsv files
#selectInput(): This is used to create dropdown menu for selecting specific column or row

#UI section
ui <- fluidPage(
  titlePanel("Final Project App"),
  tabsetPanel(
    # Tab 1: Sample Tab
    tabPanel(
      "Sample",
      fluidRow(
        column(
          width = 3,
          wellPanel(
            fileInput("file", "Upload Metadata File", accept = c(".csv",".tsv")),
            textOutput("fileValidationMessageSample"), # Validation message
            actionButton("process", "Submit")
          )
        ),
        column(
          width = 9,
          tabsetPanel(
            tabPanel("Summary", tableOutput("summaryTable")),
            tabPanel("Data Table", DTOutput("dataTable")),
            tabPanel(
              "Plots",
              fluidRow(
                column(
                  width = 4,
                  selectInput("plotVar", "Select a Variable for Plotting", choices = NULL),
                  selectInput("groupVar", "Select a Grouping Variable", choices = c("Diagnosis", "vonsattel_grade")),
                  selectInput("plotType", "Select Plot Type", choices = c("Histogram", "Density Plot", "Violin Plot")),
                  actionButton("generatePlot", "Generate Plot")
                ),
                column(width = 8, plotOutput("variablePlot"))
              )
            )
          )
        )
      )
    ),
    # Tab 2: Counts Tab
    tabPanel(
      "Counts",
      fluidRow(
        column(
          width = 3,
          wellPanel(
            fileInput("countsFile", "Upload Normalized Counts File", accept = c(".csv", ".tsv")),
            textOutput("fileValidationMessageCounts"), # Validation message
            sliderInput("varPercentile", "Select percentile of variance", min = 0, max = 100, value = 50, step = 1),
            sliderInput("nonZeroSamples", "Select non-zero samples", min = 0, max = 100, value = 50, step = 1),
            actionButton("filterGenes", "Apply Filter")
          )
        ),
        column(
          width = 9,
          tabsetPanel(
            tabPanel("Summary", tableOutput("filterSummary")),
            tabPanel("Filtered Table", DTOutput("filteredCountsTable")),
            tabPanel(
              "Diagnostics",
              fluidRow(
                column(width = 6, plotOutput("medianVariancePlot")),
                column(width = 6, plotOutput("medianZerosPlot"))
              )
            ),
            tabPanel(
              "Clustered Heatmap",
              fluidRow(
                column(width = 3, checkboxInput("logTransform", "Log-transform counts for visualization", value = FALSE)),
                column(width = 9, plotOutput("clusteredHeatmap", height = "600px"))
              )
            ),
            tabPanel(
              "PCA Scatter Plot",
              fluidRow(
                column(
                  width = 3,
                  wellPanel(
                    radioButtons(
                      "pcaPlotType", "Select Plot Type:",
                      choices = c("Scatter Plot (PC1 vs PC2)" = "scatter", "Beeswarm Plot (Top N PCs)" = "beeswarm"),
                      selected = "scatter"
                    ),
                    conditionalPanel(
                      condition = "input.pcaPlotType == 'scatter'",
                      selectInput("xPC", "X-axis Principal Component:", choices = NULL),
                      selectInput("yPC", "Y-axis Principal Component:", choices = NULL)
                    ),
                    conditionalPanel(
                      condition = "input.pcaPlotType == 'beeswarm'",
                      numericInput("topN", "Number of Top PCs to Plot:", value = 5, min = 2, step = 1)
                    ),
                    actionButton("runPCA", "Generate PCA Plot")
                  )
                ),
                column(width = 9, plotOutput("pcaPlot", height = "600px"))
              )
            )
          )
        )
      )
    ),
    # Tab 3: Differential Expression (DE) Tab
    tabPanel(
      "DE",  
      fluidRow(
        column(
          width = 3,
          wellPanel(
            fileInput("deFile", "Upload Differential Expression Results File", accept = c(".csv", ".tsv")),
            textOutput("fileValidationMessageDE"), # Validation message
            actionButton("submitDE", "Submit") # Submit button
          )
        ),
        column(
          width = 9,
          tabsetPanel(
            tabPanel(
              "DGEA Table", # Tab for displaying the table
              DTOutput("deTable") # Table Output for Differential Expression
            ),
            tabPanel(
              "Volcano Plot", # New tab for Volcano Plot
              plotOutput("volcanoPlot", height = "600px") # Output for the Volcano Plot
            ),
            tabPanel(
              "Histogram p-value", # New Tab for Histogram
              plotOutput("histogramPlot", height = "600px") # Output for the Histogram
            ),
            tabPanel(
              "Histogram Log2FC", # New Tab for Histogram of Log2FoldChange
              fluidRow(
                column(
                  width = 3,
                  wellPanel(
                    sliderInput(
                      "padjThreshold", 
                      "Adjusted p-value (padj) Threshold:", 
                      min = 0, 
                      max = 1, 
                      value = 0.10, 
                      step = 0.01
                    ) # Slider to dynamically adjust padj threshold
                  )
                ),
                column(
                  width = 9,
                  plotOutput("log2fcHistogram", height = "600px") # Output for the Log2FoldChange Histogram
                )
              )
            ),
            tabPanel(
              "Jitter Plot", # New Tab for Jitter Plot
              plotOutput("jitterPlot", height = "600px") # Output for the Jitter Plot
            )
          )
        )
      )
    ),
    # Tab 4: GSEA Tab
    tabPanel(
      "GSEA",
      fluidRow(
        column(
          width = 3,
          wellPanel(
            # File input for GSEA
            fileInput("fgseaFile", "Upload fgsea file", accept = c(".csv",".tsv")),
            textOutput("fileValidationMessageGSEA"), # Validation message
            # Sidebar content will change based on the selected subtab
            conditionalPanel(
              condition = "input.gseaSubTab === 'fgsea Barplot'",
              sliderInput("topPathways", "Select Number of Top Pathways:", 
                          min = 5, max = 50, value = 10)
            ),
            conditionalPanel(
              condition = "input.gseaSubTab === 'fgsea Table'",
              sliderInput("pvalThreshold", "Adjusted p-value Threshold:", 
                          min = 0, max = 1, value = 0.05, step = 0.01),
              radioButtons("nesDirection", "Select NES Direction:",
                           choices = c("All" = "All", "Positive" = "Positive", "Negative" = "Negative"),
                           selected = "All"),
              downloadButton("downloadFgsea", "Download Filtered Results")
            ),
            conditionalPanel(
              condition = "input.gseaSubTab === 'Scatter Plot'",
              sliderInput("scatterPvalThreshold", "Adjusted p-value Threshold:", 
                          min = 0, max = 1, value = 0.05, step = 0.01)
            )
          )
        ),
        column(
          width = 9,
          tabsetPanel(
            id = "gseaSubTab",
            tabPanel(
              "fgsea Barplot", # Subtab for GSEA Barplot
              plotOutput("gseaBarplot", height = "600px")
            ),
            tabPanel(
              "fgsea Table", # Subtab for GSEA Table
              DTOutput("fgseaTable")
            ),
            tabPanel(
              "Scatter Plot", # Subtab for Scatter Plot
              plotOutput("scatterPlot", height = "600px")
            )
          )
        )
      )
    )
  )
)


# Server interface  section 
# Server function defines the logic of the app and how the app behaves
server <- function(input, output, session) {
  # Validation logic for file uploads
  # validate_file function validates whether the uploaded file has an extension .csv or .tsv. If file input is invalid, give an error message 
  validate_file <- function(file) {
    req(file)
    ext <- tools::file_ext(file$datapath)
    if (!ext %in% c("csv", "tsv")) {
      stop("Invalid file format! Please upload a valid CSV or TSV file.")
    }
    return(TRUE)
  }
  
  # Tab 1: Sample
  # Reactive to read and process the metadata file
  # This reactive function validates the uploaded file using validate_file and reads it as a CSV into a data frame
  metadata <- reactive({
    validate_file(input$file) # Validate the file
    read.csv(input$file$datapath, stringsAsFactors = FALSE)
  })
  
  # Reactive to read and process the counts file
  counts <- reactive({
    validate_file(input$countsFile) # Validate the file
    read.csv(input$countsFile$datapath, row.names = 1, stringsAsFactors = FALSE)
  })
  
  # Update variable choices for plotting
  # This function observes changes in the metadata reactive and updates the plotVar select input with the names of numeric columns from the uploaded metadata file.
  observeEvent(metadata(), {
    numeric_vars <- names(metadata())[sapply(metadata(), is.numeric)]
    updateSelectInput(session, "plotVar", choices = numeric_vars)
  })
  
  # Reactive to create the summary table
  # Here it creates a reactive summary table by extracting and processing metadata from the uploaded file 
  summary_table <- reactive({
    req(metadata()) # Ensure metadata is loaded
    data <- metadata()
    
    # Here it iterates over each column of the dataframe to calculate summary statistics or identify distinct values based on the column type
    summary <- lapply(colnames(data), function(colname) {
      col <- data[[colname]]
      # Checks if the column is numeric 
      if (is.numeric(col)) {
        # Calculates the mean of the column values, ignoring NA values and rounds it to 2 decimal places
        # paste0 (): Combines the mean and standard deviation into a single string in the format Mean (+/- SD)
        mean_sd <- paste0(
          round(mean(col, na.rm = TRUE), 2), 
          " (+/-", 
          # Calculates the standard deviation of the column values
          round(sd(col, na.rm = TRUE), 2), 
          ")"
        )
        # Checks if the column is of type character or factor
        mean_sd
      } else if (is.character(col) || is.factor(col)) {
        # Identifies all unique values in the column using unique(col) and combines them into a single string, separated by commas
        paste(unique(col), collapse = ", ")
        # For any column type that is neither numeric nor character/factor, assigns "N/A" as the summary value
      } else {
        "N/A"
      }
    })
    
    # This code section creates a data.frame summarizing the structure and descriptive statistics of a dataset data
    data.frame(
       # This retrieves all the column names from the data object and includes them in the data frame as a list under this column
      `Column Name` = colnames(data),
       # This applies the class function to each column of the dataset data and determines its data type 
      `Type` = sapply(data, class),
      `Mean(sd) or Distinct Values` = unlist(summary),
      stringsAsFactors = FALSE
    )
  })
  
  # Subcomponent 1: Output of the summary table for Sample tab
  # This function is the output and renders the summary_table as a table in the UI
  output$summaryTable <- renderTable({
    req(input$process)
    summary_table()
  })
  
  # Subcomponent 2: Output of the data table for Sample tab
  # renderDT function is used in Shiny applications to render a DataTable (interactive table) on the user interface (UI)
  output$dataTable <- renderDT({
    req(input$process)
    datatable(metadata(), options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE))
  })
  
  # Code for selected variable and plot type
  plot_data <- reactive({
    req(metadata(), input$plotVar, input$groupVar)
    metadata() %>%
      # using filter function to remove NA values
      filter(!is.na(.data[[input$plotVar]])) %>%
      # select function to select only the grouping variable and the plotting variable for use in creating plots
      select(.data[[input$groupVar]], .data[[input$plotVar]])
  })
  
  # Subcomponent 3: Generating and displaying the selected plot
  output$variablePlot <- renderPlot({
    req(plot_data(), input$plotType)
    data <- plot_data()
    variable <- input$plotVar
    group <- input$groupVar
    
    # generates a plot based on the selected plot type (Histogram, Density Plot, or Violin Plot)
    # This line of code checks if the user has selected "Histogram" as the plot type in the app, if yes then it uses geom_histogram to creates the histogram
    if (input$plotType == "Histogram") {
      ggplot(data, aes(x = .data[[variable]], fill = as.factor(.data[[group]]))) +
        geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
        labs(x = variable, y = "Count", fill = group) +
        theme_minimal()
    } else if (input$plotType == "Density Plot") {
      ggplot(data, aes(x = .data[[variable]], fill = as.factor(.data[[group]]))) +
        geom_density(alpha = 0.7) +
        labs(x = variable, y = "Density", fill = group) +
        theme_minimal()
    } else if (input$plotType == "Violin Plot") {
      ggplot(data, aes(x = as.factor(.data[[group]]), y = .data[[variable]], fill = as.factor(.data[[group]]))) +
        geom_violin(alpha = 0.7) +
        labs(x = group, y = variable, fill = group) +
        theme_minimal()
    }
  })
  
  
  
  # Tab 2: Counts 
  # Function to filter genes
  filtered_counts <- reactive({
    req(counts(), input$varPercentile, input$nonZeroSamples)
    data <- counts()
    # Calculates the variance for each gene (row-wise)
    gene_variance <- apply(data, 1, var, na.rm = TRUE)
    # Calculates the number of non-zero samples for each gene
    non_zero_counts <- rowSums(data > 0)
    # Computes the variance threshold based on the user-defined percentile using quantile()
    var_threshold <- quantile(gene_variance, probs = input$varPercentile / 100, na.rm = TRUE)
    # Filters genes that have variance greater than or equal to the threshold
    # Filters genes that have at least X non-zero samples
    filtered_data <- data[gene_variance >= var_threshold & non_zero_counts >= input$nonZeroSamples, ]
    # Ensures that genes with all zero values after filtering are removed.
    filtered_data <- filtered_data[rowSums(filtered_data, na.rm = TRUE) > 0, ]
    
    # This returns a list containing, the filtered dataset based on the sliders, original dataset for reference, calculated variance for each gene and number of non-zero samples for each gene
    list(
      filtered = filtered_data,
      original = data,
      variance = gene_variance,
      nonZero = non_zero_counts
    )
  })
  
  # Output of the filtered data: filter summary
  output$filterSummary <- renderTable({
    req(filtered_counts())
    filtered <- filtered_counts()
    # This calculates the number of columns in the original data, which corresponds to the number of sample
    total_genes <- nrow(filtered$original)
    # The number of genes passing is the number of rows in the filtered dataset
    passing_genes <- nrow(filtered$filtered)
    failing_genes <- total_genes - passing_genes
    
    # Creating a dataframe to assign all calculated values for filter summary table
    data.frame(
      Metric = c("Total Samples", "Total Genes", "Genes Passing Filter", "Genes Failing Filter"),
      Count = c(ncol(filtered$original), total_genes, passing_genes, failing_genes),
      # the percentage is calculated by dividing the passing genes by the total genes and multiplying by 100 and same for failing genes
      Percentage = c(NA, 100, round((passing_genes / total_genes) * 100, 2), round((failing_genes / total_genes) * 100, 2))
    )
  })
  
  # Subcomponent 1: Output of the filtered counts table
  # Render a data table (using DT) for filtered counts after applying the gene filtering criteria
  output$filteredCountsTable <- renderDT({
    req(filtered_counts())
    datatable(filtered_counts()$filtered, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Subcomponent 2: Diagnostic Scatter Plot: Median Count vs Variance (Log Scale)
  # This code initializes the Median Count vs Variance (Log Scale) plot by ensuring the filtered counts data is available and processed
  output$medianVariancePlot <- renderPlot({
    req(filtered_counts())
    filtered <- filtered_counts()
    
    # Initiate the creation of a data.frame named plot_data that will store the data for scatter plot
    plot_data <- data.frame(
      # Calculates the median value for each row (gene) in the filtered$original dataset
      # apply(..., 1, ...) this function is used to apply the median function across rows (indicated by 1)
      Median = apply(filtered$original, 1, median, na.rm = TRUE),
      # Calculates the variance for each row (gene) in the filtered$original dataset
      Variance = apply(filtered$original, 1, var, na.rm = TRUE),
      # Create a logical vector which checks whether each row (gene) in the original dataset is present in the filtered dataset.
      Passed = rownames(filtered$original) %in% rownames(filtered$filtered)
    )
    
    # Function ggplot(), combined with geom_point(), scale_*, labs, and theme for scatter plotting
    ggplot(plot_data, aes(x = log10(Median + 1), y = log10(Variance + 1), color = Passed)) +
      geom_point(alpha = 0.7, size = 2) +
      # Adjusts the x-axis and y-axis scale to be logarithmic (log10)
      scale_x_log10(labels = scales::scientific, name = "Median Count (Log Scale)") +  
      scale_y_log10(labels = scales::scientific, name = "Variance (Log Scale)") +      
      scale_color_manual(values = c("FALSE" = "lightgray", "TRUE" = "darkblue")) +
      theme_minimal() +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank()  # Simplify grid for clarity
      )
  })
  
  # Subcomponent 2: Diagnostic Scatter Plot: Median Count vs Number of Zeros
  output$medianZerosPlot <- renderPlot({
    req(filtered_counts())
    filtered <- filtered_counts()
    
    # Logic is same as above explained, only difference is instead of having variance here we calculates the number of samples (columns) where the count is zero for each gene (row) 
    plot_data <- data.frame(
      Median = apply(filtered$original, 1, median, na.rm = TRUE),
      # ncol(filtered$original) - rowSums(...): Subtracts the number of non-zero counts from the total number of samples to get the number of zero counts for each gene
      Zeros = ncol(filtered$original) - rowSums(filtered$original > 0),
      Passed = rownames(filtered$original) %in% rownames(filtered$filtered)
    )
    
    # Function ggplot(), combined with geom_point(), scale_*, labs, and theme for scatter plotting
    ggplot(plot_data, aes(x = log10(Median + 1), y = Zeros, color = Passed)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "lightgray", "TRUE" = "darkblue")) +
      labs(title = "Median Count vs Number of Zeros",
           x = "Median Count",
           y = "Number of Zeros",
           color = "Passed Filter") +
      theme_minimal()
  })
  
  
  # Subcomponent 3: Generate Clustered Heatmap
  output$clusteredHeatmap <- renderPlot({
    req(filtered_counts())  # Ensure filtered counts are available
    filtered <- filtered_counts()$filtered
    
    # Optionally log-transform the data for visualization
    if (input$logTransform) {
      heatmap_data <- log10(filtered + 1)
    } else {
      heatmap_data <- filtered
    }
    
    # Create the heatmap with clustering
    # Generates a clustered heatmap using the pheatmap function from the pheatmap package
    pheatmap::pheatmap(
      mat = heatmap_data,
      # TRUE to cluster the rows and columns hierarchically
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      scale = "row",  
      # Color palette
      color = colorRampPalette(c("navy", "white", "firebrick"))(50),  
      main = "Clustered Heatmap of Filtered Counts",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize = 10,  
      fontsize_row = 8,  
      fontsize_col = 8)
  })
  
  # Subcomponent 4: PCA Plot
  # PCA: Reactive to perform PCA
  pca_result <- reactive({
    req(filtered_counts())
    data <- filtered_counts()$filtered
    data <- log2( data + 1) 
    pca <- prcomp(t(data), scale. = TRUE)
    
    # Calculate variance explained
    variance_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100
    list(
      pca = pca,
      variance_explained = variance_explained
    )
  })
  
  # Update PCA component choices based on results
  observeEvent(pca_result(), {
    pcs <- paste0("PC", seq_along(pca_result()$variance_explained))
    updateSelectInput(session, "xPC", choices = pcs, selected = pcs[1])
    updateSelectInput(session, "yPC", choices = pcs, selected = pcs[2])
  })
  
  # Render PCA plot
  output$pcaPlot <- renderPlot({
    req(pca_result(), input$runPCA)
    pca <- pca_result()$pca
    variance_explained <- pca_result()$variance_explained
    pc_data <- as.data.frame(pca$x)  # PCA projections
    pc_data$Sample <- rownames(pc_data)
    
    # Scatter plot for selected PCs
    if (input$pcaPlotType == "scatter") {
      req(input$xPC, input$yPC)
      pc_x_num <- as.numeric(sub("PC", "", input$xPC))
      pc_y_num <- as.numeric(sub("PC", "", input$yPC))
      
      pca_data <- data.frame(
        PCX = pca$x[, pc_x_num],
        PCY = pca$x[, pc_y_num],
        Sample = rownames(pca$x)
      )
      
      ggplot(pca_data, aes(x = PCX, y = PCY)) +
        geom_point(size = 3, alpha = 0.7, color = "#1E90FF") +
        theme_minimal() +
        labs(
          title = "PCA Scatter Plot",
          x = sprintf("%s (%.2f%% Variance)", input$xPC, variance_explained[pc_x_num]),
          y = sprintf("%s (%.2f%% Variance)", input$yPC, variance_explained[pc_y_num])
        ) +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.line = element_line(color = "black")
        )
    }
    
    # Beeswarm plot for top PCs
    else if (input$pcaPlotType == "beeswarm") {
      req(input$topN)
      top_n <- min(input$topN, ncol(pc_data))  # Ensure topN does not exceed available PCs
      top_pcs <- pc_data[, 1:top_n]
      top_pcs$Sample <- rownames(pc_data)
      
      melted_pcs <- reshape2::melt(top_pcs, id.vars = "Sample")
      melted_pcs$variable <- factor(melted_pcs$variable, levels = paste0("PC", 1:top_n))  # Maintain order
      
      ggplot(melted_pcs, aes(x = variable, y = value, color = variable)) +
        geom_quasirandom(alpha = 0.7, size = 2) +
        theme_minimal() +
        scale_color_viridis_d() +
        labs(
          title = "Beeswarm Plot of Top Principal Components",
          x = "Principal Component",
          y = "Projection Value",
          color = "PC"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })
  
  
  
  # Tab 3: Differential Gene Expression Analysis 
  # Code to read and process the DE file
  de_data <- reactive({
    # validate_file function ensures that the uploaded file (input$deFile) is either a .csv or .tsv file
    validate_file(input$deFile) # Validate the file
    read.csv(input$deFile$datapath, stringsAsFactors = FALSE)
  })
  
  
  # Output code for the table in DGE Results tab
  # It assigns a reactive data table (DT) to the output object deTable
  output$deTable <- renderDT({
    req(de_data())
    # datatable function, by the DT package, creates an interactive and customizable data table from an R data frame
    datatable(
      # Display the data exactly as it is in the input file
      de_data(), 
      options = list(
        pageLength = 10, # Show 10 rows per page
        scrollX = TRUE,  # Allow horizontal scrolling
        autoWidth = TRUE # Adjust column width automatically
      )
    )
  })
  
  # Generate the Volcano Plot
  # renderPlot function ensures that the volcano plot dynamically updates whenever the reactive data de_data() changes or updates
  output$volcanoPlot <- renderPlot({
    # This assigns the plot to the output object volcanoPlot, which will display the plot in the corresponding UI element
    req(de_data())
    data <- de_data()
    # Classifying genes as "Significant" if their adjusted p-value (padj) is less than 0.05 and the absolute value of their log2 fold change (log2FoldChange) is greater than 1; otherwise, they are classified as "Not Significant."
    data$Significance <- ifelse(data$padj < 0.05 & abs(data$log2FoldChange) > 1, "Significant", "Not Significant")
    # Initializes the plot using data (reactive data from de_data()) using ggplot function
    ggplot(data, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
      geom_point(alpha = 0.8, size = 2) +
      scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
      labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(P-value)") +
      theme_minimal()
  })
  
  # Generate the Histogram of p-values
  output$histogramPlot <- renderPlot({
    req(de_data())
    data <- de_data()
    # mapping the pvalue to the x-axis
    ggplot(data, aes(x = pvalue)) +
      geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
      labs(title = "Histogram of Raw P-values", x = "p-value", y = "count") +
      theme_minimal()
  })
  
  # Generate the Histogram of Log2FoldChange
  output$log2fcHistogram <- renderPlot({
    req(de_data(), input$padjThreshold)  # Require data and the padj threshold slider
    data <- de_data()
    
    # Dynamically filter genes based on the padj threshold from the slider
    filtered_genes <- data %>% filter(padj < input$padjThreshold)
    
    # Plot using ggplot() and geom_histogram()
    ggplot(filtered_genes, aes(x = log2FoldChange)) +
      geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black", alpha = 0.7) +
      labs(title = paste0("Histogram of Log2FoldChange for Significant Genes (padj < ", input$padjThreshold, ")"), 
           x = "Log2 Fold Change", y = "Count") +
      theme_minimal()
  })
  
  # Generate the Jitter Plot for Top 10 Genes
  output$jitterPlot <- renderPlot({
    req(de_data())  # Ensure the DESeq2 results file is loaded
    de_data <- de_data()  # Fetch the DESeq2 results
    
    # Extract relevant columns and filter to the top 10 genes by smallest padj
    selected_data <- de_data %>%
      select(Gene_ID, Control.mean, HD.mean, padj) %>%  # Select relevant columns
      filter(!is.na(Control.mean), !is.na(HD.mean), !is.na(Gene_ID)) %>%  # Ensure no missing values
      arrange(padj) %>%  # Sort by padj in ascending order
      head(10)  # Select top 10 rows
    
    
    # Prepare the data for plotting
    plotData <- data.frame(
      Gene = rep(selected_data$Gene_ID, each = 2),  # Replicate each Gene ID twice for Control and HD
      SampleType = rep(c("Control", "HD"), times = nrow(selected_data)),  # Define SampleType
      NormCount = c(selected_data$Control.mean, selected_data$HD.mean)  # Combine normalized counts
    )
    
    # Create the Jitter Plot using ggplot2
    ggplot(plotData, aes(x = Gene, y = log10(NormCount + 1), color = SampleType)) +
      geom_point(size = 2, position = position_jitter(width = 0.2, height = 0)) +  # Add jitter for better visualization
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for clarity
      labs(
        title = "Jitter Plot of Log10(Normalized Counts) for Top 10 Genes",
        x = "Gene ID",
        y = "Log10(Normalized Counts)",
        color = "Sample Type"
      )
  }, height = 500, width = 900)
  
  
  # Tab 4: GSEA
  # Reactive to read and process the fgsea file
  fgsea_data <- reactive({
    validate_file(input$fgseaFile) # Validate the file
    read.csv(input$fgseaFile$datapath, stringsAsFactors = FALSE)
  })
  
  # Generate the GSEA Barplot
  output$gseaBarplot <- renderPlot({
    req(fgsea_data())  # Ensure data is available
    data <- fgsea_data()
    
    # Select top pathways based on adjusted p-value
    top_pathways <- data %>%
      arrange(padj) %>%
      # Select top N pathways based on slider input
      slice(1:input$topPathways)  
    
    # Create the barplot
    ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(
        title = "Top Pathways by NES",
        x = "Pathway",
        y = "Adjusted p_value"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  # Generate the Table Entry for Selected Pathway
  output$pathwayDetails <- renderTable({
    req(input$barplotClick, fgsea_data())  # Ensure click and data are available
    data <- fgsea_data()  # Load the fgsea data
    
    # Get the index of the clicked bar
    clicked_index <- round(input$barplotClick$y)
    
    # Check if the clicked index is within bounds
    if (clicked_index > 0 && clicked_index <= nrow(data)) {
      # Get the pathway corresponding to the clicked bar
      clicked_pathway <- data %>%
        arrange(padj) %>%
        slice(clicked_index) %>%
        select(pathway) %>%
        pull()
      
      # Filter the data for the clicked pathway
      data %>%
        filter(pathway == clicked_pathway)
    } else {
      NULL  # Return NULL if the click is outside bounds
    }
  })
  
  # Generate the filtered fgsea data table for Tab 2
  filtered_fgsea_data <- reactive({
    req(fgsea_data())
    data <- fgsea_data()
    
    # Filter by adjusted p-value
    filtered_data <- data %>%
      filter(padj <= input$pvalThreshold)
    
    # Filter by NES direction
    if (input$nesDirection == "Positive") {
      filtered_data <- filtered_data %>% filter(NES > 0)
    } else if (input$nesDirection == "Negative") {
      filtered_data <- filtered_data %>% filter(NES < 0)
    }
    
    filtered_data
  })
  
  # Generate the sortable data table with compact rows
  output$fgseaTable <- renderDT({
    req(filtered_fgsea_data())
    datatable(
      filtered_fgsea_data(),
      options = list(
        pageLength = 10,               # Show 10 rows per page
        scrollX = TRUE,                # Allow horizontal scrolling
        autoWidth = TRUE,              # Automatically adjust column widths
        dom = 'lrtip',                 # Simplified table controls
        class = 'display compact'      # Compact table styling
      ),
      rownames = FALSE                # Hide row names
    ) %>%
      formatStyle(
        columns = colnames(filtered_fgsea_data()), 
        `white-space` = "nowrap"      # Prevent long text wrapping in cells
      )
  })
  
  # Download handler for filtered data
  output$downloadFgsea <- downloadHandler(
    filename = function() {
      paste("filtered_fgsea_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_fgsea_data(), file, row.names = FALSE)
    }
  )
  
  # Generate the Scatter Plot
  output$scatterPlot <- renderPlot({
    req(fgsea_data())  # Ensure data is available
    data <- fgsea_data() # Get the uploaded data
    
    # Filter data based on the adjusted p-value threshold from the slider
    filtered_data <- data %>%
      mutate(Color = ifelse(padj <= input$scatterPvalThreshold, "Above threshold", "Below threshold"))
    
    # 
    ggplot(filtered_data, aes(x = NES, y = -log10(padj), color = Color)) +
      geom_point(alpha = 0.7, size = 3) +
      scale_color_manual(values = c("Above threshold" = "blue", "Below threshold" = "grey")) +
      labs(
        title = "Scatter Plot of NES vs -log10 Adjusted p-value",
        x = "Normalized Enrichment Score (NES)",
        y = "-log10 Adjusted p-value",
        color = "Gene Sets"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  })
}


# Run the app
shinyApp(ui, server)

