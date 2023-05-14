library(shiny)
library(DT)
library(ggplot2)
library(DESeq2)
library(plotly)
library(reshape2)
library(plotly)
library(pheatmap)
library(EnhancedVolcano)

function(input, output, session) {
  
  theme_set(theme_bw(base_size = 16)) # This scales up elements of all plots
  
  # Reading files for RNA and checking if samples match
  
  ## Reading count data and filtering
  count_rna <- reactive({
    req(input$raw_rna)
    read.csv(input$raw_rna$datapath, header = T, row.names = 1)
  })
  
  ## Showing count data to the user
  output$count_rna <- renderDataTable({
    datatable(count_rna(), rownames = T, options = list(pageLength = 10), caption = "RNAseq counts")
    })
  
  ## Checking if data is numeric
  
  output$count_rna_error <- renderText({ 
    if (any(!apply(count_rna(),2,is.numeric))) {
    "Error in RNA counts: not all numeric values. Upload only counts!"
    }
    else "Raw count data: OK"
  })
  ## Reading metadata data (factors)
  metadata_rna <- reactive({
    req(input$design_rna)
    read.csv(input$design_rna$datapath, header = T, row.names = 1)
  })
  
  ## Showing metadata data to the user
  output$metadata_rna <- renderDataTable({
    datatable(metadata_rna(), rownames = T, options = list(pageLength = 10), caption = "RNAseq metadata")
  })
  
  ## Checking if metadata and count data have same samples (number and name)
  
  output$metadata_rna_error <- renderText({
    req(metadata_rna())
    if (any(!colnames(count_rna()) %in% row.names(metadata_rna())) & # Checking if all samples are present
        ncol(count_rna()) == nrow(metadata_rna())) {
      "Error in metadata: samples do not match in count table and metadata. Make sure same samples are present and have the same label"
    }
    else "Metadata: OK"
  })
  
  # Preprocessing RNAseq (filtering, and preexploration with PCA/Pheatmap for outliers)
  
  ## Getting total counts per gene
  counts_by_gene <- reactive({
    data.frame(counts=rowSums(count_rna()),
                genes=row.names(count_rna())
                )
  })
  
  ## Filtering genes by min_count
  expressed_genes <- reactive({
    req(input$min_count)
    count_rna()[counts_by_gene()[,1] >= input$min_count,]
  })
  
  ## Outputing number of genes filtered
  output$genes_filtered <- renderText({
    paste(input$min_count,"counts threshold is filtering",
          nrow(count_rna())-nrow(expressed_genes()),"genes", sep = " ")
  })
  
  ## Plotting low expressed genes (less than 1000 counts) to improve visualization of the filtering
  output$low_expressed_plot <- renderPlot({
    ggplot(data=counts_by_gene()[counts_by_gene()[,1] <= 1000,], #This is just to improve visualization
           aes(x=counts)) +
      geom_histogram(bins=100, fill="coral1") +
      geom_vline(xintercept = input$min_count, linetype="dashed", color="blue") +
      labs(title="Number of genes with X counts (up to 1000)",
           x = "Counts",
           y = "Number of genes")
    
  })
  
  ## Creating data frame with library size per sample
  library_size <- reactive({
    data.frame(samples=colnames(count_rna()),
               library_size=colSums(count_rna())
    )
  })
  
  ## Filtering by library size (min_library)
  filtered_counts <- reactive({
    req(input$min_library)
    expressed_genes()[library_size()[,2] >= input$min_library*10^4,]
  })
  
  metadata_rna_filtered <- reactive({
    req(filtered_counts)
    metadata_rna()[colnames(filtered_counts()),] # This line is here in case a sample is filtered to be also removed from the metadata. 
  })
  
  ## Outputing number of samples filtered by library size
  output$samples_lowlib <- renderText({
    paste(input$min_library*10^4,"library size threshold is filtering",
          (nrow(count_rna())-nrow(filtered_counts())),
          "samples", sep = " ")
  })
  
  ## Barplot of library size
  output$library_size_plot <- renderPlot({
    ggplot(data=library_size(), aes(x=samples,y=library_size)) +  
      geom_bar(stat = "identity", fill="coral1") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title="RNAseq library size per sample",
           x = "Samples",
           y = "Total counts") +
      geom_hline(yintercept = input$min_library*10^4)
  })
  
  ## RNAseq PCA as outlier visualization.

  ### Creating a vector where samples to be removed will be saved
  ### Check for metadata_rna object. If modified, update the input options for PCA.
  
  ### Removing samples when submit button is clicked by updating vector (done in previous block)
  
  outliers <- reactiveValues(samples=c()) # Starts value
  
  observeEvent(input$rna_submit_outlier, {
    outliers$samples <- c(outliers$samples ,input$rna_remove_outlier) # Keeps track of the samples to remove 
  })
  
  observeEvent(input$recover_all_samples, { # Restarts in case user clciks on restart button
    outliers$samples <- c()
  })
  
  samples_to_keep <- reactive({ # This var keeps sample name to save 
    colnames(filtered_counts())[!(colnames(filtered_counts()) %in% outliers$samples)]
  })
  
  final_counts <- reactive({ 
    filtered_counts()[,samples_to_keep()]
  })
  
  final_metadata_rna <- reactive({
    metadata_rna_filtered()[samples_to_keep(),]
  })
  
  observe({
    req(final_metadata_rna())
    updateSelectInput(session, "outlier_pca_colour", choices = colnames(final_metadata_rna()))
    updateSelectInput(session, "outlier_pca_shape", choices = colnames(final_metadata_rna()))
    updateSelectInput(session, "rna_remove_outlier", choices = row.names(final_metadata_rna()))
    })

  ### Creating interactive plot with ggplot2 + plotly
  output$outlier_pca <- renderPlotly({
    req(final_counts())
    req(final_metadata_rna())
    req(input$outlier_pca_colour)
    req(input$outlier_pca_shape)
    vst_counts <- vst(round(as.matrix(final_counts()))) # VST is applied (data is not normalized, so some normalization is needed)
    preprocessing_pca <- prcomp(vst_counts)
    pca_frame <- data.frame("PC1" = preprocessing_pca$rotation[,1],
                            "PC2" = preprocessing_pca$rotation[,2]
                            )
    pca_frame <- cbind(pca_frame,
                       as.factor(final_metadata_rna()[,input$outlier_pca_colour]),
                       as.factor(final_metadata_rna()[,input$outlier_pca_shape]),
                       row.names(final_metadata_rna())
                       )
    if (input$outlier_pca_colour != input$outlier_pca_shape) { #This condition is to avoid the error if the same 2 variables are selected by the user.
      colnames(pca_frame)[c(3,4,5)] <- c(input$outlier_pca_colour,input$outlier_pca_shape,"Sample") 
    }
    else colnames(pca_frame)[c(3,4,5)] <- c(input$outlier_pca_colour,"same_variable","Sample") 
    

    outlier_pca_plot <- ggplot(data=pca_frame, 
                               aes_string(x="PC1",
                                          y="PC2",
                                          label = "Sample",
                                          color=as.character(input$outlier_pca_colour),
                                          shape=as.character(input$outlier_pca_shape))) +
      geom_point(size = 3) +
      labs(title="PCA with all genes")
    
    ggplotly(outlier_pca_plot, tooltip = "Sample") # transforms ggplot into a plotly object
  })
  
  
  ## Normalizing data using DESeq2
  
  ### Creating DESeq2 object to extract normalized data
  dds <- reactive({
   DESeq(DESeqDataSetFromMatrix(countData = round(final_counts()), 
                          colData = final_metadata_rna(), 
                          design = as.formula(paste("~", input$outlier_pca_colour)) #This is a tmp fix
   ))
  })
  
  ### Extracting normalized counts
  normalized_rna <- reactive({
    req(dds())
    counts(dds(), normalized=T)
  })
  
  ### Plot before normalization
  output$rna_before_norm <- renderPlot({
    boxplot_frame <- melt(log2(t(final_counts())+1)) #long format for ggplot
    ggplot(data=boxplot_frame, aes(x=Var1,y=value)) +
      geom_boxplot(fill="coral") +  
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title="Boxplots of log(counts+1) per sample before normalization",
           x = "Samples",
           y = "log (counts + 1)") +
      geom_hline(yintercept=median(as.matrix(log2(t(final_counts())+1))), linetype="dashed", color="blue")
  })
  
  output$rna_after_norm<- renderPlot({
    boxplot_frame <- melt(log2(t(normalized_rna())+1)) #long format for ggplot
    ggplot(data=boxplot_frame, aes(x=Var1,y=value)) +
      geom_boxplot(fill="coral") +  
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title="Boxplots of log(counts+1) per sample after normalization",
           x = "Samples",
           y = "log (counts + 1)") +
      geom_hline(yintercept=median(as.matrix(log2(t(normalized_rna())+1))), linetype="dashed", color="blue")
    
  })
  
  ### Downloading normalized data and metadata
  
  output$download_counts <- downloadHandler(
    filename = function() {
      paste("normalized_",input$raw_rna, sep="")
    },
    content = function(file) {
      write.csv(normalized_rna(), file, row.names=T)
    }
  )
  
  output$download_rna_metadata <- downloadHandler(
    filename = function() {
      paste("filtered_metadata_",input$design_rna, sep="")
    },
    content = function(file) {
      write.csv(final_metadata_rna(), file, row.names=T)
    }
  )
    
    ## Differential Expression Analysis plots
  
    ### Transforming gene id into gene symbol if user chooses so
  
  ### TRYING TO FIX THIS ####
  genes_symbols <- reactiveValues(symbols="")
  observeEvent(input$gene_symbol_submit, {
    req(filtered_counts())
    ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl") # dataset is the same as species is gonna be to be selected by the user.
    gene_symbol_frame <- getBM(attributes=c('ensembl_gene_id', #This line will also depend if this is entrez ID
                                            'external_gene_name'),
                               filters = 'ensembl_gene_id', #This line will also depend if this is entrez ID
                               values = rownames(filtered_counts()),
                               mart = ensembl)
    genes_id <- rownames(filtered_counts())
    genes_to_rename <- genes_id %in% gene_symbol_frame[,1]
    genes_id[genes_to_rename] <- genes_symbol_frame[genes_to_rename,2]
    genes_symbols$symbols <- genes_id
  })

  output$tmp <- renderText({print(genes_symbols$symbols)}) ### CHECKING OUTPUT 
    
    ###Updating options to let user chose comparison
  
  observe({
    req(final_metadata_rna())
    updateSelectInput(session, "de_factor", choices = colnames(final_metadata_rna()))
    updateSelectInput(session, "de_level1", choices = c("",colnames(final_metadata_rna())))
    updateSelectInput(session, "de_level2", choices = row.names(final_metadata_rna()))
  })
    
  de_results <- reactive({
    req(final_counts())
    deseq_object <- DESeq(DESeqDataSetFromMatrix(countData = round(final_counts()),
                                                 colData = final_metadata_rna(),
                                                 design = as.formula(paste("~", input$de_factor))
                                                 
      ))
    de_results <- results(deseq_object)
    if (genes_symbols$symbols!="") {
      de_results <- cbind(de_results, ID=genes_symbols$symbols)
    }
    else de_results <- cbind(de_results, ID=rownames(de_results))
    
    ### This may be added later if 2 factor analysis is performed
    # }
    # else {
    #   deseq_object <- DESeq(DESeqDataSetFromMatrix(countData = round(final_counts()),
    #                                                colData = final_metadata_rna(),
    #                                                design = as.formula(paste("~", input$de_factor1, "+", input$de_factor2))
    #   ))
    # }  
  })
  
    output$volcano_rna <- renderPlotly({
      volcano <- EnhancedVolcano(de_results(),
                                       lab = rownames(de_results()),
                                       x = "log2FoldChange",
                                       y = "pvalue",
                                 pCutoff = input$volcano_pvalue,
                                 FCcutoff = input$volcano_fc,
                                       drawConnectors = TRUE,
                                 pointSize=1)
    ggplotly(volcano + aes(x= log2FoldChange, y= -log10(pvalue), label=pvalue, text=ID), tooltip = c("log2FoldChange","pvalue","text"))
    })
  
}