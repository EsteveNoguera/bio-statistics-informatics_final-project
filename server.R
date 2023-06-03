library(shiny)
library(DT)
library(ggplot2)
library(DESeq2)
library(plotly)
library(reshape2)
library(plotly)
library(pheatmap)
library(EnhancedVolcano)
library(biomaRt)
library(mixOmics)

function(input, output, session) {
  
  theme_set(theme_bw(base_size = 16)) # This scales up elements of all plots
  
  # Reading files for RNA and checking if samples match
  
  ## Reading count data and filtering
  
  count_rna <- reactive({
    req(input$raw_rna)
    count_rna <- read.csv(input$raw_rna$datapath, header = T, row.names = 1)
    count_rna <- as.data.frame(count_rna)
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
    metadata_rna <- read.csv(input$design_rna$datapath, header = T, row.names = 1)
    metadata_rna <- apply(metadata_rna,2,as.factor)
    metadata_rna <- as.data.frame(metadata_rna)
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
  
  
  # Normalizing data using DESeq2
  
  ## Creating DESeq2 object to extract normalized data
  
  dds <- reactive({
   DESeq(DESeqDataSetFromMatrix(countData = round(final_counts()), 
                          colData = final_metadata_rna(), 
                          design = as.formula(paste("~", input$outlier_pca_colour)) #This is a tmp fix
   ))
  })
  
  ## Extracting normalized counts
  
  normalized_rna <- reactive({
    req(dds())
    counts(dds(), normalized=T)
  })
  
  ## Plot before normalization
  
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
  
  ## Downloading normalized data and metadata
  
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
    
  # Differential Expression Analysis plots
  
  ## Translate ID into gene symbol if user chooses so
  
  gene_symbols <- reactiveValues(symbols="")
  
  translate_symbol <- function(gene_id) {
    ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    gene_translation <- getBM(
      attributes = c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = gene_id,
      mart = ensembl
    )
    return(gene_translation)
  }
  
  observeEvent(input$gene_symbol_submit, {
    req(final_counts())
    genes_id <- row.names(final_counts())
    gene_symbol_frame <- translate_symbol(genes_id)
    gene_symbol_frame[row.names(filtered_coun)]
    genes_to_rename <- genes_id %in% gene_symbol_frame[,1]
    genes_id[genes_to_rename] <- gene_symbol_frame[genes_to_rename,2]
    gene_symbols$symbols <- genes_id
  })
  
  ## Updating options to let user chose comparison
  
  observe({
    req(final_metadata_rna())
    updateSelectInput(session, "de_factor", choices = colnames(final_metadata_rna()))
  })

  observeEvent(input$de_factor, {
    req(final_metadata_rna())
    factor_selected <- input$de_factor
    if (!is.null(factor_selected) && factor_selected %in% names(final_metadata_rna())) {
      factor_variable <- as.character(final_metadata_rna()[[factor_selected]])
      levels_data <- unique(factor_variable)
      updateSelectInput(session, "de_level1", choices = levels_data)
      updateSelectInput(session, "de_level2", choices = levels_data)
    } else {
      updateSelectInput(session, "de_level1", choices = NULL)
      updateSelectInput(session, "de_level2", choices = NULL)
    }
  })
  
  deseq_objects <- reactiveValues(data="",results="")  
  observeEvent(input$de_submit, {
    req(final_counts())
    deseq_objects$data <- DESeq(DESeqDataSetFromMatrix(countData = round(final_counts()),
                                                 colData = final_metadata_rna(),
                                                 design = as.formula(paste("~", input$de_factor))                                 
      ))
    de_results <- results(deseq_objects$data, contrast=c(input$de_factor, input$de_level1, input$de_level2))
    if (length(gene_symbols$symbols) == 1) {
      de_results <- cbind(de_results,ID=rownames(de_results))
    }
    else {
      de_results <- cbind(de_results,ID=gene_symbols$symbols)
    }
    
    deseq_objects$results <- de_results
    
    ## Volcano plot
    
    output$volcano_rna <- renderPlotly({
      volcano <- EnhancedVolcano(de_results,
                                 lab = rownames(de_results),
                                 x = "log2FoldChange",
                                 y = "pvalue",
                                 pCutoff = input$volcano_pvalue,
                                 FCcutoff = input$volcano_fc,
                                 drawConnectors = TRUE,
                                 pointSize=1)
      ggplotly(volcano + aes(x= log2FoldChange, y= -log10(pvalue), label=pvalue, text=ID), tooltip = c("text","log2FoldChange","pvalue"))
    })
  })
  
  ## Download results 
  
  output$download_de <- downloadHandler(
    filename = function() {
      paste("de_",input$de_level1,"_vs_",input$de_level2, sep="")
    },
    content = function(file) {
      write.csv(deseq_objects$results, file, row.names=T)
    }
  )
  
  # Biomarker identification
  
  ## Selecting normalized counts for only those DE genes
  
  de_counts <- reactive({
    if (length(deseq_objects$results) > 1) {
      de_genes <- rownames(deseq_objects$results)[deseq_objects$results$padj <= 0.05]
      de_counts <- normalized_rna()[de_genes,]
      }
    })
  
  ## sPCA
  
  spca_model <- reactive({
    req(de_counts())
    spca_model <- spca(t(de_counts()), keepX = c(input$grid1,input$grid2))
    })
  
  output$rna_spca <- renderPlot({
    req(spca_model())
    plotIndiv(spca_model(), group=final_metadata_rna()[,input$de_factor], legend=T)
    })
  
  output$selected_spc1 <- renderDataTable({
    req(spca_model())
    selectVar(spca_model(), comp = 1)$value
    })
  
  output$selected_spc2 <- renderDataTable({
    req(spca_model())
    selectVar(spca_model(), comp = 2)$value
    })
    
    ## sPLS-DA
    
    rna_splsda <- reactiveValues(model="", performance="", tuned="")
    
    observeEvent(input$splsda_submit_rna, {
      withProgress(message = 'Building sPLS-DA model', value = 0, {
        incProgress(0.1, detail = "Creating initial model")
        splsda_model <- splsda(t(de_counts()), 
                               as.factor(final_metadata_rna()[,input$de_factor]),
                               ncomp = 10)
        
        nfold <- trunc(ncol(de_counts())/5)
        if (nfold > 10) nfold <- 10
        incProgress(0.4, detail = "Determining number of components to use for the model")
        perf_splsda_model <- perf(splsda_model, validation = "Mfold", 
                                  folds = nfold, nrepeat = 50, # function will only use this line if there Mfold is used
                                  progressBar = FALSE, auc = TRUE) 
        selected_ncomp <- perf_splsda_model$choice.ncomp[2,] # Selecting recommended ncomp by BER
        error_rate <- perf_splsda_model$error.rate.class # Getting error rates for each metric
        error_rate <- sapply(1:3,function(x,ncomp,error) error[[x]][,ncomp[x]], ncomp=selected_ncomp, error=error_rate) # Getting error rate according to ncomp
        best_metric <- which.min(apply(error_rate,2,mean)) # Which metric has the least mean BER across classes
        list.keepX <- seq(input$min_splsda_features_rna,
                          input$max_splsda_features_rna,
                          input$step_tuning_model)
        incProgress(0.4, detail = "Determining number of features to use for the model")
        tune_splsda_model <- tune.splsda(X=t(de_counts()), #data
                                         Y=final_metadata_rna()[,input$de_factor], #response vector
                                         ncomp = selected_ncomp[best_metric], #ncomp found by best_metric
                                         dist = names(selected_ncomp[best_metric]), #metric found by best_metric
                                         test.keepX = list.keepX,
                                         folds=nfold, # determined before to have at least 5 test samples
                                         nrepeat=50) # more repeats, more computing time. 
        rna_splsda$tuned <- tune_splsda_model
        optimal_ncomp <- tune_splsda_model$choice.ncomp$ncomp #final number of components
        if (optimal_ncomp==1) optimal_ncomp <- 2 # if only one component is found, force to use 2 to have plots
        optimal_keepX <- tune_splsda_model$choice.keepX[1:optimal_ncomp] #final number of features for each component
        incProgress(0.2, detail = "Building final model")
        rna_splsda$model <- splsda(X=t(de_counts()), #data
                                   Y=as.factor(final_metadata_rna()[,input$de_factor]),
                                   ncomp=optimal_ncomp,
                                   keepX = optimal_keepX)
      })
    })
    
    observe({
      if (length(rna_splsda$model) != 1) {
        if (rna_splsda$tuned$choice.ncomp$ncomp > 2) {
          updateSelectInput(session, "splsda_ncomp_rna", choices=c(1:rna_splsda$tuned$choice.ncomp$ncomp))
          updateSelectInput(session, "x_plsda_rna", choices=c(1:rna_splsda$tuned$choice.ncomp$ncomp))
          updateSelectInput(session, "y_plsda_rna", choices=c(1:rna_splsda$tuned$choice.ncomp$ncomp))
        }
      }
    })
    
    output$rna_splsda_plot <- renderPlot({
      if (length(rna_splsda$model) != 1) {
        plotIndiv(rna_splsda$model, legend=T, ellipse=T, title="sPLS-DA",
                  comp= as.numeric(c(input$x_plsda_rna,input$y_plsda_rna)))
        }
      })
      
    output$tuned_results_plot <- renderPlot({
      if (length(rna_splsda$model) != 1) {
        plot(rna_splsda$tuned)
        }
      })
    
    output$splsda_rna_signature <- renderDataTable({
      if (length(rna_splsda$model) != 1) {
        selectVar(rna_splsda$model, comp = as.numeric(input$splsda_ncomp_rna))$value
        }
      })
    
    # Metabolomics
    
    ## Reading metabolomic data and checking integrity (sample name match)
    
    intensities <- reactive({
      req(input$raw_metabo)
      intesities <- read.csv(input$raw_metabo$datapath, header = T, row.names = 1)
      intesities <- as.data.frame(intesities)
    })
    
    ## Showing intensity data to the user
    
    output$intensities_table <- renderDataTable({
      datatable(intensities(), rownames = T, options = list(pageLength = 10), caption = "Intensities")
    })
    
    ## Checking if data is numeric
    
    output$intensities_error <- renderText({ 
      if (any(!apply(intensities(),2,is.numeric))) {
        "Error in intensities: not all numeric values. Upload only intensities!"
      }
      else "Raw intensity data: OK"
    })
    
    ## Reading metadata data (factors)
    
    metadata_metabo <- reactive({
      req(input$design_metabo)
      metadata_metabo <- read.csv(input$design_metabo$datapath, header = T, row.names = 1)
      metadata_metabo <- apply(metadata_metabo,2,as.factor)
      metadata_metabo <- as.data.frame(metadata_metabo)
    })
    
    ## Showing metadata data to the user
    
    output$metadata_metabo <- renderDataTable({
      datatable(metadata_metabo(), rownames = T, options = list(pageLength = 10), caption = "Metabolomic metadata")
    })
    
    ## Checking if metadata and intensity data have same samples (number and name)
    
    output$metadata_metabo_error <- renderText({
      req(metadata_metabo())
      if (any(!colnames(intensities()) %in% row.names(metadata_metabo())) & # Checking if all samples are present
          ncol(intensities()) == nrow(metadata_metabo())) {
        "Error in metadata: samples do not match in intensity table and metadata. Make sure same samples are present and have the same label"
      }
      else "Metadata: OK"
    })
    
    # Preprocessing metabolomic data
    
    ## Getting QC and Blank samples
    
    ### Letting user chose what factor contains blank and qc 
    
    observe({
      req(metadata_metabo())
      updateSelectInput(session, "qc_blank_factor", choices = colnames(metadata_metabo()))
    })
    
    observeEvent(input$qc_blank_factor, {
      req(metadata_metabo())
      factor_selected <- input$qc_blank_factor
      if (!is.null(factor_selected) && factor_selected %in% names(metadata_metabo())) {
        factor_variable <- as.character(metadata_metabo()[[factor_selected]])
        levels_data <- unique(factor_variable)
        updateSelectInput(session, "blank_level", choices = levels_data)
        updateSelectInput(session, "qc_level", choices = levels_data)
      } else {
        updateSelectInput(session, "blank_level", choices = NULL)
        updateSelectInput(session, "qc_level", choices = NULL)
      }
    })
    
    ### Getting blank and qc samples
    
    type_sample <- reactiveValues(bio="",qc="",blank="")
    observeEvent(input$qc_level, {
      qc_samples <- metadata_metabo()[input$qc_blank_factor,] == input$qc_level
      qc_samples <- colnames(metadata_metabo())[qc_samples]
      type_sample$qc <- qc_samples
    })
    observeEvent(input$blank_level, {
      blank_samples <-  metadata_metabo()[input$qc_blank_factor,] == input$blank_level
      blank_samples <- colnames(metadata_metabo())[blank_samples]
      type_sample$blank <- blank_samples
    })
    
    ### Adding 0s to NA values to make easier computation
    
    intensities_filtered <- reactiveValues(values="")
    observe({
      req(intensities())
      intensities_filtered$values <- intensities()
      intensities_filtered$values[is.na(intensities_filtered$values)] <- 0 
    })
    
    ### CV in QC
    
    renderPlotly({
      qc_intensities <- intensities_filtered$values[type_sample$qc,]
      qc_cv <- apply(qc_intensities,1,function(x) mean(x)/sd(mean)*100)
    })
}