library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(plotly)

ui <- dashboardPage(
  dashboardHeader(title = "EasyOmics"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction", tabName = "introduction", icon = icon("door-open")),
      menuItem("RNAseq", tabName = "rnaseq", icon = icon("dna"),
               menuSubItem("Data", tabName = "rnaseq_data", icon = icon("database")),
               menuSubItem("Preprocessing", tabName = "rnaseq_preprocessing", icon = icon("filter")),
               menuSubItem("Normalization", tabName = "rnaseq_normalization", icon = icon("chart-simple")),
               menuSubItem("Differential expression", tabName = "rnaseq_de", icon = icon("arrows-left-right")),
               menuSubItem("Biomarker identification", tabName = "rnaseq_biomarker_identification", icon = icon("marker"))
               ),
      menuItem("Metabolomics", tabName = "metabolomics", icon = icon("flask"),
               menuSubItem("Data", tabName = "metabo_data", icon = icon("database")),
               menuSubItem("Preprocessing", tabName = "metabo_preprocessing", icon = icon("filter")),
               menuSubItem("Normalization", tabName = "metabo_normalization", icon = icon("chart-simple")),
               menuSubItem("Biomarker identification", tabName = "metabo_biomarker_identification", icon = icon("marker"))
               ),
      menuItem("Omics integration", tabName = "omics_integration", icon = icon("microscope"))
    )
  ),
  dashboardBody(
    
    shinyDashboardThemes(
      theme = "blue_gradient"
    ),
    
    tabItems(
      
      # RNAseq Data section
      
      tabItem(
        tabName = "rnaseq_data",
        fluidRow(
          column(
            width = 12,
            h2("Upload your data"),
            br(),
            p("Browse to your count data and metadata by clicking on the", strong("Browse..."), "button. 
              You can check the data once uploaded. Note that the data you upload has to fulfill the following conditions:"
              ),
            tags$div(
              HTML("<ul>
                  <li>Files have to be in CSV format. That means values are separated by commas (,). 
                  You can use any spreadsheet processor to transform your table into CSV format.</li>
                  <li>Your count table should have samples as columns and genes as rows.</li>
                  <li>Your metadata should have variables as columns and samples as rows.</li>
                  <li>Both count and metadata should have same samples and they should have same labels.</li>
                  </ul>"),
              p("Once you get OK messages for both files you can proceed to the",strong("Preprocessing"),
                "tab under RNAseq. These messages will appear under the following menus.")
              ),
            ),
            column(
              width=6,
              box(
                width = 12,
                solidHeader = T,
                background = "teal",
                fileInput("raw_rna", "Raw counts RNAseq"),
                orientation = "horizontal",
                br()
                ),
              ),
            column(
              width=6,
              box(
                width = 12,
                solidHeader = T,
                background = "teal",
                fileInput("design_rna", "Metadata RNAseq"),
                orientation = "horizontal",
                br()
                ),
              ),
            column(
              width = 6,
              verbatimTextOutput("count_rna_error"),
              br()
              ),
            column(
              width = 6,
              verbatimTextOutput("metadata_rna_error"),
              br()
              ),
          br(),
          column(
            width = 12,
            dataTableOutput("count_rna"),
            dataTableOutput("metadata_rna")
            )
          )
        ),
      
      # RNAseq preprocessing (filtering, PCA to check outliers)
      
      tabItem(
        tabName = "rnaseq_preprocessing",
        fluidRow(
          column(
            width = 12,
            h2("Filtering the data"),
            h3("Gene count and library size"),
            br(),
            p("In this section you are going to filter the data in order to remove low expressed genes or low quality samples.
              For that purpose, you can control in the", strong('"Counts per gene"'), "and", strong('"Library size (x10k)"'), "boxes the minimum amount of counts per gene and
              the library size. Note that each unit in the second slider it is equivalent to 10 thousand counts
              (e.g., 10 in the slider would be a filter of hundred thousand counts).You need to chose these parameters depending
              on the characteristics of your data, however the defaults values are good start."),
            p("You will find messages reporting the amount of genes/samples filtered with your chosen parameters."),
            ),
          ),
        br(),
        fluidRow(
          column(
            width=6,
            box(
              width = 12,
              solidHeader = F,
              background = "teal",
              sliderInput("min_count", "Counts per gene", min=0, max=100, value = 10),
              orientation = "horizontal"
              ),
            ),
            column(
              width=6,
              box(
                width = 12,
                solidHeader = F,
                background = "teal",
                sliderInput("min_library", "Library size (x10k)", min=0, max=100, value=20),
                orientation = "horizontal"
              ),
            ),
            column(
              width = 6,
              verbatimTextOutput("genes_filtered"),
              ),
            column(
              width = 6,
              verbatimTextOutput("samples_lowlib"),
              )
            ),
        br(),
        fluidRow(
          column(
            width = 6,
            plotOutput("low_expressed_plot")
            ),
          column(
            width = 6,
            plotOutput("library_size_plot")
            )
          ),
        fluidRow(
          column(
            width = 12,
            h3("PCA and outlier detection"),
            br(),
            p("In this section a Principal Component Analysis (PCA) using all genes is performed. 
            What PCA does is to reduce all variables (genes in this case) into 2 variables named Principal Components (PC) that summarizes the information in the dataset.
            Since all genes are used, this makes difficult to interpret or rather extract the biological meaning of these new variables, 
            but it is useful to detect any possible sample that is behaving in a unexpected manner."
              ),
            p("In order to inspect the data, 
              you can colour the samples according to any factor available in the metadata you uploaded."          
              ),
            p("Finally in case you find any sample that does not fit your expectations, 
              you can remove it using the menu in the yellow box labeled", strong("Sample"), "and click", strong("Remove"), ". In case you deleted the wrong sample you can", 
              strong("Restart"), "over by clicking over the button."
              ),
            )
          ),
        br(),
        fluidRow(
          column(
            width=12,
            box(
              width=4,
              solidHeader = T,
              background = "teal",
              selectInput("outlier_pca_colour",
                label = "Colour",
                choices = NULL
              )
            ),
            box(
              width=4,
              solidHeader = T,
              background = "teal",
              selectInput("outlier_pca_shape",
                          label = "Shape",
                          choices = NULL
              )
            ),
            box(
              width=4,
              solidHeader = T,
              background = "yellow",
              column(width = 10,
                selectInput("rna_remove_outlier",
                            label = "Sample",
                            choices = NULL
                            ),
              ),
              column(width = 2,
                     fluidRow(
                     actionButton("rna_submit_outlier", "Remove")
                     ),
                     fluidRow(
                       br(),
                       actionButton("recover_all_samples", "Restart")
                     )
              )
            ),
          )
        ),
        plotlyOutput("outlier_pca", height=700),
      ),
      
      # RNA normalization
      
      tabItem(tabName = "rnaseq_normalization",
              fluidRow(
                column(
                  width=12,
                  h3("Normalization and download preprocessed data"),
                  br(),
                  p("Normalization is an essential step in order to proceed with your statistical analysis.
                  This process is in charge of removing differences between samples due to technical variability, thus latter applied methods will be capturing biological variability.
                  This tool uses DESeq2 library in R to perform library size factor normalization. In short, what it does is to correct by sequencing depth (or library size) in order to
                  be able to compare the expression between samples."
                  ),
                  p("The following plots show distribution of counts before and after the normalization is applied. You expect that medians of each sample distribuion are located around the
                  overall median. In case there is a sample that does not scale properly, consider reviewing the preprocessing tab and check if more strict filterings may be applied or if
                  the sample could be considered an outlier in the PCA."
                  ),
                  p(strong("Download this data","for further analysis, like omics integration using mixOmics, or any other tool that fulfils your needs.")),
                  br()
                ),
              ),
              fluidRow(
                column(
                  width=4, offset=2,
                  downloadButton("download_counts","Download expression data")
                  ),
                column(
                  width=4, offset=2,
                  downloadButton("download_rna_metadata","Download expression metadadata")
                  )
                ),
              fluidRow(
                br(),
                column(
                  width=6,
                  plotOutput("rna_before_norm", height = 600)
                ),
                column(
                  width=6,
                  plotOutput("rna_after_norm", height = 600)
                )
              )
            ),
      
      #DE RNAseq
      
      tabItem(tabName = "rnaseq_de",
              fluidRow(
                column(
                  width=12,
                  h3("Differential expression analysis"),
                  br(),
                  p("Here you will perform differential expression analysis between the factor of your choice. Also, you will be able to translate to gene ID into gene symbol
                    to make easier the results to interpret"),
                  br()
                )
              ),
              fluidRow(
                column(
                  width=12,
                  box(
                    width=12,
                    solidHeader = T,
                    background = "teal",
                    column(
                      width=6,
                      h4("Comparison parameters"),
                      selectInput("de_factor", "Factor", choices=NULL),
                      selectInput("de_level1", "Reference level", choices=NULL),
                      selectInput("de_level2", "Comparison level", choices=NULL),
                      actionButton("de_submit", "Run")
                    ),
                    column(
                      width=6,
                      h4("Gene symbol translation"),
                      actionButton("gene_symbol_submit", "Translate")
                    )
                  )
                )
              ),
              fluidRow(
                column(
                  width=4,
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  br(),
                  box(
                    width=12,
                    title = "Volcano plot thresholds",
                    solidHeader = T,
                    background = "teal",
                    sliderInput("volcano_pvalue", "p-value cutoff", min=0, max=0.1, value=0.05, step=0.001),
                    sliderInput("volcano_fc", "Fold change cutoff", min=0, max=10, step=0.1, value=1),
                  ),
                  br(),
                  br(),
                  downloadButton("download_de","Download raw differential expression results"),
                ),
                column(
                  width=8,
                  br(),
                  br(),
                  plotlyOutput("volcano_rna", height = 600)
                )
              )
            ),
      
      # sPCA and sPLS-DA
      
      tabItem(tabName = "rnaseq_biomarker_identification",
              fluidRow(
                column(
                  width=12,
                  h3("Biomarker identification"),
                  br(),
                  p("In the previous section we extracted differentially expressed features from the dataset. However, extracting knowledge from all these variables makes it hard because of the amount of variables
                    to interpret. That is also know as the curse of dimensionality. In this section two different machine learning implementations offered by the library mixOmics are used, sparse Principal Component Analysis (sPCA) and sparse Partial Least Squares Discriminant Analysis (sPLS-DA)
                    Note that sparse means that these algorithms selects a certain number of features in order the create new variables that summarizes the dataset. So these algorithms can use to further reduce the differential expressed genes found in the previous section. 
                    The differences are the first one is unspervised, meaning, it is clueless to the factor you want to study (example, test vs control). The second one (sPLS-DA) uses that information in order to discriminate between groups and can outperform unspervised techniques, 
                    however with the risk of overfitting your data and possibly making it harder to generalize your results.
                    "),
                  br()
                )),
              fluidRow(
                column(
                  width=12,
                  h3("sPCA"),
                  br(),
                  p("This technique (as the following one) needs to be parametrized. That means, some sort of parameters have to be tuned. In this case, the number of features to select for each new component. 
                    Using the sliders"),
                  br()
                )),
              fluidRow(
                column(
                  width=3,
                  br(),
                  br(),
                  box(
                    width=12,
                    title = "Number of features to select",
                    solidHeader = T,
                    background = "teal",
                    sliderInput("grid1", "Features in PC1", min=0, max=100, step=1, value=20),
                    sliderInput("grid2", "Features in PC2", min=0, max=100, step=1, value=20)
                  ),
                ),
                column(
                  width=9,
                  plotOutput("rna_spca", height = 600)
                )
              ),
              fluidRow(
                br(),
                br(),
                column(
                  width=6,
                  h4("Selected features in PC1 (sPCA)"),
                  br(),
                  dataTableOutput("selected_spc1")
                ),
                column(
                  width=6,
                  h4("Selected features in PC2 (sPCA)"),
                  br(),
                  dataTableOutput("selected_spc2")
                )
              ),
              fluidRow(
                br(),
                br(),
                column(
                  width=12,
                  h3("sPLS-DA"),
                  p("sPLS-DA is computationally costly. If it takes a lot of time you can tune the parameters to reduce the time, like being more 
                    stringent with your differential expression analysis (increasing fold change and decreasing p-value filtering) and reducing the number of features to select.
                    This trade off will result probably in a less performant model, however it will reduce computation time"),
                  br()
                )
              ),
              fluidRow(
                br(),
                column(
                  width = 3,
                  box(
                    width=12,
                    title="Parameters sPLS-DA",
                    background = "teal",
                    sliderInput("min_splsda_features_rna", "Minimum features per component", min=1, max=300, value=1),
                    sliderInput("max_splsda_features_rna", "Maximum features per component", min=1, max=500, value=50),
                    sliderInput("step_tuning_model", "Step features", min=1, max=20, value=10),
                    sliderInput("nrpeats_tuning_model", "Number of repeats cross validation", min=1, max=500, value=50),
                    actionButton("splsda_submit_rna", "Run")
                  ),
                  box(
                    width = 12,
                    title="sPLS-DA plot",
                    background = "teal",
                    selectInput("x_plsda_rna", "PC in X-variate", choices=c(1,2), selected=1),
                    selectInput("y_plsda_rna", "PC in X-variate", choices=c(1,2), selected=2)
                  )
                ),
                column(
                  width = 8,
                  br(),
                  plotOutput("rna_splsda_plot", height=700)
                )
              ),
              fluidRow(
                br(),
                br(),
                column(
                  width=6,
                  plotOutput("tuned_results_plot", height=600)
                ),
                column(
                  width=6,
                  selectInput("splsda_ncomp_rna", "Signature in PC", choices=c(1,2), selected = 1),
                  dataTableOutput("splsda_rna_signature")
                  )
              )
      ),
    
      # Metabolics section
      
      tabItem(
        tabName = "metabo_data",
        fluidRow(
          column(
            width = 12,
            h2("Upload your data"),
            br(),
            p("Browse to your count data and metadata by clicking on the", strong("Browse..."), "button. 
              You can check the data once uploaded. Note that the data you upload has to fulfill the following conditions:"
            ),
            tags$div(
              HTML("<ul>
                  <li>Files have to be in CSV format. That means values are separated by commas (,). 
                  You can use any spreadsheet processor to transform your table into CSV format.</li>
                  <li>Your count table should have samples as columns and genes as rows.</li>
                  <li>Your metadata should have variables as columns and samples as rows.</li>
                  <li>Both count and metadata should have same samples and they should have same labels.</li>
                  </ul>"),
              p("Once you get OK messages for both files you can proceed to the",strong("Preprocessing"),
                "tab under RNAseq. These messages will appear under the following menus.")
            ),
          ),
          column(
            width=6,
            box(
              width = 12,
              solidHeader = T,
              background = "teal",
              fileInput("raw_metabo", "Raw intensities"),
              orientation = "horizontal",
              br()
            ),
          ),
          column(
            width=6,
            box(
              width = 12,
              solidHeader = T,
              background = "teal",
              fileInput("design_metabo", "Metadata"),
              orientation = "horizontal",
              br()
            ),
          ),
          column(
            width = 6,
            verbatimTextOutput("intensities_error"),
            br()
          ),
          column(
            width = 6,
            verbatimTextOutput("metadata_metabo_error"),
            br()
          ),
          br(),
          column(
            width = 12,
            dataTableOutput("intensities_table"),
            dataTableOutput("metadata_metabo")
          )
        )
      ),
      
      # Metaboloics preprocessing

      tabItem(
        tabName = "metabo_preprocessing",
        fluidRow(
          column(
            width = 12,
            h2("Filtering the data"),
            h3("Metabolite count and library size"),
            br(),
            p("In this section you are going to filter the data in order to remove low expressed metabolites or low quality samples.
              For that purpose, you can control in the", strong('"intensities per metabolite"'), "and", strong('"Library size (x10k)"'), "boxes the minimum amount of intensities per metabolite and
              the library size. Note that each unit in the second slider it is equivalent to 10 thousand intensities
              (e.g., 10 in the slider would be a filter of hundred thousand intensities).You need to chose these parameters depending
              on the characteristics of your data, however the defaults values are good start."),
            p("You will find messages reporting the amount of metabolites/samples filtered with your chosen parameters."),
          ),
        ),
        br(),
        fluidRow(
          column(
            width=6,
            box(
              width = 12,
              solidHeader = F,
              background = "teal",
              # sliderInput("min_count", "Intensities per metabolite", min=0, max=100, value = 10),
              selectInput("qc_blank_factor", "QC/Blank Factor", choices=NULL),
              selectInput("blank_level", "Blank level", choices=NULL),
              selectInput("qc_level", "QC level", choices=NULL),
              orientation = "horizontal"
            ),
          ),
          column(
            width=6,
            box(
              width = 12,
              solidHeader = F,
              background = "teal",
              sliderInput("min_library", "Library size (x10k)", min=0, max=100, value=20),
              orientation = "horizontal"
            ),
          ),
          column(
            width = 6,
            verbatimTextOutput("metabolites_filtered"),
          ),
          column(
            width = 6,
            verbatimTextOutput("samples_low_intensities"),
          )
        ),
        br(),
        fluidRow(
          column(
            width = 6,
            plotOutput("low_intensities_plot")
          ),
          column(
            width = 6
            # plotOutput("library_size_plot")
          )
        ),
        fluidRow(
          column(
            width = 12,
            h3("PCA and outlier detection"),
            br(),
            p("In this section a Principal Component Analysis (PCA) using all metabolites is performed.
            What PCA does is to reduce all variables (metabolites in this case) into 2 variables named Principal Components (PC) that summarizes the information in the dataset.
            Since all metabolites are used, this makes difficult to interpret or rather extract the biological meaning of these new variables,
            but it is useful to detect any possible sample that is behaving in a unexpected manner."
            ),
            p("In order to inspect the data,
              you can colour the samples according to any factor available in the metadata you uploaded."
            ),
            p("Finally in case you find any sample that does not fit your expectations,
              you can remove it using the menu in the yellow box labeled", strong("Sample"), "and click", strong("Remove"), ". In case you deleted the wrong sample you can",
              strong("Restart"), "over by clicking over the button."
            ),
          )
        ),
        br(),
        fluidRow(
          column(
            width=12,
            box(
              width=4,
              solidHeader = T,
              background = "teal",
              selectInput("metabo_outlier_pca_colour",
                          label = "Colour",
                          choices = NULL
              )
            ),
            box(
              width=4,
              solidHeader = T,
              background = "teal",
              selectInput("metabo_outlier_pca_shape",
                          label = "Shape",
                          choices = NULL
              )
            ),
            box(
              width=4,
              solidHeader = T,
              background = "yellow",
              column(width = 10,
                     selectInput("metabo_remove_outlier",
                                 label = "Sample",
                                 choices = NULL
                     ),
              ),
              column(width = 2,
                     fluidRow(
                       actionButton("metabo_submit_outlier", "Remove")
                     ),
                     fluidRow(
                       br(),
                       actionButton("metabo_recover_all_samples", "Restart")
                     )
              )
            ),
          )
        ),
        plotlyOutput("metabo_outlier_pca", height=700)
        ),
      
      # Metabolomics normalization

      tabItem(
        tabName = "metabo_normalization"
        ),
      
      # Metabolomics sPCA and sPLS-DA

      tabItem(
        tabName = "metabo_biomarker_identification"
        ),

      tabItem(
        tabName = "omics_integration"
        )
    )
  )
)