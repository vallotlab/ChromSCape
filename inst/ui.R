library(shiny)
library(shinyjs)
# 
shinyUI(shinydashboard::dashboardPage(skin='green',
                            shinydashboard::dashboardHeader(title = "ChromSCape"),
                            shinydashboard::dashboardSidebar(
                              shinydashboard::sidebarUserPanel("Institut Curie - Vallot Lab",
                                               subtitle = a(href = "#", icon("circle", class = "text-success"), "Online"),
                                               image = "curie.jpg"
                              ),
                              shinydashboard::sidebarMenu(id="tabs", style = "position: fixed; overflow: visible;",
                                          shinydashboard::menuItem("Select & Import", tabName = "select_import", icon=icon("upload")),
                                          shinydashboard::menuItem("Filter & Normalize", tabName = "filter_normalize", icon=icon("fas fa-filter")),
                                          shinydashboard::menuItem("Vizualizing Cells", tabName = "vizualize_dim_red", icon=icon("fas fa-image")),
                                          shinydashboard::menuItem("Correlation clustering", tabName = "cor_clustering", icon=icon("sitemap")),
                                          shinydashboard::menuItem("Consensus clustering", tabName = "cons_clustering", icon=icon("th")),
                                          shinydashboard::menuItem("Peak calling", tabName = "peak_calling", icon=icon("chart-area")), #mountain
                                          shinydashboard::menuItem("Differential analysis", tabName = "diff_analysis", icon=icon("chart-bar")),
                                          shinydashboard::menuItem("Enrichment analysis", tabName = "enrich_analysis", icon=icon("code-branch")),
                                          shinydashboard::menuItem("Close App & Save Analysis", tabName = "close_and_save", icon=icon("close"))
                              )
                            ),
                            shinydashboard::dashboardBody(
                              shinyjs::useShinyjs(),
                              tags$head(includeCSS('www/style.css')),
                              tags$head(includeCSS('www/introjs.min.css')),
                              tags$head(includeCSS('www/app.css')),
                              tags$style(type="text/css",
                                         ".shiny-output-error { visibility: hidden; }",
                                         ".shiny-output-error:before { visibility: hidden; }"
                              ),
                              tags$head(tags$style(HTML(
                                '.myClass {
                          font-size: 20px;
                          line-height: 50px;
                          text-align: left;
                          font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
                          padding: 0 15px;
                          overflow: hidden;
                          color: white;
                        }
                        '))),
                              tags$script(HTML('
                           $(document).ready(function() {
                           $("header").find("nav").append(\'<div id="pageHeader" class="myClass"></div>\');
                           })
                           ')),
                              shinydashboard::tabItems(
                                
                                
                                ###############################################################
                                # 1. Select Analysis & Import dataset
                                ###############################################################
                                
                                shinydashboard::tabItem(tabName = "select_import",
                                        
                                fluidPage(
                                  # Add ressources for UX
                                    shiny::includeScript(file.path(system.file(package="ChromSCape"),"js.cookie.js")),
                                    shiny::includeScript(file.path(system.file(package="ChromSCape"),"intro.min.js")),
                                    shiny::includeScript(file.path(system.file(package="ChromSCape"),"app.js")),
                                    #Load shinyJS added functions
                                    shinyjs::extendShinyjs(script = file.path(system.file(package="ChromSCape"),"shiny_js_functions.js"),
                                                           functions = c("init_directory","save_cookie","disableTab","enableTab")),
                                    #Center Panel
                                    column(width=6,align="left",
                                           shinydashboard::box(title="Select output directory", width = NULL, status="warning",
                                                                solidHeader=T,
                                                               column(12, align="left",
                                                                      shinyFiles::shinyDirButton("data_folder", "Input directory", "Upload") %>%
                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                              content = "output_directory"),
                                                                      verbatimTextOutput("directory", placeholder = TRUE)),
                                                               column(12, align="left", textOutput("data_folder_info"))),
                                           shinydashboard::box(title="Create new analysis & import raw data", width = NULL, status="warning", solidHeader=T,
                                                               column(12, align="left",
                                                                      textInput("new_analysis_name", "Enter a name for the new analysis :", value = "Analysis_1"),
                                                                      selectInput("annotation","Select reference genome:", choices=c("hg38", "mm10")),
                                                                      textOutput("data_matrices_info") %>%
                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                              content = "datamatrix_input"),
                                                                      fileInput("datafile_matrix", "Upload all data matrices (.txt or .tsv) :",
                                                                                multiple=TRUE, accept=c("text", "text/plain", ".txt", ".tsv")),
                                                                      actionButton("create_analysis", "Create analysis"))),
                                    ),
                                    column(width=6, align="right",
                                           actionButton("startHelp","Help", icon = icon("question-circle"))),
                                    column(width=5, align="right",
                                           column(width=12,br(),br(),br(),br(),
                                                  br(),br(),br()),
                                           column(width=12,
                                           shinydashboard::box(title="Select Analysis", width = NULL, status="success", solidHeader=T,
                                                               column(12, align="left",
                                                                      uiOutput("selected_analysis"),
                                                               ))
                                           )
                                    )
                                )
                        ),
                        
                        ###############################################################
                        # 2. Filter and Normalize dataset
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "filter_normalize",
                                                
                                                fluidPage(
                                                    #Right Panel
                                                    column(width=6,
                                                           shinydashboard::box(title="Filter & Normalize panel", width = NULL, status="warning", solidHeader=T,
                                                                               column(4, align="left"),
                                                                               column(8, align="left",
                                                                                      sliderInput("coverage_bins", "modify bin size of histogram :",
                                                                                                  min=5, max=100, value=50, step=5)),
                                                                               column(12, align="left",
                                                                                      plotly::plotlyOutput("cell_coverage", height=250),
                                                                                      br(), hr(),
                                                                                      uiOutput("table_QC_filt_box"),
                                                                                      br(), hr(),
                                                                                      sliderInput("min_coverage_cell", "Select minimum number of reads per cell :",
                                                                                                  min=100, max=5000, value=1600, step=100) %>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "filtering_parameters"),
                                                                                      sliderInput("quant_removal", "Select the upper percentile of cells to remove (potential doublets):", min=80, max=100, value=95, step=1),
                                                                                      sliderInput("min_cells_window", "Select minimum percentage of cells to support a window :", min=0, max=20, value=1, step=0.25),
                                                                                      checkboxInput("run_tsne", "Run T-SNE", value= TRUE) %>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "run_tsne"),
                                                                                      checkboxInput("do_subsample", "Perform subsampling", value=FALSE)%>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "subsampling"),
                                                                                      htmlOutput("do_subsample"),
                                                                                      checkboxInput("exclude_regions", "exclude specific genomic regions", value=FALSE)  %>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "exclude_region"),
                                                                                      htmlOutput("exclude_file"),
                                                                                      checkboxInput("do_batch_corr", "perform batch correction", value=FALSE) %>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "batch_correction"),
                                                                                      htmlOutput("num_batches")),
                                                                               column(5, align="left",
                                                                                      htmlOutput("batch_names")),
                                                                               column(7, align="left",
                                                                                      htmlOutput("batch_sel")),
                                                                               column(12, align="left",
                                                                                      hr(),
                                                                                      actionButton("filter_normalize_reduce", "Filter, Normalize & Reduce")))),
                                                    column(width = 6,
                                                           shinydashboard::box(title="Select filtered & normalized dataset", width = NULL, status="success", solidHeader=T,
                                                                               column(12, align="left",
                                                                                      htmlOutput("selected_reduced_dataset"),
                                                                                      textOutput("red_data_selection_info"),
                                                                                      textOutput("red_data_selection_format"))))
                                                    )
                        ),
                        
                        
                        
                        ###############################################################
                        # 2. Vizualising cells in reduced dimension
                        ###############################################################
                        
                        shinydashboard::tabItem(
                          tabName = "vizualize_dim_red",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title="PCA visualization", width = NULL, status="success", solidHeader=T,
                                             column(6, align="left", htmlOutput("color_by")),
                                             column(12, align="left", plotly::plotlyOutput("pca_plot") %>%
                                                        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                            content = "pca_plot")),
                                             column(3, align="left", htmlOutput("pc_select_x")),
                                             column(3, align="left", htmlOutput("pc_select_y"))),
                                         uiOutput("color_box")),
                                  column(width=6,
                                         uiOutput("tsne_box"),
                                         shinydashboard::box(title="UMAP visualization", width = NULL, status="success", solidHeader=T,
                                                             column(12, align="left", plotly::plotlyOutput("umap_plot") %>%
                                                                        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                            content = "umap_plot"))),
                                         )
                                )
                        ),
                        
                        ###############################################################
                        # 3. Correlation clustering
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "cor_clustering",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title="Correlation Heatmap & clustering table", width=NULL, status="success", solidHeader=T,
                                             column(12, align="center", plotOutput("corr_clust_pca_plot", height=500, width=500) %>%
                                                        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                            content = "correlation_clustering")),
                                             column(12, align="left", 
                                                    downloadButton("download_cor_clust_plot", "Download image"), br()),
                                             column(12, align="left", 
                                                    tableOutput('num_cell_before_cor_filt'))
                                             )),
                                
                                column(width=6,
                                         shinydashboard::box(title="Data filtering based on inter-cell correlation", width=NULL, status="warning", solidHeader=T,
                                             column(12, align="left", plotOutput("cell_cor_hist_plot", height=300, width=500) %>%
                                                        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                            content = "filter_correlation_distrib"),
                                                    downloadButton("download_cor_clust_hist_plot", "Download image"),
                                                    hr(),
                                                    sliderInput("corr_threshold", "correlation threshold quantile:", min=75, max=99, value=99, step=1),
                                                    sliderInput("percent_correlation", "min percent correlation of cells to others in data set:", min=0, max=15, value=1, step=0.25)),
                                             column(3, align="left", br(),
                                                    actionButton("filter_corr_cells", "Filter & save")),
                                             ),
                                       uiOutput("corr_filtered_hist")
                                       ))
                        ),
                        
                        ###############################################################
                        # 4. Consensus clustering on correlated cells
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "cons_clustering",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title="Selection of filtered dataset", width=NULL, status="warning",
                                                             solidHeader=T,
                                             column(12, align="left", htmlOutput("selected_filtered_dataset"),
                                                    textOutput("filtered_data_selection_format"))),
                                         shinydashboard::box(title="Clustering results", width=NULL, status="warning", solidHeader=T,
                                             column(4, align="left", actionButton("do_cons_clust", "Perform clustering"),
                                                    br()),
                                             column(12, align="left", textOutput("cluster_consensus_info"),
                                                    br())
                                             ),
                                         shinydashboard::box(title="Clustering results", width=NULL, status="success", solidHeader=T,
                                             column(12, align="left", uiOutput("cluster_consensus_png")),
                                             column(12, align="left", uiOutput("cons_clust_pdf")))),
                                  column(width=6,
                                         shinydashboard::box(title="Choose number of cluster", width=NULL, status="warning", solidHeader=T,
                                             column(5, align="left", selectInput("nclust", "Select number of clusters:", choices=c(2:10))),
                                             column(3, align="left", br(),
                                                    actionButton("choose_cluster", "Select & visualize")),
                                             column(12, align="left", textOutput("nclust_selection_info"))),
                                         uiOutput("anno_cc_box"),
                                         uiOutput("anno_corc_box"),
                                         uiOutput("tsne_box_cf"),
                                         uiOutput("umap_box_cf"),
                                         uiOutput("color_box_cf")))
                        ),
                        
                        ###############################################################
                        # 5. Peak calling [optional]
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "peak_calling",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title="Peak calling", width=NULL, status="warning", solidHeader=T,
                                             column(12, align="left", textOutput("peak_calling_info"), hr(),
                                                    htmlOutput("peak_calling_system"), hr(),
                                                    sliderInput("peak_distance_to_merge", "Select distance of peaks to merge:", min=0, max=50000, value=5000, step=1000),
                                                    uiOutput("bam_upload")),
                                             column(4, align="left", textOutput("pc_k_selection"),
                                                    selectInput("pc_stat","Select statistic for cutoff:", choices=c("p.value", "q.value"), selected="p.value")),
                                             column(8, align="left", br(), br(), br(), br(),
                                                    sliderInput("pc_stat_value", "select significance threshold:", min=0, max=0.25, value=0.05, step=0.01)),
                                             column(12, align="left", hr(), actionButton("do_pc", "Start"))))
                                  # ,column(width=6, uiOutput("pc_plot_box"))
                                  )
                        ),
                        
                        ###############################################################
                        # 6. Differential analysis
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "diff_analysis",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title="Differential Analysis parameters", width=NULL, status="warning", solidHeader=T,
                                             column(12, align="left", textOutput("diff_analysis_info"), br()),
                                             column(5, align="left", textOutput("selected_k")),
                                             column(5, align="left", selectInput("de_type", "Select type of cluster comparison:", choices=c("one_vs_rest","pairwise"))),
                                             column(12, align="left", sliderInput("qval.th", "adjusted p-value to select significant locations:", min=0.01, max=0.4, value=0.01, step=0.01),
                                                    sliderInput("cdiff.th", "Minimum log-fold change to select significant locations:", min=0, max=3, value=1, step=0.01),
                                                    checkboxInput("only_contrib_cells", "only use cells contributing most to the clustering", value=FALSE)),
                                             column(1, align="left", br()),
                                             column(11, align="left", uiOutput("contrib_thresh")),
                                             column(1, align="left", br()),
                                             column(8, align="left", uiOutput("contrib_hist")),
                                             column(3, align="left", br(), uiOutput("contrib_info")),
                                             column(12, align="left", hr(), actionButton("do_wilcox", "Start analysis"))),
                                         uiOutput("da_summary_box")),
                                  column(width=6,
                                         uiOutput("da_visu_box")))
                        ),
                        
                        ###############################################################
                        # 7. Enrichment analysis
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "enrich_analysis",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title="Enrichment analysis", width=NULL, status="warning", solidHeader=T,
                                             column(12, align="left", textOutput("enr_info"), br(),
                                                    uiOutput("use_peaks"),
                                                    actionButton("do_enrich", "Start enrichment analysis"))),
                                         shinydashboard::box(title="Enriched gene sets in differential regions across clusters", width=NULL, status="warning", solidHeader=T,
                                             column(4, align="left", uiOutput("enr_clust_sel"), br()),
                                             column(8, align="left", uiOutput("enr_class_sel"), br()),
                                             column(12, align="left",
                                                         mainPanel(tabsetPanel(id='enr_tables',
                                                                          tabPanel("In differential loci", div(style = 'overflow-x: scroll', DT::dataTableOutput('all_enrich_table'))),
                                                                          tabPanel("In enriched loci", div(style = 'overflow-x: scroll', DT::dataTableOutput('over_enrich_table'))),
                                                                          tabPanel("In depleted loci", div(style = 'overflow-x: scroll', DT::dataTableOutput('under_enrich_table')))), width=12),
                                                    br(), br(), br(), br(), br(),
                                                    downloadButton("download_enr_data", "Download tables")))),
                                  column(width=6,
                                         shinydashboard::box(title="Binding strength near TSS", width=NULL, status="success", solidHeader=T,
                                             column(4, align="left", uiOutput("gene_sel")),
                                             column(8, align="left", uiOutput("region_sel")),
                                             column(12, align="left", plotly::plotlyOutput("gene_tsne_plot")))))
                        ),
                        
                        ###############################################################
                        # 8. Close app
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName="close_and_save",
                                fluidPage(shinyjs::useShinyjs(),
                                          # extendShinyjs(text = jscode, functions = c("closeWindow")),
                                          column(width=6,
                                                 shinydashboard::box(title='Close App & Save Analysis', solidHeader=T, status='danger', width=NULL,
                                                     column(12, actionButton("close_and_save", "Close App & Save Analysis")))
                                          ),
                                          column(width=6,
                                                 shinydashboard::box(title="Delete analysis", width = NULL, status="success", solidHeader=T,
                                                                     column(9, align="left", uiOutput("selected_delete_analysis")),
                                                                     column(3, align="left", br(), actionButton("delete_analysis", "Delete")),
                                                                     column(12, align="left", textOutput("analysis_deletion_info"))))
                                )
                        )
                        
                      )))
)