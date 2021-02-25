library(shiny)
library(shinyjs)
library(dplyr)

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
                                          shinydashboard::menuItem("Visualize Cells", tabName = "vizualize_dim_red", icon=icon("fas fa-image")),
                                          shinydashboard::menuItem("Cluster Cells", tabName = "cons_clustering", icon=icon("th")),
                                          shinydashboard::menuItem("Peak Calling & Coverage", tabName = "peak_calling", icon=icon("chart-area")), #mountain
                                          shinydashboard::menuItem("Differential Analysis", tabName = "diff_analysis", icon=icon("chart-bar")),
                                          shinydashboard::menuItem("Pathway Enrichment Analysis", tabName = "enrich_analysis", icon=icon("code-branch")),
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
                                    shiny::includeScript(file.path(system.file(package="ChromSCape"),"www","shiny_js_functions.js")),

                                    #Load shinyJS added functions
                                    shinyjs::extendShinyjs(script = file.path(system.file(package="ChromSCape"),"www","shiny_js_functions.js"),
                                                           functions = c("init_directory","save_cookie","disableTab","enableTab")),
                                    #Center Panel
                                    column(width=6,align="left",
                                           shinydashboard::box(title="Select output directory", width = NULL, status="success",
                                                                solidHeader=TRUE,
                                                               column(12, align="left",
                                                                      shinyFiles::shinyDirButton("data_folder", "Browse", "Upload") %>%
                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                              content = "output_directory"),
                                                                      verbatimTextOutput("directory", placeholder = TRUE)),
                                                               column(12, align="left", textOutput("data_folder_info"))),
                                           shinydashboard::box(title="Import raw data & Create new analysis", width = NULL, status="success", solidHeader=TRUE,
                                                               column(12, align="left",
                                                                      column(12, align="left",textInput("new_analysis_name", "Enter a name for the new analysis :", value = "Analysis_1"),
                                                                      selectInput("annotation","Select reference genome:", choices=c("hg38", "mm10")),
                                                                      textOutput("data_matrices_info") %>%
                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                              content = "datamatrix_input")),
                                                                      column(6,
                                                                      radioButtons("data_choice_box", label = "Input data type",
                                                                                   choices = list("Count matrix(ces)"="count_mat",
                                                                                                  "Index, Peak & Barcode files" = "Index_Peak_Barcode",
                                                                                                  "Single-cell BAM files" = "BAM",
                                                                                                  "Single-cell BED files" = "BED"
                                                                                                  ))),
                                                                      column(6,
                                                                      uiOutput("input_data_ui")),
                                                                      column(12,
                                                                      uiOutput("advanced_data_input")),
                                                                      column(12,actionButton("create_analysis", "Create analysis")))),

                                                         shinydashboard::box(title="Load previous analysis", width = NULL, status="success", solidHeader=TRUE,
                                                                             column(12, align="left",
                                                                                    uiOutput("selected_analysis")
                                                                             )
                                                         )
                                    ),
                                    column(width=6, align="right",
                                           actionButton("startHelp","Help", icon = icon("question-circle")), br(), br()),
                                    column(width = 6, uiOutput("rename_file_box")),
                                )
                        ),
                        
                        ###############################################################
                        # 2. Filter and Normalize dataset
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "filter_normalize",
                                                
                                                fluidPage(
                                                    #Right Panel
                                                    column(width=6,
                                                           shinydashboard::box(title="Filter & Normalize panel", width = NULL, status="success", solidHeader=TRUE,
                                                                               column(6, align="left",
                                                                                      column(12, align = "middle", h4("Cell coverage")),
                                                                                      br(), br(),
                                                                                      plotly::plotlyOutput("cell_coverage", height=250) %>%
                                                                                          shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75),
                                                                                      br(), hr()),
                                                                               column(6, align="left",
                                                                                      column(12, align = "middle", h4("Feature coverage")),
                                                                                      br(), br(),
                                                                                      plotly::plotlyOutput("feature_coverage", height=250) %>%
                                                                                          shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75),
                                                                                      br(), hr()),
                                                                               column(12, align="left",
                                                                                      uiOutput("table_QC_filt_box"),
                                                                                      br(), hr(),
                                                                                      sliderInput("min_coverage_cell", shiny::HTML("<p><span style='color: green'>Select minimum number of reads per cell :</span></p>"),
                                                                                                  min=100, max=5000, value=1600, step=100) %>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "filtering_parameters"),
                                                                                      sliderInput("quant_removal", shiny::HTML("<p><span style='color: red'>Select the upper percentile of cells to remove (potential doublets):</span></p>"),
                                                                                                   min=80, max=100, value=95, step=1),
                                                                                      sliderInput("min_cells_window", shiny::HTML("<p><span style='color: #9C36CF'>Select minimum percentage of cells to support a window :</span></p>"), min=0, max=5, value=1, step=0.05),
                                                                                      checkboxInput("run_tsne", "Run T-SNE", value= FALSE) %>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "run_tsne"),
                                                                                      checkboxInput("do_subsample", "Perform subsampling", value=FALSE)%>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "subsampling"),
                                                                                      uiOutput("do_subsample"),
                                                                                      checkboxInput("exclude_regions", "Exclude specific genomic regions", value=FALSE)  %>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "exclude_region"),
                                                                                      uiOutput("exclude_file"),
                                                                                      checkboxInput("do_batch_corr", "Perform batch correction", value=FALSE) %>%
                                                                                          shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                                              content = "batch_correction"),
                                                                                      uiOutput("num_batches")),
                                                                               column(5, align="left",
                                                                                      uiOutput("batch_names")),
                                                                               column(7, align="left",
                                                                                      htmlOutput("batch_sel")),
                                                                               column(12, align="left",
                                                                                      hr(),
                                                                                      actionButton("filter_normalize_reduce", "Filter, Normalize & Reduce")))),
                                                    column(width = 6,
                                                           shinydashboard::box(title="Select filtered & normalized dataset", width = NULL, status="success", solidHeader=TRUE,
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
                                         shinydashboard::box(title="PCA", width = NULL, status="success", solidHeader=TRUE,
                                             column(6, align="left", uiOutput("feature_color")),
                                             column(12, align="left",plotOutput("pca_plot") %>%
                                                        shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75) %>%
                                                        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                            content = "pca_plot")
                                                        ),
                                             column(3, align="left", uiOutput("pc_select_x")),
                                             column(3, align="left", uiOutput("pc_select_y"))),
                                         uiOutput("contribution_to_pca_UI")),
                                  column(width=6,
                                         shinydashboard::box(title="UMAP", width = NULL, status="success", solidHeader=TRUE,
                                                             column(12, align="left", plotOutput("UMAP_plot") %>%
                                                                        shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75) %>%
                                                                        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                            content = "umap_plot")
                                                                    )),
                                         uiOutput("tsne_box")
                                         ),
                                  column(width=6, uiOutput("color_box"))
                                )
                        ),
                        
                        
                        ###############################################################
                        # 4. Consensus clustering on correlated cells
                        ###############################################################
                               
                        shinydashboard::tabItem(tabName = "cons_clustering",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title=tagList(shiny::icon("cube"),"Hierarchical Clustering"),
                                                             width=NULL, status="success", solidHeader=TRUE,
                                                             column(12, align ="center", plotOutput("hc_heatmap_plot", height=500, width=500) %>%
                                                                        shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75) %>%
                                                                        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                            content = "correlation_clustering"))
                                         ),
                                         shinydashboard::box(title = tagList(shiny::icon("cubes"),"Consensus Hierarchical Clustering"),
                                                             width=NULL, status="success", solidHeader=TRUE,
                                                             collapsible = TRUE, collapsed = TRUE,
                                                             column(12, align = "left", br()),
                                                             column(12, offset = 0,
                                                                    div(style="display: inline-block;vertical-align:top; width: 200px;", sliderInput("maxK", "Max cluster :", min=2, max=20, value=10, step=1)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 25px;",HTML("<br>")),
                                                                    div(style="display: inline-block;vertical-align:top; width: 200px;", sliderInput("consclust_iter", "Number of iterations:", min=10, max=1000, value=100, step=10)),
                                                                    div(style="display: inline-block;vertical-align:top; width: 25px;",HTML("<br>")),
                                                                    div(style="display: inline-block;vertical-align:top; width: 200px;", selectInput("clusterAlg", "Cluster Algorithm:", choices = c("Hierarchical","Partitioning Medoids","K-means"),selected = "hc")),
                                                                    br()),
                                                             column(3, align="left", actionButton("do_cons_clust", "Launch Consensus Clustering"),
                                                                    br()),
                                                             column(12, align = "left",
                                                                    hr()),
                                                             column(12, align ="center", uiOutput("cons_corr_clust_pca_UI")),
                                                             column(12, align="left", uiOutput("cluster_consensus_plot")),
                                                             column(12,hr()),
                                                             column(12, align="left", textOutput("cluster_consensus_info")),
                                                             column(12, align="left", uiOutput("cons_clust_pdf"))
                                         ),
                                         shinydashboard::box(title=tagList(shiny::icon("fas fa-filter"),"Filter lowly correlated cells"),
                                                             width=NULL, status="success", solidHeader=TRUE,
                                                             collapsible = TRUE, collapsed = TRUE,
                                                             column(12, align="center",
                                                                    plotOutput("cell_cor_hist_plot", height=300, width=500) %>%
                                                                        shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75) %>%
                                                                        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                                                            content = "filter_correlation_distrib")),
                                                             column(12, align = "left",
                                                                    hr(),
                                                                    sliderInput("corr_threshold", "Correlation threshold quantile:", min=75, max=99, value=99, step=1),
                                                                    sliderInput("percent_correlation", "Minimum percentage of cells to correlate with (%)", min=0, max=15, value=1, step=0.25)),
                                                             column(3, align="left", br(),
                                                                    actionButton("filter_corr_cells", "Filter & save")),
                                                             column(12,  uiOutput("table_cor_filtered")
                                                             )
                                         )
                                  ),
   
                                  column(width=6,
                                         shinydashboard::box(title="Select Cluster Number", width=NULL, status="success", solidHeader=TRUE,
                                                             align="left",
                                                             column(width=4,uiOutput("nclust_UI")),
                                                             column(width=3, br(),br(), checkboxInput("cluster_type",label = shiny::HTML("<b>Use consensus</b>"),value = FALSE)),
                                                             column(width=3,br(),br(),actionButton(inputId = "choose_cluster", label = "Choose Cluster"))
                                         ),
                                         uiOutput("contingency_table_cluster"),
                                         uiOutput("UMAP_box"),
                                         uiOutput("tsne_box_cf"),
                                         uiOutput("color_box_cf"),
                                         shinydashboard::box(title="Correlation heatmap with cluster annotation", width=NULL,
                                                             status="success", solidHeader=TRUE, align="left", 
                                                             column(width=4,
                                                                    actionButton(inputId = "do_annotated_heatmap_plot", label = "Plot Clustered Heatmap"), br()),
                                                             column(width=12, uiOutput("annotated_heatmap_UI"))
                                         ),
                                         shinydashboard::box(title=tagList(shiny::icon("fas fa-filter"), "Intra / Inter correlation"), width=NULL,
                                                             status="success", solidHeader=TRUE, align="left",
                                                             collapsible = TRUE, collapsed = TRUE,
                                                             column(12, align="left",
                                                                    column(width=4,uiOutput("violin_color")),
                                                                    column(width=4,checkboxInput("add_jitter", shiny::HTML("<b>Add single-cells</b>"), value= FALSE)),
                                                                    column(width=4,uiOutput("jitter_color")),  br(),br(),br(),br(),br(),
                                                                    mainPanel(tabsetPanel(id='inter_intra_cor',
                                                                                          
                                                                                          tabPanel("Intracorrelation", column(width=12, uiOutput("intra_corr_UI"))),
                                                                                          tabPanel("Intercorrelation", column(width=12, column(10, uiOutput("reference_group")),
                                                                                                                              br(), br(), br(), br(), br(),
                                                                                                                              uiOutput("inter_corr_UI")))
                                                                                          ))),
                                                             
                                         )
                                  )
                                )
                        ),
                        
                        ###############################################################
                        # 5. Peak calling [optional]
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "peak_calling",
                                fluidPage(
                                  column(width=4,
                                         shinydashboard::box(title="Peak calling", width=NULL, status="success", solidHeader=TRUE,
                                             column(12, align="left", textOutput("peak_calling_info"), hr()),
                                             tags$style(HTML(".large_icon { font-size: 70px; }")),
                                             column(5, align="left",
                                                    htmlOutput("peak_calling_system"), hr()),
                                             column(5,offset = 1,align="center", htmlOutput("peak_calling_icon")),
                                             column(12, align="left",
                                                    sliderInput("peak_distance_to_merge", "Select distance of peaks to merge:", min=0, max=50000, value=5000, step=100),
                                                    shinyFiles::shinyDirButton("pc_folder", "Browse folder of BAM / scBED files" ,
                                                                   title = "Please select a folder:",
                                                                   buttonType = "default", class = NULL), br(),
                                                    
                                                    textOutput("pc_dir"), br(),
                                                    uiOutput("pc_upload")),
                                             column(4, align="left", textOutput("pc_k_selection"),
                                                    selectInput("pc_stat","Select statistic for cutoff:", choices=c("p.value", "q.value"), selected="p.value")),
                                             column(12, align="left", br(), br(),
                                                    sliderInput("pc_stat_value", "Select significance threshold:", min=0, max=0.25, value=0.05, step=0.01)),
                                             column(12, align="left", hr(), actionButton("do_pc", "Start")))),
                                  column(width=8, uiOutput("coverage_UI"),
                                         uiOutput("coverage_plot_UI"))
                                  )
                        ),
                        
                        ###############################################################
                        # 6. Differential analysis
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "diff_analysis",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title="Differential Analysis parameters", width=NULL, status="success", solidHeader=TRUE,
                                             column(12, align="left", textOutput("diff_analysis_info"), br(),
                                                    htmlOutput("selected_DA_GSA_dataset")),
                                             column(8, align="left", uiOutput("selected_k"), br(), br()),
                                             column(5, align="left", selectInput("de_type", "Select type of cluster comparison:",
                                                                                 choices=c("one_vs_rest","pairwise","custom"))),
                                             column(2, offset = 0,align="left"),
                                             column(4, align="left", selectInput("da_method", "Select type of DA method:",
                                                                                 choices=list("Wilcoxon" = "wilcox",
                                                                                              "edgeR GLM" = "neg.binomial"))),
                                             column(6, align="left", offset = 0, uiOutput("name_group")),
                                             column(6, align="left", offset = 0, uiOutput("name_ref")),
                                             column(2, align="left", offset = 0, uiOutput("custom_da_group")),
                                             column(3, align="left", offset = 0, uiOutput("group_choice")),
                                             column(1, align="left", offset = 0, uiOutput("text_vs")),
                                             column(2, align="left", offset = 0, uiOutput("custom_da_ref")),
                                             column(3, align="left", offset = 0, uiOutput("ref_choice")),
                                             column(12, align="left", sliderInput("qval.th", "Adjusted p-value to select significant features:", min=0.01, max=0.4, value=0.01, step=0.01),
                                                    sliderInput("cdiff.th", "Minimum log-fold change to select significant locations:", min=0, max=3, value=1, step=0.01)),
                                                    # plotOutput("distrib_norm_hist") %>%  shinyhelper::helper(type = 'markdown',
                                                    #                                                          icon ="info-circle",
                                                    #                                                          content = "distrib_norm_hist"
                                                    #                                                          )),
                                                    # uiOutput("only_contrib_cell_ui")),
                                             # column(1, align="left", br()),
                                             # column(11, align="left", uiOutput("contrib_thresh")),
                                             # column(1, align="left", br()),
                                             # column(8, align="left", uiOutput("contrib_hist")),
                                             # column(3, align="left", br(), uiOutput("contrib_info")),
                                             column(12, align="left", hr(), actionButton("run_DA", "Start analysis"))),
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
                                         shinydashboard::box(title="Enrichment analysis", width=NULL, status="success", solidHeader=TRUE,
                                             column(12, align="left", htmlOutput("enr_info"), br(),
                                                    uiOutput("use_peaks"),
                                                    actionButton("do_enrich", "Start enrichment analysis"))),
                                         shinydashboard::box(title="Enriched gene sets in differential regions across clusters", width=NULL, status="success", solidHeader=TRUE,
                                             column(4, align="left", uiOutput("GSA_group_sel"), br()),
                                             column(8, align="left", uiOutput("enr_class_sel"), br()),
                                             column(12, align="left",
                                                         mainPanel(tabsetPanel(id='enr_tables',
                                                                          tabPanel("In differential features", div(style = 'overflow-x: scroll', DT::dataTableOutput('all_enrich_table'))),
                                                                          tabPanel("In enriched features", div(style = 'overflow-x: scroll', DT::dataTableOutput('over_enrich_table'))),
                                                                          tabPanel("In depleted features", div(style = 'overflow-x: scroll', DT::dataTableOutput('under_enrich_table')))), width=12),
                                                    br(), br(), br(), br(), br(),
                                                    downloadButton("download_enr_data", "Download tables")))),
                                  column(width=6,
                                         shinydashboard::box(title="Enrichment near TSS", width=NULL, status="success", solidHeader=TRUE,
                                             column(4, align="left", uiOutput("gene_sel")),
                                             column(8, align="left", uiOutput("region_sel")),
                                             column(12, align="left", uiOutput("gene_umap_UI")))))
                        ),
                        
                        ###############################################################
                        # 8. Close app
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName="close_and_save",
                                fluidPage(shinyjs::useShinyjs(),
                                          # extendShinyjs(text = jscode, functions = c("closeWindow")),
                                          column(width=6,
                                                 shinydashboard::box(title='Close App & Save Analysis', solidHeader=TRUE, status='danger', width=NULL,
                                                     column(12, actionButton("close_and_save", "Close App & Save Analysis")))
                                          ),
                                          column(width=6,
                                                 shinydashboard::box(title="Delete analysis", width = NULL, status="success", solidHeader=TRUE,
                                                                     column(9, align="left", uiOutput("selected_delete_analysis")),
                                                                     column(3, align="left", br(), actionButton("delete_analysis", "Delete")),
                                                                     column(12, align="left", textOutput("analysis_deletion_info"))))
                                )
                        )
                        
                      )))
)
