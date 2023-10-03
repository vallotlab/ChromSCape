library(shiny)
library(shinyjs)
library(dplyr)

mytheme <- fresh::create_theme(
    fresh::adminlte_color(
        light_blue = "#434C5E",
        red = "#D81B60",
        green = "#434C5E",
        yellow = "#f39c12",
        aqua = "#39CCCC"
    ),
    fresh::adminlte_sidebar(
        dark_bg = "#D8DEE9",
        dark_hover_bg = "#81A1C1",
        dark_color = "#2E3440"),
    fresh::adminlte_global(
        content_bg = "#FFF",
    )
)

shinyUI(shinydashboardPlus::dashboardPage(freshTheme = mytheme,
                            header = shinydashboardPlus::dashboardHeader(
                                title = tagList(
                                    span(class = "logo-lg", "ChromSCape"), 
                                    icon("ellipsis-h")),
                                leftUi = tagList(uiOutput("dropdown_feature_select_ui"),
                                                 br(),
                                                 br(),
                                                 uiOutput("download_scExp_UI"))
                                ),
                            sidebar = shinydashboardPlus::dashboardSidebar(
                                  shinydashboard::sidebarUserPanel(
                                      "Vallot Lab",
                                               subtitle = a(href = "https://github.com/vallotlab/ChromSCape/", icon("github"), "GitHub"),
                                               image = "LogoCurie.png"
                              ),
                              shinydashboard::sidebarMenu(id="tabs", 
                                          shinydashboard::menuItem("Select & Import", tabName = "select_import", icon=icon("upload")),
                                          shinydashboard::menuItem("Filter & Normalize", tabName = "filter_normalize", icon=icon("fas fa-filter", verify_fa = FALSE)),
                                          shinydashboard::menuItem("Visualize Cells", tabName = "vizualize_dim_red", icon=icon("fas fa-image")),
                                          shinydashboard::menuItem("Cluster Cells", tabName = "cons_clustering", icon=icon("th")),
                                          shinydashboard::menuItem("Coverage & Peak Calling", tabName = "coverage", icon=icon("chart-area")),
                                          shinydashboard::menuItem("Differential Analysis", tabName = "diff_analysis", icon=icon("sort-amount-up")),
                                          shinydashboard::menuItem("Gene Set Analysis", tabName = "enrich_analysis", icon=icon("code-branch")),
                                          shinydashboard::menuItem("TF analysis", tabName = "TF_analysis", icon=icon("bezier-curve")),
                                          shinydashboard::menuItem("Generate Report & Close App", tabName = "close_and_save", icon=icon("close", verify_fa = FALSE)),
                                          shinydashboard::menuItem(
                                            badgeColor = "green",
                                            expandedName = "hello",
                                            startExpanded = TRUE, 
                                            " ", column(width = 12, actionButton("popupUMAP", width = "70%", class = "btn-default",
                                                                                 label = "UMAP", icon = icon("fas fa-image")))
                                          )
                              )
                            ),
                            controlbar = shinydashboardPlus::dashboardControlbar(
                                id = "controlbar",
                                shinydashboardPlus::controlbarMenu(
                                    id = "menu", 
                                    shinydashboardPlus::controlbarItem(
                                        title = "Visualisation Settings",
                                        icon = icon("chart-bar"),
                                        h4("UMAP"), br(),
                                        sliderInput("options.dotplot_downsampling", "Downsampled cells",
                                                    min = 500, max = 100000, step = 10, value = getOption("ChromSCape_options")$dotplot_downsample),
                                        sliderInput("options.dotplot_size", "Point size",
                                                    min = 0.01, max = 5, step = 0.01, value = getOption("ChromSCape_options")$dotplot_size),
                                        sliderInput("options.dotplot_transparency", "Point transparency",
                                                    min = 0.01, max = 1, step = 0.01, value = getOption("ChromSCape_options")$dotplot_transparency),
                                        sliderInput("options.dotplot_max_distanceToTSS", "Max. distance to TSS",
                                                    min = 0, max = 100000, step = 100, value = getOption("ChromSCape_options")$dotplot_max_distanceToTSS),
                                        sliderInput("options.dotplot_min_quantile", "Remove lower range outlier",
                                                    min = 0, max = 1, step = 0.01, value = getOption("ChromSCape_options")$dotplot_min_quantile),
                                        sliderInput("options.dotplot_max_quantile", "Remove higher range outlier",
                                                    min = 0, max = 1, step = 0.01, value = getOption("ChromSCape_options")$dotplot_max_quantile),
                                        h4("Heatmap"), br(),
                                        sliderInput("options.heatmap_downsampling", "Downsampled cells",
                                                    min = 500, max = 20000, step = 10, value = getOption("ChromSCape_options")$heatmap_downsample)
                                    ),
                                    shinydashboardPlus::controlbarItem(
                                        title = "Parallel Settings",
                                        icon = icon("align-left"),
                                    selectInput("options.nb_workers","Select number of workers for parallel computation",
                                                choices= seq_len(future::availableCores()), selected = BiocParallel::bpworkers(BiocParallel::bpparam())),
                                    selectInput("options.bpparam_class","Select parallel computation type",
                                                choices = names(BiocParallel::registered()), selected = 1)
                                    )
                                )
                            ),
                            body = shinydashboard::dashboardBody(
                              shinyjs::useShinyjs(),
                              tags$head(tags$link(rel="shortcut icon", href="ChromSCape_logo.png")),
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
                                    shinyWidgets::chooseSliderSkin("Flat", color = "#112446"),
                                    #Load shinyJS added functions
                                    shinyjs::extendShinyjs(script = file.path(system.file(package="ChromSCape"),"www","shiny_js_functions.js"),
                                                           functions = c("init_directory","save_cookie","disableTab","enableTab")),
                                    #Center Panel
                                    column(width=6,align="left",
                                           shinydashboard::box(title=tagList(shiny::icon("folder"), " Select output directory"), width = NULL, status="success",
                                                                solidHeader=TRUE,
                                                               column(12, align="left",
                                                                      shinyFiles::shinyDirButton("data_folder",  icon = icon("folder-open"), label = "Browse", title = "Upload") %>%
                                                                          shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                                                                              content = "output_directory"), br(),
                                                                      verbatimTextOutput("directory", placeholder = TRUE)),
                                                               column(12, align="left", textOutput("data_folder_info"))),
                                           shinydashboard::box(title=tagList(shiny::icon("upload"), " Import raw data & Create new analysis"), width = NULL, status="success", solidHeader=TRUE,
                                                               column(8, br()) %>%
                                                                 shinyhelper::helper(type = 'markdown', colour = "#E6C717", icon = "question-circle",
                                                                                     content = "HistoneMark_recommandations"), 
                                                               column(12, align="left",
                                                                      column(12, align="left",textInput("new_analysis_name", "Enter a name for the new analysis :", value = "Analysis_1"),
                                                                      selectInput("annotation","Select reference genome:", choices=c("hg38", "mm10", "ce11")),
                                                                      textOutput("data_matrices_info") %>%
                                                                          shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                                                                              content = "datamatrix_input")),
                                                                      column(6,
                                                                      radioButtons("data_choice_box", label = "Input data type",
                                                                                   choices = list("Dense Matrix"="DenseMatrix",
                                                                                                  "Sparse Matrix" = "SparseMatrix",
                                                                                                  "Fragment File" = "FragmentFile",
                                                                                                  "Single-cell BAM files" = "scBAM",
                                                                                                  "Single-cell BED files" = "scBED"
                                                                                                  ))),
                                                                      column(6,
                                                                      uiOutput("input_data_ui")),
                                                                      column(12,
                                                                      uiOutput("advanced_data_input")),
                                                                      column(3, uiOutput("rebin_matrices_checkbox_ui") %>%
                                                                                 shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                                                                                     content = "rebin_matrices")) ,
                                                                      column(12,uiOutput("rebin_matrices_ui")),
                                                                      column(3, actionButton("create_analysis", "Create analysis")),
                                                                      column(4, uiOutput("add_to_current_analysis_checkbox_UI"))
                                                                      )) 
                                    ),
                                    column(6,
                                           column(11, align="left",
                                                  shinydashboard::box(title=tagList(shiny::icon("file"), " Current analysis"), width = NULL, status="success", solidHeader=TRUE,
                                                                      column(12, align="left",
                                                                             uiOutput("selected_analysis")
                                                                      )),
                                                  uiOutput("rename_file_box")),
                                           column(1, actionButton("startHelp","Help", icon = icon("question-circle")), br()),
                                    )
                                )
                        ),
                        
                        ###############################################################
                        # 2. Filter and Normalize dataset
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "filter_normalize",
                                                
                                                fluidPage(
                                                    #Right Panel
                                                    column(width=7,
                                                           shinydashboard::box(title = tagList(shiny::icon("fas fa-filter", verify_fa = FALSE ), " Filter & Normalize"),
                                                                               width = NULL, status="success", solidHeader=TRUE,
                                                                               column(6, align="left",
                                                                                      column(12, align = "middle", h4("Cell coverage")),
                                                                                      br(), br(),
                                                                                      plotly::plotlyOutput("cell_coverage", height=250) %>%
                                                                                          shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75),
                                                                                      br(), hr()),
                                                                               column(6, align="left",
                                                                                      column(12, align = "middle", h4("Feature coverage")),
                                                                                      br(), br(),
                                                                                      plotly::plotlyOutput("feature_coverage", height=250) %>%
                                                                                          shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75),
                                                                                      br(), hr()),
                                                                               column(12, align="left",
                                                                                      uiOutput("table_QC_filt_box"),
                                                                                      br(), hr(),
                                                                                      uiOutput('min_coverage_cell_ui'),
                                                                                      uiOutput('quant_removal_ui'),
                                                                                      uiOutput('n_feature_UI'),
                                                                                      selectInput("norm_type", "Normalization method: ", choices = c("TFIDF", "CPM"), width = "20%"),
                                                                                      uiOutput("remove_PC_UI"),
                                                                                      checkboxInput("run_tsne", "Run T-SNE", value= FALSE) %>%
                                                                                          shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                                                                                              content = "run_tsne"),
                                                                                      checkboxInput("do_subsample", "Perform subsampling", value=FALSE)%>%
                                                                                          shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                                                                                              content = "subsampling"),
                                                                                      uiOutput("do_subsample"),
                                                                                      checkboxInput("exclude_regions", "Exclude specific genomic regions", value=FALSE)  %>%
                                                                                          shinyhelper::helper(type = 'markdown',  colour = "#434C5E", icon ="info-circle",
                                                                                                              content = "exclude_region"),
                                                                                      uiOutput("exclude_file"),
                                                                                      checkboxInput("do_batch_corr", "Perform batch correction", value=FALSE) %>%
                                                                                          shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                                                                                              content = "batch_correction"),
                                                                                      uiOutput("num_batches")),
                                                                               column(5, align="left",
                                                                                      uiOutput("batch_names")),
                                                                               column(7, align="left",
                                                                                      htmlOutput("batch_sel")),
                                                                               column(12, align="left",
                                                                                      hr(),
                                                                                      actionButton("filter_normalize_reduce", "Filter, Normalize & Reduce")))),
                                                    column(width = 5,
                                                           shinydashboard::box(title=tagList(shiny::icon("file"), " Select filtered & normalized dataset"), width = NULL, status="success", solidHeader=TRUE,
                                                                               column(12, align="left",
                                                                                      htmlOutput("selected_reduced_dataset"),
                                                                                      textOutput("red_data_selection_info"),
                                                                                      textOutput("red_data_selection_format")))
                                                           # uiOutput("cell_feature_coverage_after")
                                                           )
                                                    )
                        ),
                        
                        
                        
                        ###############################################################
                        # 3. Vizualising cells in reduced dimension
                        ###############################################################
                        
                        shinydashboard::tabItem(
                          tabName = "vizualize_dim_red",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title = tagList(shiny::icon("fas fa-image"), " PCA"), width = NULL, status="success", solidHeader=TRUE,
                                             column(6, align="left", uiOutput("feature_color")),
                                             column(12, align="left",plotOutput("pca_plot") %>%
                                                        shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75) %>%
                                                        shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                                                            content = "pca_plot")
                                                        ),
                                             column(3, align="left", uiOutput("pc_select_x")),
                                             column(3, align="left", uiOutput("pc_select_y"))),
                                         uiOutput("contribution_to_pca_UI"),
                                         uiOutput("correlation_to_pca_UI")),
                                  column(width=6,
                                         shinydashboard::box(title=tagList(shiny::icon("fas fa-image"), " UMAP"), width = NULL, status="success", solidHeader=TRUE,
                                                             column(12, align="left", plotOutput("UMAP_plot") %>%
                                                                        shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75) %>%
                                                                        shinyhelper::helper(type = 'markdown',  colour = "#434C5E", icon ="info-circle",
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
                                         shinydashboard::box(title=tagList(shiny::icon("cube")," Hierarchical Clustering"),
                                                             width=NULL, status="success", solidHeader=TRUE,
                                                             column(12, align ="center", plotOutput("hc_heatmap_plot", height=500, width=500) %>%
                                                                        shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75) %>%
                                                                        shinyhelper::helper(type = 'markdown',  colour = "#434C5E", icon ="info-circle",
                                                                                            content = "correlation_clustering"))
                                         ),
                                         uiOutput("consensus_clustering_box"),
                                         uiOutput("correlation_filtering_box")
                                  ),
   
                                  column(width=6,
                                         shinydashboard::box(title=tagList(shiny::icon("project-diagram"), "Cluster Cells"), width=NULL, status="success", solidHeader=TRUE,
                                                             align="left",
                                                             column(width=4, uiOutput("clustering_method_UI")),
                                                             column(width=4, uiOutput("clustering_UI")),
                                                             column(width=3, br(),br(), uiOutput("consensus_clustering_UI")),
                                                             column(width=3, br(),br(),actionButton(inputId = "choose_cluster", label = "Launch clustering"))
                                         ),
                                         uiOutput("contingency_table_cluster"),
                                         uiOutput("UMAP_box"),
                                         uiOutput("IDclust_box"),
                                         uiOutput("tsne_box_cf"),
                                         uiOutput("color_box_cf"),
                                         shinydashboard::box(title=tagList(shiny::icon("cube")," Hierarchical Clustering with cluster annotation"), width=NULL,
                                                             status="success", solidHeader=TRUE, align="left", 
                                                             column(width=4,
                                                                    actionButton(inputId = "do_annotated_heatmap_plot", label = "Plot Clustered Heatmap"), br()),
                                                             column(width=12, uiOutput("annotated_heatmap_UI"))
                                         ),
                                         shinydashboard::box(title=tagList(shiny::icon("ruler-combined"), " Intra / Inter correlation"), width=NULL,
                                                             status="success", solidHeader=TRUE, align="left",
                                                             collapsible = TRUE, collapsed = TRUE,
                                                             column(12, align="left",
                                                                    column(width=4,uiOutput("violin_color")),
                                                                    column(width=4,checkboxInput("add_jitter", shiny::HTML("<b>Add single-cells</b>"), value= FALSE),
                                                                    column(3, align = "left", actionButton(inputId = "save_plots_violins",
                                                                                                                  label = "Save HQ plots",
                                                                                                                  icon = icon("fas fa-image", verify_fa = FALSE))) %>%
                                                                               shinyhelper::helper(type = 'markdown',  colour = "#434C5E", icon ="info-circle",
                                                                                                   content = "intra_inter_correlation")),
                                                                    column(width=4,uiOutput("jitter_color")),  br(),br(),br(),br(),br(),
                                                                    mainPanel(tabsetPanel(id='inter_intra_cor',
                                                                                          
                                                                                          tabPanel("Intracorrelation", column(width=12, uiOutput("intra_corr_UI") 
                                                                                                                              )),
                                                                                          tabPanel("Intercorrelation", column(width=12, column(10, uiOutput("reference_group")),
                                                                                                                              br(), br(), br(), br(), br(),
                                                                                                                              uiOutput("inter_corr_UI")))
                                                                                          ))),
                                                             
                                         )
                                  )
                                )
                        ),
                        ###############################################################
                        # 5. Coverages [optional]
                        ###############################################################
                        shinydashboard::tabItem(tabName = "coverage",
                                                fluidPage(
                                                    column(width=4,
                                                           shinydashboard::box(title=tagList(shiny::icon("chart-area"), " Generate Coverage Plots"), width=NULL, status="success", solidHeader=TRUE,
                                                                               column(12, align="left", htmlOutput("coverage_info"), hr()),
                                                                               tags$style(HTML(".large_icon { font-size: 70px; }")),
                                                                               
                                                                               column(12, align="left",
                                                                                      uiOutput("coverage_folder_ui"),
                                                                                      br(),hr(),br(),
                                                                                      uiOutput("coverage_upload")),
                                                                               column(12, align="left", sliderInput(
                                                                                   inputId = "quantile_for_peak_calling", label = "Quantile of signal to define peak",
                                                                                   min = 0, max = 1, value = 0.85, step = 0.005)),
                                                                               column(4, align="left", textInput(
                                                                                   inputId = "coverage_bin_size", label = "Resolution", value = 150)),
                                                                               column(4, align="left", textInput(
                                                                                   inputId = "coverage_n_smoothBin", label = "Smoothing", value = 5)),
                                                                               column(12, align="left", actionButton("do_coverage", "Create coverage & Call Peaks")))),
                                                    column(width=8, uiOutput("coverage_UI"),
                                                           uiOutput("coverage_plot_UI"))
                                                )
                        ),
                        
                        ###############################################################
                        # 6. Peak calling [optional]
                        ###############################################################
                        
                        # shinydashboard::tabItem(tabName = "peak_calling",
                        #         fluidPage(
                        #           column(width=6,
                        #                  shinydashboard::box(title=tagList(shiny::icon("fab fa-mountain", verify_fa = FALSE), " Peak calling"), width=NULL, status="success", solidHeader=TRUE,
                        #                      column(12, align="left", textOutput("peak_calling_info"), hr()),
                        #                      tags$style(HTML(".large_icon { font-size: 70px; }")),
                        #                      column(5, align="left",
                        #                             htmlOutput("peak_calling_system"), hr()),
                        #                      column(5,offset = 1,align="center", htmlOutput("peak_calling_icon")),
                        #                      column(12, align="left",
                        #                             sliderInput("peak_distance_to_merge", "Select distance of peaks to merge:", min=0, max=50000, value=5000, step=100),
                        #                             shinyFiles::shinyDirButton("pc_folder", "Browse folder of BAM / scBED files" ,
                        #                                                        icon = icon("folder-open"),
                        #                                                        title = "Please select a folder:",
                        #                                                        buttonType = "default", class = NULL),
                        #                             uiOutput("pc_upload")),
                        #                      br(), br(),
                        #                      column(4, align="left", textOutput("pc_k_selection"),
                        #                             selectInput("pc_stat","Select statistic for cutoff:", choices=c("p.value", "q.value"), selected="p.value")),
                        #                      column(12, align="left", br(), br(),
                        #                             sliderInput("pc_stat_value", "Select significance threshold:", min=0, max=0.25, value=0.05, step=0.01)),
                        #                      column(12, align="left", hr(), actionButton("do_pc", "Call peaks / Create Consensus Peaks"))))
                        #           )
                        # ),
                        
                        ###############################################################
                        # 7. Differential analysis
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "diff_analysis",
                                fluidPage(
                                  column(width=5,
                                         shinydashboard::box(title= tagList(shiny::icon("sort-amount-up"), " Run Differential Analysis"), width=NULL, status="success", solidHeader=TRUE,
                                             column(12, align="left", textOutput("diff_analysis_info") %>%
                                                        shinyhelper::helper(type = 'markdown',  colour = "#434C5E", icon ="info-circle",
                                                                            content = "differential_analysis"), br()),
                                             column(8, align="left", uiOutput("selected_k"), br(), br()),
                                             column(5, align="left", selectInput("de_type", "Select type of comparison:",
                                                                                 choices=list("One vs Rest (activation)" = "one_vs_rest_fast",
                                                                                              "One vs Rest (classic)" = "one_vs_rest",
                                                                                               "Pairwise"  = "pairwise",
                                                                                               "Custom" = "custom"))),
                                             column(2, offset = 0,align="left"),
                                             column(4, align="left", selectInput("da_method", "Select method:",
                                                                                 choices=list("Wilcoxon" = "wilcox",
                                                                                              "edgeR GLM" = "neg.binomial"))),
                                             column(6, align="left", offset = 0, uiOutput("name_group")),
                                             column(6, align="left", offset = 0, uiOutput("name_ref")),
                                             column(3, align="left", offset = 0, uiOutput("da_group")),
                                             column(3, align="left", offset = 0, uiOutput("group_choice")),
                                             column(1, align="left", offset = 0, uiOutput("text_vs")),
                                             column(2, align="left", offset = 0, uiOutput("custom_da_ref")),
                                             column(3, align="left", offset = 0, uiOutput("ref_choice")),
                                             column(12, align="left", hr(), actionButton("run_DA", "Start analysis"))),
                                         uiOutput("da_summary_box")),
                                  column(width=7,
                                         shinydashboard::box(title= tagList(shiny::icon("bullseye"), " Find Significantly Differential Features"), width=NULL, status="success", solidHeader=TRUE,
                                                             column(8, align="left", htmlOutput("selected_DA_GSA_dataset")),
                                                             column(8, align="left",
                                                shinyWidgets::sliderTextInput(inputId = "qval.th", label = "Adjusted p-value to select significant features:",
                                                                                                choices =  c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.04, 0.05, 0.1, 0.15, 0.20, 0.25),
                                                                                                selected = 0.01, grid = TRUE), 
                                                sliderInput("logFC.th", "Minimum log-fold change to select significant locations:", min=0, max=3, value=1, step=0.01),
                                                sliderInput("min.percent", "Minimum fraction of active cells in feature:", min=0, max=0.5, value=0.01, step=0.005)),
                                                column(4, align="center", br(), br(), br(),br(), br(), br(),
                                                       shiny::actionButton("apply_DA_filters", label = "Apply filters"))),
                                         uiOutput("da_visu_box")))
                        ),
                        
                        ###############################################################
                        # 8. Enrichment analysis
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "enrich_analysis",
                                fluidPage(
                                  column(width=6,
                                         shinydashboard::box(title=tagList(shiny::icon("code-branch"), " Gene Set Analysis"), width=NULL, status="success", solidHeader=TRUE,
                                             column(12, align="left", htmlOutput("enr_info"), br(),
                                                    uiOutput("use_peaks"),
                                                    actionButton("do_enrich", "Start enrichment analysis"))),
                                         shinydashboard::box(title=tagList(icon("th"), " Enriched Gene sets in differential features"), width=NULL, status="success", solidHeader=TRUE,
                                             column(4, align="left", uiOutput("GSA_group_sel"), br()),
                                             column(8, align="left", uiOutput("enr_class_sel"), br()),
                                             column(8, align="left", sliderInput("pathway_size_GSA", label = "Select pathway size",
                                                                                 min = 1, max = 2500, value = c(25, 1000), step = 1), br()),
                                             column(12, align="left",
                                                         mainPanel(tabsetPanel(id='enr_tables',
                                                                          tabPanel("In differential features", div(style = 'overflow-x: scroll', DT::dataTableOutput('all_enrich_table'))),
                                                                          tabPanel("In enriched features", div(style = 'overflow-x: scroll', DT::dataTableOutput('over_enrich_table'))),
                                                                          tabPanel("In depleted features", div(style = 'overflow-x: scroll', DT::dataTableOutput('under_enrich_table')))), width=12),
                                                    br(), br(), br(), br(), br(),
                                                    downloadButton("download_enr_data", "Download tables")))),
                                  column(width=6,
                                         shinydashboard::box(title=tagList(icon("image"), " UMAP of enrichment near TSS"), width=NULL, status="success", solidHeader=TRUE,
                                             column(3, align="left", uiOutput("gene_sel")),
                                             column(3, align="left", checkboxInput("label_cluster_umap_GSA", "Label cluster", value = T) ),
                                             column(3, align="left", selectInput("color_by_violin_GSA", "Color by", choices = c("cell_cluster","sample_id")) ),
                                             column(3, align="left", actionButton("save_plot_GSA", "Save HQ plot")),
                                             column(12, align="left", uiOutput("gene_umap_UI")), br(),
                                             column(12, align="left", uiOutput("gene_violin_UI"))
                                             ),
                                         shinydashboard::box(title=tagList(icon("image"), " Enrichment in Gene Sets"), width=NULL, status="success", solidHeader=TRUE,
                                                             column(4, align="left", actionButton("plot_pathways", "Plot pathways")), 
                                                             column(8, align="left", uiOutput("pathways_sel")), 
                                                             column(12, align="left",uiOutput("pathways_umap_UI"))
                                         )
                                         ))
                        ),
                        
                        ###############################################################
                        # 8. TF analysis
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName = "TF_analysis",
                                                fluidPage(
                                                  column(width=6,
                                                         shinydashboard::box(title=tagList(shiny::icon("bezier-curve"), " Transcription Factor Analysis using ChEA3"), width=NULL, status="success", solidHeader=TRUE,
                                                                             column(12, align="left", htmlOutput("TF_info"), br(),
                                                                                    uiOutput("use_peaks_TF"),
                                                                                    actionButton("do_enrich_TF", "Start TF enrichment analysis"))),
                                                         shinydashboard::box(title=tagList(icon("th"), " Enriched TF in differential features"), width=NULL, status="success", solidHeader=TRUE,
                                                                             column(4, align="left", uiOutput("TF_group_sel"), br()),
                                                                             column(12, align="left",
                                                                                    mainPanel(tabsetPanel(id='enr_tables_TF',
                                                                                                          tabPanel("In differential features", div(style = 'overflow-x: scroll', DT::dataTableOutput('all_enrich_table_TF'))),
                                                                                                          tabPanel("In enriched features", div(style = 'overflow-x: scroll', DT::dataTableOutput('over_enrich_table_TF'))),
                                                                                                          tabPanel("In depleted features", div(style = 'overflow-x: scroll', DT::dataTableOutput('under_enrich_table_TF')))), width=12),
                                                                                    br(), br(), br(), br(), br(),
                                                                                    downloadButton("download_enr_data_TF", "Download tables"))),
                                                         ),
                                                  column(width=6, uiOutput("TF_results_UI"))
                                                  
                                                  )
                        ),
                        
                        ###############################################################
                        # 9. Close app
                        ###############################################################
                        
                        shinydashboard::tabItem(tabName="close_and_save",
                                fluidPage(shinyjs::useShinyjs(),
                                          # extendShinyjs(text = jscode, functions = c("closeWindow")),
                                          column(width=6,
                                                 shinydashboard::box(title='Close App & Save Analysis', solidHeader=TRUE, status='danger', width=NULL,
                                                     column(12, actionButton("generate_report", "Generate HTML Report")), br(), br(),
                                                     column(12, actionButton("close_and_save", "Close App & Save Analysis"))
                                                     )
                                          ),
                                          column(width=6,
                                                 shinydashboard::box(title="Delete analysis", width = NULL, status="success", solidHeader=TRUE,
                                                                     column(9, align="left", uiOutput("selected_delete_analysis")),
                                                                     column(3, align="left", br(), actionButton("delete_analysis", "Delete")),
                                                                     column(12, align="left", textOutput("analysis_deletion_info")))),
                                )
                        )
                        
                      )))
)
