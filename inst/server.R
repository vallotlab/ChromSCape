
shinyServer(function(input, output, session) {
  
  library(ggplot2)
  library(dplyr)
  
  options(shiny.maxRequestSize = 5000*1024^2) # allow upload of files with max 5GB
  
  ###############################################################
  # 0. Global variables and functions
  ###############################################################
  
  #Initializating user experience functions
  js$init_directory() #Getting cookie for the directory
  volumes = c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
shinyhelper::observe_helpers(help_dir = "www/helpfiles",withMathJax = TRUE)
  
  # addResourcePath("www", system.file("www", package="ChromSCape"))
  tab_vector = c("filter_normalize",
                 "vizualize_dim_red",
                 "cons_clustering",
                 "peak_calling",
                 "diff_analysis",
                 "enrich_analysis") #list of all lockable tabs
  unlocked = reactiveValues(list = list(selected_analysis = FALSE,
                                        selected_reduced_dataset = FALSE,
                                        pca = FALSE,
                                        affectation = FALSE,
                                        diff_my_res = FALSE)) #list of all required items to unlock a tab
  for(tab in tab_vector){
    js$disableTab(tab) #Disabling all tabs but the first one
  }
  
  observeEvent(input$startHelp,{
    print("Started help")
    # on click, send custom message to start help
    session$sendCustomMessage(type = 'startHelp', message = list(""))
  })
  
  #Global reactives values
  scExp = reactiveVal(NULL)
  scExp_cf = reactiveVal(NULL)
  
  cell_cov_df <- reactive ({data.frame(coverage = sort(unname(Matrix::colSums(init$datamatrix)))) })  # used for plotting cell coverage on first page
  analysis_name <- reactive({ input$selected_analysis })
  annotation_id_norm <- reactive({ read.table(file.path(init$data_folder, 'ChromSCape_analyses', input$selected_analysis, 'annotation.txt'), header = FALSE, stringsAsFactors = FALSE)[[1]] })
  annotation_id <- reactive({ read.table(file.path(init$data_folder, 'ChromSCape_analyses', analysis_name(), 'annotation.txt'), header = FALSE, stringsAsFactors = FALSE)[[1]] })
  
  #Global Functions
  init <- reactiveValues(data_folder =  getwd(), datamatrix = data.frame(), annot_raw = data.frame(), available_raw_datasets = NULL,
                         available_reduced_datasets = NULL)
  reduced_datasets <- reactive({ 
    if (is.null(init$available_reduced_datasets)) c() else gsub('.{6}$', '', basename(init$available_reduced_datasets)) })
  
  observeEvent({analysis_name()},{
    init$available_reduced_datasets = get.available.reduced.datasets(analysis_name())
  })
  annotCol <- reactive({ c("sample_id","total_counts","batch_id") })
  
  
  observeEvent(analysis_name(), { # application header (tells you which data set is selected)
    req(analysis_name())
    header <- paste0('<b>Analysis : ', analysis_name(), ' </b>')
    shinyjs::html("pageHeader", header)
  })
  
  get.available.reduced.datasets <- function(selected_analysis){
    list.files(path = file.path(init$data_folder, "ChromSCape_analyses", selected_analysis,"Filtering_Normalize_Reduce"), full.names = FALSE, recursive = TRUE,
               pattern="[[:print:]]+_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+_(uncorrected|batchCorrected).RData")
  }
  
  get.available.filtered.datasets <- function(name, preproc){
    list.files(path = file.path(init$data_folder, "ChromSCape_analyses", name, "correlation_clustering"), full.names = FALSE, recursive = FALSE, pattern = paste0(preproc, "_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?.RData"))
  }
  
  able_disable_tab <- function(variables_to_check, tab_id) {

    able_or_disable = c()
    for(var in variables_to_check){
      if (unlocked$list[[var]]==T) {
        able_or_disable = c(able_or_disable,TRUE)
      } else{
        able_or_disable = c(able_or_disable,FALSE)
      }}
    for (tab in tab_id) {
      if (FALSE %in% able_or_disable) {
        js$disableTab(tab_id)
      }
      else{
        js$enableTab(tab)
      }}
  }
  
  batchUsed <- reactive({ grepl("batchCorrected", input$selected_reduced_dataset) })
  
  ###############################################################
  # 1. Select Analysis & Import dataset
  ###############################################################
  
  output$selected_analysis <- renderUI({ 
    selectInput("selected_analysis", "Select an Analysis:",
                choices = init$available_raw_datasets, multiple = F) 
    })
  
  output$data_folder_info <- renderText({
    "All your analyses will be saved in this folder."
    })

  output$data_matrices_info <- renderText({"The filename of each selected matrix must be the sample name, e.g. T23_K4me3.txt (Must contain only alpha-numeric charachter and underscore !)"})
  
  
  shinyFiles::shinyDirChoose(
    input, "data_folder", roots = volumes,
    session = session, restrictions = system.file(package = "base")
  )
    
  directory <- reactive(input$data_folder)
  output$directory <- renderText({
    init$data_folder
  })
  
  #Look for existing cookie
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$path_cookie
    },
    handlerExpr = {
       if(.Platform$OS.type != "windows"){
         if ( (input$path_cookie != "[null]") && !is.null(input$path_cookie) && !is.na(input$path_cookie)) {
           #Uploading the name displayed in Data Folder
           
           init$data_folder <- gsub(pattern = "\"|\\[|\\]|\\\\", "",
                                    as.character(input$path_cookie))
           
           init$available_raw_datasets <- list.dirs(path = file.path(init$data_folder, "ChromSCape_analyses"), full.names = FALSE, recursive = FALSE)
           init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())
         }
       }
      })
  
  #Selecting a working directory using shinyDirectoryInput::readDirectoryInput(input$data_folder) and saving cookie
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$data_folder  
    },
    handlerExpr = {
      if (!"path" %in% names(directory())) return()

      init$data_folder <- shinyFiles::parseDirPath(volumes, directory())
      
      init$available_raw_datasets <- list.dirs(
        path = file.path(init$data_folder, "ChromSCape_analyses"),
        full.names = FALSE, recursive = FALSE)
      init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())

      if(.Platform$OS.type != "windows"){
        js$save_cookie(init$data_folder)
      }
    }
  )
  
  observeEvent(input$create_analysis, {  # save new dataset
    req(input$new_analysis_name, input$annotation, input$datafile_matrix)
    datamatrix <- NULL
    annot_raw <- NULL
    withProgress(message='Compiling new data set...',value = 0, {
      dir.create(file.path(init$data_folder, "ChromSCape_analyses"), showWarnings = FALSE)
      dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name))
      dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "Filtering_Normalize_Reduce"))
      dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "correlation_clustering"))
      dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "correlation_clustering","Plots"))
      dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "diff_analysis_GSEA"))
      write.table(input$annotation, file.path(init$data_folder, 'ChromSCape_analyses', input$new_analysis_name, 'annotation.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
      incProgress(0.3, detail="reading data matrices")
  
      tmp_list = import_scExp(file_names = input$datafile_matrix$name,
                   path_to_matrix = input$datafile_matrix$datapath)
      datamatrix = tmp_list$datamatrix
      annot_raw = tmp_list$annot_raw
      incProgress(0.7, detail="saving raw data")
      save(datamatrix, annot_raw, file = file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "scChIP_raw.RData"))
      init$available_raw_datasets <- list.dirs(path = file.path(init$data_folder, "ChromSCape_analyses"), full.names = FALSE, recursive = FALSE)
      init$datamatrix <- datamatrix
      init$annot_raw <- annot_raw
      
      updateActionButton(session, "create_analysis", label="Added successfully", icon = icon("check-circle"))
    })
  })
  
  observeEvent(input$selected_analysis, {  # load precompiled dataset and update coverage plot
    if(!is.null(input$selected_analysis) & input$selected_analysis != ""){
      myData = new.env()
      load(file.path(init$data_folder,"ChromSCape_analyses", input$selected_analysis, "scChIP_raw.RData"), envir = myData)
      init$datamatrix <- myData$datamatrix
      init$annot_raw <- myData$annot_raw
          
      }
  })
  
  # observeEvent(input$selected_analysis,{
  #   req(input$selected_analysis)
  #   
  #   if(!is.null(input$selected_reduced_dataset)){
  #     delay(1500, {
  #       if(gsub(pattern ="_\\d*_\\d*_\\d*_\\w*","",input$selected_reduced_dataset) != input$selected_analysis){
  #         
  #         showNotification(paste0("Warning : Selected Analysis '",input$selected_analysis,
  #                                 "' is different from selected reduced dataset '", input$selected_reduced_dataset,"'"),
  #                          duration = 5, closeButton = TRUE, type="warning")
  #       }
  #     })
  #   }
  # })
  
  observeEvent(input$new_analysis_name, {  # reset label on actionButtion when new dataset is added
    updateActionButton(session, "create_analysis", label="Create Analysis", icon = character(0))
  })
  
  observeEvent(input$selected_analysis,{if(nchar(input$selected_analysis)>0){unlocked$list$selected_analysis=T}else{for(i in names(unlocked$list)){unlocked$list[[i]]=F}}})
  observeEvent(unlocked$list,{able_disable_tab("selected_analysis","filter_normalize")}) 
  
  ###############################################################
  # 2. Filter and Normalize dataset
  ###############################################################
  
  output$selected_reduced_dataset <- renderUI({ 
    selectInput("selected_reduced_dataset", "Select filtered & normalized set :",
                choices = reduced_datasets()) 
    })
  output$red_data_selection_info <- renderText({"The selected data set is automatically loaded and will be used for all subsequent analysis. 
    If you try different filtering parameters for one analysis, you can select the results using each parameter set here."})
  output$red_data_selection_format <- renderText({"The name of the filtered & normalized 
    dataset is composed of the following information: analysis name, upper percentile
    of cells to remove (potential doublets), min percentage of cells to support a window,
    quantile of cell read counts to keep and batch correction type."})
  
  output$exclude_file <- renderUI({ if(input$exclude_regions){
    fileInput("exclude_file", ".bed file containing the regions to exclude from data set:", multiple = FALSE, accept = c(".bed",".txt"))
  }})
  output$num_batches <- renderUI({ if(input$do_batch_corr){
    selectInput("num_batches","Select number of batches (each can have one or multiple samples):", choices=c(2:10))
  }})
  output$batch_names <- renderUI({ if(input$do_batch_corr & !is.null(input$num_batches)){
    lapply(1:input$num_batches, function(i){
      textInput(paste0('batch_name_', i), paste0('Batch name ', i, ':'), value = paste0('batch', i))
    })
  }})
  batch_choice <- reactive({ unique(init$annot_raw$sample_id) })
  output$batch_sel <- renderUI({ if(input$do_batch_corr & dim(init$annot_raw)[1] > 0 & !is.null(input$num_batches)){
    lapply(1:input$num_batches, function(i){
      selectInput(paste0('batch_sel_', i), paste0('Select samples for batch ', i, ':'), choices=batch_choice(), multiple=TRUE)
    })
  }})
  output$do_subsample <- renderUI({ if(input$do_subsample){
    sliderInput("subsample_n", "Select number of cells to subsample for each sample:", min=100, max=5000, value=500, step=10) }})
  
  observeEvent(input$filter_normalize_reduce, {  # perform QC filtering and dim. reduction
    num_batches <- if(is.null(input$num_batches)) 0 else input$num_batches
    batch_names <- if(is.null(input$num_batches)) c() else sapply(1:input$num_batches, function(i){ input[[paste0('batch_name_', i)]] })
    batch_sels <- if(is.null(input$num_batches)) list() else lapply(1:input$num_batches, function(i){ input[[paste0('batch_sel_', i)]] })
    names(batch_sels) = batch_names
    subsample_n <- if(input$do_subsample){input$subsample_n}
    
    annotationId <- annotation_id_norm()
    exclude_regions <- if(input$exclude_regions) setNames(read.table(
      input$exclude_file$datapath, header = FALSE, stringsAsFactors = FALSE), c("chr", "start", "stop")) else NULL
  
    callModule(Module_preprocessing_filtering_and_reduction, "Module_preprocessing_filtering_and_reduction", reactive({input$selected_analysis}), reactive({input$min_coverage_cell}),
               reactive({input$min_cells_window}), reactive({input$quant_removal}), reactive({init$datamatrix}), reactive({init$annot_raw}),
               reactive({init$data_folder}),reactive({annotationId}), reactive({exclude_regions}), reactive({annotCol()}),reactive({input$do_batch_corr}),
                reactive({batch_sels}), reactive({input$run_tsne}), reactive({subsample_n}))
    init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())
    updateActionButton(session, "filter_normalize_reduce", label="Processed and saved successfully", icon = icon("check-circle"))
    
  })

  observeEvent({input$selected_analysis  # reset label on actionButtion when new filtering should be filtered
   input$min_coverage_cell
   input$quant_removal
   input$min_cells_window
   input$do_batch_corr}, {
   updateActionButton(session, "filter_normalize_reduce", label="Filter, Normalize & Reduce", icon = character(0))
  })
  
  reduced_dataset <- eventReactive(input$selected_reduced_dataset, { # load reduced data set to work with on next pages
    req(input$selected_reduced_dataset)
    if(file.exists(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots"))) 
      addResourcePath('Plots', file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots"))

    file_index <- match(c(input$selected_reduced_dataset), reduced_datasets())
    filename_sel <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),"Filtering_Normalize_Reduce",init$available_reduced_datasets[file_index])

    myData = new.env()
    load(filename_sel, envir = myData)
    if(is.reactive(myData$scExp)) myData$scExp = isolate(myData$scExp())
    myData
    
  })

  observeEvent(input$selected_reduced_dataset, {  # load scExp, add colors, add correlation
    req(input$selected_reduced_dataset, reduced_dataset())
    # Retrieve the scExp filtered
    scExp. = reduced_dataset()$scExp # retrieve filtered scExp
    scExp. = correlation_and_hierarchical_clust_scExp(scExp.)
    scExp(scExp.)
    rm(scExp.)
  })

  quantile_threshold =  reactive({
    index = round(ncol(init$datamatrix) * as.numeric(input$quant_removal) * 0.01)
    q = cell_cov_df()$coverage[index]
    q
  })
  
  cell_cov_plot <- reactive({
    ggplot(cell_cov_df(), aes(x = coverage)) + 
      geom_histogram(color="black", fill="steelblue", bins = input$coverage_bins) +
      labs(x="Log10(Reads per cell)")  + 
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), 
            panel.background=element_blank(), axis.line=element_line(colour="black"),
            panel.border=element_rect(colour="black", fill=NA)) +
      geom_vline(xintercept = as.numeric(input$min_coverage_cell), color = "#22AD18") + 
      geom_vline(xintercept = quantile_threshold(), color = "#D61111") +
      scale_x_log10()
      
    })
  
  output$cell_coverage <- plotly::renderPlotly( plotly::ggplotly(cell_cov_plot(), 
                                                                 tooltip="Sample", dynamicTicks=T) )
  
  
  
  output$num_cell_after_QC_filt <- function(){
    req(reduced_dataset(), scExp())
    num_cell_after_QC_filt_scExp(scExp(),init$annot_raw)
  }
  
  output$table_QC_filt_box <- renderUI({
    if(!is.null(reduced_dataset())){
          column(12, align="left", tableOutput("num_cell_after_QC_filt"))
    }
  })
  
  observeEvent(input$selected_reduced_dataset,{
    if(suppressWarnings(nchar(input$selected_reduced_dataset)>0)){
      unlocked$list$selected_reduced_dataset=T
      }else{
        for(i in names(unlocked$list)){unlocked$list[[i]]=F}
        }
    })
  observeEvent(unlocked$list,{able_disable_tab(c("selected_reduced_dataset"),"vizualize_dim_red")}) 
  
  
  ###############################################################
  # 2. PCA
  ###############################################################
  
  output$pc_select_x <- renderUI({ selectInput("pc_select_x", "X",choices=paste0("Component_", c(1:15)), selected="Component_1") })
  output$pc_select_y <- renderUI({ selectInput("pc_select_y", "Y",choices=paste0("Component_", c(1:15)), selected="Component_2") })
  output$color_by <- renderUI({selectInput("color_by", "Color by", choices=annotCol()) })

  pca_plot <- reactive({
    req(scExp(), input$pc_select_x,input$pc_select_y,  input$color_by)

    p = plot_reduced_dim_scExp(scExp(),input$color_by, "PCA",
                               select_x = input$pc_select_x,
                               select_y = input$pc_select_y
    )
    unlocked$list$pca=T
    p
  })
  output$pca_plot <- plotly::renderPlotly( plotly::ggplotly(pca_plot(), tooltip="Sample", dynamicTicks=T) )
  
  output$tsne_box <- renderUI({
    req(scExp(), input$color_by)
    if("TSNE" %in% SingleCellExperiment::reducedDimNames(scExp())){
      p = plot_reduced_dim_scExp(scExp(),input$color_by, "TSNE")
      output$tsne_plot = plotly::renderPlotly( plotly::ggplotly(p, tooltip="Sample", dynamicTicks=T) )
      shinydashboard::box(title="tSNE visualization", width = NULL, status="success", solidHeader=T,
                          column(12, align="left", plotly::plotlyOutput("tsne_plot") %>%
                                   shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                       content = "tsne_plot")))
    }
  })
  
  umap_plot <- reactive({
    req(scExp(), input$color_by)
    p = plot_reduced_dim_scExp(scExp(),input$color_by, "UMAP"
    )
    p
  })
  output$umap_plot <- plotly::renderPlotly( plotly::ggplotly(umap_plot(), tooltip="Sample", dynamicTicks=T) )
  
  output$color_box <- renderUI({
    req(input$color_by)
    if(input$color_by != 'total_counts'){
      shinydashboard::box(title = tagList("Color settings ",shiny::icon("palette")),
                          width = NULL, status = "warning", solidHeader = T,
          column(6, htmlOutput("color_picker")),
          column(6 , br(), actionButton("col_reset", "Default colours", icon = icon("undo")),
                 br(), br(), actionButton("save_color", "Save colors & apply to all", icon = icon("save"))))
    }
  })
  
  observeEvent(input$col_reset, {
    cols <- gg_fill_hue(length(levels_selected()))
    lapply(seq_along(levels_selected()), function(i) {
      do.call(what="updateColourInput", args=list(session=session, inputId=paste0("color_", levels_selected()[i]), value=cols[i]))
    })
  })
  
  output$color_picker <- renderUI({
    #Color picker
    colsModif <- SummarizedExperiment::colData(scExp())[,c(input$color_by,paste0(input$color_by,"_color"))] %>% unique()
    lapply(seq_along(levels_selected()), function(i) {
      colourpicker::colourInput(inputId=paste0("color_", levels_selected()[i]),
                                label=paste0("Choose colour for ", levels_selected()[i]),
                                value=colsModif[i,paste0(input$color_by,"_color")], returnName=TRUE) ## Add ", palette = "limited"" to get selectable color panel       
    })
  })
  
  observeEvent(input$save_color, {  
    req(scExp(), input$color_by)

    color_df = get_color_dataframe_from_input(input,levels_selected(),input$color_by)

    scExp. = colors_scExp(scExp(),annotCol = input$color_by,color_by = input$color_by, color_df = color_df)
    scExp(scExp.)
    
    save("scExp", file = file.path(init$data_folder, "ChromSCape_analyses",
                                   analysis_name(), "Filtering_Normalize_Reduce",
                                   paste0(input$selected_reduced_dataset,".RData")))
    rm(scExp.)
  })
  
  levels_selected <- reactive({
    req(scExp(),input$color_by)
    levels_selected = SummarizedExperiment::colData(scExp())[,input$color_by] %>% unique() %>% as.vector()
  })

  observeEvent(unlocked$list,able_disable_tab(c("selected_reduced_dataset","pca"),"cons_clustering")) # if conditions are met, unlock tab Correlation Clustering
  
  
  ###############################################################
  # 3. Consensus clustering on correlated cells
  ###############################################################

  corColors <- grDevices::colorRampPalette(c("royalblue","white","indianred1"))(256)
  
  observeEvent(priority = 10,input$tabs, once = T,
               if(input$tabs == "cons_clustering"){
                 file = file.path(init$data_folder, "ChromSCape_analyses",
                                  analysis_name(), "correlation_clustering",
                                  paste0(selected_filtered_dataset(),".RData"))
                 if(file.exists(file)){
                   myData = new.env()
                   load(file, envir = myData)
                   scExp_cf(myData$data$scExp_cf)
                 } else {
                   scExp_cf(scExp())
                 }
               }
  )
  
  consensus_ran = reactive({
    req(scExp_cf())
    "consclust" %in% names(scExp_cf()@metadata)
  })
  
  observeEvent(priority = 10, {selected_filtered_dataset()
               input$nclust
               input$tabs
               scExp_cf()
               },{
                 req(input$nclust, scExp_cf())
                 if(input$tabs == "cons_clustering"){
                   scExp_cf(choose_cluster_scExp(scExp_cf(), nclust = as.numeric(input$nclust),
                                                 consensus = consensus_ran()))
                 }
               })
  
  hc_pca_plot <- reactive({
    req(scExp_cf())
    unlocked$list$cor_clust_plot=TRUE;
    unlocked$list$affectation=TRUE;
    plot_heatmap_scExp(scExp_cf())
  })

  output$corr_clust_pca_plot <- renderPlot(hc_pca_plot())
  output$cons_corr_clust_pca_plot <- renderPlot(plot_heatmap_scExp(scExp_cf()))
  
  output$cons_corr_clust_pca_UI <- renderUI({
    if(consensus_ran()){
      plotOutput("cons_corr_clust_pca_plot", height = 500, width = 500) %>%
        shinyhelper::helper(type = 'markdown', icon ="info-circle",
                            content = "correlation_clustering")
    }
  })
  


  
  output$download_cor_clust_plot <- downloadHandler(
    filename=function(){ paste0("correlation_clustering_", input$selected_reduced_dataset, ".png")},
    content=function(file){
      grDevices::png(file, width=1200, height=1400, res=300,pointsize = 8)
      plot_heatmap_scExp(scExp_cf())
      grDevices::dev.off()

  })
  
  correlation_values <- reactiveValues(limitC=vector(length=500))
  corChIP <- reactive({ SingleCellExperiment::reducedDim(scExp(),"Cor") })
  z <- reactive({ matrix(sample(t(SingleCellExperiment::reducedDim(scExp(),"PCA"))), nrow=ncol(SingleCellExperiment::reducedDim(scExp(),"PCA"))) })
  thresh2 <- reactive({quantile(cor(z()), probs=seq(0,1,0.01))})
  limitC <- reactive({thresh2()[input$corr_threshold+1]})
  
  
  cell_cor_hist <- reactive({
    req(scExp())
    hist(corChIP(), prob=TRUE, col=alpha("steelblue", 0.8), breaks=50, ylim=c(0,4), main="Distribution of cell to cell correlation scores", xlab="Pearson Corr Scores")
    lines(density(corChIP()), col="blue", lwd=2)
    lines(density(cor(z())), col="black", lwd=2)
    abline(v=limitC(), lwd=2, col="red", lty=2)
    legend("topleft", legend=c("dataset", "randomized data", "correlation threshold"), col=c("blue", "black", "red"), lty=c(1, 1, 2), cex=0.8)
    
  })
  
  output$cell_cor_hist_plot <- renderPlot(cell_cor_hist())
  
  output$download_cor_clust_hist_plot <- downloadHandler(
    filename=function(){ paste0("correlation_distribution_", input$selected_reduced_dataset, ".png")},
    content=function(file){
      grDevices::png(file, width=2000, height=1400, res=300)
      hist(corChIP(), prob=TRUE, col=alpha("steelblue", 0.8), breaks=50, ylim=c(0,4),cex=0.4, main="Distribution of cell to cell correlation scores", xlab="Pearson Corr Scores")
      lines(density(corChIP()), col="blue", lwd=2)
      lines(density(cor(z())), col="black", lwd=2)
      abline(v=limitC(), lwd=2, col="red", lty=2)
      legend("topleft", legend=c("dataset", "randomized data", "correlation threshold"), col=c("blue", "black", "red"), lty=c(1, 1, 2), cex=0.8)
      grDevices::dev.off()
  })
  
  observeEvent(input$filter_corr_cells, {  # retreiveing cells with low correlation score
    withProgress(message='Filtering correlated cells...', value = 0, {
      incProgress(amount=0.6, detail=paste("Filtering"))
      scExp_cf(filter_correlated_cell_scExp(scExp_cf(), random_iter = 50, corr_threshold = input$corr_threshold,
                                            percent_correlation = input$percent_correlation))
      scExp_cf(choose_cluster_scExp(scExp_cf(), nclust = as.numeric(input$nclust), consensus = consensus_ran()))
      incProgress(amount=0.2, detail=paste("Saving"))
      data = list("scExp_cf" = scExp_cf())
      save(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                              paste0(input$selected_reduced_dataset,  ".RData")))
      incProgress(amount=0.2, detail=paste("Finished"))
      
      updateActionButton(session, "filter_corr_cells", label="Saved", icon = icon("check-circle"))
  
    })
  })
  
  observeEvent({input$selected_reduced_dataset  # reset label on actionButtion when new filtering should be filtered
    input$percent_correlation}, {
      updateActionButton(session, "filter_corr_cells", label="Filter & save", icon=character(0))
    })
  
  output$num_cell_before_cor_filt <- renderTable(
  {
    req(scExp())
    num_cell_before_cor_filt_scExp(scExp())
  })
  
  output$num_cell_after_cor_filt <- renderTable(
  {
    req(scExp(),scExp_cf())
    num_cell_after_cor_filt_scExp(scExp(),scExp_cf())
  })
  
  output$table_cor_filtered = renderUI({
    if(!is.null(scExp()@metadata$limitC)){
      return(tableOutput("num_cell_before_cor_filt"))
    } else{
      return(tableOutput("num_cell_after_cor_filt_scExp"))
    }
  })
  
  selected_filtered_dataset <- reactive({input$selected_reduced_dataset})
  
  output$filtered_data_selection_format <- renderText({"The name of the filtered dataset is composed of the following information: data set name, min percentage of reads per cell, 
    min percentage of cells to support a window, quantile of cell read counts to keep, correlation threshold, percent of cell correlation. To work on a different dataset or different preprocessing state, select it on the first page."})
  output$select_n_clust_hc = renderUI({
    selectInput("nclust", "Select number of clusters:", choices=c(2:10))
    })
  output$select_n_clust_chc = renderUI({
    selectInput("nclust", "Select number of clusters:", choices=c(2:10))
  })
  
  filtered_dataset <- observeEvent(selected_filtered_dataset(), {
    if(!is.null(selected_filtered_dataset()) ){
      file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering", paste0(selected_filtered_dataset(), ".RData"))
      if(file.exists(file)){
        myData = new.env()
        load(file, envir = myData)
        scExp_cf(myData$data$scExp_cf)
      }
    }
  })
  
  plotting_directory <- reactive({
    req( selected_filtered_dataset())
    if(!dir.exists(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots")))
      dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots"))
    file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots", selected_filtered_dataset())
  })
  
  clust <- reactiveValues(clust_pdf=NULL)
  
  observeEvent(selected_filtered_dataset(), priority = 10, {
 
  if(file.exists(file.path(init$data_folder, "ChromSCape_analyses",
                           analysis_name(), "correlation_clustering",
                           "Plots", selected_filtered_dataset(), "consensus.pdf"))){
    clust$clust_pdf <- file.path("Plots", selected_filtered_dataset(), "consensus.pdf")
  }
  })
  
  observeEvent(selected_filtered_dataset(), priority = 11, {
    filename <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                          paste0(input$selected_reduced_dataset, ".RData"))
      if(file.exists(filename)){
      myData = new.env()
      load(filename, envir = myData)
      scExp_cf(myData$data$scExp_cf)
      if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
        updateSelectInput(session, inputId = "nclust", label = "Select number of clusters:", choices=c(2:10),
                            selected = dplyr::n_distinct(SummarizedExperiment::colData(scExp_cf())$cell_cluster))
      }
    } else {
      NULL
    }
  })
  
  clusterAlg <- reactive({if(input$clusterAlg=="K-means") "kmdist"
    else if(input$clusterAlg == "Partitioning Medoids") "pam"
    else "hc"})
  
  observeEvent(input$do_cons_clust, {
    withProgress(message='Performing consensus clustering...', value = 0, {
      incProgress(amount=0.4, detail=paste("part one"))

      scExp_cf(consensus_clustering_scExp(scExp_cf(), reps = as.numeric(input$consclust_iter),
                                          maxK=as.numeric(input$maxK),
                                          clusterAlg = clusterAlg(),
                                          prefix = plotting_directory()))
      data = list("scExp_cf" = scExp_cf())
      save(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                  paste0(input$selected_reduced_dataset, ".RData")))

      clust$clust_pdf <- NULL  # needed in order to update the pdf output
      clust$clust_pdf <- file.path("Plots", selected_filtered_dataset(), "consensus.pdf")
  
      incProgress(amount=0.2, detail=paste("Finished"))
    })
  })
  
 output$cons_clust_pdf <- renderUI({
    req(clust$clust_pdf)
    if(!is.null(clust$clust_pdf)) tags$iframe(style = "height:500px; width:100%", src = clust$clust_pdf)
    })

 output$icl_plot <- renderPlot(plot_cluster_consensus_scExp(scExp_cf()))
  output$cluster_consensus_plot <- renderUI({
    if(!is.null(scExp_cf())){
      if("icl" %in% names(scExp_cf()@metadata)) {
         plotOutput("icl_plot", height = 350, width = 600)
      }
    }
  })
  
  output$cluster_consensus_info <- renderText({
    if(dir.exists(plotting_directory())){
      paste0("If the PDF of consensus clustering results doesn't display correctly,
             you can take a look at the PDF saved locally at :", 
             plotting_directory())
  } else ''
  })
  output$nclust_selection_info <- renderText({"After performing the clustering
    and checking the results for different numbers of clusters, select here the
    preferred number of clusters to make additional annotated plots."})

  observeEvent(input$choose_cluster, {
    if(is.null(scExp_cf())){
      print("scExp cf is null. Run consensus clustering first.")
    }
    else{
      if(! "consclust" %in% names(scExp_cf()@metadata)){
        showNotification("Could not find clustering results. Please make sure to perform the clustering before generating the final figures.", duration=7, closeButton=TRUE, type="warning")
      } else{
        withProgress(message='Preparing final figures...', value = 0, {
          incProgress(amount=0.2, detail=paste("Choosing cluster & recalculating tsne..."))
          scExp_cf(choose_cluster_scExp(scExp_cf(), as.integer(input$nclust)))
          incProgress(amount=0.4, detail = paste("Finishing..."))
          data = list("scExp_cf" = scExp_cf())
          save(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                      paste0(input$selected_reduced_dataset,".RData")))
          incProgress(amount=0.4, detail = paste("Saved"))
        })
      }
    }
  })
  
output$cons_clust_anno_plot <- renderPlot({
  if(! is.null(scExp_cf())){
    if("ConsensusAssociation" %in% names(scExp_cf()@metadata)){
      colors <- SummarizedExperiment::colData(scExp_cf())[scExp_cf()@metadata$hc_consensus_association$order,"cell_cluster_color"]
      heatmap(SingleCellExperiment::reducedDim(scExp_cf(),"ConsensusAssociation")[scExp_cf()@metadata$hc_consensus_association$order,],
              Colv = as.dendrogram(scExp_cf()@metadata$hc_consensus_association),
              Rowv = NA, symm = FALSE, scale="none", col = grDevices::colorRampPalette(c("white", "blue"))(100),
              na.rm = TRUE, labRow = F, labCol = F, mar = c(5, 5), main = paste("consensus matrix k=", input$nclust, sep=""),
              ColSideCol = colors)
    }
  }
    })
  
output$anno_cc_box <- renderUI({
  if(! is.null(scExp_cf())){
    if("ConsensusAssociation" %in% names(scExp_cf()@metadata)){
      shinydashboard::box(title="Annotated consensus clustering", width = NULL, status="success", solidHeader = T,
          column(12, align="left", plotOutput("cons_clust_anno_plot", height = 500, width = 500),
                 downloadButton("download_anno_cc_plot", "Download image")))
    }
  }
  })
  
  output$download_anno_cc_plot <- downloadHandler(
    filename = function(){ paste0("consensus_clustering_k", input$nclust, "_", selected_filtered_dataset(), ".png")},
    content = function(file){
      grDevices::png(file, width = 1200, height = 800, res = 300)
      colors <- SummarizedExperiment::colData(scExp_cf())[scExp_cf()@metadata$hc_consensus_association$order,"cell_cluster_color"]
      heatmap(SingleCellExperiment::reducedDim(scExp_cf(),"ConsensusAssociation")[scExp_cf()@metadata$hc_consensus_association$order,],
              Colv = as.dendrogram(scExp_cf()@metadata$hc_consensus_association),
              Rowv = NA, symm = FALSE, scale="none", col = colorRampPalette(c("white", "blue"))(100),
              na.rm = TRUE, labRow = F, labCol = F, mar = c(5, 5), main = paste("consensus matrix k=", input$nclust, sep=""),
              ColSideCol = colors)
      grDevices::dev.off()
  })
  
  output$contingency_table_cluster <- renderUI({
    if(! is.null(scExp_cf())){
      if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
        shinydashboard::box(title="Samples & Cluster table", width = NULL, status="success", solidHeader = T,
            column(12, align="left", tableOutput("hcor_kable")),
            column(5,offset = 2, align="left", htmlOutput("chi_info"),br())
        )
      }
    }
  })

    
  output$hcor_kable <- function(){
    req(scExp_cf())
    if(! is.null(scExp_cf())){
      if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
        num_cell_in_cluster_scExp(scExp_cf())
      } else{
        NULL
      }
    } else NULL
  }
 
  output$tsne_box_cf <- renderUI({
    req(scExp_cf(),input$color_by)
    if(!is.null(scExp_cf())){
      if("TSNE" %in% SingleCellExperiment::reducedDimNames(scExp_cf())){
        req(scExp_cf(), input$color_by_cf)
        p = plot_reduced_dim_scExp(scExp_cf(),input$color_by_cf, "TSNE",
                                   select_x = "Component_1",
                                   select_y = "Component_2") +
          ggtitle("t-SNE")
        output$tsne_plot_cf <- plotly::renderPlotly(
          plotly::ggplotly(p, tooltip="Sample", dynamicTicks = T))
        shinydashboard::box(title="tSNE visualization", width = NULL, status="success", solidHeader=T,
                            column(12, align="left", plotly::plotlyOutput("tsne_plot_cf")))
      }
    }
  })
  
  umap_p_cf <- reactive({
    req(scExp_cf(), input$color_by_cf)
    p = plot_reduced_dim_scExp(scExp_cf(),input$color_by_cf, "UMAP",
                               select_x = "Component_1",
                               select_y = "Component_2") +
      ggtitle("UMAP")
    p
  })
  output$umap_plot_cf <- plotly::renderPlotly(
    plotly::ggplotly(umap_p_cf(), tooltip="Sample", dynamicTicks = T) )
  
  levels_selected_cf <- reactive({
    req(scExp_cf(),input$color_by_cf)
    levels_selected_cf =
      SummarizedExperiment::colData(scExp_cf())[,input$color_by_cf] %>%
      unique() %>% as.vector()
  })
  
  output$umap_box_cf <- renderUI({
    if(! is.null(scExp_cf())){
      if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf())) ){
        shinydashboard::box(title="Vizualisation in reduced dimensions", width = NULL, status="success", solidHeader = T,
                            column(6, align="left",
                                   selectInput("color_by_cf", "Color by",
                                               selected = "cell_cluster",
                                               choices = c('sample_id', 'total_counts',
                                                           'cell_cluster','batch_id'))),
                            column(12, align="left", plotly::plotlyOutput("umap_plot_cf")))
      }
    }
  })
  
  output$color_box_cf <- renderUI({
    req(input$color_by_cf)
    if(! is.null(scExp_cf())){
      if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
        if(input$color_by_cf != 'total_counts'){
          shinydashboard::box(title=tagList("Color settings ",shiny::icon("palette")), width = NULL, status = "warning", solidHeader = T,
                              column(6, htmlOutput("color_picker_cf")),
                              column(4 , br(), actionButton("col_reset_cf", "Default colours", icon = icon("undo")),
                                     br(), br(), actionButton("save_color_cf",
                                                              "Save colors & apply to all", icon = icon("save"))))
        }
      }
    }
  })
  
  observeEvent(input$save_color_cf, {  
    req(scExp_cf(), input$color_by_cf)
    color_df = get_color_dataframe_from_input(input,levels_selected_cf(), input$color_by_cf, "color_cf_")
    scExp_cf. = colors_scExp(scExp_cf(), annotCol = input$color_by_cf,
                             color_by = input$color_by_cf, color_df = color_df)
    scExp_cf(scExp_cf.)
    rm(scExp_cf.)
  })
  
  observeEvent(input$col_reset_cf, {
    cols <- gg_fill_hue(length(levels_selected_cf()))
    lapply(seq_along(levels_selected_cf()), function(i) {
      do.call(what="updateColourInput",
              args=list(session=session,
                        inputId=paste0("color_cf", levels_selected_cf()[i]), value=cols[i]))
    })
  })

  output$color_picker_cf <- renderUI({
    #Color picker
    colsModif <- as.data.frame(SummarizedExperiment::colData(scExp_cf()))[,c(input$color_by_cf,paste0(input$color_by_cf,"_color"))] %>% unique()
    lapply(seq_along(levels_selected_cf()), function(i) {
      colourpicker::colourInput(inputId=paste0("color_cf_", levels_selected_cf()[i]),
                                label=paste0("Choose colour for ", levels_selected_cf()[i]),
                                value=colsModif[i,paste0(input$color_by_cf,"_color")],
                                returnName = TRUE) ## Add ", palette = "limited"" to get selectable color panel       
    })
  })
  
  # observeEvent({
  #   selected_filtered_dataset()
  #   scExp_cf()
  # },
  # {
  #   if(length(selected_filtered_dataset())>0 & !is.null(scExp_cf())){
  #     if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
  #       unlocked$list$affectation = T
  #     } else {
  #       unlocked$list$affectation = F
  #     }
  #   } else unlocked$list$affectation = F
  # }
  # )
  
  observeEvent(unlocked$list, {able_disable_tab(c("selected_analysis","selected_reduced_dataset","affectation"),c("peak_calling","diff_analysis"))}) 
  
  ###############################################################
  # 5. Peak calling [optional]
  ###############################################################
  
  output$peak_calling_info <- renderText({"This module is optional, but recommended in order to obtain more precise results for enrichment analysis. 
    Based on the BAM files for the samples in your project, peaks will be called using MACS2 [only work on unix systems] so that the counts can be mapped to the gene TSS more specifically.
    "})
  
  can_run = reactiveVal({FALSE})
  
  output$peak_calling_system <- renderText({
    platform = as.character(.Platform[1])
    if(length(grep("nix",platform,ignore.case = T)) ){
      macs2=""
      samtools=""
      try({
        macs2 = system2("which",args = "macs2",stdout = TRUE)
      })
      try({
        samtools = system2("which",args = "samtools",stdout = TRUE)
      })
      if(length(macs2)>0 & length(samtools)>0){
        can_run(TRUE)
        return(paste0("<b>You are running on an ", platform, " OS.<br>Found MACS2 at ", macs2, ".<br>Found samtools at ", samtools,"</b>"))
      }
      if(length(macs2)>0 & length(samtools)==0)
        return(paste0("<b>You are running on an ", platform, " OS.<br>Found MACS2 at ", macs2,".<br>Didn't find samtools ! Please install samtools or skip this step.</b>"))
      if(length(macs2)==0 & length(samtools)>0)
        return(paste0("<b>You are running on an ", platform, " OS.<br>Didn't find MACS2, please install MACS2 or skip this step.<br>Found samtools at ", samtools,"</b>"))
      if(length(macs2)==0 & length(samtools)==0)
        return(paste0("<b>You are running on an ", platform, " OS.<br>Didn't find MACS2, please install MACS2 or skip this step.<br>Didn't find samtools ! Please install samtools or skip this step.</b>"))
      
      } else {
      return(paste0("<b>You are running on a non unix system, peak calling is not available, you can move directly to differential analysis.</b> "))
    }
    })
  
  output$peak_calling_icon = renderText({
    if(can_run()) {
      return( as.character(icon("check-circle", class = "large_icon")))}
    else{
      return( as.character(icon("times-circle", class = "large_icon")))
    }
  })
      
  shinyFiles::shinyDirChoose(input, "bam_folder", roots = volumes, session = 
                   session)
  list_bams = reactive({
    if(!is.null(bam_folder())){
      list.files(bam_folder(), full.names = T, pattern = "*.bam$")
    }
  })
  
  bam_folder = reactiveVal(NULL)
  
  observeEvent(input$bam_folder,
               {
               if(!is.null(input$bam_folder)){
                 # browser()
                 bam_folder(parseDirPath(volumes, input$bam_folder))
                 output$bam_dir <- renderText(bam_folder())
               }
})

  output$bam_upload <- renderUI({
    req(scExp_cf(), list_bams())
    if(!is.null(list_bams()))
      selectInput("bam_selection", label = "Selected BAM files",
                  choices = basename(list_bams()), multiple = T,
                  selected = basename(list_bams()) )
  })
  
  observeEvent(input$do_pc, {
    nclust = length(unique(scExp_cf()$cell_cluster))  
    withProgress(message='Performing enrichment analysis...', value = 0, {
        
        incProgress(amount = 0.1, detail = paste("Starting Peak Calling..."))
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks"), showWarnings = FALSE)
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", nclust)), showWarnings = FALSE)
        odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", nclust))
        sample_ids <- unique(SummarizedExperiment::colData(scExp_cf())$sample_id)
        # inputBams <- as.character(unlist(sapply(sample_ids, function(x){ input[[paste0('bam_', x)]] })))
        inputBams <- as.character(list_bams())

        checkBams <- sapply(inputBams, function(x){ if(file.exists(x)){ 0 } else { 
          showNotification(paste0("Could not find file ", x, ". Please make sure to give a full path including the file name."),
                           duration = 7, closeButton = TRUE, type="warning"); 1}  
        })
        incProgress(amount = 0.3, detail = paste("Running Peak Calling..."))
        if(sum(checkBams)==0){
          scExp_cf(subset_bam_call_peaks(scExp_cf(), odir, inputBams, as.numeric(input$pc_stat_value), annotation_id(), input$peak_distance_to_merge))
          data = list("scExp_cf" = scExp_cf())
          save(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                      paste0(input$selected_reduced_dataset, ".RData")))
          pc$new <- Sys.time()
          updateActionButton(session, "do_pc", label="Finished successfully", icon = icon("check-circle"))
        }
        incProgress(amount = 0.3, detail = paste("Finished Peak Calling..."))
      })
    
  })
  
  pc <- reactiveValues(new="")
  
  observeEvent({selected_filtered_dataset()  # reset label on actionButtion when new peak calling should be performed
    input$pc_stat
    input$pc_stat_value}, {
      pc$available_pc <- has_available_pc()
      updateActionButton(session, "do_pc", label="Start", icon = character(0))
  })
  
  has_available_pc <- reactive({
    if(!is.null(scExp_cf())){
      if("refined_annotation" %in% names(scExp_cf()@metadata) ){
        return(TRUE)
      } else return(FALSE)
    } else{
      return(FALSE)
    }
  })
  # available_pc_plots <- reactive({
  #   print(paste0("new peaks done: ", pc$new))
  #   fe <- sapply(c(1:input$pc_k_selection), function(i){file.exists(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", input$pc_k_selection), paste0("C", i, "_model.r")))})
  #   which(fe==TRUE)
  # })

  # output$pc_plot_box <- renderUI({
  #   if(has_available_pc()){
  #     shinydashboard::box(title="Peak calling visualization", width = NULL, status="success", solidHeader = T,
  #         column(8, align="left", selectInput("pc_cluster","Select cluster (only those shown for which plots are available):", choices = paste0("C", available_pc_plots()))),
  #         column(12, align="left",
  #                plotOutput("peak_model_plot", height = 500, width = 500),
  #                plotOutput("cross_corr_plot", height = 500, width = 500)))
  #   }
  # })
  # 
  # pc_data_p <- reactive({ req(has_available_pc()); as.numeric(unlist(read.table(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", input$pc_k_selection), paste0(input$pc_cluster, "_data_p.txt"))))) })
  # pc_data_m <- reactive({ req(has_available_pc()); as.numeric(unlist(read.table(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", input$pc_k_selection), paste0(input$pc_cluster, "_data_m.txt"))))) })
  # pc_data_xcorr <- reactive({ req(has_available_pc()); as.numeric(unlist(read.table(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", input$pc_k_selection), paste0(input$pc_cluster, "_data_xcorr.txt"))))) })
  # pc_data_ycorr <- reactive({ req(has_available_pc()); as.numeric(unlist(read.table(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", input$pc_k_selection), paste0(input$pc_cluster, "_data_ycorr.txt"))))) })
  # peak_model_p <- reactive({
  #   req(has_available_pc())
  #   x <- seq.int((length(pc_data_p())-1)/2*-1,(length(pc_data_p())-1)/2)
  #   plot(x,pc_data_p(),type='l',col = c('red'),main='Peak Model',xlab='Distance to the middle',ylab='Percentage')
  #   lines(x,pc_data_m(),col = c('blue'))
  #   legend('topleft',c('forward tags','reverse tags'),lty = c(1,1,1),col = c('red','blue'))
  # })
  # output$peak_model_plot <- renderPlot(peak_model_p())
  # cross_corr_p <- reactive({
  #   req(has_available_pc())
  #   altd  <- c(297)
  #   plot(pc_data_xcorr(),pc_data_ycorr(),type='l',col = c('black'),main='Cross-Correlation',xlab='Lag between + and - tags',ylab='Correlation')
  #   abline(v = altd,lty = 2,col = c('red'))
  #   legend('topleft','alternative lag(s)',lty = 2,col='red')
  #   legend('topright','alt lag(s) : 297',bty='n')
  # })
  # output$cross_corr_plot <- renderPlot(cross_corr_p())
  
  
  ###############################################################
  # 6. Differential analysis
  ###############################################################
  
  output$diff_analysis_info <- renderText({"Differential analysis will be performed using the cluster assignment obtained on the Consensus clustering page. To use a different number of clusters, 
    go to this page and first perform the clustering, then select the preferred number of clusters in the box on the right in order to display and save the data. 
    It will then appear here for selection."})
  output$selected_k <- renderText({
    paste0( "Number of clusters selected  = ", dplyr:::n_distinct(SummarizedExperiment::colData(scExp_cf())$cell_cluster))
  })
  output$contrib_hist <- renderUI({ if(input$only_contrib_cells){ plotOutput("contrib_hist_p", height = 250, width = 500) }})
  output$contrib_hist_p <- renderPlot(contrib_hist_plot())
  contrib_hist_plot <- reactive({
    maxCons <- tapply(scExp_cf()@metadata$icl$itemConsensus$itemConsensus[scExp_cf()@metadata$icl$itemConsensus$k==length(unique(scExp_cf()$cell_cluster))],
                      scExp_cf()@metadata$icl$itemConsensus$item[scExp_cf()@metadata$icl$itemConsensus$k==length(unique(scExp_cf()$cell_cluster))], max)
    hist(maxCons, col="steelblue", breaks = 80, main="Max cluster contribution per cell", xlab="", ylab="number of cells")
    abline(v = input$contrib_thresh, lwd = 2, col="red", lty = 2)
    legend("topleft", legend = c("cluster contribution threshold"), col = c("red"), lty = c(2), cex = 0.8)
  })
  output$contrib_thresh <- renderUI({ if(input$only_contrib_cells){ sliderInput("contrib_thresh", "Select minimum cluster contribution for cells:", min = 0.6, max = 1, value = 0.9, step = 0.01) }})
  output$contrib_info <- renderUI({ if(input$only_contrib_cells){ textOutput("contrib_info_text") }})
  output$contrib_info_text <- renderText({
    total_cells <- length(unique(scExp_cf()@metadata$icl$itemConsensus$item[scExp_cf()@metadata$icl$itemConsensus$k==length(unique(scExp_cf()$cell_cluster))]))
    sel_cells <- length(unique(scExp_cf()@metadata$icl$itemConsensus$item[scExp_cf()@metadata$icl$itemConsensus$k==length(unique(scExp_cf()$cell_cluster)) & scExp_cf()@metadata$icl$itemConsensus$itemConsensus >= input$contrib_thresh]))
    paste("Selected top", sel_cells, "cells out of", total_cells)
  })
  
  observeEvent(c(input$qval.th, input$tabs, input$cdiff.th, input$de_type, selected_filtered_dataset()), priority = 10,{
    if(input$tabs == "diff_analysis"){
      if(!is.null(selected_filtered_dataset()) && !is.null(input$qval.th) && !is.null(input$cdiff.th)){
        filename <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "diff_analysis_GSEA",
                              paste0(selected_filtered_dataset(), "_", input$nclust,
                                     "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".RData"))
        
        if(file.exists(filename)){
          myData = new.env()
          load(filename, envir = myData)
          scExp_cf(myData$data$scExp_cf)
        } else {
          NULL
        }
      }
    }
    })
  
  observeEvent(input$do_wilcox, {  # perform differential analysis based on wilcoxon test
    withProgress(message='Performing differential analysis...', value = 0, {
      incProgress(amount = 0.2, detail = paste("Initializing DA"))
      if(batchUsed()) block = T else block = F
      scExp_cf(differential_analysis_scExp(scExp_cf(), de_type = input$de_type,
                                           cdiff.th = input$cdiff.th, qval.th = input$qval.th, block)) 
      incProgress(amount = 0.6, detail = paste("Finishing DA..."))
      data = list("scExp_cf" = scExp_cf())
      save(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "diff_analysis_GSEA",
                                  paste0(selected_filtered_dataset(), "_", length(unique(scExp_cf()$cell_cluster)),
                                         "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".RData")))
  
      incProgress(amount = 0.2, detail = paste("Saving DA"))
    })
  })
  
  output$da_summary_box <- renderUI({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        shinydashboard::box(title="Number of differentially bound regions", width = NULL, status="success", solidHeader = T,
            column(5, align="left", br(), tableOutput("da_summary_kable")),
            column(7, align="left", plotOutput("da_barplot", height = 270, width = 250)),
            column(12, align="left", div(style = 'overflow-x: scroll', DT::dataTableOutput('da_table')), br()),
            column(4, align="left", downloadButton("download_da_barplot", "Download barplot")),
            column(4, align="left", downloadButton("download_da_table", "Download table")))
      }
    }
  })
  
  output$da_summary_kable <- function(){
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        scExp_cf()@metadata$diff$summary %>%
          kableExtra::kable(escape = F, align="c") %>%
          kableExtra::kable_styling(c("striped", "condensed"), full_width = F)
      }
    }
  }
  output$da_barplot <- renderPlot({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        plot_differential_summary_scExp(scExp_cf())
      }
    }
  })
  
  output$download_da_barplot <- downloadHandler(
    filename = function(){ paste0("diffAnalysis_numRegions_barplot_", selected_filtered_dataset(), "_", length(unique(scExp_cf()$cell_cluster)), "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".png")},
    content = function(file){
      grDevices::png(file, width = 800, height = 600, res = 150)
      plot_differential_summary_scExp(scExp_cf())
      grDevices::dev.off()
    })
  
  output$da_table <- DT::renderDataTable({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        table <- scExp_cf()@metadata$diff$res[, -c(1)]
        rownames(table) <- NULL
        for(i in seq(from = 0, to=(dim(table)[2]-8)/5, by = 1)){
          table[, (5*i+5)] <- round(table[, (5*i+5)], 3) #counts
          table[, (5*i+6)] <- round(table[, (5*i+6)], 3) #cdiff
        }
        table <- table[order(table$Rank.C1),]
        DT::datatable(table, options = list(dom='tpi'))
      }
    }
  })
  
  output$download_da_table <- downloadHandler(
    filename = function(){ paste0("diffAnalysis_data_", selected_filtered_dataset(),
                                  "_", length(unique(scExp_cf()$cell_cluster)), "_",
                                  input$qval.th, "_", input$cdiff.th, "_",
                                  input$de_type, ".csv")},
    content = function(file){
      write.table(scExp_cf()@metadata$diff$res, file, row.names = F, quote = F, sep=",")
    })
  
  output$da_visu_box <- renderUI({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        shinydashboard::box(title="Visualization", width = NULL, status="success", solidHeader = T,
            column(4, align="left", selectInput("gpsamp", "Select cluster:", choices = scExp_cf()@metadata$diff$groups)),
            column(12, align="left", plotOutput("h1_prop", height = 300, width = 500),
                   plotOutput("da_volcano", height = 500, width = 500)),
            column(4, align="left", downloadButton("download_h1_plot", "Download histogram")),
            column(4, align="left", downloadButton("download_da_volcano", "Download volcano plot")))
      }
    }
  })
  
  output$h1_prop <- renderPlot({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        plot_differential_H1_scExp(scExp_cf(), input$gpsamp)
      }  
    }
  })
  
  output$download_h1_plot <- downloadHandler(
    filename = function(){ paste0("diffAnalysis_numRegions_barplot_", selected_filtered_dataset(),
                                  "_", length(unique(scExp_cf()$cell_cluster)), "_", input$qval.th,
                                  "_", input$cdiff.th, "_", input$de_type, "_", input$gpsamp, ".png")},
    content = function(file){
      grDevices::png(file, width = 1000, height = 600, res = 300)
      plot_differential_H1_scExp(scExp_cf(), input$gpsamp)
      grDevices::dev.off()
  })
  
  output$da_volcano <- renderPlot({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        plot_differential_volcano_scExp(scExp_cf(),cell_cluster = input$gpsamp,
                                        cdiff.th = input$cdiff.th, qval.th = input$qval.th)
      }
    }
  })
  
  output$download_da_volcano <- downloadHandler(
    filename = function(){ paste0("diffAnalysis_numRegions_barplot_", selected_filtered_dataset(),
                                  "_", length(unique(scExp_cf()$cell_cluster)), "_",
                                  input$qval.th, "_", input$cdiff.th, "_", 
                                  input$de_type, "_", input$gpsamp, ".png")},
        content = function(file){
        grDevices::png(file, width = 900, height = 900, res = 300)
          plot_differential_volcano_scExp(scExp_cf(),
                                          cell_cluster = input$gpsamp,
                                          cdiff.th = input$cdiff.th,
                                          qval.th = input$qval.th)
        grDevices::dev.off()
    })
  
    observeEvent(unlocked$list, {
      able_disable_tab(c("selected_analysis","selected_reduced_dataset",
                         "diff_my_res"),"enrich_analysis")
    }) 
    
    observeEvent(scExp_cf(), {
      if(!is.null(scExp_cf()@metadata$diff)){
        unlocked$list$diff_my_res = TRUE 
        } else{
          unlocked$list$diff_my_res = FALSE
        }
    }) 
    
  
  ###############################################################
  # 7. Enrichment analysis
  ###############################################################
  
  MSIG.classes <- reactive({
    myData = new.env()
    eval(parse(text = paste0("data(", annotation_id(), ".MSIG.gs)")))
    eval(parse(text = paste0("classes = ",annotation_id(), ".MSIG.gs$Class")))
    unique(classes)
  })

  annotFeat_long <- reactive({
    af = as.data.frame(SummarizedExperiment::rowData(scExp_cf()))
    af = tidyr::separate_rows(af, Gene,sep = ", ")
    af
  })
  
  output$enr_info <- renderText({"Enrichment will be performed based on the significant genes per cluster that were computed on the previous page. 
    Please make sure that you have run the differential analysis on the clustering that you prefer before running the enrichment analysis."})
  
  canUsePeaks <- reactive({
    print("Can usepeak calling exist ?")
    if(!is.null(scExp_cf())){
      if("refined_annotation" %in% names(scExp_cf()@metadata) ){
        return(TRUE)
      } else return(FALSE)
    } else{
      return(FALSE)
    }
    })
  
  output$use_peaks <- renderUI({
    if(canUsePeaks()){ checkboxInput("use_peaks", "use peak calling results to refine analysis", value = TRUE) }
  })
  
  enr <- reactiveValues(Both = NULL, Overexpressed = NULL, Underexpressed = NULL)
  
  GencodeGenes <- reactive({
    myData = new.env()
    eval(parse(text = paste0("data(",annotation_id(),".GeneTSS, envir = myData)")))
    as.character(unique(
      eval(parse(text = paste0("myData$",annotation_id(),".GeneTSS$gene")))
    ))
  })
  
  observeEvent(input$do_enrich, {
    withProgress(message='Performing enrichment analysis...', value = 0, {

      incProgress(amount = 0.3, detail = paste("Running GSEA..."))
      scExp_cf(gene_set_enrichment_analysis_scExp(scExp_cf(), enrichment_qval = 0.01, qval.th = input$qval.th,
                                                  ref = annotation_id(), cdiff.th = input$cdiff.th,
                                                  peak_distance = 1000, use_peaks = input$use_peaks))
      incProgress(amount = 0.6, detail = paste("Finishing GSEA..."))
      data = list("scExp_cf" = scExp_cf() )
      save(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "diff_analysis_GSEA",
                                  paste0(selected_filtered_dataset(), "_", length(unique(scExp_cf()$cell_cluster)),
                                         "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".RData")))
      incProgress(amount = 0.6, detail = paste("Saving GSEA"))
      
    })
  })
  
  output$enr_clust_sel <- renderUI({ 
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
      selectInput("enr_clust_sel", "Select cluster:", choices = scExp_cf()@metadata$diff$groups) 
      }
    }
    })
  
  output$enr_class_sel <- renderUI({shiny::checkboxGroupInput(
    inputId = "enr_class_sel", inline = T,
    label =  "Select classes to display:",
    selected = MSIG.classes(), choiceNames = MSIG.classes(),
    choiceValues = MSIG.classes())
  })
  
  output$all_enrich_table <- DT::renderDataTable({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$enr) && !is.null(input$enr_clust_sel) && !is.null(input$enr_class_sel)){
        table_enriched_genes_scExp(scExp_cf(),set = "Both", cell_cluster = input$enr_clust_sel, input$enr_class_sel)
      }
    }
  })
  
  output$over_enrich_table <- DT::renderDataTable({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$enr)){
        table_enriched_genes_scExp(scExp_cf(), set = "Overexpressed", cell_cluster = input$enr_clust_sel, input$enr_class_sel)
      }
    }
  })
  
  output$under_enrich_table <- DT::renderDataTable({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$enr)){
        table_enriched_genes_scExp(scExp_cf(), set = "Underexpressed", cell_cluster = input$enr_clust_sel, input$enr_class_sel)
      }
    }
  })
  
  output$download_enr_data <- downloadHandler(
    filename = paste(selected_filtered_dataset(), length(unique(scExp_cf()$cell_cluster)), input$qval.th, input$cdiff.th, "enrichment_tables.zip", sep="_"),
    content = function(fname){
      fis <- c()
      for(i in 1:length(scExp_cf()@metadata$diff$groups)){
        if(!is.null(scExp_cf()@metadata$enr$Both[[i]])){
          filename <- paste0(scExp_cf()@metadata$diff$groups[i], "_significant_gene_sets.csv")
          fis <- c(fis, filename)
          write.table(scExp_cf()@metadata$enr$Both[[i]], file = filename, quote = FALSE, row.names = FALSE, sep=",")
        }
        if(!is.null(scExp_cf()@metadata$enr$Overexpressed[[i]])){
          filename <- paste0(scExp_cf()@metadata$diff$groups[i], "_enriched_gene_sets.csv")
          fis <- c(fis, filename)
          write.table(scExp_cf()@metadata$enr$Overexpressed[[i]], file = filename, quote = FALSE, row.names = FALSE, sep=",")
        }
        if(!is.null(scExp_cf()@metadata$enr$Underexpressed[[i]])){
          filename <- paste0(scExp_cf()@metadata$diff$groups[i], "_depleted_gene_sets.csv")
          fis <- c(fis, filename)
          write.table(scExp_cf()@metadata$enr$Underexpressed[[i]], file = filename, quote = FALSE, row.names = FALSE, sep=",")
        }
      }
      zip(zipfile = fname, files = fis)},
    contentType = "application/zip"
  )
  
  output$gene_sel <- renderUI({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$enr)){
        most_diff = scExp_cf()@metadata$diff$res %>% dplyr::select(ID,starts_with("qval."))
        most_diff[,"qval"] = Matrix::rowMeans(as.matrix(most_diff[,-1]))
        most_diff = dplyr::left_join(most_diff[order(most_diff$qval),], annotFeat_long(),by = c("ID"))
        most_diff = most_diff %>% dplyr::filter(!is.na(Gene)) 
        genes = base::intersect(most_diff$Gene,unique(GencodeGenes()))

        selectizeInput(inputId = "gene_sel", "Select gene:",options= list(maxOptions = 250),genes)
      }
    }
  })
  
  output$region_sel <- renderUI({
    req(input$gene_sel)
    subset <- annotFeat_long()[which(annotFeat_long()$Gene==input$gene_sel), ]
    subset <- subset[order(subset$distance),]
    regions <- paste0(subset$ID, " (distance to gene TSS: ", subset$distance, ")")
    selectInput("region_sel", "Select associated genomic region:", choices = regions)
  })
  

  gene_tsne_p <- reactive({
    req(input$gene_sel, input$region_sel)
    region <- strsplit(input$region_sel, " ")[[1]][1]
    if(region %in% rownames(scExp_cf())){
      p <- ggplot(as.data.frame(SingleCellExperiment::reducedDim(scExp_cf(), "TSNE")),
                  aes(x = Component_1, y = Component_2)) +
        geom_point(alpha = 0.5, aes(color = SingleCellExperiment::normcounts(scExp_cf())[region, ],
                                    shape = SummarizedExperiment::colData(scExp_cf())$cell_cluster)) +
        labs(color="norm. count for region", shape="Cluster", x="t-SNE 1", y="t-SNE 2") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour="black"),
              panel.border = element_rect(colour="black", fill = NA)) +
        viridis::scale_color_viridis(direction=-1)
    }
  })
  
  output$gene_tsne_plot <- plotly::renderPlotly({
      req(input$gene_sel, input$region_sel)
      region <- strsplit(input$region_sel, " ")[[1]][1]
      if(region %in% rownames(scExp_cf())){
        plotly::ggplotly(gene_tsne_p(), tooltip="Sample", dynamicTicks = T)
      }
  })
  
  
  ###############################################################
  # 8. Close app
  ###############################################################
  
  output$analysis_deletion_info <- renderText({"The selected data set will be fully deleted from the computer, including all reduced data versions that have been produced so far for this set."})
  output$selected_delete_analysis <- renderUI({ selectInput("selected_delete_analysis", "Select data set :", choices = init$available_raw_datasets) })
  
  observeEvent(input$delete_analysis, {  # delete selected dataset
    withProgress(message='Deleting data set', value = 0, {
      incProgress(amount=0.5, detail=paste("..."))
      unlink(file.path(init$data_folder, "ChromSCape_analyses", input$selected_delete_analysis), recursive=TRUE)
      init$available_raw_datasets <- list.dirs(path=file.path(init$data_folder, "ChromSCape_analyses"), full.names=FALSE, recursive=FALSE)
      init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())
      incProgress(amount=0.5, detail=paste("... finished"))
    })
    showNotification("Data set successfully deleted.", duration=5, closeButton=TRUE, type="warning")
  })
  
  observeEvent(input$close_and_save, {
    unlink(file.path("www", "images", "*"))
    unlink(file.path(".", "*.csv"))
    
    if(!is.null(scExp()) & !is.null(input$selected_reduced_dataset)){
      scExp = isolate(scExp())
      save(scExp, file = file.path(init$data_folder, "ChromSCape_analyses",
                                   analysis_name(), "Filtering_Normalize_Reduce",
                                   paste0(input$selected_reduced_dataset,".RData")))
      
    }
    if(!is.null(scExp_cf()) & !is.null(selected_filtered_dataset())){
      data = list()
      data$scExp_cf = isolate(scExp_cf())
      
      if("consclust" %in% names(scExp_cf()@metadata)){
        if("diff" %in% names(scExp_cf()@metadata)){
          dir =file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), 
                         "diff_analysis_GSEA", paste0(
                           selected_filtered_dataset(), "_",
                           length(unique(scExp_cf()$cell_cluster)),
                           "_", input$qval.th, "_", input$cdiff.th, "_",
                           input$de_type, ".RData"))
        } else{
          dir = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),
                          "correlation_clustering",
                         paste0(selected_filtered_dataset(), 
                                 ".RData"))
        }
        save(scExp, file = dir)
      } 
      
    }
    
    lapply(names(resourcePaths()), removeResourcePath)
    
    print("Thank you & see you next time !")
    stopApp()
  })
  
})
