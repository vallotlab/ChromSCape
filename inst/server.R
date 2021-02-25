
shinyServer(function(input, output, session) {
  
  library(ggplot2)
  library(dplyr)
  
  options(shiny.maxRequestSize = 50000*1024^2) # allow upload of files with max 15GB
  
  ###############################################################
  # 0. Global variables and functions
  ###############################################################
  
  #Initializating user experience functions
  js$init_directory() #Getting cookie for the directory
  volumes = c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
shinyhelper::observe_helpers(help_dir = "www/helpfiles",withMathJax = TRUE)
  
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
  
  analysis_name <- reactive({ input$selected_analysis })
  annotation_id_norm <- reactive({ read.table(file.path(init$data_folder, 'ChromSCape_analyses', input$selected_analysis, 'annotation.txt'), header = FALSE, stringsAsFactors = FALSE)[[1]] })
  annotation_id <- reactive({ read.table(file.path(init$data_folder, 'ChromSCape_analyses', analysis_name(), 'annotation.txt'), header = FALSE, stringsAsFactors = FALSE)[[1]] })
  
  #Global Functions
  init <- reactiveValues(data_folder =  getwd(), datamatrix = data.frame(), annot_raw = data.frame(),
                         available_analyses = list.dirs(path = file.path(getwd(), "ChromSCape_analyses"), full.names = FALSE, recursive = FALSE),
                         available_reduced_datasets = NULL, available_DA_GSA_datasets = NULL)
  reduced_datasets <- reactive({
    if (is.null(init$available_reduced_datasets)) c() else gsub('.{3}$', '', basename(init$available_reduced_datasets)) })
  
  observeEvent({analysis_name()},{
    init$available_reduced_datasets = get.available.reduced.datasets(analysis_name())
    print("observeEvent({analysis_name()}")
    print(init$available_reduced_datasets)
  })
  annotCol <- reactive({
    if("batch_name" %in% colnames(SummarizedExperiment::colData(scExp()))){
      c("sample_id","total_counts","batch_name")} else{
        c("sample_id","total_counts")
    }
    })
  output$feature_color <- renderUI({selectInput("color_by", "Color by", choices=annotCol())})
  
  observeEvent({
    analysis_name()
    input$selected_DA_GSA_dataset
    input$selected_reduced_dataset
    }, { # application header (tells you which data set is selected)
    req(analysis_name())
      header <- paste0('<b>Analysis : ', analysis_name(), ' </b>')
      if(!is.null(input$selected_DA_GSA_dataset)){
        if(input$tabs %in% c("diff_analysis","enrich_analysis")){
          pattern = paste0(input$selected_reduced_dataset, "_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+_")
          
          suffix = gsub(pattern,"", input$selected_DA_GSA_dataset)
          header <- paste0('<b>Analysis : ', analysis_name(), ' - ', suffix, ' </b>')
        }
      } 
      shinyjs::html("pageHeader", header) 
  })
  
  
  get.available.reduced.datasets <- function(selected_analysis){
    print("get.available.reduced.datasets")
    print(list.files(path = file.path(init$data_folder, "ChromSCape_analyses", selected_analysis,"Filtering_Normalize_Reduce"), full.names = FALSE, recursive = TRUE,
                     pattern="[[:print:]]+_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+_(uncorrected|batchCorrected).qs"))
    list.files(path = file.path(init$data_folder, "ChromSCape_analyses", selected_analysis,"Filtering_Normalize_Reduce"), full.names = FALSE, recursive = TRUE,
               pattern="[[:print:]]+_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+_(uncorrected|batchCorrected).qs")
    
  }
  
  able_disable_tab <- function(variables_to_check, tab_id) {

    able_or_disable = c()
    for(var in variables_to_check){
      if (unlocked$list[[var]]==TRUE) {
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
    selectInput("selected_analysis", "Choose analysis:",
                choices = init$available_analyses, multiple = FALSE) 
    })
  
  output$data_folder_info <- renderText({
    "All your analyses will be saved in this folder."
    })

  output$data_matrices_info <- renderText({"The file(s) name(s) for each matrix must be the sample name and must contain only alpha-numeric character and underscores."})
  
  
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
           
           init$available_analyses <- list.dirs(path = file.path(init$data_folder, "ChromSCape_analyses"), full.names = FALSE, recursive = FALSE)
           init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())
           print("observeEvent({ input$path_cookie}")
           print(init$available_reduced_datasets)
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
      
      init$available_analyses <- list.dirs(
        path = file.path(init$data_folder, "ChromSCape_analyses"),
        full.names = FALSE, recursive = FALSE)
      init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())
      print("observeEvent({ input$data_folder}")
      print(init$available_reduced_datasets)
      if(.Platform$OS.type != "windows"){
        js$save_cookie(init$data_folder)
      }
    }
  )
  
  
  output$input_data_ui <- renderUI({
    if(input$data_choice_box== "count_mat"){
      column(12, br(),fileInput("datafile_matrix", "Upload all data matrices (.txt, .tsv or .csv) :",
                multiple=TRUE, accept=c("text", "text/plain", ".txt", ".tsv", ".csv", ".gz")),
             checkboxInput("is_combined_mat", "Single Multi-sample count matrix ?",value = FALSE),
             uiOutput("nb_samples_mat")
             )
      
    }
    else{
      column(12,
             br(),
             HTML(paste0("<b>Upload folder (", input$data_choice_box,")</b><br>")),
             shinyFiles::shinyDirButton(id = "datafile_folder", label = "Select folder containing raw files",
                                        title =  paste0("Select a directory containing your ",
                                                        input$data_choice_box," files."),
                                        icon = icon("folder-open")),
             selectInput(inputId = "nb_samples_to_find",label = "Number of samples:",
                         choices = 1:100,selected = 1,multiple = FALSE)
             )
    }
    
  })
  
  output$advanced_data_input <- renderUI({
   if(input$data_choice_box != "count_mat"){
    column(12,
           shinydashboard::box(title="Counting parameters", width = NULL, status="success", solidHeader = TRUE,
                               column(6, 
                                      radioButtons("count_on_box", label = "Select a count method",
                                            choices = list("Count on bins (width)"="bin_width",
                                                           "Count on bins (number of bins)" = "n_bins",
                                                           "Count on peaks (must provide a .bed file)" = "peak_file",
                                                           "Count around gene TSS" = "geneTSS")),
                                      
                               ),
                               column(6,
                                      uiOutput("bin_width"),
                                      uiOutput("n_bins"),
                                      uiOutput("peak_file"),
                                      uiOutput("aroundTSS"))
           )
    )
   }
  })
  
  output$nb_samples_mat <- renderUI({ if(input$is_combined_mat == TRUE){
    selectInput(inputId = "nb_samples_to_find",label = "Number of samples:",
                choices = 1:100,selected = 1,multiple = FALSE)
  }})
  
  output$bin_width <- renderUI({ if(input$count_on_box == "bin_width"){
    textInput("bin_width", label = "Width of bins to count on (in bp) :",value = 50000)
  }})
  output$n_bins <- renderUI({ if(input$count_on_box == "n_bins" ){
    textInput("n_bins", label = "Number of bins to count on :", value = 10000)
  }})
  output$peak_file <- renderUI({ if(input$count_on_box == "peak_file"){
    fileInput("peak_file", ".bed file containing the peaks to count on:", multiple = FALSE, accept = c(".bed",".txt"))
  }})
  output$aroundTSS <- renderUI({ if(input$count_on_box == "geneTSS" ){
    textInput("aroundTSS", label = "Distance Up/Downstream of TSS(bp):", value = 2500)
  }})

  shinyFiles::shinyDirChoose(input, "datafile_folder", roots = volumes, session = 
                               session)
  
  observeEvent(input$create_analysis, {  # save new dataset
    req(input$new_analysis_name, input$annotation)
    if(is.null(input$datafile_folder) & is.null(input$datafile_matrix)) return()
    
    datamatrix <- NULL
    annot_raw <- NULL
    type_file = as.character(input$data_choice_box)
    if(dir.exists(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name))){
      showNotification(paste0("Warning : The name : '",input$new_analysis_name,
                              " is already taken by a preexisting analysis. Please
                              choose another name for your analysis."),
                              duration = 5, closeButton = TRUE, type="warning")
    }else{
      withProgress(message='Creating new data set...',value = 0, {
        
        if(type_file == "count_mat" & !is.null(input$datafile_matrix)){
          incProgress(0.3, detail="Reading count matrices")
          if(input$is_combined_mat == TRUE){
            if(length(input$datafile_matrix$name)>1){
              showNotification(paste0("Warning : When checking the 
                                      'The matrix contains multiple samples ?' button,
                                      you have to input a single count matrix."),
                               duration = 5, closeButton = TRUE, type="warning")
              return()
            }
          }
          tmp_list = import_scExp(file_names = input$datafile_matrix$name,
                                  path_to_matrix = input$datafile_matrix$datapath)
          datamatrix = tmp_list$datamatrix
          if(input$is_combined_mat == TRUE) {
            samples_ids = detect_samples(colnames(datamatrix),
                                         nb_samples = as.numeric(input$nb_samples_to_find))
            annot_raw = data.frame(barcode = colnames(datamatrix),
                                   cell_id = colnames(datamatrix),
                                   sample_id = samples_ids,
                                   batch_id = factor(rep(1, ncol(datamatrix))))
          } else{ annot_raw = tmp_list$annot_raw }
        }
        else if(type_file %in% c("BAM","BED","Index_Peak_Barcode") & !is.null(input$datafile_folder)) {
          datafile_folder = shinyFiles::parseDirPath(volumes, input$datafile_folder)
          send_warning = FALSE
          if(type_file == "BAM") if(length(list.files(datafile_folder,pattern = "*.bam$"))==0) send_warning = TRUE
          if(type_file == "BED") if(length(list.files(datafile_folder,pattern = "*.bed$|.*.bed.gz"))==0) send_warning = TRUE
          if(type_file == "Index_Peak_Barcode") 
            if(length(list.files(datafile_folder,pattern = "*.index.txt|.*.barcodes.txt|.*.peak.bed"))==0) send_warning = TRUE
          
          if(send_warning) {
            showNotification(paste0("Warning : Can't find any specified file types in the upload folder. 
                                    Select another upload folder or another data type."),
                             duration = 5, closeButton = TRUE, type="warning")
            return()
          }
          incProgress(0.2, detail=paste0("Reading ",type_file," files to create matrix. This might take a while."))
          
          if(input$count_on_box == "bin_width") datamatrix = ChromSCape:::raw_counts_to_feature_count_files(
            files_dir = datafile_folder,
            file_type = type_file,
            bin_width = as.numeric(input$bin_width),
            ref = input$annotation)
          
          if(input$count_on_box == "n_bins") datamatrix = ChromSCape:::raw_counts_to_feature_count_files(
            files_dir = datafile_folder,
            file_type = type_file,
            n_bins = as.numeric(input$n_bins),
            ref = input$annotation)
          
          
          if(input$count_on_box == "peak_file") datamatrix = ChromSCape:::raw_counts_to_feature_count_files(
            files_dir = datafile_folder,
            file_type = type_file,
            peak_file = as.character(input$peak_file$datapath),
            ref = input$annotation)
          
          if(input$count_on_box == "geneTSS") datamatrix = ChromSCape:::raw_counts_to_feature_count_files(
            files_dir = datafile_folder,
            file_type = type_file,
            geneTSS = TRUE,
            aroundTSS = as.numeric(input$aroundTSS),
            ref = input$annotation)
          
          incProgress(0.3, detail=paste0("Finished creating matrix, assigning sample labels to cells heuristically."))
          samples_ids = detect_samples(colnames(datamatrix), nb_samples = as.numeric(input$nb_samples_to_find))
          annot_raw = data.frame(barcode = colnames(datamatrix),
                                 cell_id = colnames(datamatrix),
                                 sample_id = samples_ids,
                                 batch_id = factor(rep(1, ncol(datamatrix)))
          )
          
                  
        } else {
          stop("No data folder or data files selected.")
        }
        incProgress(0.4, detail="Saving matrix & annotation...")
        dir.create(file.path(init$data_folder, "ChromSCape_analyses"), showWarnings = FALSE)
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name))
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "Filtering_Normalize_Reduce"))
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "correlation_clustering"))
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "correlation_clustering","Plots"))
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "Diff_Analysis_Gene_Sets"))
        write.table(input$annotation, file.path(init$data_folder, 'ChromSCape_analyses', input$new_analysis_name, 'annotation.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        qs::qsave(datamatrix, file = file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "datamatrix.qs"))
        qs::qsave(annot_raw, file = file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "annot_raw.qs"))
        
        init$available_analyses <- list.dirs(path = file.path(init$data_folder, "ChromSCape_analyses"), full.names = FALSE, recursive = FALSE)

        updateSelectInput(session = session, inputId = "selected_analysis",
                          label =  "Select an Analysis:",
                          choices = init$available_analyses,
                          selected =  input$new_analysis_name)

        init$datamatrix <- datamatrix
        init$annot_raw <- annot_raw
        incProgress(0.1, detail="Import successfully finished! ")
        updateActionButton(session, "create_analysis", label="Added successfully", icon = icon("check-circle"))
      })
    }
  })
  
  observeEvent(input$selected_analysis, {  # load precompiled dataset and update coverage plot
    if(!is.null(input$selected_analysis) & input$selected_analysis != ""){
      init$datamatrix <- qs::qread(file.path(init$data_folder, "ChromSCape_analyses", input$selected_analysis, "datamatrix.qs"))
      init$annot_raw <-  qs::qread(file.path(init$data_folder, "ChromSCape_analyses", input$selected_analysis, "annot_raw.qs"))
    }
  })
  
  output$rename_file_box <- renderUI({
    shinydashboard::box(title = tagList("Rename Samples", shiny::icon("signature")),
                        width = NULL, status = "success", solidHeader = TRUE, 
                        collapsible = TRUE, collapsed = TRUE,
                        column(8, htmlOutput("rename_file_UI")),
                        br(),
                        column(6, actionButton("save_rename", "Rename", icon = icon("save")))
                        )
  })
  
  output$rename_file_UI <- renderUI({
    req(analysis_name())
    print("Inside rename_file_UI")
    if(nrow(init$datamatrix) > 0 & nrow(init$annot_raw) > 0 ){
      sample_names = unique(init$annot_raw$sample_id)
      print(sample_names)
      lapply(sample_names, function(i) {
        shiny::textInput(inputId = paste0("sample_", i),
                         placeholder = i,
                         label = i,
                         value = i
        )
      })
    }
  })
  
  observeEvent(input$save_rename, {  
    req(analysis_name(), init$data_folder, init$datamatrix, init$annot_raw)
    
    if(nrow(init$datamatrix) > 0 & nrow(init$annot_raw) > 0 ){
      annot_raw = init$annot_raw
      datamatrix = init$datamatrix
      sample_names = unique(annot_raw$sample_id)
      for(i in sample_names) {
        expr <- paste0("re_name = input$sample_", i)
        eval(parse(text = expr))

        if(re_name != i){
          cell_idxs = which( annot_raw$sample_id == i)
          annot_raw[cell_idxs, "sample_id"] = re_name
          annot_raw[cell_idxs, "cell_id"] = gsub(i,"",annot_raw[cell_idxs, "cell_id"])
          annot_raw[cell_idxs, "cell_id"] = paste0(re_name,annot_raw[cell_idxs, "cell_id"])
          rownames(annot_raw) = annot_raw$cell_id

          colnames(datamatrix)[cell_idxs] = annot_raw$cell_id[cell_idxs]
        }
      }
     
      init$annot_raw <- annot_raw
      init$datamatrix <- datamatrix

      qs::qsave(datamatrix, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "datamatrix.qs"))
      qs::qsave(annot_raw, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "annot_raw.qs"))
      
      if(length(list.files(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering"), pattern = ".qs"))>0){
        showNotification(paste0("Please re-run analysis from filtering in order to rename downstream analysis."),
                         duration = 15, closeButton = TRUE, type="warning")
        
      }
      
      print("Renaming file done !")
      
    }
  })

  observeEvent(input$new_analysis_name, {  # reset label on actionButtion when new dataset is added
    updateActionButton(session, "create_analysis", label="Create Analysis", icon = character(0))
  })
  
  observeEvent(input$selected_analysis,
               {
                 if(nchar(input$selected_analysis)>0){
                   unlocked$list$selected_analysis=TRUE
                   }
                 else{
                   for(i in names(unlocked$list)){
                     unlocked$list[[i]]=FALSE
                   }
                 }
               })
  observeEvent(input$selected_analysis,
          {
    unlocked$list

    able_disable_tab("selected_analysis","filter_normalize")}
    ) 
  
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
    fileInput("exclude_file", ".bed / .bed.gz / .txt file containing the feature to exclude from data set:", multiple = FALSE, accept = c(".bed",".txt", ".bed.gz"))
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
    print(input$exclude_file)
    if(input$exclude_regions) {
      if(!is.null(input$exclude_file) && file.exists(as.character(input$exclude_file$datapath))){
        tmp_file = tempfile(fileext = ".bed.gz")
        file.copy(input$exclude_file$datapath, tmp_file, overwrite = T)
        exclude_regions <- rtracklayer::import(tmp_file)
      } else {
        warning("The BED file you specified doesn't exist. No specific features will be
                removed during filtering.")
        exclude_regions <- NULL
      }
    }
    else {
      exclude_regions <- NULL
    }

    callModule(Module_preprocessing_filtering_and_reduction, "Module_preprocessing_filtering_and_reduction", reactive({input$selected_analysis}), reactive({input$min_coverage_cell}),
               reactive({input$min_cells_window}), reactive({input$quant_removal}), reactive({init$datamatrix}), reactive({init$annot_raw}),
               reactive({init$data_folder}),reactive({annotationId}), reactive({exclude_regions}) ,reactive({input$do_batch_corr}),
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
  
  observeEvent(input$selected_reduced_dataset, { # load reduced data set to work with on next pages
    req(input$selected_reduced_dataset)
    
    if(file.exists(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots"))) 
      addResourcePath('Plots', file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots"))

    file_index <- match(c(input$selected_reduced_dataset), reduced_datasets())
    filename_sel <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),"Filtering_Normalize_Reduce",init$available_reduced_datasets[file_index])
    
    
    t1 = system.time({
    scExp. = qs::qread(filename_sel)
    if(is.reactive(scExp.)) {
      scExp. = isolate(scExp.())
    }
    scExp(scExp.) # retrieve filtered scExp
    rm(scExp.)
    gc()
    })
    cat("Loaded reduced data in ",t1[3]," secs\n")
  })


  cell_cov_df <- reactive ({
    req(init$datamatrix)
    df = data.frame(coverage = sort(unname(Matrix::colSums(init$datamatrix)))) 
    df
    })  # used for plotting cell coverage on first page
  
  quantile_threshold =  reactive({
    req(init$datamatrix)
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
                                                                 tooltip="Sample", dynamicTicks=TRUE) )
  
  
  output$num_cell <- function(){
    req(init$annot_raw)
    tab = num_cell_scExp(init$annot_raw)
    tab
  }
  
  output$num_cell_after_QC_filt <- function(){
    req(input$selected_reduced_dataset,scExp())
    tab = num_cell_after_QC_filt_scExp(scExp(),init$annot_raw)
    tab
  }
  
  output$table_QC_filt_box <- renderUI({
    if(!is.null(input$selected_reduced_dataset) && 
       input$selected_reduced_dataset != ""){
          column(12, align="left", tableOutput("num_cell_after_QC_filt"))
    } else {
      column(12, align="left", tableOutput("num_cell"))
    }
  })
  
  observeEvent(input$selected_reduced_dataset,{

    if(suppressWarnings(nchar(input$selected_reduced_dataset)>0)){
      unlocked$list$selected_reduced_dataset=TRUE
      }else{
        for(i in names(unlocked$list)){unlocked$list[[i]]=FALSE}
        }
    })
  observeEvent(unlocked$list,{
    
    able_disable_tab(c("selected_reduced_dataset"),"vizualize_dim_red")}) 
  
  
  ###############################################################
  # 2. PCA
  ###############################################################
  
  output$pc_select_x <- renderUI({ selectInput("pc_select_x", "X",choices=paste0("Component_", c(1:15)), selected="Component_1") })
  output$pc_select_y <- renderUI({ selectInput("pc_select_y", "Y",choices=paste0("Component_", c(1:15)), selected="Component_2") })

  pca_plot <- reactive({
    req(scExp(), annotCol(), input$pc_select_x,input$pc_select_y,  input$color_by)
    if(input$color_by %in% colnames(SingleCellExperiment::colData(scExp())) ){
    p = plot_reduced_dim_scExp(scExp(),input$color_by, "PCA",
                               select_x = input$pc_select_x,
                               select_y = input$pc_select_y
    )
    unlocked$list$pca=TRUE
    p
    }
  })
  output$pca_plot <- renderPlot(pca_plot())
  
  contrib_features_plot <- reactive({
    req(scExp(),  input$pc_select_x)
    if(input$pc_select_x %in% colnames(SingleCellExperiment::reducedDim(scExp(), "PCA"))){
      p = plot_most_contributing_features(scExp(), component = input$pc_select_x,
                                          n_top_bot = as.numeric(input$n_top_features))
      p
    }
  })
  
  output$contrib_features_plot <- renderPlotly(contrib_features_plot())
  
  contrib_chr_plot <- reactive({
    req(scExp(),  input$pc_select_x)
    if(input$pc_select_x %in% colnames(SingleCellExperiment::reducedDim(scExp(), "PCA"))){
      p = plot_pie_most_contributing_chr(scExp(), component = input$pc_select_x)
      p
    }
  })
  output$contrib_chr_plot <- renderPlot(contrib_chr_plot())
  
  output$contribution_to_pca_UI <- renderUI({
    req(scExp(), input$pc_select_x)
    if("PCA" %in% SingleCellExperiment::reducedDimNames(scExp())){
      if(input$pc_select_x %in% colnames(SingleCellExperiment::reducedDim(scExp(), "PCA"))){
        shinydashboard::box(title="Contribution to PCA", width = NULL, status="success", solidHeader=TRUE,
                            column(12, align="left",
                                   selectInput("n_top_features", label =  "Top features",
                                               choices = 5:100,
                                               multiple = FALSE, selected = 10),
                                   h3(paste0("Most contributing features to '", input$pc_select_x,"' .")),
                                   plotlyOutput("contrib_features_plot") %>% 
                                     shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75)
                            ),
                            column(12, align="left", 
                                   h3("Contribution of chromosomes in top 100 features."),
                                   plotOutput("contrib_chr_plot") %>% 
                                     shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75)
                                   ))
      }
    }
  })
  
  output$tsne_box <- renderUI({
    req(scExp(), annotCol(), input$color_by)
    if("TSNE" %in% SingleCellExperiment::reducedDimNames(scExp())){
      if(input$color_by %in% colnames(SingleCellExperiment::colData(scExp())) ){
      p = plot_reduced_dim_scExp(scExp(),input$color_by, "TSNE")
      output$tsne_plot = renderPlot(p)
      shinydashboard::box(title="t-SNE vizualisation 1", width = NULL, status="success", solidHeader=TRUE,
                          column(12, align="left", plotOutput("tsne_plot") %>% 
                                   shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75) %>%
                                   shinyhelper::helper(type = 'markdown', icon ="info-circle",
                                                       content = "tsne_plot")
                                 ))
      }
    }
  })

  output$UMAP_plot <- renderPlot({
    req(scExp(), annotCol(), input$color_by)
    if(input$color_by %in% colnames(SingleCellExperiment::colData(scExp())) ){
      p = plot_reduced_dim_scExp(scExp(), input$color_by, "UMAP")
      p
    }
     })
  
  output$color_box <- renderUI({
    req(input$color_by)
    if(input$color_by != 'total_counts'){
      shinydashboard::box(title = tagList("Color settings ",shiny::icon("palette")),
                          width = NULL, status = "success", solidHeader = TRUE,
          column(6, htmlOutput("color_picker")),
          column(6 , br(), actionButton("col_reset", "Default colours", icon = icon("undo")),
                 br(), br(), actionButton("save_color", "Save colors & apply to all", icon = icon("save"))))
    }
  })
  
  observeEvent(input$col_reset, {
    cols <- ChromSCape:::gg_fill_hue(length(levels_selected()))
    for(i in seq_along(levels_selected())){
      colourpicker::updateColourInput(session=session, inputId=paste0("color_", levels_selected()[i]),
                                      value=cols[i])
    }
  })
  
  output$color_picker <- renderUI({
    #Color picker
    if(input$color_by != "total_counts"){
      if(input$color_by %in% colnames(SingleCellExperiment::colData(scExp())) ){
      colsModif <- SummarizedExperiment::colData(scExp())[,c(input$color_by,paste0(input$color_by,"_color"))] %>% unique()
      lapply(seq_along(levels_selected()), function(i) {
        colourpicker::colourInput(inputId=paste0("color_", levels_selected()[i]),
                                  label=paste0("Choose colour for ", levels_selected()[i]),
                                  value=colsModif[i,paste0(input$color_by,"_color")], returnName=TRUE) ## Add ", palette = "limited"" to get selectable color panel       
      })
      }
    }
  })
  
  observeEvent(input$save_color, {  
    req(scExp(), input$color_by)
  
    color_df = ChromSCape:::get_color_dataframe_from_input(input,levels_selected(),input$color_by)

    scExp. = colors_scExp(scExp(),annotCol = input$color_by,color_by = input$color_by, color_df = color_df)
    scExp(scExp.)
    
    qs::qsave(scExp, file = file.path(init$data_folder, "ChromSCape_analyses",
                                   analysis_name(), "Filtering_Normalize_Reduce",
                                   paste0(input$selected_reduced_dataset,".qs")))
    
    rm(scExp.)
    rm(color_df)
    gc()
  })
  
  levels_selected <- reactive({
    req(scExp(),input$color_by)
    if(input$color_by != "total_counts") 
      levels_selected = SummarizedExperiment::colData(scExp())[,input$color_by] %>% unique() %>% as.vector()
    else NULL
  })

  observeEvent(unlocked$list,able_disable_tab(c("pca"),"cons_clustering")) # if conditions are met, unlock tab Correlation Clustering
  
  
  ###############################################################
  # 3. Consensus clustering on correlated cells
  ###############################################################

  corColors <- grDevices::colorRampPalette(c("royalblue","white","indianred1"))(256)
  
  selected_filtered_dataset <- reactive({input$selected_reduced_dataset})
  
  annotCol_cf <- reactive({
    req(scExp_cf())
    cols = c("sample_id")
    if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
      cols = c(cols,"cell_cluster")
      } 
    if("batch_name" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
      cols = c(cols,"batch_name")
      } 
    cols = c(cols, "total_counts")
    cols
  })
  
  observeEvent(input$tabs, 
               {
               if(input$tabs == "cons_clustering"){
                 file = file.path(init$data_folder, "ChromSCape_analyses",
                                  analysis_name(), "correlation_clustering",
                                  paste0(selected_filtered_dataset(),".qs"))
                 if(file.exists(file)){
                   data = qs::qread(file)
                   scExp_cf(data$scExp_cf)
                   rm(data)
                   gc()

                 } else {
                   scExp_cf(scExp())
                   gc()
                 }
                 output$hc_heatmap_plot <- renderPlot({
                   plot_heatmap_scExp(scExp(), color_by = annotCol())
                 })
               }
               })
  
  cluster_type = reactive({
    if(!is.null(scExp_cf())){
      if("consclust" %in% names(scExp_cf()@metadata)) {
        input$cluster_type
      } else {
        FALSE
      }
    } else{
      showNotification("Run Consensus Hiearchical Clustering first..",type="warning")
      updateCheckboxInput(session,"cluster_type",value = FALSE)
      FALSE
    }
    
  })
  output$nclust_UI = renderUI({
    selectInput("nclust", br("Number of Clusters:"), choices=c(2:input$maxK))
  })
  
  observeEvent({input$choose_cluster},{
    
                 req(input$nclust, scExp_cf())
                 if(input$nclust != ""){
                   scExp_cf(choose_cluster_scExp(scExp_cf(), nclust = as.numeric(input$nclust),
                                                 consensus = cluster_type()))
                   unlocked$list$cor_clust_plot=TRUE;
                   unlocked$list$affectation=TRUE;
                   gc()
                   file = file.path(init$data_folder, "ChromSCape_analyses",
                                    analysis_name(), "correlation_clustering",
                                    paste0(selected_filtered_dataset(),".qs"))
                   if(!file.exists(file)){
                     data = list("scExp_cf" = scExp_cf())
                     qs::qsave(data,file=file )
                     rm(data)
                     gc()

                   }
                 }
               })
  
  

  # output$corr_clust_pca_plot <- renderPlot(hc_pca_plot())
  
  # output$cons_corr_clust_pca_plot <- renderPlot(plot_heatmap_scExp(scExp_cf()))
  # 
  # output$cons_corr_clust_pca_UI <- renderUI({
  #   if(consensus_ran()){
  #     plotOutput("cons_corr_clust_pca_plot", height = 500, width = 500) %>%
  #       shinyhelper::helper(type = 'markdown', icon ="info-circle",
  #                           content = "correlation_clustering")
  #   }
  # })

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
      scExp_cf(choose_cluster_scExp(scExp_cf(), nclust = as.numeric(input$nclust), consensus = cluster_type()))
      gc()
      incProgress(amount=0.2, detail=paste("Saving"))
      data = list("scExp_cf" = scExp_cf())
      qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                              paste0(input$selected_reduced_dataset, ".qs")))
      incProgress(amount=0.2, detail=paste("Finished"))
      rm(data)
      gc()
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
    tab = num_cell_before_cor_filt_scExp(scExp())
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
  
  
  output$filtered_data_selection_format <- renderText({"The name of the filtered dataset is composed of the following information: data set name, min percentage of reads per cell, 
    min percentage of cells to support a window, quantile of cell read counts to keep, correlation threshold, percent of cell correlation. To work on a different dataset or different preprocessing state, select it on the first page."})
  output$select_n_clust_hc = renderUI({
    selectInput("nclust", "Select number of clusters:", choices=c(2:10))
    })
  output$select_n_clust_chc = renderUI({
    selectInput("nclust", "Select number of clusters:", choices=c(2:10))
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
      gc()
      data = list("scExp_cf" = scExp_cf())
      qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                  paste0(input$selected_reduced_dataset, ".qs")))
      rm(data)
      gc()
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
  
  observeEvent(input$do_annotated_heatmap_plot,
               {
                 if(input$nclust != ""){
                   output$annotated_heatmap_UI <- renderUI({
                     
                     output$annotated_heatmap_plot = renderPlot(plot_heatmap_scExp(isolate(scExp_cf())))
                     plotOutput("annotated_heatmap_plot",width = 500,height = 500) %>%
                       shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75)
                 })
               }
               }
  )
  
  output$violin_color <- renderUI({
    req(annotCol_cf())
    selectInput("violin_color","Color by", choices=annotCol_cf()[-which(annotCol_cf() == "total_counts")])
  })
  
  output$jitter_color <- renderUI({
    req(annotCol_cf(), input$add_jitter)
    if(input$add_jitter) selectInput("jitter_color","Color cells", choices=annotCol_cf())
  })
  
  output$intra_corr_UI <- renderUI({
    req(annotCol_cf(), input$violin_color)
   
    if(input$add_jitter) jitter_col = input$jitter_color else jitter_col = NULL
     
   output$intra_corr_plot = renderPlotly(
     plot_intra_correlation_scExp(scExp_cf(), by = input$violin_color,
                                  jitter_by = jitter_col, seed = 47))
    plotly::plotlyOutput("intra_corr_plot",width = 500,height = 500) %>%
      shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75)
    
    })
  
  output$reference_group <- renderUI({
    req(annotCol_cf(), scExp_cf(), input$violin_color)
    selectInput("reference_group","Correlate with",
                choices = unique(scExp_cf()[[input$violin_color]]) )
  })
  

  output$inter_corr_UI <- renderUI({
    req(annotCol_cf(), input$violin_color)
    
    if(input$add_jitter) jitter_col = input$jitter_color else jitter_col = NULL
    output$inter_corr_plot = renderPlotly(
      plot_inter_correlation_scExp(scExp_cf(), by = input$violin_color,
                                   jitter_by = jitter_col,
                                   reference_group = input$reference_group,
                                   seed = 47))
    
    plotly::plotlyOutput("inter_corr_plot",width = 500,height = 500) %>%
      shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75)
    
  })
  
#   
# output$cons_clust_anno_plot <- renderPlot({
#   if(! is.null(scExp_cf())){
#     if("ConsensusAssociation" %in% names(scExp_cf()@metadata)){
#       colors <- SummarizedExperiment::colData(scExp_cf())[scExp_cf()@metadata$hc_consensus_association$order,"cell_cluster_color"]
#       heatmap(SingleCellExperiment::reducedDim(scExp_cf(),"ConsensusAssociation")[scExp_cf()@metadata$hc_consensus_association$order,],
#               Colv = as.dendrogram(scExp_cf()@metadata$hc_consensus_association),
#               Rowv = NA, symm = FALSE, scale="none", col = grDevices::colorRampPalette(c("white", "blue"))(100),
#               na.rm = TRUE, labRow = FALSE, labCol = FALSE, mar = c(5, 5), main = paste("consensus matrix k=", input$nclust, sep=""),
#               ColSideCol = colors)
#     }
#   }
#     })
#   
# output$anno_cc_box <- renderUI({
#   if(! is.null(scExp_cf())){
#     if("ConsensusAssociation" %in% names(scExp_cf()@metadata)){
#       shinydashboard::box(title="Annotated consensus clustering", width = NULL, status="success", solidHeader = TRUE,
#           column(12, align="left", plotOutput("cons_clust_anno_plot", height = 500, width = 500),
#                  downloadButton("download_anno_cc_plot", "Download image")))
#     }
#   }
#   })
#   
#   output$download_anno_cc_plot <- downloadHandler(
#     filename = function(){ paste0("consensus_clustering_k", input$nclust, "_", selected_filtered_dataset(), ".png")},
#     content = function(file){
#       grDevices::png(file, width = 1200, height = 800, res = 300)
#       colors <- SummarizedExperiment::colData(scExp_cf())[scExp_cf()@metadata$hc_consensus_association$order,"cell_cluster_color"]
#       heatmap(SingleCellExperiment::reducedDim(scExp_cf(),"ConsensusAssociation")[scExp_cf()@metadata$hc_consensus_association$order,],
#               Colv = as.dendrogram(scExp_cf()@metadata$hc_consensus_association),
#               Rowv = NA, symm = FALSE, scale="none", col = colorRampPalette(c("white", "blue"))(100),
#               na.rm = TRUE, labRow = FALSE, labCol = FALSE, mar = c(5, 5), main = paste("consensus matrix k=", input$nclust, sep=""),
#               ColSideCol = colors)
#       grDevices::dev.off()
#   })
#   
  output$contingency_table_cluster <- renderUI({
    if(! is.null(scExp_cf())){
      if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
        shinydashboard::box(title="Samples & Cluster association table", width = NULL, status="success", solidHeader = TRUE,
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
    req(scExp_cf(), annotCol(), input$color_by)
    if(!is.null(scExp_cf())){
      if("TSNE" %in% SingleCellExperiment::reducedDimNames(scExp_cf())){
        req(scExp_cf(), input$color_by_cf)
        if(input$color_by_cf %in% colnames(SingleCellExperiment::colData(scExp_cf())) ){
        p = plot_reduced_dim_scExp(scExp_cf(),input$color_by_cf, "TSNE",
                                   select_x = "Component_1",
                                   select_y = "Component_2") +
          ggtitle("t-SNE")
        output$tsne_plot_cf <- renderPlot(p)
        shinydashboard::box(title="t-SNE", width = NULL, status="success", solidHeader=TRUE,
                            column(12, align="left", plotOutput("tsne_plot_cf")))
        }
      }
    }
  })
  
  umap_p_cf <- reactive({
    req(scExp_cf(), annotCol(), input$color_by_cf)
    if(input$color_by_cf %in% colnames(SingleCellExperiment::colData(scExp_cf())) ){
      p = plot_reduced_dim_scExp(scExp_cf(),input$color_by_cf, "UMAP",
                                 select_x = "Component_1",
                                 select_y = "Component_2")
      p
    }
  })

  output$plot_CF_UMAP <- renderPlot({
    umap_p_cf()
    })

  
  levels_selected_cf <- reactive({
    req(scExp_cf(),input$color_by_cf)
    if(input$color_by_cf != "total_counts") levels_selected_cf =
      SummarizedExperiment::colData(scExp_cf())[,input$color_by_cf] %>%
      unique() %>% as.vector() else NULL
  })
  
  output$UMAP_box <- renderUI({
    if(! is.null(scExp_cf())){
      if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf())) ){
        shinydashboard::box(title="UMAP vizualisation 2", width = NULL, status="success", solidHeader = TRUE,
                            column(6, align="left",
                                   selectInput("color_by_cf", "Color by",
                                               selected = "cell_cluster",
                                               choices = c(annotCol(),'cell_cluster'))),
                            column(12, align="left", plotOutput("plot_CF_UMAP")))
      }
    }
  })
  
  output$color_box_cf <- renderUI({
    req(input$color_by_cf)
    if(! is.null(scExp_cf())){
      if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
        if(input$color_by_cf != 'total_counts'){
          shinydashboard::box(title=tagList("Color settings ",shiny::icon("palette")), width = NULL, status = "success", solidHeader = TRUE,
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
    color_df = ChromSCape:::get_color_dataframe_from_input(input,levels_selected_cf(), input$color_by_cf, "color_cf_")
    scExp_cf. = colors_scExp(scExp_cf(), annotCol = input$color_by_cf,
                             color_by = input$color_by_cf, color_df = color_df)
    scExp_cf(scExp_cf.)
    rm(scExp_cf.)
  })
  
  observeEvent(input$col_reset_cf, {
    cols <- ChromSCape:::gg_fill_hue(length(levels_selected_cf()))
    for(i in seq_along(levels_selected_cf())){
      colourpicker::updateColourInput(session=session, inputId=paste0("color_cf_", levels_selected_cf()[i]),
                                      value=cols[i])
    }
  })

  output$color_picker_cf <- renderUI({
    #Color picker
    if( (input$color_by_cf %in% colnames(SummarizedExperiment::colData(scExp_cf()))) &
        (paste0(input$color_by_cf,"_color") %in% colnames(SummarizedExperiment::colData(scExp_cf()))) ) {
    colsModif <- as.data.frame(SummarizedExperiment::colData(scExp_cf()))[,c(input$color_by_cf,paste0(input$color_by_cf,"_color"))] %>% unique()
    lapply(seq_along(levels_selected_cf()), function(i) {
      colourpicker::colourInput(inputId=paste0("color_cf_", levels_selected_cf()[i]),
                                label=paste0("Choose colour for ", levels_selected_cf()[i]),
                                value=colsModif[i,paste0(input$color_by_cf,"_color")],
                                returnName = TRUE) ## Add ", palette = "limited"" to get selectable color panel       
    })
    }
  })
  
  # observeEvent({
  #   selected_filtered_dataset()
  #   scExp_cf()
  # },
  # {
  #   if(length(selected_filtered_dataset())>0 & !is.null(scExp_cf())){
  #     if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
  #       unlocked$list$affectation = TRUE
  #     } else {
  #       unlocked$list$affectation = FALSE
  #     }
  #   } else unlocked$list$affectation = FALSE
  # }
  # )
  
  observeEvent(unlocked$list, {able_disable_tab(c("selected_reduced_dataset","affectation"),c("peak_calling","diff_analysis"))}) 
  
  ###############################################################
  # 5. Peak calling [optional]
  ###############################################################
  
  output$peak_calling_info <- renderText({"This module is optional, but recommended in order to obtain the most meaningful results for pathway enrichment analysis. Peaks will be called from the BAM files of the samples selected in your project, using MACS2 [only works on unix systems] so that counts can be assigned more specifically to genes TSS . If you have MACS2 installed but ChromSCape can’t find these softwares, try relaunching R from the terminal and start ChromSCape again."})
  
  can_run = reactiveVal({FALSE})
  
  output$peak_calling_system <- renderText({
    platform = as.character(.Platform[1])
    if(length(grep("nix",platform,ignore.case = TRUE)) ){
      macs2=""
      try({
        macs2 = system2("which",args = "macs2",stdout = TRUE)
      })
      if(length(macs2)>0){
        can_run(TRUE)
        return(paste0("<b>You are running on an ", platform, " OS.<br>Found MACS2 at ", macs2))
      }
      if(length(macs2)==0)
        return(paste0("<b>You are running on an ", platform, " OS.<br>Didn't find MACS2, please install MACS2 or skip this step."))
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
   
  output$select_raw_files <- renderUI({
    req(scExp_cf(), list_bams())
    
  })  
  
  shinyFiles::shinyDirChoose(input, "pc_folder", roots = volumes, session = 
                               session)
  
  pc_folder = reactiveVal(NULL)
  
  list_files_pc = reactive({
    req(pc_folder())
    if(!is.null(pc_folder())){
      list.files(pc_folder(), full.names = TRUE, pattern = "*.bam$|*.bed|*.bed.gz")
    }
  })
  
  bam_or_bed <- reactive({
    req(list_files_pc())
    nBAMs = length(grep(".*.bam$", list_files_pc()))
    nBEDs = length(grep(".*.bed.gz|.*.bed", list_files_pc()))
    print("is bam or bed ?")
    print(nBAMs)
    print(nBEDs)
    ifelse(nBAMs>nBEDs,"BAM","BED")
  })
  
  observeEvent(input$pc_folder,
               {
                 req(input$pc_folder)
                 if(!is.null(input$pc_folder)){
                   # browser()
                   pc_folder(shinyFiles::parseDirPath(volumes, input$pc_folder))
                   output$pc_dir <- renderText(pc_folder())
                 }
               })
  
  output$pc_upload <- renderUI({
    req(list_files_pc(), bam_or_bed())
    print("Inside output$pc_upload UI")
    if(!is.null(list_files_pc()) & !is.null(bam_or_bed())){
      if(bam_or_bed() == "BAM") {
        files = list_files_pc()[grep(".bam$", list_files_pc())]
        print("Inside output$pc_upload UI... BAM")
        selectInput("bam_selection", label = "Selected BAM files",
                    choices = basename(files), multiple = TRUE,
                    selected = basename(files) )
      } else {
        print("Inside output$pc_upload UI... BED")
        NULL
      }
    }
  })
  
  observeEvent(input$do_pc, {
    req(bam_or_bed(), list_files_pc(), scExp_cf())
    print("Doing peak call")
    print(bam_or_bed())
    if(bam_or_bed() == "BAM") {
    input_files_pc <- as.character(file.path(dirname(list_files_pc()),input$bam_selection))
    } else {
      input_files_pc = as.character(list_files_pc()[grep(".bed|.bed.gz", list_files_pc())])
    }
    if(length(input_files_pc)==0){
      warning("Can't find any input BAM / BED files.")
    } else{
      
      nclust = length(unique(scExp_cf()$cell_cluster))  
      withProgress(message='Performing enrichment analysis...', value = 0, {
        
        incProgress(amount = 0.1, detail = paste("Starting Peak Calling..."))
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks"), showWarnings = FALSE)
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", nclust)), showWarnings = FALSE)
        odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", nclust))
        sample_ids <- unique(SummarizedExperiment::colData(scExp_cf())$sample_id)
        # inputBams <- as.character(unlist(sapply(sample_ids, function(x){ input[[paste0('bam_', x)]] })))

        
        checkFiles <- sapply(input_files_pc, function(x){ if(file.exists(x)){ 0 } else { 
          showNotification(paste0("Could not find file ", x, ". Please make sure to give a full path including the file name."),
                           duration = 7, closeButton = TRUE, type="warning"); 1}  
        })
        incProgress(amount = 0.3, detail = paste("Running Peak Calling..."))
        if(sum(checkFiles)==0){
          scExp_cf(subset_bam_call_peaks(scExp_cf(), odir, input_files_pc, format = bam_or_bed(),
                                         as.numeric(input$pc_stat_value), annotation_id(),
                                         input$peak_distance_to_merge))
          data = list("scExp_cf" = scExp_cf())
          qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                      paste0(input$selected_reduced_dataset, ".qs")))
          pc$new <- Sys.time()
          updateActionButton(session, "do_pc", label="Finished successfully", icon = icon("check-circle"))
        }
        incProgress(amount = 0.3, detail = paste("Finished Peak Calling..."))
      })
    }
  })
  
  pc <- reactiveValues(new="")
  
  has_available_pc <- reactive({
    if(!is.null(scExp_cf())){
      if("refined_annotation" %in% names(scExp_cf()@metadata) ){
        return(TRUE)
      } else return(FALSE)
    } else{
      return(FALSE)
    }
  })
  
  observeEvent({selected_filtered_dataset()  # reset label on actionButtion when new peak calling should be performed
    input$pc_stat
    input$pc_stat_value}, {
      pc$available_pc <- has_available_pc()
      updateActionButton(session, "do_pc", label="Start", icon = character(0))
  })
  
  
  output$coverage_UI <- renderUI({
    req(GenePool(),has_available_pc())
    print("Inside output$coverage_UI ... ")
      if(has_available_pc()){
        print("Inside output$coverage_UI ... 2 ")
        s <- shinydashboard::box(title="Coverage visualization", width = NULL, status="success", solidHeader = TRUE,
                                 column(3, align="left", selectizeInput(inputId = "select_cov_gene", "Select gene:", choices = NULL, selected = 1)),
                                 column(2, align="left", selectInput("cov_chr","Select chromosome:", choices = unique(rowRanges(scExp_cf())$chr))),
                                 column(2, align="left", textInput("cov_start","Start", value = 15000000)),
                                 column(2, align="left", textInput("cov_end","End", 16000000)),
                                 column(2, align="left", actionButton("make_plot_coverage", "Plot Coverage", icon = icon("chart-area"))))
        print("Inside output$coverage_UI ... 3")
        updateSelectizeInput(session = session,
                             inputId = 'select_cov_gene',
                             choices = GenePool(),
                             server = TRUE)
        return(s)
      }
    })
  
  output$coverage_plot_UI <- renderUI({
    req(GenePool(),has_available_pc(), display_coverage_plot())
    print("Inside output$coverage_UI ... ")
    if(has_available_pc()){
      if(display_coverage_plot()){
      shinydashboard::box(title="Coverage visualization", width = NULL, status="success", solidHeader = TRUE,
                          column(12, align="left", plotOutput("coverage_region_plot") %>%
                                   shinycssloaders::withSpinner(type=8,color="#0F9D58",size = 0.75))) 
      }
    }
  })
  
  
  GenePool <- reactive({
    req(annotation_id())
    print("Inside GenePool reactive ...")
    eval(parse(text = paste0("data(", annotation_id(), ".GeneTSS)")))
    GenePool = unique(eval(parse(text = paste0("", annotation_id(), ".GeneTSS$gene"))))
    print("Inside GenePool reactive2 ...")
    c("Enter Gene...", GenePool)
  })

  
  coverages <- reactive({
    req(has_available_pc(), scExp_cf())
    print("Inside coverage reactive ...")
    display_coverage_plot(TRUE)
    if(has_available_pc()){
      nclust = length(unique(scExp_cf()$cell_cluster))
      odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),
                        "peaks", paste0(selected_filtered_dataset(), "_k", nclust))
      BigWig_filenames <- list.files(odir,pattern = "*.bw", full.names = TRUE)
      print("BigWig_filenames")
      print(BigWig_filenames)
      if(length(BigWig_filenames)>0){
        print("coverages = sapply(BigWig_filenames, rtracklayer::import)")
        coverages = sapply(BigWig_filenames, rtracklayer::import)
        print(length(coverages))
      } else {
        coverages = NULL
      }
      coverages
    }
  })
  
  display_coverage_plot <- reactiveVal(FALSE)
  
  observeEvent(input$make_plot_coverage, { # load reduced data set to work with on next pages
      req(input$cov_chr, has_available_pc(),
          scExp_cf(), coverages(), annotation_id())
    if(input$select_cov_gene != "Enter Gene..."){
      print("Uploading Gene selection")
      
      eval(parse(text = paste0("data(", annotation_id(), ".GeneTSS)")))
      gene_annot = eval(parse(text = paste0("", annotation_id(), ".GeneTSS")))
      print(input$select_cov_gene )
      print(head(gene_annot))
      print(head(gene_annot$gene))
      print(which(gene_annot$gene == input$select_cov_gene))
      updateSelectInput(session, "cov_chr", selected = gene_annot$chr[which(gene_annot$gene == input$select_cov_gene)])
      updateSelectInput(session, "cov_start", selected = gene_annot$start[which(gene_annot$gene == input$select_cov_gene)] - 25000)
      updateSelectInput(session, "cov_end", selected = gene_annot$end[which(gene_annot$gene == input$select_cov_gene)] + 25000)
      updateSelectInput(session, "select_cov_gene", selected = "Enter Gene...")
    }
    label_color_list = setNames(unique(scExp_cf()$cell_cluster_color), unique(scExp_cf()$cell_cluster))
        print("input$make_plot_coverage")
        output$coverage_region_plot <- renderPlot({
          plot_coverage_BigWig(
            coverages(), 
            label_color_list,
            chrom = input$cov_chr,
            start = as.numeric(input$cov_start),
            end =  as.numeric(input$cov_end),
            ref = annotation_id())
        })
        
  })
  
  

 
  # available_pc_plots <- reactive({
  #   fe <- sapply(c(1:input$pc_k_selection), function(i){file.exists(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", input$pc_k_selection), paste0("C", i, "_model.r")))})
  #   which(fe==TRUE)
  # })

  # output$pc_plot_box <- renderUI({
  #   if(has_available_pc()){
  #     shinydashboard::box(title="Peak calling visualization", width = NULL, status="success", solidHeader = TRUE,
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
  
  get.available.DA_GSEA.datasets <- function(name, preproc){
    list.files(path = file.path(init$data_folder, "ChromSCape_analyses", name, "Diff_Analysis_Gene_Sets"),
               full.names = FALSE, recursive = FALSE, pattern = paste0(preproc, "_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+_"))
  }
  
  observeEvent(c(input$qval.th, input$tabs, input$cdiff.th, input$de_type, selected_filtered_dataset()), priority = 10,{
    if(input$tabs == "diff_analysis"){
    init$available_DA_GSA_datasets = get.available.DA_GSEA.datasets(analysis_name(), input$selected_reduced_dataset)
    }
  })
  
  DA_GSA_datasets <- reactive({
  if (is.null(init$available_DA_GSA_datasets)) c() else gsub('.{3}$', '', basename(init$available_DA_GSA_datasets)) })
  
  output$selected_DA_GSA_dataset <- renderUI({ 
    selectInput("selected_DA_GSA_dataset", "Select set with Differential Analysis:",
                choices = DA_GSA_datasets()) 
  })
  
  observeEvent(input$selected_DA_GSA_dataset, { # load reduced data set to work with on next pages
    req(input$selected_DA_GSA_dataset)
  
    file_index <- match(c(input$selected_DA_GSA_dataset), DA_GSA_datasets())
    filename_sel <- file.path(init$data_folder, "ChromSCape_analyses", 
                              analysis_name(),"Diff_Analysis_Gene_Sets",
                              init$available_DA_GSA_datasets[file_index])
    t1 = system.time({
      scExp_cf. = qs::qread(filename_sel)
      if(is.reactive(scExp_cf.)) {
        scExp_cf. = isolate(scExp_cf.())
      }
      scExp_cf(scExp_cf.$scExp_cf) # retrieve filtered scExp
      rm(scExp_cf.)
    })
    cat("Loaded reduced data in ",t1[3]," secs\n")
  })

  output$diff_analysis_info <- renderText({"Differential analysis is performed using
    the cluster assignment obtained with the Cluster Cells module. To change the
    partition of the datasets - number of k clusters - go back to this module and
    select the preferred number of clusters in the box in the upper right panel.
    User can choose either a non-parametric test (Wilcoxon) or a parametric test
    (here edgeR GLM) depending on the observed distribution of the reads."})
  output$selected_k <- renderUI({
    req(scExp_cf())
    HTML(paste0("<h3><b>Number of clusters selected  = ", 
                                         dplyr:::n_distinct(SummarizedExperiment::colData(scExp_cf())$cell_cluster),"</b></h3>"))
  })
  # output$only_contrib_cell_ui <- renderUI({
  #   if("icl" %in% names(scExp_cf()@metadata) && !is.null(scExp_cf()@metadata$icl)){
  #     checkboxInput("only_contrib_cells", "Only use cells contributing most to the clustering ?", value=FALSE)
  #   }
  # })
  # output$contrib_hist <- renderUI({ if(!is.null(input$only_contrib_cells) && input$only_contrib_cells){ plotOutput("contrib_hist_p", height = 250, width = 500) }})
  # output$contrib_hist_p <- renderPlot(contrib_hist_plot())
  # contrib_hist_plot <- reactive({
  #   maxCons <- tapply(scExp_cf()@metadata$icl$itemConsensus$itemConsensus[scExp_cf()@metadata$icl$itemConsensus$k==length(unique(scExp_cf()$cell_cluster))],
  #                     scExp_cf()@metadata$icl$itemConsensus$item[scExp_cf()@metadata$icl$itemConsensus$k==length(unique(scExp_cf()$cell_cluster))], max)
  #   hist(maxCons, col="steelblue", breaks = 80, main="Max cluster contribution per cell", xlab="", ylab="number of cells")
  #   abline(v = input$contrib_thresh, lwd = 2, col="red", lty = 2)
  #   legend("topleft", legend = c("cluster contribution threshold"), col = c("red"), lty = c(2), cex = 0.8)
  # })
  # output$distrib_norm_hist <- renderPlot({
  #   signal = log10(as.numeric(counts(scExp_cf()))+1)
  #   h = hist(signal,
  #   col="steelblue", breaks = 120, main="Distribution of raw signal",
  #        xlab="log10(Raw Signal)", ylab="frequency")
  #   })
  # output$contrib_thresh <- renderUI({ if(!is.null(input$only_contrib_cells) && input$only_contrib_cells){ sliderInput("contrib_thresh", "Select minimum cluster contribution for cells:", min = 0.6, max = 1, value = 0.9, step = 0.01) }})
  # output$contrib_info <- renderUI({ if(!is.null(input$only_contrib_cells) && input$only_contrib_cells){ textOutput("contrib_info_text") }})
  # output$contrib_info_text <- renderText({
  #   total_cells <- length(unique(scExp_cf()@metadata$icl$itemConsensus$item[scExp_cf()@metadata$icl$itemConsensus$k==length(unique(scExp_cf()$cell_cluster))]))
  #   sel_cells <- length(unique(scExp_cf()@metadata$icl$itemConsensus$item[scExp_cf()@metadata$icl$itemConsensus$k==length(unique(scExp_cf()$cell_cluster)) & scExp_cf()@metadata$icl$itemConsensus$itemConsensus >= input$contrib_thresh]))
  #   paste("Selected top", sel_cells, "cells out of", total_cells)
  # })
  # 
  observeEvent(c(input$qval.th, input$tabs, input$cdiff.th, input$de_type, selected_filtered_dataset()), priority = 10,{
    if(input$tabs == "diff_analysis"){
      if(!is.null(selected_filtered_dataset()) && !is.null(input$qval.th) && !is.null(input$cdiff.th)){
        DA_GSA_suffix = input$de_type
        if(input$de_type == "custom") DA_GSA_suffix = paste0(gsub("[^[:alnum:]|_]","",input$name_group),"_vs_",
                                                             gsub("[^[:alnum:]|_]","",input$name_ref))
        filename <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "Diff_Analysis_Gene_Sets",
                              paste0(selected_filtered_dataset(), "_", input$nclust,
                                     "_", input$qval.th, "_", input$cdiff.th, "_", DA_GSA_suffix, ".qs"))
        
        if(file.exists(filename)){
          data = qs::qread(filename)
          scExp_cf(data$scExp_cf)
          rm(data)
          gc()
        } else {
          NULL
        }
      }
    }
    })
  output$custom_da_ref <- renderUI({
    req(input$de_type)
    if(input$de_type == "custom"){
      shiny::radioButtons("ref_type","Reference type", choices = c("sample_id","cell_cluster"))
    }
  })
  output$custom_da_group <- renderUI({
    req(input$de_type)
    if(input$de_type == "custom"){
      shiny::radioButtons("group_type","Group type", choices= c("sample_id","cell_cluster"))
    }
  })
  
  output$ref_choice <- renderUI({
    req(scExp_cf(), input$ref_type, input$de_type)
    if(input$de_type == "custom"){
      if(input$ref_type == "sample_id"){
        choices = unique(scExp_cf()$sample_id)
        selectInput("ref_choice","Reference sample(s)", choices= choices, 
                    multiple = TRUE)
      } else if(input$ref_type == "cell_cluster"){
        choices = unique(scExp_cf()$cell_cluster)
        selectInput("ref_choice","Reference cluster(s)", choices=choices,
                    multiple = TRUE)
      } else{NULL}
    }
  })
  
  output$group_choice <- renderUI({
    req(scExp_cf(), input$group_type,input$de_type)
    if(input$de_type == "custom"){
      if(input$group_type == "sample_id"){
        choices = unique(scExp_cf()$sample_id)
        selectInput("group_choice","Group sample(s)", choices= choices, 
                    multiple = TRUE)
      } else if(input$group_type == "cell_cluster"){
        choices = unique(scExp_cf()$cell_cluster)
        selectInput("group_choice","Group cluster(s)", choices=choices,
                    multiple = TRUE)
      } else{ NULL}
    }
  })
  
  output$text_vs <- renderUI({
    req(input$de_type)
    if(input$de_type == "custom") h3(" VS ")
  })
  
  output$name_group <- renderUI({
    req(input$de_type)
    if(input$de_type == "custom") textInput("name_group", label = "Name Group:", value = "group")
  })
  
  output$name_ref <- renderUI({
    req(input$de_type)
    if(input$de_type == "custom") textInput("name_ref", label = "Name Reference:", value = "ref")
  })
  
  observeEvent(input$run_DA, {  # perform differential analysis based on wilcoxon test
  print("RUN DA")
  print(input$group_choice)
  print(input$ref_choice)
    if(input$de_type == "custom"){
      
        if(TRUE %in% all.equal(input$group_choice, input$ref_choice)){
          showNotification(paste0("Warning : Please select different group and references."),
                           duration = 7, closeButton = TRUE, type="warning")
          return(0)
        }
      
  }
    withProgress(message='Performing differential analysis...', value = 0, {
      incProgress(amount = 0.2, detail = paste("Initializing DA"))
      if(batchUsed()) block = TRUE else block = FALSE
      gc()
      group = ""
      ref = ""
      if(input$de_type == "custom"){
        group = data.frame(input$group_choice)
        ref = data.frame(input$ref_choice)
        colnames(group) = gsub("[^[:alnum:]|_]","",input$name_group)
        colnames(ref) = gsub("[^[:alnum:]|_]","",input$name_ref)
      }

      scExp_cf(differential_analysis_scExp(scExp = scExp_cf(),
                                           method= input$da_method,
                                           de_type = input$de_type,
                                           cdiff.th = input$cdiff.th,
                                           qval.th = input$qval.th,
                                           block = block,
                                           group = group,
                                           ref = ref)) 
      gc()
      incProgress(amount = 0.6, detail = paste("Finishing DA..."))
      data = list("scExp_cf" = scExp_cf())
      DA_GSA_suffix = input$de_type
      if(input$de_type == "custom") DA_GSA_suffix = paste0(gsub("[^[:alnum:]|_]","",input$name_group),"_vs_",
                                                           gsub("[^[:alnum:]|_]","",input$name_ref))
      
      suffix = paste0(selected_filtered_dataset(), "_", length(unique(scExp_cf()$cell_cluster)),
             "_", input$qval.th, "_", input$cdiff.th, "_", DA_GSA_suffix, ".qs")
      
      qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "Diff_Analysis_Gene_Sets",
                                       suffix))
      rm(data)
      gc()
      init$available_DA_GSA_datasets = get.available.DA_GSEA.datasets(analysis_name(), input$selected_reduced_dataset)

      updateSelectInput(session = session, inputId = "selected_DA_GSA_dataset",
                        label =  "Select set with Differential Analysis:",
                        choices = DA_GSA_datasets(),
                        selected =  gsub(".qs$","",suffix))
      incProgress(amount = 0.2, detail = paste("Saving DA"))
    })
  })
  
  output$da_summary_box <- renderUI({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        shinydashboard::box(title="Number of differentially enriched regions", width = NULL, status="success", solidHeader = TRUE,
            column(5, align="left", br(), tableOutput("da_summary_kable")),
            column(7, align="left", plotOutput("da_barplot", height = 270, width = 250)),
            column(4, align="left", downloadButton("download_da_barplot", "Download barplot"))
        )
      }
    }
  })
  
  output$da_summary_kable <- function(){
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        scExp_cf()@metadata$diff$summary %>%
          kableExtra::kable(escape = FALSE, align="c") %>%
          kableExtra::kable_styling(c("striped", "condensed"), full_width = FALSE)
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
    req(input$gpsamp, scExp_cf())
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        if(input$gpsamp %in% scExp_cf()@metadata$diff$groups){
        diff = scExp_cf()@metadata$diff$res[,-grep("Rank|pval|ID",colnames(scExp_cf()@metadata$diff$res))]
        if(!is.null(SummarizedExperiment::rowRanges(scExp_cf())$Gene)){
          diff = cbind(SummarizedExperiment::rowRanges(scExp_cf())$distanceToTSS, diff)
          diff = cbind(SummarizedExperiment::rowRanges(scExp_cf())$Gene, diff)
          colnames(diff)[1:2] = c("Gene","distanceToTSS")
        }
        rownames(diff) <- NULL
        # for(i in seq(from = 0, to=(dim(diff)[2]-8)/5, by = 1)){
        #   diff[, (5*i+5)] <- round(diff[, (5*i+5)], 3) #counts
        #   diff[, (5*i+6)] <- round(diff[, (5*i+6)], 3) #cdiff
        # }
        
        diff = diff[,c("Gene","distanceToTSS", "chr","start","end",
                       colnames(diff)[grep(input$gpsamp,colnames(diff))] )]
        diff <- diff[order(diff[,paste0("cdiff.",input$gpsamp)]),]
        DT::datatable(diff, options = list(dom='tpi'), class = "display",
                      rownames = FALSE, autoHideNavigation = TRUE,
                      )
      }
      }
    }
  })
  
  output$download_da_table <- downloadHandler(
    filename = function(){ paste0("diffAnalysis_data_", selected_filtered_dataset(),
                                  "_", length(unique(scExp_cf()$cell_cluster)), "_",
                                  input$qval.th, "_", input$cdiff.th, "_",
                                  input$de_type, ".csv")},
    content = function(file){
      write.table(scExp_cf()@metadata$diff$res, file, row.names = FALSE, quote = FALSE, sep=",")
    })
  
  output$da_visu_box <- renderUI({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
        shinydashboard::box(title="Detailed differential analysis per cluster", width = NULL, status="success", solidHeader = TRUE,
            column(4, align="left", selectInput("gpsamp", "Select group:", choices = scExp_cf()@metadata$diff$groups)),
            column(4, align="left", downloadButton("download_da_table", "Download table")),
            column(12, align="left", div(style = 'overflow-x: scroll', DT::dataTableOutput('da_table')), br()),
            column(12, align="left", plotOutput("h1_prop", height = 300, width = 500),
                   plotOutput("da_volcano", height = 500, width = 500)),
            column(4, align="left", downloadButton("download_h1_plot", "Download histogram")),
            column(4, align="left", downloadButton("download_da_volcano", "Download volcano plot")))
      }
    }
  })
  
  output$h1_prop <- renderPlot({
    req(scExp_cf())
    if(!is.null(scExp_cf())){
      if(input$gpsamp %in% scExp_cf()@metadata$diff$groups){
        if(!is.null(scExp_cf()@metadata$diff)){
          plot_differential_H1_scExp(scExp_cf(), input$gpsamp)
        }  
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
    req(scExp_cf())
    if(!is.null(scExp_cf())){
      if(input$gpsamp %in% scExp_cf()@metadata$diff$groups){
        if(!is.null(scExp_cf()@metadata$diff)){
          plot_differential_volcano_scExp(scExp_cf(),cell_cluster = input$gpsamp,
                                          cdiff.th = input$cdiff.th, qval.th = input$qval.th)
        }
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
      able_disable_tab(c("diff_my_res"),"enrich_analysis")
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
    c("c1_positional","c2_curated","c3_motif","c4_computational",
      "c5_GO","c6_oncogenic","c7_immunologic","hallmark")
  })

  annotFeat_long <- reactive({
    if(input$tabs == "enrich_analysis"){
      af = as.data.frame(SummarizedExperiment::rowData(scExp_cf()))
      af = tidyr::separate_rows(af, Gene,sep = ", ")
      af
    }
  })
  
  url <- a("MSigDB homepage", href="https://www.gsea-msigdb.org/gsea/msigdb/index.jsp",
           target='"_blank">gsea-msigdb.org/gsea/msigdb/index.jsp')
  output$enr_info <- renderUI({tagList("Enrichment will be performed based on the
                                         significant genes per cluster that were computed on the previous page. 
                                         Genes in vincinity of differential features are tested using hypergeometric test
                                         against MSigDB pathway lists (", url,").")})
  
  canUsePeaks <- reactive({
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
    if(input$tabs == "enrich_analysis"){
      myData = new.env()
      eval(parse(text = paste0("data(",annotation_id(),".GeneTSS, envir = myData)")))
      as.character(unique(
        eval(parse(text = paste0("myData$",annotation_id(),".GeneTSS$gene")))
      ))
    }
  })
  
  observeEvent(input$do_enrich, {
    withProgress(message='Running Pathway Enrichment Analysis...', value = 0, {

      incProgress(amount = 0.3, detail = paste("Running Hypergeometric Enrichment Testing against MSigDB..."))
      gc()
      scExp_cf(gene_set_enrichment_analysis_scExp(scExp_cf(), enrichment_qval = 0.01, qval.th = input$qval.th,
                                                  ref = annotation_id(), cdiff.th = input$cdiff.th,
                                                  peak_distance = 1000, use_peaks = input$use_peaks,
                                                  GeneSetClasses = MSIG.classes()))
      gc()
      incProgress(amount = 0.6, detail = paste("Finishing Pathway Enrichment Analysis..."))
      data = list("scExp_cf" = scExp_cf() )
      DA_GSA_suffix = input$de_type
      if(input$de_type == "custom") DA_GSA_suffix = paste0(gsub("[^[:alnum:]|_]","",input$name_group),"_vs_",
                                                           gsub("[^[:alnum:]|_]","",input$name_ref))
      qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "Diff_Analysis_Gene_Sets",
                                  paste0(selected_filtered_dataset(), "_", length(unique(scExp_cf()$cell_cluster)),
                                         "_", input$qval.th, "_", input$cdiff.th, "_", DA_GSA_suffix, ".qs")))
      rm(data)
      gc()
      incProgress(amount = 0.6, detail = paste("Saving Pathway Enrichment Analysis"))
      
    })
  })
  
  output$GSA_group_sel <- renderUI({ 
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$diff)){
      selectInput("GSA_group_sel", "Select group:", choices = scExp_cf()@metadata$diff$groups) 
      }
    }
    })
  
  output$enr_class_sel <- renderUI({
    req(MSIG.classes())
    shiny::checkboxGroupInput(
    inputId = "enr_class_sel", inline = TRUE,
    label =  "Select classes to display:",
    selected = MSIG.classes(), choiceNames = MSIG.classes(),
    choiceValues = MSIG.classes())
  })
  
  output$all_enrich_table <- DT::renderDataTable({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$enr) && !is.null(input$GSA_group_sel) && !is.null(input$enr_class_sel)){
        if(input$GSA_group_sel %in% scExp_cf()@metadata$diff$groups)
          table_enriched_genes_scExp(scExp_cf(),set = "Both", group = input$GSA_group_sel, input$enr_class_sel)
      }
    }
  })
  
  output$over_enrich_table <- DT::renderDataTable({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$enr)){
        if(input$GSA_group_sel %in% scExp_cf()@metadata$diff$groups)
          table_enriched_genes_scExp(scExp_cf(), set = "Overexpressed", group = input$GSA_group_sel, input$enr_class_sel)
      }
    }
  })
  
  output$under_enrich_table <- DT::renderDataTable({
    if(!is.null(scExp_cf())){
      if(!is.null(scExp_cf()@metadata$enr)){
        if(input$GSA_group_sel %in% scExp_cf()@metadata$diff$groups) 
          table_enriched_genes_scExp(scExp_cf(), set = "Underexpressed", group = input$GSA_group_sel, input$enr_class_sel)
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
    req(annotFeat_long())
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
    req(input$gene_sel, annotFeat_long())
    subset <- annotFeat_long()[which(annotFeat_long()$Gene==input$gene_sel), ]
    if(!is.null(subset)){
      subset <- subset[order(subset$distanceToTSS),]
      regions <- paste0(subset$ID, " (distanceToTSS to gene TSS: ", subset$distanceToTSS, ")")
      selectInput("region_sel", "Select associated genomic region:", choices = regions)
    }
  })
  

  gene_umap_p <- reactive({
    req(input$gene_sel, input$region_sel)
    region <- strsplit(input$region_sel, " ")[[1]][1]
    if(region %in% rownames(scExp_cf())){
      p <- ggplot(as.data.frame(SingleCellExperiment::reducedDim(scExp_cf(), "UMAP")),
                  aes(x = Component_1, y = Component_2)) +
        geom_point(alpha = 0.5, aes(color = SingleCellExperiment::normcounts(scExp_cf())[region, ],
                                    shape = SummarizedExperiment::colData(scExp_cf())$cell_cluster)) +
        labs(color="norm. count for region", title = "UMAP", shape="Cluster", x="Component 1", y="Component 2") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour="black"),
              panel.border = element_rect(colour="black", fill = NA)) +
        viridis::scale_color_viridis(direction=-1)
    }
  })
  
  output$gene_umap_UI <- renderUI({
    req(input$gene_sel, input$region_sel)
    region <- strsplit(input$region_sel, " ")[[1]][1]
    if(region %in% rownames(scExp_cf())){
      output$gene_umap_plot <- plotly::renderPlotly({
          plotly::ggplotly(gene_umap_p(), tooltip="Sample", dynamicTicks = TRUE)
      })
      plotly::plotlyOutput("gene_umap_plot")
    }
  })
  
  
  ###############################################################
  # 8. Close app
  ###############################################################
  
  output$analysis_deletion_info <- renderText({"The selected data set will be fully deleted from the computer, including all reduced data versions that have been produced so far for this set."})
  output$selected_delete_analysis <- renderUI({ selectInput("selected_delete_analysis", "Select data set :", choices = init$available_analyses) })
  
  observeEvent(input$delete_analysis, {  # delete selected dataset
    withProgress(message='Deleting data set', value = 0, {
      incProgress(amount=0.5, detail=paste("..."))
      unlink(file.path(init$data_folder, "ChromSCape_analyses", input$selected_delete_analysis), recursive=TRUE)
      init$available_analyses <- list.dirs(path=file.path(init$data_folder, "ChromSCape_analyses"), full.names=FALSE, recursive=FALSE)
      init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())
      print("observeEvent({ input$delete_analysis}")
      print(init$available_reduced_datasets)
      incProgress(amount=0.5, detail=paste("... finished"))
    })
    showNotification("Data set successfully deleted.", duration=5, closeButton=TRUE, type="warning")
  })
  
  observeEvent(input$close_and_save, {
    unlink(file.path("www", "images", "*"))
    unlink(file.path(".", "*.csv"))
    
    if(!is.null(scExp()) & !is.null(input$selected_reduced_dataset)){
      scExp = isolate(scExp())
      qs::qsave(scExp, file = file.path(init$data_folder, "ChromSCape_analyses",
                                   analysis_name(), "Filtering_Normalize_Reduce",
                                   paste0(input$selected_reduced_dataset,".qs")))
      
    }
    if(!is.null(scExp_cf()) & !is.null(selected_filtered_dataset())){
      data = list()
      data$scExp_cf = isolate(scExp_cf())
      
      if("consclust" %in% names(scExp_cf()@metadata)){
        if("diff" %in% names(scExp_cf()@metadata)){
          dir =file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), 
                         "Diff_Analysis_Gene_Sets", paste0(
                           selected_filtered_dataset(), "_",
                           length(unique(scExp_cf()$cell_cluster)),
                           "_", input$qval.th, "_", input$cdiff.th, "_",
                           input$de_type, ".qs"))
        } else{
          dir = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),
                          "correlation_clustering",
                         paste0(selected_filtered_dataset(), 
                                 ".qs"))
        }
        qs::qsave(scExp, file = dir)
      } 
      
    }
    
    lapply(names(resourcePaths()), removeResourcePath)
    
    print("Thank you & see you next time !")
    stopApp()
  })
  
})
