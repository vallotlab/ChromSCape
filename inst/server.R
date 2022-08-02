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
    shinyhelper::observe_helpers(help_dir = "www/helpfiles", withMathJax = TRUE)
    
    tab_vector = c("filter_normalize",
                   "vizualize_dim_red",
                   "cons_clustering",
                   "coverage",
                   "diff_analysis",
                   "enrich_analysis",
                   "TF_analysis") #list of all lockable tabs
    unlocked = reactiveValues(list = list(selected_analysis = FALSE,
                                          selected_reduced_dataset = FALSE,
                                          pca = FALSE,
                                          affectation = FALSE,
                                          diff_my_res = FALSE)) #list of all required items to unlock a tab
    for(tab in tab_vector){
        js$disableTab(tab) #Disabling all tabs but the first one
    }
    
    observeEvent(input$startHelp,{
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
    
    CS_options.BPPARAM = reactive({
        req(input$options.bpparam_class, input$options.nb_workers)
        BPPARAM = BiocParallel::bpparam(input$options.bpparam_class)
        BiocParallel::bpworkers(BPPARAM) = input$options.nb_workers
        BPPARAM
    })
    
    reduced_datasets <- reactive({
        if (is.null(init$available_reduced_datasets)) c() else gsub('.{3}$', '', basename(init$available_reduced_datasets)) })
    
    observeEvent({analysis_name()},{
        init$available_reduced_datasets = get.available.reduced.datasets(analysis_name())
    })
    annotCol <- reactive({
        req(scExp())
        counts_cols = paste0("counts_", getExperimentNames(scExp()))
        if("batch_name" %in% colnames(SummarizedExperiment::colData(scExp())))
        {
            c("sample_id",counts_cols,"batch_name")
        } else{
            c("sample_id",counts_cols)
        }
    })
    output$feature_color <- renderUI({selectInput("color_by", "Color by", choices=annotCol())})
    
    observeEvent({
        analysis_name()
        input$feature_select
    }, { # application header (tells you which data set is selected)
        req(analysis_name())
        header <- paste0('<b>Analysis : ', analysis_name(), ' - ',input$feature_select,' </b>')
        shinyjs::html("pageHeader", header) 
    })
    
    
    get.available.reduced.datasets <- function(selected_analysis){
        list.files(path = file.path(init$data_folder, "ChromSCape_analyses", selected_analysis,"Filtering_Normalize_Reduce"), full.names = FALSE, recursive = TRUE,
                   pattern="[[:print:]]+_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+(.[[:digit:]]+)?_(uncorrected|batchCorrected).qs")
        
    }
    get.available.alternative.datasets <- function(selected_analysis){
        gsub(".qs$","",gsub("^datamatrix_","",list.files(path = file.path(init$data_folder, "ChromSCape_analyses", selected_analysis),
                                                         full.names = FALSE, recursive = FALSE,
                                                         pattern="datamatrix_.*.qs")
        ))
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
            if(.Platform$OS.type != "windows"){
                js$save_cookie(init$data_folder)
            }
        }
    )
    
    
    output$input_data_ui <- renderUI({
        if(input$data_choice_box == "DenseMatrix"){
            column(12, br(),fileInput("datafile_matrix", "Upload all data matrices (.tsv / .gz) :",
                                      multiple=TRUE, accept=c("text", "text/plain", ".txt", ".tsv", ".csv", ".gz")),
                   checkboxInput("is_combined_mat", "Single Multi-sample dense matrix ?",value = FALSE),
                   uiOutput("nb_samples_mat")
            )
            
        }
        else if (input$data_choice_box == "SparseMatrix"){
            column(12,
                   br(),
                   HTML(paste0("<b>Upload folder (", input$data_choice_box,")</b><br>")),
                   shinyFiles::shinyDirButton(id = "datafile_folder", label = "Select folder containing raw files",
                                              title =  paste0("Select a directory containing your ",
                                                              input$data_choice_box," files."),
                                              icon = icon("folder-open")),
                   checkboxInput("is_combined_mat", "Single Multi-sample dense matrix ?",value = FALSE),
                   uiOutput("nb_samples_mat")
            )
        }else if(input$data_choice_box == "FragmentFile" | input$data_choice_box == "scBAM" | input$data_choice_box == "scBED"){
            column(12,
                   br(),
                   HTML(paste0("<b>Upload folder (", input$data_choice_box,")</b><br>")),
                   shinyFiles::shinyDirButton(id = "datafile_folder", label = "Select folder containing raw files",
                                              title =  paste0("Select a directory containing your ",
                                                              input$data_choice_box," files."),
                                              icon = icon("folder-open"))
                   
            )
        }
        
    })
    
    input_sc <- function(){
        modalDialog(
            title = "How to select my samples:",
            HTML("
      <p> Before starting uploading your data, make sure you have placed all your 
      datasets in a given <u>directory</u>. Each sample input data, no matter it's format,
      should be place in its own <b>folder</b> in this <u>directory</u>. Labels will be given to
      each sample based on the name of the <b>folders</b>. </p>
      <br>
      <p> Here, simply pick the <u>directory</u> containing the samples folders, and
      ChromSCape will automatically detect your <b>samples</b>. </p>
      <br>
      Example of folder structure (for scBED, but applies in a similar manner for SparseMatrix, Fragment File & scBAM) : <br>
      <u>scChIPseq</u><br>
          * <b>Sample_1</b><br>
                 - cell_1.bed<br>
                 - cell_2.bed<br>
                 - cell_3.bed<br>
                 - cell_4.bed<br>
          * <b>Sample_2</b><br>
                 - cell_1.bed<br>
                 - cell_2.bed<br>
                 - cell_3.bed<br>
                 - cell_4.bed<br>

      <br>
      You're almost there !
      <br>
      ")
            ,
            easyClose = TRUE
        )
    }
    
    showHelpSC = reactiveVal(TRUE)
    observeEvent(input$datafile_folder, priority = 10000, {
        if(showHelpSC()) {
            showModal(input_sc())
            showHelpSC(FALSE)
        } else { showHelpSC(TRUE)}
    })
    
    
    output$rawcount_data_input <- renderUI({
        if(input$data_choice_box != "DenseMatrix"){
            column(12,
                   shinydashboard::box(title="Counting parameters", width = NULL, status="success", solidHeader = TRUE,
                                       column(6, 
                                              radioButtons("count_on_box", label = "Select a count method",
                                                           choices = list("Count on bins (width)"="bin_width",
                                                                          "Count on peaks (must provide a .bed file)" = "peak_file",
                                                                          "Count on genes (body + promoter)" = "genebody"))
                                       ),
                                       column(6,
                                              uiOutput("bin_width"),
                                              uiOutput("peak_file"),
                                              uiOutput("extendPromoter"))
                   )
            )
        }
    })
    
    output$advanced_data_input <- renderUI({
        if(input$data_choice_box != "DenseMatrix"){
            if(input$data_choice_box == "SparseMatrix"){
                column(12, uiOutput("datafile_folder_upload_UI"))
            } else{
                column(12, uiOutput("rawcount_data_input"),
                       uiOutput("datafile_folder_upload_UI"))
            }
        }
    })
    
    
    output$nb_samples_mat <- renderUI({ if(input$is_combined_mat == TRUE){
        selectInput(inputId = "nb_samples_to_find",label = "Number of samples:",
                    choices = 1:100, selected = 1, multiple = FALSE)
    }})
    
    output$bin_width <- renderUI({ if(input$count_on_box == "bin_width"){
        textInput("bin_width", label = "Width of bins to count on (in bp) :",value = 50000)
    }})
    output$peak_file <- renderUI({ if(input$count_on_box == "peak_file"){
        fileInput("peak_file", ".bed file containing the peaks to count on:", multiple = FALSE, accept = c(".bed",".txt"))
    }})
    output$extendPromoter <- renderUI({ if(input$count_on_box == "genebody" ){
        textInput("extendPromoter", label = "Extend promoter by", value = 2500)
    }})
    
    shinyFiles::shinyDirChoose(input, "datafile_folder", roots = volumes, session = 
                                   session)
    
    files_dir <- reactiveVal(NULL)
    
    observeEvent(input$datafile_folder,
                 {
                     req(input$datafile_folder)
                     if(!is.null(input$datafile_folder)){
                         datafile_folder = shinyFiles::parseDirPath(volumes, input$datafile_folder)
                         files_dir_tmp = list.dirs(datafile_folder, recursive = FALSE, full.names = TRUE)
                         names(files_dir_tmp) = basename(files_dir_tmp)
                         files_dir(files_dir_tmp)
                     }
                 })
    
    output$datafile_folder_upload_UI <- renderUI({
        req(files_dir())
        if(!is.null(files_dir()) & length(files_dir()) > 0 ){
            selectInput("sample_selection", label = "Selected samples",
                        choices = names(files_dir()), multiple = TRUE,
                        selected = names(files_dir()) )
        }
    })
    
    output$rebin_matrices_checkbox_ui <- renderUI({
        req(input$data_choice_box)
        if(input$data_choice_box %in% c("DenseMatrix", "SparseMatrix")) 
            checkboxInput("rebin_matrices", label = "Rebin matrices", value = FALSE)
    })
    
    output$rebin_matrices_ui <- renderUI({
        req(input$rebin_matrices)
        if(isTRUE(input$rebin_matrices)){
            column(6, textInput("rebin_bin_size", label = "Re-count on genomic bins (bp)", value = 50000),
                   fileInput("rebin_custom_annotation", label = "Re-count on features (BED)"),  multiple = FALSE, accept = c(".bed",".txt", ".bed.gz"),
                   textInput("minoverlap", label = "Minimum Overlap (bp)", value = 500))
        }
    })
    
    output$add_to_current_analysis_checkbox_UI <- renderUI({
        req(input$feature_select)
        column(12, 
               shinyWidgets::materialSwitch(inputId = "add_to_current_analysis", value = FALSE,
                                            label = "Add to current analysis",width = "80%", inline = F, status = "primary") %>%
                   shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                       content = "add_to_current_analysis", size = "l")
        )
    })
    
    add_to_current_analysis_Modal <- function(failed = FALSE){
        modalDialog(title = "Additional feature layer on the same cells",
                    size = "l",
                    HTML(paste0('You are about to add another layer of features in the analysis <b>"', analysis_name(),
                                '"</b> that contains <b>',ncol(init$datamatrix), '</b> cells (barcodes: <b>',
                                paste0(head(init$annot_raw$barcode,2),collapse = ", "),
                                '...)</b> and <b>',length(unique(init$annot_raw$sample_id)),'</b> samples (samples: <b>',
                                paste0(head(init$annot_raw$sample_id,2),collapse = ", "), '</b>...).'),
                         '<br> The data you have selected <b> must have the same cell barcodes an sample names</b> as in the current',
                         ' analysis. '),
                    br(), br(),
                    span('Choose a name corresponding to the feature you are counting on, ',
                         'e.g. "bins_50k", "peaks", "genebodies"...'),
                    br(),br(),
                    textInput("alt_name", "Name of the features",
                              placeholder = 'Try "bins_50k" or "peaks"'
                    ),
                    
                    if (failed)
                        div(tags$b("Name must not contains non-alphanumerical characters except '_', be different than 'main' and must not be empty and be shorter than 25 characters.", style = "color: red;")),
                    
                    footer = tagList(
                        actionButton("add_to_current_analysis_Modal_cancel", "Cancel"),
                        actionButton("add_to_current_analysis_Modal_ok", "OK")
                    )
        )
    }
    
    alt_name <- reactiveVal("")
    
    observeEvent(input$add_to_current_analysis_Modal_ok, {
        # Check that data object exists and is data frame.
        if (!is.null(input$alt_name) && nchar(input$alt_name) > 0 && nchar(input$alt_name) < 15
            && !grepl("[[:punct:]]", gsub("_","",input$alt_name)) && input$alt_name!="main") {
            alt_name(input$alt_name)
            updateActionButton(session, "create_analysis", label= "Create Analysis", icon = character(0))
            removeModal()
        } else {
            showModal(add_to_current_analysis_Modal(failed = TRUE))
        }
    })
    
    observeEvent(input$add_to_current_analysis_Modal_cancel, {
        # Check that data object exists and is data frame.
        alt_name("")
        shinyWidgets::updateMaterialSwitch(session = session, inputId = "add_to_current_analysis", value = FALSE)
        removeModal()
    })
    
    observeEvent(input$add_to_current_analysis, {
        if(!is.null(input$selected_analysis) & nchar(input$selected_analysis) > 1){
            if(input$add_to_current_analysis){
                showModal(add_to_current_analysis_Modal())
            }
        } else{
            if(input$add_to_current_analysis) 
                showNotification(paste0("Warning : Please create an analysis before",
                                        " adding new features..."),
                                 duration = 10, closeButton = TRUE, type="warning")
            shinyWidgets::updateMaterialSwitch(session = session, inputId = "add_to_current_analysis", value = FALSE)
        }
    })
    
    observeEvent(input$create_analysis, {  # save new dataset
        req(input$new_analysis_name, input$annotation)
        raw_mat = NULL
        
        if(is.null(input$datafile_folder) && is.null(input$datafile_matrix)) return()
        if(!is.null(input$datafile_folder) && is.null(input$datafile_matrix)){
            if(is.null(files_dir())) return()
        }
        
        datamatrix <- NULL
        annot_raw <- NULL
        type_file = as.character(input$data_choice_box)
        if(!input$add_to_current_analysis && 
           dir.exists(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name))){
            showNotification(paste0("Warning : The name : '",input$new_analysis_name,
                                    " is already taken by a preexisting analysis. Please
                              choose another name for your analysis."),
                             duration = 5, closeButton = TRUE, type="warning")
        } else{
            
            progress <- shiny::Progress$new(session, min=0, max=1)
            on.exit(progress$close())
            if(input$add_to_current_analysis) 
                progress$set(message=paste0('Adding new feature : ', alt_name(), "..."), value = 0.1) else
                    progress$set(message='Creating new data set..', value = 0.1)
            
            if(type_file == "DenseMatrix" & !is.null(input$datafile_matrix)){
                
                progress$inc(detail="Reading Dense Matrices", amount = 0.3)
                if(input$is_combined_mat == TRUE){
                    if(length(input$datafile_matrix$name)>1){
                        showNotification(paste0("Warning : When checking the 
                                      'The matrix contains multiple samples ?' button,
                                      you have to input a single Dense Matrix."),
                                         duration = 5, closeButton = TRUE, type="warning")
                        return()
                    }
                }
                tmp_list = import_scExp(file_paths = input$datafile_matrix$name,
                                        temp_path = input$datafile_matrix$datapath)
                datamatrix = tmp_list$datamatrix
                if(input$is_combined_mat == TRUE) {
                    samples_ids = detect_samples(colnames(datamatrix),
                                                 nb_samples = as.numeric(input$nb_samples_to_find))
                    annot_raw = data.frame(barcode = colnames(datamatrix),
                                           cell_id = colnames(datamatrix),
                                           sample_id = samples_ids,
                                           batch_id = factor(rep(1, ncol(datamatrix))))
                } else{ annot_raw = tmp_list$annot_raw }
            } else {
                selected_sample_folders = files_dir()[input$sample_selection]
                send_warning = FALSE
                if(type_file == "scBAM") if(length(list.files(selected_sample_folders[1],pattern = "*.bam$"))==0) send_warning = TRUE
                if(type_file == "scBED") if(length(list.files(selected_sample_folders[1],pattern = "*.bed$|.*.bed.gz"))==0) send_warning = TRUE
                if(type_file == "FragmentFile") if(length(list.files(selected_sample_folders[1],pattern = "*.tsv|.*.tsv.gz"))==0) send_warning = TRUE
                if(type_file == "FragmentFile") if(!requireNamespace("Signac", quietly=TRUE)){
                    showNotification(paste0("Warning : In order to read in Fragment Files, you must install Signac Package first.",
                                            "Run install.packages('Signac') in console. "), duration = 20, closeButton = TRUE, type="error")
                }
                if(type_file == "SparseMatrix") {
                    combin = expand.grid(c(".*features", ".*barcodes", ".*matrix"), c(".mtx",".tsv",".txt",".bed",".*.gz"))[-c(1,2,6,9,11,12),]
                    pattern = paste(combin$Var1,combin$Var2, sep="", collapse = "|")
                    if(length(list.files(selected_sample_folders[1],
                                         pattern = pattern))!=3) send_warning = TRUE
                }
                
                if(send_warning) {
                    showNotification(paste0("Warning : Can't find any specified file types in the upload folder. 
                                    Select another upload folder or another data type."),
                                     duration = 5, closeButton = TRUE, type="warning")
                    return()
                }
                progress$inc(detail=paste0("Reading ",type_file,
                                           " files to create matrix. This might take a while."), amount = 0.1)
                if(type_file == "SparseMatrix"){
                    out =  ChromSCape:::read_sparse_matrix(
                        files_dir_list = selected_sample_folders,
                        ref = input$annotation)
                } else if(type_file %in% c("FragmentFile","scBAM","scBED") & !is.null(input$datafile_folder)) {
                    
                    if(input$count_on_box == "bin_width") out = ChromSCape:::raw_counts_to_sparse_matrix(
                        files_dir_list = selected_sample_folders,
                        file_type = type_file,
                        bin_width = as.numeric(input$bin_width),
                        ref = input$annotation,
                        progress = progress,
                        BPPARAM = CS_options.BPPARAM())
                    
                    if(input$count_on_box == "peak_file") out = ChromSCape:::raw_counts_to_sparse_matrix(
                        files_dir_list = selected_sample_folders,
                        file_type = type_file,
                        peak_file = as.character(input$peak_file$datapath),
                        ref = input$annotation,
                        progress = progress,
                        BPPARAM = CS_options.BPPARAM())
                    
                    if(input$count_on_box == "genebody")  out = ChromSCape:::raw_counts_to_sparse_matrix(
                        files_dir_list = selected_sample_folders,
                        file_type = type_file,
                        genebody = TRUE,
                        extendPromoter = as.numeric(input$extendPromoter),
                        ref = input$annotation,
                        progress = progress,
                        BPPARAM = CS_options.BPPARAM())
                } else {
                    stop("No data folder or data files selected.")
                }
                
                datamatrix = out$datamatrix
                annot_raw = out$annot_raw
                
                scExp. = create_scExp(datamatrix[1:min(nrow(datamatrix),100),1:10], annot_raw[1:10,], FALSE, FALSE, FALSE, FALSE,verbose = FALSE)
                original_bin_size = mean(GenomicRanges::width(get_genomic_coordinates(scExp.)))
                rm(scExp.)
                
                if(original_bin_size < 300 & input$rebin_matrices == TRUE){
                    print("Saving raw matrix as average bin size is lesser than 300bp, for later use (coverage)...")
                    raw_mat = datamatrix
                }
                if(input$rebin_matrices == TRUE){
                    
                    if(!is.na(as.numeric(input$rebin_bin_size)) & !is.na(as.numeric(input$minoverlap))){
                        if(as.numeric(input$minoverlap) > as.numeric(input$rebin_bin_size) | as.numeric(input$minoverlap) > original_bin_size) {
                            showNotification(paste0("Warning : To rebin the matrices, 
                                      'The minimum overlap of original matrices should be smaller than both
                                      the original and new bin sizes."),
                                             duration = 10, closeButton = TRUE, type="warning")
                            return()
                        } 
                        if(as.numeric(input$rebin_bin_size) < original_bin_size){
                            showNotification(paste0("Warning : To rebin the matrices, 
                                      'The new bin size must be larger than the
                                      the original bin size."),
                                             duration = 10, closeButton = TRUE, type="warning")
                            return()
                        }
                        if(!is.null(input$rebin_custom_annotation) && file.exists(as.character(input$rebin_custom_annotation$datapath))){
                            tmp_file = tempfile(fileext = ".bed.gz")
                            file.copy(input$rebin_custom_annotation$datapath, tmp_file, overwrite = TRUE)
                            rebin_custom_annotation <- rtracklayer::import(tmp_file)
                        } else {
                            rebin_custom_annotation = NULL
                        }
                        progress$set(message='Rebinning data set into new bins..', detail = "", value = 0.3)
                        print(input$rebin_bin_size)
                        datamatrix = rebin_matrix(mat = datamatrix,
                                                  bin_width = as.numeric(input$rebin_bin_size),
                                                  custom_annotation = rebin_custom_annotation,
                                                  minoverlap = as.numeric(input$minoverlap),
                                                  ref = input$annotation,
                                                  verbose = T)
                        progress$set(message='Finished rebinning data set..', value = 0.8)
                    }
                }
                
            }
            progress$inc(detail=paste0("Finished creating matrix. Saving..."), amount = 0.1)
            
            dir.create(file.path(init$data_folder, "ChromSCape_analyses"), showWarnings = FALSE)
            if(!input$add_to_current_analysis){
                dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name))
                dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "Filtering_Normalize_Reduce"))
                dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "correlation_clustering"))
                dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "correlation_clustering","Plots"))
                dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "Diff_Analysis_Gene_Sets"))
                dir.create(file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "Plots"))
                write.table(input$annotation, file.path(init$data_folder, 'ChromSCape_analyses', input$new_analysis_name, 'annotation.txt'),
                            row.names = FALSE, col.names = FALSE, quote = FALSE)
            }
            
            if(input$add_to_current_analysis){
                qs::qsave(datamatrix, file = file.path(init$data_folder, 
                                                       "ChromSCape_analyses", input$selected_analysis,
                                                       paste0("datamatrix_",alt_name(),".qs")), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                updateRadioButtons(
                    inputId = "feature_select",
                    choices = c("main", get.available.alternative.datasets(input$selected_analysis))
                )
                alt_name("")
                shinyWidgets::updateMaterialSwitch(session = session, inputId = "add_to_current_analysis", value = FALSE)
            } else{
                qs::qsave(datamatrix, file = file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "datamatrix.qs"), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                qs::qsave(annot_raw, file = file.path(init$data_folder, "ChromSCape_analyses", input$new_analysis_name, "annot_raw.qs"), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                init$available_analyses <- list.dirs(path = file.path(init$data_folder, "ChromSCape_analyses"), full.names = FALSE, recursive = FALSE)
                updateSelectInput(session = session, inputId = "selected_analysis",
                                  label =  "Select an Analysis:",
                                  choices = init$available_analyses,
                                  selected =  input$new_analysis_name)
                
                init$datamatrix <- datamatrix
                init$annot_raw <- annot_raw
                
                if(!is.null(raw_mat)) qs::qsave(raw_mat, file = file.path(init$data_folder, 
                                                                          "ChromSCape_analyses", input$new_analysis_name,
                                                                          paste0("raw_mat.qs")), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                rm(raw_mat)
                gc()
            }
            
            progress$inc(detail=paste0("Import successfully finished! "), amount = 0.1)
            updateActionButton(session, "create_analysis", label="Added successfully", icon = icon("check-circle"))
        }
    })
    
    output$download_scExp_UI <- renderUI({
        if(!is.null(scExp()) | !is.null(scExp_cf())){
           
            shinyWidgets::downloadBttn(outputId = "download_scExp", size = "sm",
                                          label = "Download Object", color = "default",
                                       style = "minimal")
        } 
    })
    
    output$download_scExp <- downloadHandler(
        filename = function() {
            paste('scExp_', gsub(":| |-", "_", gsub(":.. ", "",Sys.time())), '.qs', sep='')
        },
        content = function(con) {
            if(!is.null(scExp_cf())){
                scExp. = isolate(scExp_cf())
            } else {
                scExp. = isolate(scExp())
            }
            qs::qsave(scExp., con)
        }
    )
    
    output$dropdown_feature_select_ui <- renderUI({
        if(!is.null(input$selected_analysis) && input$selected_analysis != ""){
            shinydashboardPlus::dropdownBlock(
                id = "feature_select_dropdown",
                badgeStatus = NULL,
                title = shiny::HTML(paste0("<span style='color: white'><h4> <i class='far fa-caret-square-down'",
                                           "role='presentation' aria-label='caret-square-down icon'></i> &nbsp; Features</h4></span>")),
                icon = NULL, # shiny::HTML(paste0("<span style='color: white'>",icon("circle"),"</span>"))
                radioButtons(
                    inputId = "feature_select",
                    label = "Select a feature",
                    choices = c("main", get.available.alternative.datasets(input$selected_analysis))
                ))
        } else {
            shinydashboardPlus::dropdownBlock(
                id = "feature_select_dropdown",
                badgeStatus = NULL,
                title = shiny::HTML(paste0("<span style='color: white'><h4> <i class='far fa-caret-square-down'",
                                           "role='presentation' aria-label='caret-square-down icon'></i> &nbsp; Features</h4></span>")),
                icon = NULL, # shiny::HTML(paste0("<span style='color: white'>",icon("circle"),"</span>"))
                radioButtons(
                    inputId = "feature_select",
                    label = "Select a feature",
                    choices = c("main")
                ))
        }
    })
    
    is_main <- reactive({
        req(input$feature_select)
        ifelse(input$feature_select == "main", TRUE, FALSE)
    })
    
    observeEvent({input$selected_analysis
        input$feature_select}, {  # load precompiled dataset and update coverage plot
            req(input$selected_analysis)
            if(!is.null(input$selected_analysis) && input$selected_analysis != ""){
                if(input$feature_select != "main"){
                    if(input$feature_select %in% get.available.alternative.datasets(input$selected_analysis)) {
                        init$datamatrix <- qs::qread(file.path(init$data_folder, "ChromSCape_analyses", input$selected_analysis,
                                                               paste0("datamatrix_",input$feature_select,".qs")), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                    } else{
                        updateRadioButtons(
                            inputId = "feature_select",
                            choices = c("main", get.available.alternative.datasets(input$selected_analysis))
                        )
                    }
                } else{
                    init$datamatrix <- qs::qread(file.path(init$data_folder, "ChromSCape_analyses", input$selected_analysis, "datamatrix.qs"), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                    init$annot_raw <-  qs::qread(file.path(init$data_folder, "ChromSCape_analyses", input$selected_analysis, "annot_raw.qs"), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                }
                
            }
        })
    
    # When switching analysis, set scExp & scExp_cf to NULL
    observeEvent({input$selected_analysis},
                 {
                     scExp(NULL)
                     scExp_cf(NULL)
                     gc()
                 })
    
    
    output$rename_file_box <- renderUI({
        shinydashboard::box(title = tagList(shiny::icon("signature"), " Rename Samples"),
                            width = NULL, status = "success", solidHeader = TRUE, 
                            collapsible = TRUE, collapsed = TRUE,
                            column(8, htmlOutput("rename_file_UI")),
                            br(),
                            column(6, actionButton("save_rename", "Rename", icon = icon("save")))
        )
    })
    
    output$rename_file_UI <- renderUI({
        req(analysis_name())
        if(nrow(init$datamatrix) > 0 & nrow(init$annot_raw) > 0 ){
            sample_names = unique(init$annot_raw$sample_id)
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
            
            qs::qsave(datamatrix, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "datamatrix.qs"), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
            qs::qsave(annot_raw, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "annot_raw.qs"), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
            
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
    
    output$min_coverage_cell_ui <- renderUI({
        req(input$feature_select)
        req(init$datamatrix)
        if(input$feature_select == "main"){
            sliderInput("min_coverage_cell", shiny::HTML("<p><span style='color: green'>Select minimum number of reads per cell :</span></p>"),
                        min=min(ncol(init$datamatrix),50), max=ncol(init$datamatrix), value=min(200,ncol(init$datamatrix)), step=50) %>%
                shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                    content = "filtering_parameters")
        }
    })
    
    output$quant_removal_ui <- renderUI({
        req(input$feature_select)
        if(input$feature_select == "main"){
            sliderInput("quant_removal", shiny::HTML("<p><span style='color: red'>Select the upper percentile of cells to remove (potential doublets):</span></p>"),
                        min=90, max=100, value=99, step=0.01) 
        } else{
            column(12, align = "middle", br(),h4("Filtering features only, keeping same cells as in 'main' features."),br())
        }
    }) 
    
    output$n_feature_UI <- renderUI({
        req(init$datamatrix)
        sliderInput("n_top_features", shiny::HTML("<p><span style='color: #A020F0'>Select number of top covered features to keep:</span></p>"), 
                    min=min(500,nrow(init$datamatrix)), max=nrow(init$datamatrix), value=nrow(init$datamatrix), step=100)
    })
    
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
    
    output$remove_PC_UI <- renderUI({ 
        if(input$norm_type == "TFIDF"){
            checkboxGroupInput(inputId = "remove_PC",
                               label = "Remove First Components:",
                               choices = paste0("Component_",1:3),
                               selected = "Component_1", 
                               inline = F )
        } else{
            checkboxGroupInput(inputId = "remove_PC",
                               label = "Remove First Components:",
                               choices = paste0("Component_",1:3),
                               inline = F )
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
        print(annotationId)
        if(input$exclude_regions) {
            if(!is.null(input$exclude_file) && file.exists(as.character(input$exclude_file$datapath))){
                tmp_file = tempfile(fileext = ".bed.gz")
                file.copy(input$exclude_file$datapath, tmp_file, overwrite = TRUE)
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
        main_scExp = NULL
        prefix = NULL
        if(input$feature_select != "main"){
            if(!is.null(scExp())){
                main_scExp = getMainExperiment(scExp())
                prefix = input$selected_reduced_dataset
            } else{
                shiny::showNotification(
                    "Please run Filter, Normalize & Reduce on the 'main' features first. (Switch features)",
                    duration = 15, closeButton = TRUE, type="error")
                return(0)
            }
        }
        callModule(Module_preprocessing_filtering_and_reduction, "Module_preprocessing_filtering_and_reduction",
                   reactive({input$selected_analysis}),
                   reactive({input$feature_select}), reactive({main_scExp}), reactive({prefix}),
                   reactive({input$min_coverage_cell}), reactive({input$n_top_features}),
                   reactive({input$quant_removal}),
                   reactive({init$datamatrix}), reactive({init$annot_raw}),
                   reactive({init$data_folder}),reactive({annotationId}),
                   reactive({input$norm_type}), reactive({input$remove_PC}),
                   reactive({exclude_regions}),
                   reactive({input$do_batch_corr}), reactive({batch_sels}),
                   reactive({input$run_tsne}), reactive({subsample_n}))
        
        if(input$feature_select != "main" && !is.null(scExp())){
            file_index <- match(c(input$selected_reduced_dataset), reduced_datasets())
            filename_sel <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),"Filtering_Normalize_Reduce",init$available_reduced_datasets[file_index])
            
            t1 = system.time({
                scExp. = qs::qread(filename_sel, nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                if(is.reactive(scExp.)) {
                    scExp. = isolate(scExp.())
                }
                scExp(scExp.) # retrieve filtered scExp
                rm(scExp.)
                gc()
            })
            cat("Loaded reduced data in ",t1[3]," secs\n")  
            scExp(swapAltExp_sameColData(scExp(), input$feature_select))
        } else {
            init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())
        }
        updateActionButton(session, "filter_normalize_reduce", label="Processed and saved successfully", icon = icon("check-circle"))
    })
    
    observeEvent({input$selected_analysis  # reset label on actionButtion when new filtering should be filtered
        input$min_coverage_cell
        input$quant_removal
        input$n_top_features
        input$do_batch_corr}, {
            updateActionButton(session, "filter_normalize_reduce", label="Filter, Normalize & Reduce", icon = character(0))
        })
    
    observeEvent(input$feature_select,{
        req(scExp())
        if(input$feature_select %in% SingleCellExperiment::altExpNames(scExp())){
            shiny::showNotification(paste0("Switching to ", input$feature_select," !"), duration = 10, closeButton = TRUE, type="message")
            scExp(swapAltExp_sameColData(scExp(), input$feature_select))
        } else{
            shiny::showNotification(paste0("The features ", input$feature_select," were not processed yet.",
                                           " Please please run 'Filter and Normalize' first."), duration = 10, closeButton = TRUE, type="message")
        }
    })
    
    observeEvent(input$selected_reduced_dataset, { # load reduced data set to work with on next pages
        req(input$selected_reduced_dataset)
        
        if(file.exists(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots"))) 
            addResourcePath('Plots', file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering","Plots"))
        
        file_index <- match(c(input$selected_reduced_dataset), reduced_datasets())
        filename_sel <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),"Filtering_Normalize_Reduce",init$available_reduced_datasets[file_index])
        
        scExp(NULL)
        t1 = system.time({
            scExp. = qs::qread(filename_sel, nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
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
        req(init$datamatrix, input$quant_removal)
        index = round(ncol(init$datamatrix) * as.numeric(input$quant_removal) * 0.01)
        q = cell_cov_df()$coverage[index]
        q
    })
    
    cell_cov_plot <- reactive({
        req(input$feature_select)
        p = ggplot(cell_cov_df(), aes(x = coverage)) + 
            geom_histogram(color="black", fill="steelblue", bins = 75) +
            labs(x="Log10(Reads per cell)", y = "nCells")  + 
            theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), 
                  panel.background=element_blank(), axis.line=element_line(colour="black"),
                  panel.border=element_rect(colour="black", fill=NA)) + scale_x_log10()
        if(input$feature_select == "main"){
            p = p + geom_vline(xintercept = as.numeric(input$min_coverage_cell), color = "#22AD18") + 
                geom_vline(xintercept = quantile_threshold(), color = "#D61111")
            
        }
        p
        
    })
    
    output$cell_coverage <- plotly::renderPlotly( plotly::ggplotly(cell_cov_plot(), 
                                                                   tooltip="Sample", dynamicTicks=TRUE) )
    
    feature_cov_df <- reactive ({
        req(init$datamatrix)
        df = data.frame(coverage = sort(unname(Matrix::rowSums(init$datamatrix)),decreasing = TRUE))
        df = df[which(df$coverage>10),,drop=FALSE]
        df
    })  # used for plotting feature coverage on first page
    
    top_feature_min_reads =  reactive({
        req(feature_cov_df(), input$n_top_features)
        min_reads = min(feature_cov_df()$coverage[seq_len(input$n_top_features)])
        min_reads
    })
    
    feature_cov_plot <- reactive({
        ggplot(feature_cov_df(), aes(x = coverage)) + 
            geom_histogram(color="black", fill="#F2B066", bins = 75) +
            labs(x="Log10(Cells 'ON' per feature)", y = "nFeatures")  + 
            theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
                  panel.background=element_blank(), axis.line=element_line(colour="black"),
                  panel.border=element_rect(colour="black", fill=NA)) +
            geom_vline(xintercept = as.numeric(top_feature_min_reads()), color = "#9C36CF") +
            scale_x_log10()
        
    })
    
    output$feature_coverage <- plotly::renderPlotly(plotly::ggplotly(feature_cov_plot(), 
                                                                     tooltip="Feature", dynamicTicks=TRUE) )
    
    output$num_cell <- function(){
        req(init$annot_raw)
        tab = num_cell_scExp(init$annot_raw, init$datamatrix)
        tab
    }
    
    output$num_cell_after_QC_filt <- function(){
        req(input$selected_reduced_dataset,scExp())
        tab = num_cell_after_QC_filt_scExp(scExp(),init$annot_raw, init$datamatrix)
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
    
    output$pc_select_x <- renderUI({
        req(scExp())
        selectInput("pc_select_x", "X",choices=colnames(SingleCellExperiment::reducedDim(scExp(),"PCA"))) })
    output$pc_select_y <- renderUI({
        req(scExp())
        selectInput("pc_select_y", "Y",choices=colnames(SingleCellExperiment::reducedDim(scExp(),"PCA")),
                    selected = colnames(SingleCellExperiment::reducedDim(scExp(),"PCA"))[2]) })
    
    pca_plot <- reactive({
        req(scExp(), annotCol(), input$pc_select_x,input$pc_select_y,  input$color_by)
        if(input$color_by %in% colnames(SingleCellExperiment::colData(scExp())) ){
            if(input$pc_select_x  %in% colnames(SingleCellExperiment::reducedDim(scExp(), "PCA"))){
                p = plot_reduced_dim_scExp(scExp(),input$color_by, "PCA",
                                           select_x = input$pc_select_x,
                                           select_y = input$pc_select_y,
                                           transparency = input$options.dotplot_transparency,
                                           size = input$options.dotplot_size,
                                           max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                           downsample = input$options.dotplot_downsampling,
                                           min_quantile = input$options.dotplot_min_quantile,
                                           max_quantile = input$options.dotplot_max_quantile
                )
                unlocked$list$pca=TRUE
                p
            }
        }
    })
    output$pca_plot <- renderPlot(pca_plot())
    
    contrib_features_plot <- reactive({
        req(scExp(),  input$pc_select_x)
        if(input$pc_select_x %in% colnames(SingleCellExperiment::reducedDim(scExp(), "PCA"))){
            p = plot_most_contributing_features(scExp(), component = input$pc_select_x,
                                                n_top_bot = as.numeric(input$n_features_contributing))
            p
        }
    })
    
    output$contrib_features_plot <- plotly::renderPlotly(contrib_features_plot())
    
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
                shinydashboard::box(title=tagList(shiny::icon("fas fa-chart-bar"), " Contribution to PCA"), collapsible = TRUE, collapsed = TRUE, width = NULL, status="success", solidHeader=TRUE,
                                    column(12, align="left",
                                           selectInput("n_features_contributing", label =  "Top features",
                                                       choices = 5:100,
                                                       multiple = FALSE, selected = 10),
                                           h3(paste0("Most contributing features to '", input$pc_select_x,"' .")),
                                           plotly::plotlyOutput("contrib_features_plot") %>% 
                                               shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75) %>%
                                               shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                                                   content = "most_contributing_features")
                                    ),
                                    column(12, align="left", 
                                           h3("Contribution of chromosomes in top 100 features."),
                                           plotOutput("contrib_chr_plot") %>% 
                                               shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)
                                    ))
            }
        }
    })
    
    
    output$correlation_to_pca_UI <- renderUI({
        req(scExp())
        if("PCA" %in% SingleCellExperiment::reducedDimNames(scExp())){
            n = ncol(SingleCellExperiment::reducedDim(scExp(), "PCA"))
            choices_group_by = c("None",intersect(colnames(SingleCellExperiment::colData(scExp())), c("sample_id", "batch_id")))
            shinydashboard::box(title=tagList(shiny::icon("fas fa-chart-bar"), " Correlation of library size to PCA"),
                                collapsible = TRUE, collapsed = TRUE, width = NULL, status="success", solidHeader=TRUE,
                                column(12, align="left",
                                       selectInput("color_by_pca_cor", "Group by", choices = choices_group_by),
                                       selectInput("n_PCA", label =  "Top Components",
                                                   choices = 1:n, multiple = FALSE, selected = n),
                                       plotly::plotlyOutput("correlation_to_PCA_plot") %>% 
                                           shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)
                                ))
            
        }
    })
    
    correlation_to_PCA_plot <- reactive({
        req(scExp(), input$color_by_pca_cor, input$n_PCA)
        if(input$color_by_pca_cor == "None") color_by_pca_cor = NULL else color_by_pca_cor = input$color_by_pca_cor
        p = plot_correlation_PCA_scExp(scExp(),correlation_var = "total_counts", color_by = color_by_pca_cor, topPC = as.numeric(input$n_PCA) )
        p
        
    })
    
    output$correlation_to_PCA_plot <- plotly::renderPlotly(correlation_to_PCA_plot())
    
    
    
    output$tsne_box <- renderUI({
        req(scExp(), annotCol(), input$color_by)
        if("TSNE" %in% SingleCellExperiment::reducedDimNames(scExp())){
            if(input$color_by %in% colnames(SingleCellExperiment::colData(scExp())) ){
                p = plot_reduced_dim_scExp(scExp(), input$color_by, "TSNE",
                                           transparency = input$options.dotplot_transparency,
                                           size = input$options.dotplot_size,
                                           max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                           downsample = input$options.dotplot_downsampling,
                                           min_quantile = input$options.dotplot_min_quantile,
                                           max_quantile = input$options.dotplot_max_quantile)
                output$tsne_plot = renderPlot(p)
                shinydashboard::box(title="t-SNE vizualisation 1", width = NULL, status="success", solidHeader=TRUE,
                                    column(12, align="left", plotOutput("tsne_plot") %>% 
                                               shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75) %>%
                                               shinyhelper::helper(type = 'markdown',  colour = "#434C5E", icon ="info-circle",
                                                                   content = "tsne_plot")
                                    ))
            }
        }
    })
    
    output$UMAP_plot <- renderPlot({
        req(scExp(), annotCol(), input$color_by)
        if(input$color_by %in% colnames(SingleCellExperiment::colData(scExp())) ){
            p = plot_reduced_dim_scExp(scExp(), input$color_by, "UMAP",
                                       transparency = input$options.dotplot_transparency,
                                       size = input$options.dotplot_size,
                                       max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                       downsample = input$options.dotplot_downsampling,
                                       min_quantile = input$options.dotplot_min_quantile,
                                       max_quantile = input$options.dotplot_max_quantile)
            p
        }
    })
    
    output$UMAP_popup_plot <- renderPlot({
        req(scExp(), annotCol(), input$color_by_popup)
        if(input$color_by_popup %in% colnames(SingleCellExperiment::colData(scExp()))){
            p = plot_reduced_dim_scExp(scExp(), input$color_by_popup, "UMAP",
                                       transparency = input$options.dotplot_transparency,
                                       size = input$options.dotplot_size,
                                       max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                       downsample = input$options.dotplot_downsampling,
                                       min_quantile = input$options.dotplot_min_quantile,
                                       max_quantile = input$options.dotplot_max_quantile)
            p
        }
    })
    
    output$UMAP_popup_plot_cf <- renderPlot({
        req(scExp_cf(), annotCol_cf(), input$color_by_popup_cf)
        if(input$color_by_popup_cf %in% colnames(SingleCellExperiment::colData(scExp_cf()))){
            p = plot_reduced_dim_scExp(scExp_cf(), input$color_by_popup_cf, "UMAP",
                                       transparency = input$options.dotplot_transparency,
                                       size = input$options.dotplot_size,
                                       max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                       downsample = input$options.dotplot_downsampling,
                                       min_quantile = input$options.dotplot_min_quantile,
                                       max_quantile = input$options.dotplot_max_quantile,
                                       annotate_clusters = input$label_cluster_umap_popup,)
            p
        }
    })
    output$popup_UI <- renderUI({
        req(annotCol())
        
        selectInput("color_by_popup", "Color by", choices = annotCol())
    })
    
    output$popup_UI_cf <- renderUI({
        req(annotCol_cf())
        column(12,selectInput("color_by_popup_cf", "Color by", choices = annotCol_cf()),
               checkboxInput("label_cluster_umap_popup", "Label cluster", TRUE)
        )
    })
    
    observeEvent(input$popupUMAP, {
        req(scExp(), annotCol())
        if(!is.null(scExp_cf())){
            if(!is.null(annotCol_cf())){
                showModal(modalDialog(
                    shinydashboard::box(title=tagList(shiny::icon("fas fa-image"), " UMAP colored by sample"), width = NULL, status="success", solidHeader=TRUE,
                                        column(4, align = "left", uiOutput("popup_UI_cf")),
                                        column(12, align = "center", plotOutput("UMAP_popup_plot_cf", width = "800px", height = "600px") %>% 
                                                   shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75))),
                    footer = NULL,
                    easyClose = TRUE,
                    size = "l" 
                ))
            }
        } else{
            if(!is.null(scExp())){
                showModal(modalDialog(
                    shinydashboard::box(title=tagList(shiny::icon("fas fa-image"), " UMAP colored by sample"), width = NULL, status="success", solidHeader=TRUE,
                                        column(3, align = "left", uiOutput("popup_UI")),
                                        column(12, align = "center", plotOutput("UMAP_popup_plot", width = "800px", height = "600px") %>% 
                                                   shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75))),
                    footer = NULL,
                    easyClose = TRUE,
                    size = "l" 
                ))
            }
        }
        
    })
    
    output$color_box <- renderUI({
        req(input$color_by)
        if(!grepl("counts",input$color_by)){
            shinydashboard::box(title = tagList(shiny::icon("palette"), " Color settings"),
                                width = NULL, status = "success", solidHeader = TRUE,
                                column(6, htmlOutput("color_picker")),
                                column(6 , br(), actionButton("col_reset", "Default colours", icon = icon("undo")),
                                       br(), br(), actionButton("save_color", "Save colors & apply to all", icon = icon("save")),
                                       br(), br(), actionButton("save_plots_PCA", "Save HQ plots", icon = icon("fas fa-image"))))
        } else{
            shinydashboard::box(title = tagList(shiny::icon("palette"), " Color settings"),
                                width = NULL, status = "success", solidHeader = TRUE,
                                column(6 ,br(), br(), actionButton("save_plots_PCA", "Save HQ plots", icon = icon("fas fa-image"))))
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
        if(!grepl("counts",input$color_by)){
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
        
        scExp(colors_scExp(scExp(), annotCol = input$color_by, color_by = input$color_by, color_df = color_df))
        
        qs::qsave(getMainExperiment(scExp()), file = file.path(init$data_folder, "ChromSCape_analyses",
                                                               analysis_name(), "Filtering_Normalize_Reduce",
                                                               paste0(input$selected_reduced_dataset,".qs")),
                  nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
        
        rm(color_df)
        gc()
    })
    
    levels_selected <- reactive({
        req(scExp(),input$color_by)
        if(!grepl("counts",input$color_by)) 
            levels_selected = SummarizedExperiment::colData(scExp())[,input$color_by] %>% unique() %>% as.vector()
        else NULL
    })
    
    plot_dir <- reactive({
        req(input$selected_analysis)
        file.path(init$data_folder, "ChromSCape_analyses", input$selected_analysis, "Plots")
    })
    
    observeEvent(input$save_plots_PCA,{
        if(!dir.exists(plot_dir())) dir.create(plot_dir())
        pdf(file = file.path(plot_dir(), paste0("PCA_UMAP_",input$color_by,".pdf")))
        print(pca_plot())
        print(plot_reduced_dim_scExp(scExp(), input$color_by, "UMAP",
                                     transparency = input$options.dotplot_transparency,
                                     size = input$options.dotplot_size,
                                     max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                     downsample = input$options.dotplot_downsampling,
                                     min_quantile = input$options.dotplot_min_quantile,
                                     max_quantile = input$options.dotplot_max_quantile))
        print(pca_plot() + theme(text = element_blank(), legend.position = "none"))
        print(plot_reduced_dim_scExp(scExp(), input$color_by, "UMAP",
                                     transparency = input$options.dotplot_transparency,
                                     size = input$options.dotplot_size,
                                     max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                     downsample = input$options.dotplot_downsampling,
                                     min_quantile = input$options.dotplot_min_quantile,
                                     max_quantile = input$options.dotplot_max_quantile) +
                  theme(text = element_blank(), legend.position = "none"))
        if("t-SNE" %in% names(scExp()@metadata)) print(tsne_plot()  + theme(text = element_blank(),legend.position = "none"))
        dev.off()
        shiny::showNotification(paste0("Plots saved in '",plot_dir(),"' !"),
                                duration = 15, closeButton = TRUE, type="message")
        updateActionButton(session = session, inputId = "save_plots_PCA", label = "Save HQ plots", icon = icon("check-circle"))
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
        counts_cols = paste0("counts_", getExperimentNames(scExp_cf()))
        cluster_cols = paste0("cluster_", getExperimentNames(scExp_cf()))
        if(any(cluster_cols %in% colnames(SummarizedExperiment::colData(scExp_cf())))){
            cols = c(cols, cluster_cols[which(cluster_cols %in% colnames(SummarizedExperiment::colData(scExp_cf())))])
        } 
        if("batch_name" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
            cols = c(cols,"batch_name")
        } 
        cols = c(cols, counts_cols)
        cols
    })
    
    observeEvent(input$tabs, 
                 {
                     if(input$tabs == "cons_clustering"){
                         gc()
                         file = file.path(init$data_folder, "ChromSCape_analyses",
                                          analysis_name(), "correlation_clustering",
                                          paste0(selected_filtered_dataset(),".qs"))
                         print(file)
                         if(file.exists(file)){
                             if(is.null(scExp_cf())){
                                 cat("Loading scExp_cf - ", selected_filtered_dataset(),"...\n")
                                 data = qs::qread(file, nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                                 scExp_cf(data$scExp_cf)
                                 rm(data)
                                 gc()
                                 
                                 if(length(setdiff(get.available.alternative.datasets(input$selected_analysis),
                                                   getExperimentNames(scExp_cf()))) > 0 ){
                                     cat("scExp_cf - to scExp...\n")
                                     scExp_cf(scExp())
                                     gc()  
                                 }
                             }
                             unlocked$list$cor_clust_plot=TRUE;
                             unlocked$list$affectation=TRUE;
                         } else {
                             scExp_cf(scExp())
                             gc()
                         }
                         output$hc_heatmap_plot <- renderPlot({
                             plot_heatmap_scExp(scExp(),downsample = input$options.heatmap_downsampling, color_by = annotCol())
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
    
    output$clustering_method_UI <- renderUI({
        req(scExp_cf())
        if("Cor" %in% SingleCellExperiment::reducedDimNames(scExp_cf())) {
            selectInput("clustering_method", br("Clustering:"), choices=c("louvain", "hierarchical"))
            
        } else{
            selectInput("clustering_method", br("Clustering:"), choices=c("louvain"))
        }
    })
    
    output$clustering_UI = renderUI({
        req(input$clustering_method)
        if(input$clustering_method == "hierarchical") {
            selectInput("nclust", br("Number of Clusters:"), choices=c(2:30))
            
        } else{
            column(12,
                   shinyWidgets::sliderTextInput(inputId = "resolution", label = "Resolution:",
                                                 choices =  c(1e-5, 1e-4, 1e-3, 0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.20, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2),
                                                 selected = 1, grid = TRUE),
                   sliderInput("k_SNN", br("Number of neighbors:"), min = 2, max = 500, step = 1, value = 10)
            )
        }
    })
    
    output$consensus_clustering_UI = renderUI({
        req(input$clustering_method)
        if(input$clustering_method == "hierarchical") {
            checkboxInput("cluster_type",label = shiny::HTML("<b>Use consensus</b>"),value = FALSE)
            
        }
    })
    
    observeEvent(input$feature_select,{
        req(scExp_cf())
        
        if(input$feature_select %in% SingleCellExperiment::altExpNames(scExp_cf())){
            shiny::showNotification(paste0("Switching to ", input$feature_select," !"), duration = 10, closeButton = TRUE, type="message")
            scExp_cf(swapAltExp_sameColData(scExp_cf(), input$feature_select))
        } else{
            shiny::showNotification(paste0("The features ", input$feature_select," were not processed yet.",
                                           " Please please run 'Filter and Normalize' first."), duration = 10, closeButton = TRUE, type="message")
        }
    })
    
    observeEvent({input$choose_cluster},{
        
        req(scExp_cf())
        if(input$clustering_method == "louvain"){
            withProgress(message='Running Louvain...', value = 0.1, {
                incProgress(amount=0.3, detail=paste("Building SNN graph with k=", input$k_SNN))
                print(input$k_SNN)
                scExp_cf(find_clusters_louvain_scExp(scExp_cf(),
                                                     k = input$k_SNN,
                                                     resolution = input$resolution,
                                                     BPPARAM = CS_options.BPPARAM()))
                incProgress(amount=0.6, detail=paste("Finished graph with k=", input$k_SNN) )
                shiny::showNotification(paste0("Found ", length(unique(scExp_cf()$cell_cluster))," clusters with Louvain clustering."),
                                        duration = 10, closeButton = TRUE, type="message")
            })
        } else {
            if(input$nclust != ""){
                scExp_cf(choose_cluster_scExp(scExp_cf(), nclust = as.numeric(input$nclust),
                                              consensus = cluster_type()))
            } else{
                return(NULL)
            }
        }
        odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "coverage", paste0(selected_filtered_dataset(), "_k", length(unique(scExp_cf()$cell_cluster))))
        
        unlocked$list$cor_clust_plot=TRUE;
        unlocked$list$affectation=TRUE;
        gc()
        file = file.path(init$data_folder, "ChromSCape_analyses",
                         analysis_name(), "correlation_clustering",
                         paste0(selected_filtered_dataset(),".qs"))
        # if(!file.exists(file)){
            data = list("scExp_cf" = getMainExperiment(scExp_cf()))
            qs::qsave(data, file=file, nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
            rm(data)
            gc()
        # }
        
        if(file.exists(file.path(odir,"refined_annotation.qs"))){
            # Loading refined peak annotation
            scExp_cf. = scExp_cf()
            scExp_cf.@metadata[["refined_annotation"]] = qs::qread(file = file.path(odir, "refined_annotation.qs"))
            scExp_cf(scExp_cf.)
            rm(scExp_cf.)
        } else{
            
            if("refined_annotation" %in% names(scExp_cf()@metadata)){
                # Reset any refined annotation found with other nclust
                scExp_cf. = scExp_cf()
                scExp_cf.@metadata[["refined_annotation"]] <- NULL
                scExp_cf(scExp_cf.)
                rm(scExp_cf.)
            } 
        }
        cat("New number of clusters:", set_numclust(),"\n")
    })
    
    set_numclust <- reactive({
        req(scExp_cf())
        if("cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp_cf()))){
            length(unique(scExp_cf()$cell_cluster))
        }
    })
    
    observeEvent({
        input$choose_cluster
    }, { # application header (tells you which data set is selected)
        req(analysis_name(), scExp_cf(), set_numclust())
        header <- paste0('<b>Analysis:</b> ', analysis_name(), ' - ',input$feature_select, ' - k = ', set_numclust(), ' ')
        shinyjs::html("pageHeader", header) 
    })
    
    
    output$download_cor_clust_plot <- downloadHandler(
        filename=function(){ paste0("correlation_clustering_", input$selected_reduced_dataset, ".png")},
        content=function(file){
            grDevices::png(file, width=1200, height=1400, res=300,pointsize = 8)
            plot_heatmap_scExp(scExp_cf(), downsample = input$options.heatmap_downsampling)
            grDevices::dev.off()
            
        })
    
    output$consensus_clustering_box = renderUI({
        if("Cor" %in% SingleCellExperiment::reducedDimNames(scExp_cf())){
            shinydashboard::box(title = tagList(shiny::icon("cubes")," Consensus Hierarchical Clustering"),
                                width=NULL, status="success", solidHeader=TRUE,
                                collapsible = TRUE, collapsed = TRUE,
                                column(12, align = "left", br()),
                                column(12, offset = 0,
                                       div(style="display: inline-block;vertical-align:top; width: 200px;", sliderInput("maxK", "Max cluster :", min=2, max=30, value=10, step=1)),
                                       div(style="display: inline-block;vertical-align:top; width: 25px;",HTML("<br>")),
                                       div(style="display: inline-block;vertical-align:top; width: 200px;", sliderInput("consclust_iter", "Number of iterations:", min=10, max=1000, value=100, step=10)),
                                       div(style="display: inline-block;vertical-align:top; width: 25px;",HTML("<br>")),
                                       div(style="display: inline-block;vertical-align:top; width: 200px;", selectInput("clusterAlg", "Cluster Algorithm:", choices = c("Hierarchical","Partitioning Medoids"))),
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
            )
        }
    })
    
    output$correlation_filtering_box = renderUI({
        if("Cor" %in% SingleCellExperiment::reducedDimNames(scExp_cf())){
            shinydashboard::box(title=tagList(shiny::icon("fas fa-filter", verify_fa = FALSE)," Filter lowly correlated cells"),
                                width=NULL, status="success", solidHeader=TRUE,
                                collapsible = TRUE, collapsed = TRUE,
                                column(12, align="center",
                                       plotOutput("cell_cor_hist_plot", height=300, width=500) %>%
                                           shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75) %>%
                                           shinyhelper::helper(type = 'markdown',  colour = "#434C5E", icon ="info-circle",
                                                               content = "filter_correlation_distrib")),
                                column(12, align = "left",
                                       hr(),
                                       sliderInput("corr_threshold", "Correlation threshold quantile:", min=75, max=99, value=99, step = 0.01),
                                       sliderInput("percent_correlation", "Minimum percentage of cells to correlate with (%)", min=0, max=15, value=1, step=0.01)),
                                column(3, align="left", br(),
                                       actionButton("filter_corr_cells", "Filter & save")),
                                column(3, align="left", br(),
                                       actionButton("reset_corr_cells", "Reset filtering")),
                                column(12,  uiOutput("table_cor_filtered")
                                )
            )
        }
    })
    
    sample_cf <- reactive({
        sample(seq_len(ncol(scExp_cf())), min(2000,ncol(scExp_cf())), replace = FALSE)
    })
    corChIP <- reactive({
        as.matrix(SingleCellExperiment::reducedDim(scExp_cf(),"Cor")[sample_cf(),sample_cf()])
    })
    
    z <- reactive({ 
        pca = SingleCellExperiment::reducedDim(scExp_cf(),"PCA")[sample_cf(),]
        random_pca = matrix(data = sample(as.numeric(as.matrix(pca))), nrow = nrow(pca), ncol = ncol(pca))
        z = coop::pcor(t(random_pca), inplace = TRUE)
        diag(z) <- NA
        z = z[which(!is.na(z))]
        z
    })
    
    thresh2 <- reactive({
        quantile(z(), probs=seq(0,1,0.01))
    })
    
    limitC <- reactive({
        thresh2()[input$corr_threshold+1]
    })
    
    cell_cor_hist <- reactive({
        req(scExp())
        hist(corChIP(), prob=TRUE, col=alpha("steelblue", 0.8), breaks=50, ylim=c(0,4), main="Distribution of cell to cell correlation scores", xlab="Pearson Corr Scores")
        lines(density(corChIP()), col="blue", lwd=2)
        lines(density(z()), col="black", lwd=2)
        abline(v=limitC(), lwd=2, col="red", lty=2)
        legend("topleft", legend=c("dataset", "randomized data", "correlation threshold"),
               col=c("blue", "black", "red"), lty=c(1, 1, 2), cex=0.8)
    })
    
    output$cell_cor_hist_plot <- renderPlot(cell_cor_hist())
    
    output$download_cor_clust_hist_plot <- downloadHandler(
        filename=function(){ paste0("correlation_distribution_", input$selected_reduced_dataset, ".png")},
        content=function(file){
            grDevices::png(file, width=2000, height=1400, res=300)
            hist(corChIP(), prob=TRUE, col=alpha("steelblue", 0.8), breaks=50, ylim=c(0,4),cex=0.4, main="Distribution of cell to cell correlation scores", xlab="Pearson Corr Scores")
            lines(density(corChIP()), col="blue", lwd=2)
            lines(density(z()), col="black", lwd=2)
            abline(v=limitC(), lwd=2, col="red", lty=2)
            legend("topleft", legend=c("dataset", "randomized data", "correlation threshold"), col=c("blue", "black", "red"), lty=c(1, 1, 2), cex=0.8)
            grDevices::dev.off()
        })
    
    observeEvent(input$filter_corr_cells, {  # retreiveing cells with low correlation score
        withProgress(message='Filtering correlated cells...', value = 0, {
            incProgress(amount=0.6, detail=paste("Filtering"))
            print("Filtering...")
            scExp_cf(filter_correlated_cell_scExp(scExp_cf(), random_iter = 50, corr_threshold = input$corr_threshold,
                                                  percent_correlation = input$percent_correlation,
                                                  BPPARAM = CS_options.BPPARAM()))
            print("Choosing cluster...")
            print(set_numclust())
            scExp_cf(choose_cluster_scExp(scExp_cf(), nclust = 3, consensus = cluster_type()))
            print("Saving...")
            gc()
            incProgress(amount=0.2, detail=paste("Saving"))
            data = list("scExp_cf" = getMainExperiment(scExp_cf()))
            qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                             paste0(input$selected_reduced_dataset, ".qs")), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
            incProgress(amount=0.2, detail=paste("Finished"))
            rm(data)
            gc()
            updateActionButton(session, "filter_corr_cells", label="Saved", icon = icon("check-circle"))
            
        })
    })
    
    observeEvent({input$reset_corr_cells  # reset label on actionButtion when new filtering should be filtered
    }, {
        req(scExp_cf())
        
        if("Unfiltered" %in% names(scExp_cf()@metadata)){
            shiny::showNotification(paste0("Resetting original cells ..."),
                                    duration = 5, closeButton = TRUE, type="message")
            scExp_cf(scExp_cf()@metadata$Unfiltered)
            updateActionButton(session, "reset_corr_cells", label="Resetted", icon = icon("check-circle"))
            updateActionButton(session, "filter_corr_cells", "Filter & save", icon=character(0))
        } else {
            shiny::showNotification(paste0("Already original cells ..."),
                                    duration = 5, closeButton = TRUE, type="warning")
        }
        
        
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
        # withProgress(message='Performing consensus clustering...', value = 0, {
        # incProgress(amount=0.4, detail=paste("part one"))
        progress <- shiny::Progress$new(session, min=0, max=1)
        on.exit(progress$close())
        progress$set(message='Performing consensus clustering...', value = 0.1)
        scExp_cf(consensus_clustering_scExp(scExp_cf(), reps = as.numeric(input$consclust_iter),
                                            maxK=as.numeric(input$maxK),
                                            clusterAlg = clusterAlg(),
                                            prefix = plotting_directory()))
        gc()
        progress$set(value = 0.7)
        data = list("scExp_cf" = getMainExperiment(scExp_cf()))
        qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "correlation_clustering",
                                         paste0(input$selected_reduced_dataset, ".qs")), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
        rm(data)
        gc()
        clust$clust_pdf <- NULL  # needed in order to update the pdf output
        clust$clust_pdf <- file.path("Plots", selected_filtered_dataset(), "consensus.pdf")
        progress$set(message='Consensus clustering done !', value = 0.9)
        # incProgress(amount=0.2, detail=paste("Finished"))
        # })
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
                     print("Doing annotated heatmap")
                     req(input$clustering_method)
                     if("cell_cluster" %in% colnames(SingleCellExperiment::colData(scExp_cf()))){
                         output$annotated_heatmap_UI <- renderUI({
                             output$annotated_heatmap_plot = renderPlot(plot_heatmap_scExp(scExp_cf(), downsample = input$options.heatmap_downsampling))
                             plotOutput("annotated_heatmap_plot",width = 500,height = 500) %>%
                                 shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)
                         })
                     }
                 }
    )
    
    output$violin_color <- renderUI({
        req(annotCol_cf())
        selectInput("violin_color","Color by", choices=annotCol_cf()[-grep("counts", annotCol_cf())])
    })
    
    output$jitter_color <- renderUI({
        req(annotCol_cf(), input$add_jitter)
        if(input$add_jitter) selectInput("jitter_color","Color cells", choices=annotCol_cf())
    })
    
    output$intra_corr_UI <- renderUI({
        req(annotCol_cf(), input$violin_color)
        
        if(input$add_jitter) jitter_col = input$jitter_color else jitter_col = NULL
        
        output$intra_corr_plot = plotly::renderPlotly(
            plot_intra_correlation_scExp(scExp_cf(), by = input$violin_color,
                                         jitter_by = jitter_col))
        plotly::plotlyOutput("intra_corr_plot",width = 500,height = 500) %>%
            shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)
        
    })
    
    output$reference_group <- renderUI({
        req(annotCol_cf(), scExp_cf(), input$violin_color)
        selectInput("reference_group","Correlate with",
                    choices = unique(scExp_cf()[[input$violin_color]]) )
    })
    
    
    output$inter_corr_UI <- renderUI({
        req(annotCol_cf(), input$violin_color)
        
        if(input$add_jitter) jitter_col = input$jitter_color else jitter_col = NULL
        output$inter_corr_plot = plotly::renderPlotly(
            plot_inter_correlation_scExp(scExp_cf(), by = input$violin_color,
                                         jitter_by = jitter_col,
                                         reference_group = input$reference_group))
        
        plotly::plotlyOutput("inter_corr_plot",width = 500,height = 500) %>%
            shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)
        
    })
    
    observeEvent(input$save_plots_violins,{
        if(!dir.exists(plot_dir())) dir.create(plot_dir())
        if(input$add_jitter) jitter_col = input$jitter_color else jitter_col = NULL
        pdf(file = file.path(plot_dir(), paste0("violins_intra_inter_correlation.pdf")))
        print(plot_intra_correlation_scExp(scExp_cf(), by = input$violin_color,
                                           jitter_by = jitter_col))
        if(!is.null(input$reference_group)) print(plot_inter_correlation_scExp(scExp_cf(), by = input$violin_color,
                                                                               jitter_by = jitter_col,
                                                                               reference_group = input$reference_group))
        
        if(!is.null(input$reference_group))  print(plot_inter_correlation_scExp(scExp_cf(), by = input$violin_color,
                                                                                jitter_by = jitter_col,
                                                                                reference_group = input$reference_group) + theme(legend.position = "none"))
        print(plot_intra_correlation_scExp(scExp_cf(), by = input$violin_color,
                                           jitter_by = jitter_col) + theme(legend.position = "none"))
        dev.off()
        shiny::showNotification(paste0("Plots saved in '",plot_dir(),"' !"),
                                duration = 15, closeButton = TRUE, type="message")
        
        updateActionButton(session = session, inputId = "save_plots_COR", label = "Save HQ plots", icon = icon("check-circle"))
    })
    
    output$contingency_table_cluster <- renderUI({
        if(! is.null(scExp_cf())){
            if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
                shinydashboard::box(title=tagList(icon("th-large"), " Samples & Cluster association table"), width = NULL, status="success", solidHeader = TRUE,
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
                                               select_y = "Component_2",
                                               transparency = input$options.dotplot_transparency,
                                               size = input$options.dotplot_size,
                                               max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                               downsample = input$options.dotplot_downsampling,
                                               min_quantile = input$options.dotplot_min_quantile,
                                               max_quantile = input$options.dotplot_max_quantile) +
                        ggtitle("t-SNE")
                    output$tsne_plot_cf <- renderPlot(p)
                    shinydashboard::box(title="t-SNE", width = NULL, status="success", solidHeader=TRUE,
                                        column(12, align="left", plotOutput("tsne_plot_cf")))
                }
            }
        }
    })
    
    umap_p_cf <- reactive({
        req(scExp_cf(), annotCol_cf(), input$color_by_cf)
        if(input$color_by_cf %in% colnames(SingleCellExperiment::colData(scExp_cf())) ){
            p = plot_reduced_dim_scExp(scExp_cf(),input$color_by_cf, "UMAP",
                                       select_x = "Component_1",
                                       select_y = "Component_2",
                                       transparency = input$options.dotplot_transparency,
                                       size = input$options.dotplot_size,
                                       max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                       downsample = input$options.dotplot_downsampling,
                                       min_quantile = input$options.dotplot_min_quantile,
                                       max_quantile = input$options.dotplot_max_quantile,
                                       annotate_clusters = input$label_cluster_umap)
            p
        }
    })
    
    output$plot_CF_UMAP <- renderPlot({
        umap_p_cf()
    })
    
    
    levels_selected_cf <- reactive({
        req(scExp_cf(),input$color_by_cf)
        if( !grepl("counts",input$color_by_cf) ) levels_selected_cf =
                SummarizedExperiment::colData(scExp_cf())[,input$color_by_cf] %>%
                unique() %>% as.vector() else NULL
    })
    
    output$UMAP_box <- renderUI({
        if(!is.null(scExp_cf())){
            if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf())) ){
                shinydashboard::box(title= tagList(shiny::icon("fas fa-image"), " UMAP"), width = NULL, status="success", solidHeader = TRUE,
                                    column(4, align="left",
                                           selectInput("color_by_cf", "Color by",
                                                       selected = annotCol_cf()[max(grep("cluster",annotCol_cf())[1],1)],
                                                       choices = annotCol_cf())),
                                    column(4, align = "left", checkboxInput("label_cluster_umap", "Label cluster", TRUE)),
                                    column(3, align = "left", actionButton(inputId = "save_plots_COR",
                                                                           label = "Save HQ plots",
                                                                           icon = icon("fas fa-image"))),
                                    column(12, align="left", plotOutput("plot_CF_UMAP")))
            }
        }
    })
    
    observeEvent(input$save_plots_COR,{
        if(!dir.exists(plot_dir())) dir.create(plot_dir())
        pdf(file = file.path(plot_dir(), paste0("UMAP_heatmap_",input$color_by_cf,".pdf")))
        print(umap_p_cf())
        print(umap_p_cf() + theme(text = element_blank(), legend.position = "none"))
        plot_heatmap_scExp(isolate(scExp_cf()), downsample = input$options.heatmap_downsampling)
        dev.off()
        shiny::showNotification(paste0("Plots saved in '",plot_dir(),"' !"),
                                duration = 15, closeButton = TRUE, type="message")
        updateActionButton(session = session, inputId = "save_plots_COR", label = "Save HQ plots", icon = icon("check-circle"))
    })
    
    output$color_box_cf <- renderUI({
        req(input$color_by_cf)
        if(! is.null(scExp_cf())){
            if("cell_cluster" %in% colnames(SummarizedExperiment::colData(scExp_cf()))){
                if(!grepl("counts",input$color_by_cf)){
                    shinydashboard::box(title=tagList(shiny::icon("palette"), " Color settings "),collapsible = TRUE, collapsed = TRUE,
                                        width = NULL, status = "success", solidHeader = TRUE,
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
    
    observeEvent(unlocked$list, {able_disable_tab(c("selected_reduced_dataset","affectation"), c("coverage", "diff_analysis"))}) 
    
    ###############################################################
    # 5. Coverage plot [optional]
    ###############################################################
    
    output$coverage_info <- renderUI({
        HTML(paste0("In order to visualize the cumulative signal of each subpopulation (cluster), create & plot the coverage tracks of some loci of interest.",
                    "The coverage tracks are saved as bigwig files in the analysis directory (under coverage/). <br/>",
                    "The recommended way to create the coverage tracks is to start your analysis from genomic bins <300bp (with re-binning to larger bin sizes). If this",
                    " is the case, the raw matrix has been saved and the coverage tracks will be created ith this bin size for each cluster",
                    " directly from the raw matrix. <br/>",
                    "If you did not start with a small bin size you have to upload the root directory containing the folders of scBED files of each sample.",
                    " The names of the folders placed in the selected directory must exactly match the sample names and the scBED must be named as in the SingleCellExperiment object. <br/>" 
        )
        )
    })
    
    shinyFiles::shinyDirChoose(input, "coverage_folder", roots = volumes, session = session)
    
    coverage_folder = reactiveVal(NULL)
    
    observeEvent(input$coverage_folder, priority = 10000, {
        if(showHelpSC()) {
            showModal(input_sc())
            showHelpSC(FALSE)
        } else { showHelpSC(TRUE)}
    })
    
    list_files_coverage = reactive({
        req(coverage_folder())
        if(!is.null(coverage_folder())){
            list.files(coverage_folder(), full.names = TRUE, pattern = "*.bed|*.bed.gz", recursive = TRUE)
        }
    })
    
    list_dirs_coverage = reactive({
        req(coverage_folder())
        if(!is.null(coverage_folder())){
            setNames(list.dirs(coverage_folder(), full.names = TRUE, recursive = FALSE), basename(list.dirs(coverage_folder(), recursive = FALSE)))
        }
    })
    
    raw_mat_generated <- reactive({
        req(init$data_folder, analysis_name())
        odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name())
        ifelse("raw_mat.qs" %in% list.files(odir),TRUE,FALSE)
    })
    
    coverage_folder_ui <- renderUI({
        if(!raw_mat_generated()){
            shinyFiles::shinyDirButton("coverage_folder", icon = icon("folder-open"), 
                                       "Browse directory of raw signal (scBED)" ,
                                       title = "Please select a folder:",
                                       buttonType = "default", class = NULL) %>%
                shinyhelper::helper(type = 'markdown', colour = "#434C5E", icon ="info-circle",
                                    content = "coverage")
        }
        
    })
    observeEvent(input$coverage_folder,
                 {
                     req(input$coverage_folder)
                     if(!is.null(input$coverage_folder)){
                         # browser()
                         coverage_folder(shinyFiles::parseDirPath(volumes, input$coverage_folder))
                     }
                 })
    
    output$coverage_upload <- renderUI({
        req(list_files_coverage(),  list_dirs_coverage())
        if(!is.null(list_files_coverage())){
            dirs = basename(list_dirs_coverage())
            selectInput("coverage_selection", label = "Selected samples:",
                        choices = dirs, multiple = TRUE,
                        selected = dirs)
            
        }
    })
    
    has_available_coverage <- reactiveVal(FALSE)

    observe({
        req(input$tabs, scExp_cf(), set_numclust(), analysis_name(), selected_filtered_dataset())
        
        if(input$tabs == "coverage"){
            if(!is.null(scExp_cf())){
                nclust = set_numclust()
                odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "coverage", paste0(selected_filtered_dataset(), "_k", nclust))
                if(length(list.files(odir, pattern = ".bw|.bigWig|.bigwig")) > 0){
                    has_available_coverage(TRUE)
                } else has_available_coverage(FALSE)
            } else{
                has_available_coverage(FALSE)
            }
        }
        
    })
    
    observeEvent(input$do_coverage, {
        req(scExp_cf(), set_numclust())
        progress <- shiny::Progress$new(session, min=0, max=1)
        on.exit(progress$close())
        progress$set(message='Generating coverage tracks...', value = 0.0)
        nclust = set_numclust()
        
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "coverage"), showWarnings = FALSE)
        dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "coverage", paste0(selected_filtered_dataset(), "_k", nclust)), showWarnings = FALSE)
        odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "coverage", paste0(selected_filtered_dataset(), "_k", nclust))
        
        if(raw_mat_generated()){

            raw_mat =  qs::qread(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "raw_mat.qs"))
            generate_coverage_tracks(scExp_cf = scExp_cf(),
                                     input = raw_mat,
                                     odir = odir,
                                     format = "raw_mat",
                                     bin_width = as.numeric(input$coverage_bin_size),
                                     n_smoothBin = as.numeric(input$coverage_n_smoothBin),
                                     quantile_for_peak_calling =  as.numeric(input$quantile_for_peak_calling),
                                     ref_genome = annotation_id(),
                                     progress = progress)
            has_available_coverage(TRUE)
            updateActionButton(session, "do_coverage", label="Finished successfully", icon = icon("check-circle"))
        } else{
            if(is.null(input$coverage_selection) | length(input$coverage_selection) ==0 ) {
                warning("Can't find any input BED files.")
                return()
            }
            input_files_coverage = lapply(list_dirs_coverage()[input$coverage_selection], function(i) list.files(i, full.names = TRUE, pattern = ".bed|.bed.gz"))
            names(input_files_coverage) = basename(list_dirs_coverage()[input$coverage_selection])
            
            if(length(input_files_coverage)==0){
                warning("Can't find any input single-cell BED files. Please make sure you selected the root of a directory",
                        " containing one folder per sample. Each folder should contain single-cell raw reads as .bed or .bed.gz file (one file per cell).")
            } else{

                checkFiles = 0
                if(sum(checkFiles)==0){
                    generate_coverage_tracks(scExp_cf = scExp_cf(),
                                             input = input_files_coverage,
                                             odir = odir,
                                             format = "BED",
                                             bin_width = as.numeric(input$coverage_bin_size),
                                             n_smoothBin = as.numeric(input$coverage_n_smoothBin),
                                             quantile_for_peak_calling =  as.numeric(input$quantile_for_peak_calling),
                                             ref_genome = annotation_id(),
                                             progress = progress)
                    has_available_coverage(TRUE)
                    updateActionButton(session, "do_coverage", label="Finished successfully", icon = icon("check-circle"))
                }
            }
        }
        # Move peaks to peaks folder. Create the folder if not existing
        peak_files = list.files(path = odir, pattern = ".bed|.bed.gz", full.names = TRUE)
        if(length(peak_files) > 0){
            message("Creating consensus peak annotation by merging peaks called on each distinct cluster...")
            merged_peaks = sapply(peak_files, rtracklayer::import.bed)
            
            eval(parse(text = paste0("data(", annotation_id(), ".GeneTSS)")))
            geneTSS_annotation = eval(parse(text = paste0("", annotation_id(), ".GeneTSS")))
            geneTSS_annotation = as(geneTSS_annotation, "GRanges")
            scExp_cf. = scExp_cf()
            refined_annotation <- annotation_from_merged_peaks(scExp_cf(), odir, merged_peaks, geneTSS_annotation)
            
            scExp_cf.@metadata$refined_annotation = refined_annotation
            scExp_cf(scExp_cf.)
            rm(scExp_cf.)
            gc()
            qs::qsave(refined_annotation, file = file.path(odir, "refined_annotation.qs"))
            
        }
    })

    
    output$coverage_UI <- renderUI({
        req(GenePool(),has_available_coverage())
        s <- shinydashboard::box(title="Coverage visualization", width = NULL, status="success", solidHeader = TRUE,
                                 column(3, align="left", selectizeInput(inputId = "select_cov_gene", "Select gene:", choices = NULL, selected = 1)),
                                 column(2, align="left", selectInput("cov_chr","Chromosome", choices = unique(as.character(SummarizedExperiment::seqnames(SummarizedExperiment::rowRanges(scExp_cf())))))),
                                 column(2, align="left", textInput("cov_start","Start", value = 15000000)),
                                 column(2, align="left", textInput("cov_end","End", 16000000)),
                                 column(2, align="left", actionButton("make_plot_coverage", "Plot Coverage", icon = icon("chart-area"))),
                                 column(12, align="left", textOutput("top_genes")))
        updateSelectizeInput(session = session,
                             inputId = 'select_cov_gene',
                             choices = GenePool(),
                             server = TRUE)
        return(s)
        
    })
    
    output$coverage_plot_UI <- renderUI({
        req(GenePool(), has_available_coverage(), display_coverage_plot())
        if(has_available_coverage()){
            if(display_coverage_plot()){
                shinydashboard::box(title="Coverage visualization", width = NULL, status="success", solidHeader = TRUE,
                                    column(12, align="left", plotOutput("coverage_region_plot") %>%
                                               shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)),
                                    column(3, actionButton("save_plots_coverage", "Save HQ plot", icon =  icon("fas fa-image")))) 
            }
        }
    })
    
    
    GenePool <- reactive({
        req(annotation_id())
        eval(parse(text = paste0("data(", annotation_id(), ".GeneTSS)")))
        GenePool = unique(eval(parse(text = paste0("", annotation_id(), ".GeneTSS$Gene"))))
        c("Enter Gene...", GenePool)
    })
    
    
    coverages <- reactive({
        withProgress(message='Loading coverage (.bw) tracks...', value = 0, {
            incProgress(amount = 0.3)
            req(has_available_coverage(), scExp_cf())
            display_coverage_plot(TRUE)
            if(has_available_coverage()){
                nclust = set_numclust()
                odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),
                                  "coverage", paste0(selected_filtered_dataset(), "_k", nclust))
                BigWig_filenames <- list.files(odir, pattern = "*.bw", full.names = TRUE)
                if(length(BigWig_filenames)>0){
                    incProgress(amount = 0.3)
                    coverages = sapply(BigWig_filenames, rtracklayer::import)
                    incProgress(amount = 0.3, detail = paste("Done !"))
                } else {
                    coverages = NULL
                }
                coverages
            }
        })
    })
    
    display_coverage_plot <- reactiveVal(FALSE)
    top_genes <- reactive({
        req(scExp_cf())
        top_feats = ChromSCape:::retrieve_top_bot_features_pca(SingleCellExperiment::reducedDim(scExp_cf(),"PCA"),
                                                               component = colnames(SingleCellExperiment::reducedDim(scExp_cf(),"PCA"))[1],
                                                               counts = SingleCellExperiment::counts(scExp_cf()),
                                                               n_top_bot = 50, absolute = TRUE)
        
        genes = SummarizedExperiment::rowRanges(scExp_cf())$Gene[which(rownames(top_feats) %in% rownames(scExp_cf()))]
        genes = paste0(paste(genes[seq_len(min(10, length(genes)))], collapse = ", "),"...")
        genes
    })
    output$top_genes = renderText({paste0("Top contributing genes: ", top_genes())})
    
    observeEvent(input$make_plot_coverage, { # load reduced data set to work with on next pages
        req(input$cov_chr, has_available_coverage(),
            scExp_cf(), coverages(), annotation_id())
        if(input$select_cov_gene != "Enter Gene..."){
            eval(parse(text = paste0("data(", annotation_id(), ".GeneTSS)")))
            gene_annot = eval(parse(text = paste0(annotation_id(), ".GeneTSS")))
            updateSelectInput(session, "cov_chr", selected = gene_annot$chr[which(gene_annot$Gene == input$select_cov_gene)])
            updateSelectInput(session, "cov_start", selected = gene_annot$start[which(gene_annot$Gene == input$select_cov_gene)] - 25000)
            updateSelectInput(session, "cov_end", selected = gene_annot$end[which(gene_annot$Gene == input$select_cov_gene)] + 25000)
            updateSelectInput(session, "select_cov_gene", selected = "Enter Gene...")
        }
        label_color_list = setNames(unique(scExp_cf()$cell_cluster_color), unique(scExp_cf()$cell_cluster))
        nclust = set_numclust()
        odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),
                          "coverage", paste0(selected_filtered_dataset(), "_k", nclust))
        if(file.exists(file.path(odir, "consensus_peaks.bed"))){
            peaks = rtracklayer::import.bed(file.path(odir, "consensus_peaks.bed"))
        } else{peaks = NULL}
        
        output$coverage_region_plot <- renderPlot({
            plot_coverage_BigWig(
                coverages(), 
                label_color_list,
                chrom = input$cov_chr,
                start = as.numeric(input$cov_start),
                end =  as.numeric(input$cov_end),
                ref = annotation_id(),
                peaks = peaks)
        })
        
    })
    
    observeEvent(input$save_plots_coverage, {
        label_color_list = setNames(unique(scExp_cf()$cell_cluster_color), unique(scExp_cf()$cell_cluster))
        nclust = set_numclust()
        odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),
                          "coverage", paste0(selected_filtered_dataset(), "_k", nclust))
        if(file.exists(file.path(odir, "consensus_peaks.bed"))){
            peaks = rtracklayer::import.bed(file.path(odir, "consensus_peaks.bed"))
        } else{peaks = NULL}
        
        if(!dir.exists(plot_dir())) dir.create(plot_dir())
        pdf(file.path(plot_dir(), paste0("coverage_",input$cov_chr,"_",input$cov_start,"_",input$cov_end, ".pdf") ))
        
        plot_coverage_BigWig(
            coverages(), 
            label_color_list,
            chrom = input$cov_chr,
            start = as.numeric(input$cov_start),
            end =  as.numeric(input$cov_end),
            ref = annotation_id())
        dev.off()
        
        shiny::showNotification(paste0("Plots saved in '",plot_dir(),"' !"),
                                duration = 15, closeButton = TRUE, type="message")
        
        updateActionButton(session = session, inputId = "save_plots_coverage", label = "Save HQ plots", icon = icon("check-circle"))
    })
    
    ###############################################################
    # 5. Peak calling [optional]
    ###############################################################
    
    # output$peak_calling_info <- renderText({"This module is optional, but recommended in order to obtain the most meaningful results for pathway enrichment analysis. Peaks will be called from the BAM files of the samples selected in your project, using MACS2 [only works on unix systems] so that counts can be assigned more specifically to genes TSS . If you have MACS2 installed but ChromSCape can???t find these softwares, try relaunching R from the terminal and start ChromSCape again."})
    # 
    # can_run = reactiveVal({FALSE})
    # 
    # output$peak_calling_system <- renderText({
    #     platform = as.character(.Platform[1])
    #     if(length(grep("nix",platform,ignore.case = TRUE)) ){
    #         macs2=""
    #         try({
    #             macs2 = system2("which",args = "macs2",stdout = TRUE)
    #         })
    #         if(length(macs2)>0){
    #             can_run(TRUE)
    #             return(paste0("<b>You are running on an ", platform, " OS.<br>Found MACS2 at ", macs2))
    #         }
    #         if(length(macs2)==0)
    #             return(paste0("<b>You are running on an ", platform, " OS.<br>Didn't find MACS2, please install MACS2 or skip this step."))
    #     } else {
    #         return(paste0("<b>You are running on a non unix system, MACS2 peak calling is not available. If you ran coverage you",
    #         "you can move directly to differential analysis.</b> "))
    #     }
    # })
    # 
    # output$peak_calling_icon = renderText({
    #     if(can_run()) {
    #         return( as.character(icon("check-circle", class = "large_icon")))}
    #     else{
    #         return( as.character(icon("times-circle", class = "large_icon")))
    #     }
    # })
    # 
    # shinyFiles::shinyDirChoose(input, "pc_folder", roots = volumes, session = 
    #                                session)
    # 
    # pc_folder = reactiveVal(NULL)
    # 
    # observeEvent(input$pc_folder, priority = 10000, {
    #     if(showHelpSC()) {
    #         showModal(input_sc())
    #         showHelpSC(FALSE)
    #     } else { showHelpSC(TRUE)}
    # })
    # 
    # list_files_pc = reactive({
    #     req(pc_folder())
    #     if(!is.null(pc_folder())){
    #         list.files(pc_folder(), full.names = TRUE, pattern = "*.bam$|*.bed|*.bed.gz", recursive = TRUE)
    #     }
    # })
    # 
    # list_dirs_pc = reactive({
    #     req(pc_folder())
    #     if(!is.null(pc_folder())){
    #         setNames(list.dirs(pc_folder(), full.names = TRUE, recursive = FALSE), basename(list.dirs(pc_folder(), recursive = FALSE)))
    #     }
    # })
    # 
    # bam_or_bed <- reactive({
    #     req(list_files_pc())
    #     nBAMs = length(grep(".*.bam$", list_files_pc()))
    #     nBEDs = length(grep(".*.bed.gz$|.*.bed$", list_files_pc()))
    #     ifelse(nBAMs>nBEDs,"BAM","BED")
    # })
    # 
    # observeEvent(input$pc_folder,
    #              {
    #                  req(input$pc_folder)
    #                  if(!is.null(input$pc_folder)){
    #                      # browser()
    #                      pc_folder(shinyFiles::parseDirPath(volumes, input$pc_folder))
    #                      output$pc_dir <- renderText(bam_or_bed())
    #                  }
    #              })
    # 
    # output$pc_upload <- renderUI({
    #     req(list_files_pc(), bam_or_bed(), list_dirs_pc())
    #     if(!is.null(list_files_pc()) & !is.null(bam_or_bed())){
    #         if(bam_or_bed() == "BAM") {
    #             files = list_files_pc()[grep(".bam$", list_files_pc())]
    #             selectInput("bam_selection", label = "Selected samples:",
    #                         choices = basename(files), multiple = TRUE,
    #                         selected = basename(files) )
    #         } else {
    #             dirs = basename(list_dirs_pc())
    #             selectInput("bed_selection", label = "Selected samples:",
    #                         choices = dirs, multiple = TRUE,
    #                         selected = dirs)
    #         }
    #     }
    # })
    # 
    # observeEvent(input$do_pc, {
    #     req(scExp_cf())
    #     if(is.null(input$pc_folder) | length(input$pc_folder)==0 | length(input$bed_selection)==0){
    #         shiny::showNotification(paste0("Please select a valid folder for ",
    #                                        "the raw data before running peak_calling..."),
    #                                 duration = 10, closeButton = TRUE, type="warning")
    #         return()
    #     }
    #     
    #     progress <- shiny::Progress$new(session, min=0, max=1)
    #     on.exit(progress$close())
    #     progress$set(message='Performing peak calling and coverage...', value = 0.0)
    #     
    #     if(bam_or_bed() == "BAM") {
    #         input_files_pc <- as.character(file.path(dirname(list_files_pc()),input$bam_selection))
    #     } else {
    #         if(is.null(input$bed_selection) | length(input$bed_selection) ==0 ) {
    #             warning("Can't find any input BED files.")
    #             return()
    #         }
    #         input_files_pc = lapply(list_dirs_pc()[input$bed_selection], function(i) list.files(i, full.names = TRUE, pattern = ".bed|.bed.gz"))
    #         names(input_files_pc) = basename(list_dirs_pc()[input$bed_selection])
    #     }
    #     if(length(input_files_pc)==0){
    #         warning("Can't find any input BAM / BED files.")
    #     } else{
    #         nclust = set_numclust()
    #         
    #         dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks"), showWarnings = FALSE)
    #         dir.create(file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", nclust)), showWarnings = FALSE)
    #         odir <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "peaks", paste0(selected_filtered_dataset(), "_k", nclust))
    #         sample_ids <- unique(SummarizedExperiment::colData(scExp_cf())$sample_id)
    #         
    #         # checkFiles <- sapply(input_files_pc, function(x){ if(file.exists(x)){ 0 } else {
    #         #   showNotification(paste0("Could not find file ", x, ". Please make sure to give a full path including the file name."),
    #         #                    duration = 7, closeButton = TRUE, type="warning"); 1}
    #         # })
    #         checkFiles = 0
    #         if(sum(checkFiles)==0){
    #             
    #             progress$set(message='Performing peak calling ...', value = 0.1)
    #             scExp_cf(subset_bam_call_peaks(scExp_cf(), odir, input_files_pc, format = bam_or_bed(),
    #                                            as.numeric(input$pc_stat_value), annotation_id(),
    #                                            input$peak_distance_to_merge, progress = progress))
    #             progress$set(detail = "Done !", value = 0.95)
    #             
    #             # Export rowRanges as peaks
    #             refined_annotation = scExp_cf()@metadata$refined_annotation
    #             qs::qsave(refined_annotation, file = file.path(odir, "refined_annotation.qs"), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
    #             
    #             pc$new <- Sys.time()
    #             updateActionButton(session, "do_pc", label="Finished successfully", icon = icon("check-circle"))
    #         }
    #     }
    # })
    # 
    # pc <- reactiveValues(new="")
    # 
    # has_available_pc <- reactive({
    #     if(!is.null(scExp_cf())){
    #         if("refined_annotation" %in% names(scExp_cf()@metadata) ){
    #             return(TRUE)
    #         } else return(FALSE)
    #     } else{
    #         return(FALSE)
    #     }
    # })
    # 
    # observeEvent({selected_filtered_dataset()  # reset label on actionButtion when new peak calling should be performed
    #     input$pc_stat
    #     input$pc_stat_value}, {
    #         pc$available_pc <- has_available_pc()
    #         updateActionButton(session, "do_pc", label="Start", icon = character(0))
    #     })
    
    ###############################################################
    # 6. Differential analysis
    ###############################################################
    
    get.available.DA_GSEA.datasets <- function(name, preproc, numclust){
        
        l = list.files(path = file.path(init$data_folder, "ChromSCape_analyses", name, "Diff_Analysis_Gene_Sets"),
                       full.names = FALSE, recursive = FALSE, pattern = paste0(preproc, "_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+(.[[:digit:]]+)?_"))
        
        l = l[grep(paste0("orrected_", numclust),l)]
    }
    
    observeEvent(c(set_numclust(), input$qval.th, input$tabs, input$logFC.th, input$de_type, selected_filtered_dataset()), {
        if(input$tabs == "diff_analysis"){
            init$available_DA_GSA_datasets = get.available.DA_GSEA.datasets(analysis_name(), input$selected_reduced_dataset, set_numclust())
            
        }
    })
    
    DA_GSA_datasets <- reactive({
        if (is.null(init$available_DA_GSA_datasets)) c() else gsub('.{3}$', '', basename(init$available_DA_GSA_datasets)) })
    
    output$selected_DA_GSA_dataset <- renderUI({ 
        selectInput("selected_DA_GSA_dataset", "Select set with Differential Analysis:",
                    choices = DA_GSA_datasets()) 
    })
    
    observeEvent({
        input$selected_DA_GSA_dataset
    }, { # application header (tells you which data set is selected)
        req(analysis_name(), scExp_cf(), input$selected_DA_GSA_dataset)
        if(!is.null(input$selected_DA_GSA_dataset)){
            if(input$tabs %in% c("diff_analysis","enrich_analysis", "TF_analysis")){
                pattern = paste0(input$selected_reduced_dataset, "_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+(.[[:digit:]]+)?_")
                suffix = gsub(pattern,"", input$selected_DA_GSA_dataset)
                header <- paste0('<b>Analysis : ', analysis_name(), ' - ',input$feature_select, ' - k = ', length(unique(scExp_cf()$cell_cluster)), ' - ', suffix, ' </b>')
                shinyjs::html("pageHeader", header) 
            }
        }
    })
    
    observeEvent(input$selected_DA_GSA_dataset, { # load reduced data set to work with on next pages
        req(input$selected_DA_GSA_dataset)
        
        cat("Loading ", input$selected_DA_GSA_dataset)
        
        file_index <- match(c(input$selected_DA_GSA_dataset), DA_GSA_datasets())
        filename_sel <- file.path(init$data_folder, "ChromSCape_analyses",
                                  analysis_name(),"Diff_Analysis_Gene_Sets",
                                  init$available_DA_GSA_datasets[file_index])
        t1 = system.time({
            data = qs::qread(filename_sel, nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
            scExp_cf(NULL)
            if(input$feature_select %in% getExperimentNames(data$scExp_cf))
                scExp_cf. = swapAltExp_sameColData(data$scExp_cf,input$feature_select) else
                    scExp_cf. = data$scExp_cf
            
            scExp_cf(scExp_cf.)
            rm(data, scExp_cf.)
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
        req(set_numclust())
        HTML(paste0("<h5><b>Number of clusters selected  = ", set_numclust(),"</b></h5>"))
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
    # observeEvent(c(input$tabs, selected_filtered_dataset()),{
    #   if(input$tabs == "diff_analysis"){
    #     req(input$selected_DA_GSA_dataset)
    #     cat("WTF ", input$selected_DA_GSA_dataset, "\n")
    #     cat("WTF ", set_numclust(), "\n")
    #     if(!is.null(selected_filtered_dataset()) && !is.null(input$qval.th) && !is.null(input$logFC.th) &
    #        length(input$selected_DA_GSA_dataset) > 0){
    #       
    #       filename <- file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "Diff_Analysis_Gene_Sets",
    #                             paste0(input$selected_DA_GSA_dataset, ".qs"))
    #       cat(filename,"\n")
    #       if(file.exists(filename)){
    #         data = qs::qread(filename, nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
    #         
    #         if(input$feature_select %in% getExperimentNames(data$scExp_cf))
    #           scExp_cf. = swapAltExp_sameColData(data$scExp_cf,input$feature_select) else
    #             scExp_cf. = data$scExp_cf
    #         
    #         scExp_cf(NULL)
    #         scExp_cf(scExp_cf.)
    #         rm(data)
    #         gc()
    #       } else {
    #         NULL
    #       }
    #     }
    #   }
    #   })
    # 
    output$custom_da_ref <- renderUI({
        req(input$de_type)
        if(input$de_type == "custom"){
            shiny::selectInput("ref_type","Reference type", choices = comparable_variables(scExp_cf()))
        }
    })
    output$da_group <- renderUI({
            shiny::selectInput("group_type","Group type", choices = comparable_variables(scExp_cf()))
    })
    
    output$ref_choice <- renderUI({
        req(scExp_cf(), input$ref_type, input$de_type)
        if(input$de_type == "custom"){
                choices = unique(unlist(SingleCellExperiment::colData(scExp_cf())[[input$ref_type]]))
                selectInput("ref_choice","Reference sample(s)", choices= choices, 
                            multiple = TRUE)
            } else{NULL}
    })
    
    output$group_choice <- renderUI({
        req(scExp_cf(), input$group_type, input$de_type)
        if(input$de_type == "custom"){ 
            
                choices = unique(unlist(SingleCellExperiment::colData(scExp_cf())[[input$group_type]]))
                selectInput("group_choice","Reference sample(s)", choices= choices, 
                            multiple = TRUE)
        } else{
            NULL
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
        print("run_DA")
        if(input$de_type == "custom"){
          print("input$de_type == custom")
            if(TRUE %in% all.equal(input$group_choice, input$ref_choice)){
                showNotification(paste0("Warning : Please select different group and references."),
                                 duration = 7, closeButton = TRUE, type="warning")
                return(0)
            }
            
        }
        progress <- shiny::Progress$new(session, min=0, max=1)
        on.exit(progress$close())
        progress$set(message='Performing differential analysis...', value = 0.05)
        progress$set(detail='Initialization...', value = 0.15)
        
        if(batchUsed()) block = TRUE else block = FALSE
        gc()
        group = ""
        ref = ""
        by = input$group_type
        if(input$de_type == "custom"){
          print("input$de_type == custom")
            by = c(by, input$ref_type)
            group = data.frame(input$group_choice)
            ref = data.frame(input$ref_choice)
            colnames(group) = gsub("[^[:alnum:]|_]","",input$name_group)
            colnames(ref) = gsub("[^[:alnum:]|_]","",input$name_ref)
        }
        if(!is.null(scExp_cf()@metadata$enr)){
            scExp_cf. = scExp_cf()
            scExp_cf.@metadata$enr = NULL
            scExp_cf(scExp_cf.)
        } 
        scExp_cf(differential_analysis_scExp(scExp = scExp_cf(),
                                             by = by,
                                             method= input$da_method,
                                             de_type = input$de_type,
                                             block = block,
                                             group = group,
                                             ref = ref,
                                             progress = progress,
                                             BPPARAM = CS_options.BPPARAM())) 
        
        gc()
        
        scExp_cf. = scExp_cf()
        
        scExp_cf.@metadata$DA_parameters$logFC.th = input$logFC.th
        scExp_cf.@metadata$DA_parameters$qval.th = input$qval.th
        scExp_cf.@metadata$DA_parameters$min.percent = input$min.percent
        
        scExp_cf(scExp_cf.)
        data = list("scExp_cf" = getMainExperiment(scExp_cf()))
        
        DA_GSA_suffix = scExp_cf()@metadata$DA_parameters$de_type
        if(DA_GSA_suffix == "custom") DA_GSA_suffix = paste0(gsub("[^[:alnum:]|_]","",input$name_group),"_vs_",
                                                             gsub("[^[:alnum:]|_]","",input$name_ref))
        
        suffix = paste0(selected_filtered_dataset(), "_", length(unique(scExp_cf()$cell_cluster)),
                        "_", input$qval.th, "_", input$logFC.th, "_", DA_GSA_suffix, ".qs")
        
        qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "Diff_Analysis_Gene_Sets",
                                         suffix), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
        rm(data)
        gc()
        init$available_DA_GSA_datasets = get.available.DA_GSEA.datasets(analysis_name(), input$selected_reduced_dataset, set_numclust())
        
        updateSelectInput(session = session, inputId = "selected_DA_GSA_dataset",
                          label =  "Select set with Differential Analysis:",
                          choices = DA_GSA_datasets(),
                          selected =  gsub(".qs$","",suffix))
        
        updateActionButton(session = session, inputId = "run_DA", label = "Start analysis ", icon = icon("check-circle"))
        updateActionButton(session = session, inputId = "apply_DA_filters", label = "Apply filters")
        # incProgress(amount = 0.2, detail = paste("Saving DA"))
        # })
    })
    
    observeEvent(input$apply_DA_filters,
                 {
                     
                     # Add Differential Analysis parameter to the DA 
                     scExp_cf. = scExp_cf()
                     
                     scExp_cf.@metadata$DA_parameters$logFC.th = input$logFC.th
                     scExp_cf.@metadata$DA_parameters$qval.th = input$qval.th
                     scExp_cf.@metadata$DA_parameters$min.percent = input$min.percent
                     
                     scExp_cf(scExp_cf.)
                     data = list("scExp_cf" = getMainExperiment(scExp_cf()))
                     
                     DA_GSA_suffix = scExp_cf()@metadata$DA_parameters$de_type
                     if(DA_GSA_suffix == "custom") DA_GSA_suffix = paste0(gsub("[^[:alnum:]|_]","",input$name_group),"_vs_",
                                                                          gsub("[^[:alnum:]|_]","",input$name_ref))
                     
                     suffix = paste0(selected_filtered_dataset(), "_", length(unique(scExp_cf()$cell_cluster)),
                                     "_", input$qval.th, "_", input$logFC.th, "_", DA_GSA_suffix, ".qs")
                     
                     qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "Diff_Analysis_Gene_Sets",
                                                      suffix), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
                     rm(data)
                     gc()
                     init$available_DA_GSA_datasets = get.available.DA_GSEA.datasets(analysis_name(), input$selected_reduced_dataset, set_numclust())
                     
                     updateSelectInput(session = session, inputId = "selected_DA_GSA_dataset",
                                       label =  "Select set with Differential Analysis:",
                                       choices = DA_GSA_datasets(),
                                       selected =  gsub(".qs$","",suffix))
                     updateActionButton(session = session, inputId = "apply_DA_filters", label = "Saved ", icon = icon("check-circle"))
                     
                 })
    
    observeEvent({input$logFC.th; input$qval.th; input$min.percent},{
        updateActionButton(session = session, inputId = "apply_DA_filters", label = "Apply filters ")
        
    })
    
    output$da_summary_box <- renderUI({
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$DA_parameters)){
                input$logFC.th; input$qval.th; input$min.percent
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
            if(!is.null(scExp_cf()@metadata$DA_parameters)){
                summary_DA(scExp_cf(), qval.th = input$qval.th, 
                           logFC.th = input$logFC.th, min.percent = input$min.percent) %>%
                    kableExtra::kable(escape = FALSE, align="c") %>%
                    kableExtra::kable_styling(c("striped", "condensed"), full_width = FALSE)
            }
        }
    }
    
    output$da_barplot <- renderPlot({
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$DA_parameters)){
                plot_differential_summary_scExp(scExp_cf(), qval.th = input$qval.th, 
                                                logFC.th = input$logFC.th, min.percent = input$min.percent)
            }
        }
    })
    
    output$download_da_barplot <- downloadHandler(
        filename = function(){ paste0("diffAnalysis_numRegions_barplot_", selected_filtered_dataset(), "_", length(unique(scExp_cf()$cell_cluster)), "_", input$qval.th, "_", input$logFC.th, "_", input$de_type, ".png")},
        content = function(file){
            grDevices::png(file, width = 800, height = 600, res = 150)
            plot_differential_summary_scExp(scExp_cf(), qval.th = input$qval.th, 
                                            logFC.th = input$logFC.th, min.percent = input$min.percent)
            grDevices::dev.off()
        })
    
    output$da_table <- DT::renderDataTable({
        req(input$gpsamp, scExp_cf(), input$logFC.th, input$qval.th, input$min.percent)
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$DA_parameters)){
                rowD = as.data.frame(SingleCellExperiment::rowData(scExp_cf()))
                if(any(grepl(input$gpsamp,colnames(rowD))) ){
                    diff = rowD[,-grep("ID",colnames(rowD))]
                    if(!is.null(SummarizedExperiment::rowRanges(scExp_cf())$Gene)){
                        rowD = SummarizedExperiment::rowRanges(scExp_cf())
                        rowD = rowD[which(SummarizedExperiment::rowData(scExp_cf())$top_feature)]
                    }
                    rownames(diff) <- NULL
                    
                    diff = diff[,c("Gene", 
                                   colnames(diff)[grep(input$gpsamp, colnames(diff))], "distanceToTSS")]
                    diff <- diff[order(diff[,paste0("logFC.",input$gpsamp)]),]
                    diff[,paste0("logFC.",input$gpsamp)] = round(diff[,paste0("logFC.",input$gpsamp)],3)
                    diff[,paste0("group_activation.",input$gpsamp)] = round(diff[,paste0("group_activation.",input$gpsamp)],3)
                    diff[,paste0("reference_activation.",input$gpsamp)] = round(diff[,paste0("reference_activation.",input$gpsamp)],3)
                    diff_up = diff %>% 
                        filter(
                            .data[[paste0("logFC.",input$gpsamp)]] >  input$logFC.th &
                                .data[[paste0("qval.",input$gpsamp)]]  <  input$qval.th &
                                .data[[paste0("group_activation.",input$gpsamp)]] > input$min.percent
                        )
                    diff_down = diff %>% 
                        filter(
                            .data[[paste0("logFC.",input$gpsamp)]] <  -input$logFC.th &
                                .data[[paste0("qval.",input$gpsamp)]]  <  input$qval.th &
                                .data[[paste0("reference_activation.",input$gpsamp)]] > input$min.percent
                        )
                    diff = rbind(diff_down, diff_up)
                    DT::datatable(diff, options = list(dom='tpi'), class = "display",
                                  rownames = FALSE)
                }
            }
        }
    })
    
    output$download_da_table <- downloadHandler(
        filename = function(){ paste0("diffAnalysis_data_", selected_filtered_dataset(),
                                      "_", length(unique(scExp_cf()$cell_cluster)), "_",
                                      input$qval.th, "_", input$logFC.th, "_",
                                      input$de_type, ".tsv")},
        content = function(file){
            write.table(SingleCellExperiment::rowData(scExp_cf()), file, row.names = FALSE, quote = FALSE, sep="\t")
        })
    
    output$da_visu_box <- renderUI({
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$DA_parameters)){
                shinydashboard::box(title="Detailed differential analysis per cluster", width = NULL, status="success", solidHeader = TRUE,
                                    column(4, align="left", selectInput("gpsamp", "Select group:", choices = DA_groups())),
                                    column(4, align="left", downloadButton("download_da_table", "Download table")),
                                    column(12, align="left", div(style = 'overflow-x: scroll', DT::dataTableOutput('da_table')), br()),
                                    column(12, align="left", plotOutput("da_volcano", height = "500px", width = "100%")),
                                    column(4, align="left", actionButton("save_plots_DA", "Save HQ plots", icon = icon("fas fa-image"))
                                    ))
            }
        }
    })
    
    observeEvent(input$save_plots_DA,{
        if(!dir.exists(plot_dir())) dir.create(plot_dir())
        pdf(file.path(plot_dir(), paste0("differential_analysis_", input$gpsamp,".pdf")))
        plot_differential_volcano_scExp(scExp_cf(), group = input$gpsamp,
                                        logFC.th = input$logFC.th, qval.th = input$qval.th,
                                        min.percent = input$min.percent)
        dev.off()
        shiny::showNotification(paste0("Plots saved in '",plot_dir(),"' !"),
                                duration = 15, closeButton = TRUE, type="message")
        updateActionButton(session = session, inputId = "save_plots_DA", label = "Save HQ plots", icon = icon("check-circle"))
    })
    
    output$da_volcano <- renderPlot({
        req(scExp_cf())
        if(!is.null(scExp_cf())){
            if(input$gpsamp %in% DA_groups()){
                if(!is.null(scExp_cf()@metadata$DA_parameters)){
                    plot_differential_volcano_scExp(scExp_cf(), group = input$gpsamp,
                                                    logFC.th = input$logFC.th, 
                                                    qval.th = input$qval.th, 
                                                    min.percent = input$min.percent)
                }
            }
        }
    })
    
    output$download_da_volcano <- downloadHandler(
        filename = function(){ paste0("diffAnalysis_numRegions_barplot_", selected_filtered_dataset(),
                                      "_", length(unique(scExp_cf()$cell_cluster)), "_",
                                      input$qval.th, "_", input$logFC.th, "_", 
                                      input$de_type, "_", input$gpsamp, ".png")},
        content = function(file){
            grDevices::png(file, width = 900, height = 900, res = 300)
            plot_differential_volcano_scExp(scExp_cf(),
                                            group = input$gpsamp,
                                            logFC.th = input$logFC.th,
                                            qval.th = input$qval.th,
                                            min.percent = input$min.percent)
            grDevices::dev.off()
        })
    
    DA_groups <- reactive({
        gsub(".*\\.","", grep("qval",colnames(SingleCellExperiment::rowData(scExp_cf())), value = TRUE))
    })
    
    observeEvent(unlocked$list, {
        able_disable_tab(c("diff_my_res"),"enrich_analysis")
        able_disable_tab(c("diff_my_res"),"TF_analysis")
    }) 
    
    observeEvent(scExp_cf(), {
        if(any(grepl("qval",colnames(SingleCellExperiment::rowData(scExp_cf()))))){
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
                eval(parse(text = paste0("myData$",annotation_id(),".GeneTSS$Gene")))
            ))
        }
    })
    
    observeEvent(input$do_enrich, {
        
        progress <- shiny::Progress$new(session, min=0, max=1)
        on.exit(progress$close())
        progress$set(message='Performing Gene Set Analysis...', value = 0.05)
        progress$set(detail='Initialization...', value = 0.15)
        scExp_cf(gene_set_enrichment_analysis_scExp(scExp_cf(), enrichment_qval = 0.01, qval.th = input$qval.th,
                                                    ref = annotation_id(), logFC.th = input$logFC.th,
                                                    min.percent = input$min.percent,
                                                    peak_distance = 1000, use_peaks = input$use_peaks,
                                                    GeneSetClasses = MSIG.classes(), progress = progress))
        gc()
        progress$inc(detail='Finished Gene Set Analysis - Saving...', amount = 0.0)
        data = list("scExp_cf" = getMainExperiment(scExp_cf()) )
        
        qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "Diff_Analysis_Gene_Sets",
                                         paste0(input$selected_DA_GSA_dataset, ".qs")), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
        rm(data)
        gc()
        progress$inc(detail='Done !', amount = 0.5)
        
    })
    
    output$GSA_group_sel <- renderUI({ 
        if(!is.null(scExp_cf())){
            if(!is.null(DA_groups())){
                selectInput("GSA_group_sel", "Select group:",
                            choices = DA_groups()) 
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
            if(!is.null(scExp_cf()@metadata$enr) && !is.null(input$GSA_group_sel) && !is.null(input$enr_class_sel) &&
               !is.null(input$pathway_size_GSA)){
                if(input$GSA_group_sel %in% DA_groups()){
                    if(length(scExp_cf()@metadata$enr$Both)>0){
                        DT::datatable(table_enriched_genes_scExp(
                            scExp_cf(), set = 'Both',  group = input$GSA_group_sel,
                            enr_class_sel = input$enr_class_sel) %>%
                                dplyr::filter(as.numeric(gsub(".*/","",Num_deregulated_genes)) >
                                                  input$pathway_size_GSA[1] &
                                                  as.numeric(gsub(".*/","",Num_deregulated_genes)) < 
                                                  input$pathway_size_GSA[2]),
                            options = list(pageLength = 10,  dom = 'tpi'), class = 'display', rownames = FALSE)
                    }
                }
            }
        }
    })
    
    output$over_enrich_table <- DT::renderDataTable({
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$enr)){
                if(input$GSA_group_sel %in% DA_groups()){
                    if(length(scExp_cf()@metadata$enr$Overexpressed)>0){
                        DT::datatable(table_enriched_genes_scExp(
                            scExp_cf(), set = 'Overexpressed',  group = input$GSA_group_sel,
                            enr_class_sel = input$enr_class_sel) %>%
                                dplyr::filter(as.numeric(gsub(".*/","",Num_deregulated_genes)) >
                                                  input$pathway_size_GSA[1] &
                                                  as.numeric(gsub(".*/","",Num_deregulated_genes)) < 
                                                  input$pathway_size_GSA[2]),
                            options = list(pageLength = 10,  dom = 'tpi'), class = 'display', rownames = FALSE)
                    }
                }
                
            }
        }
    })
    
    output$under_enrich_table <- DT::renderDataTable({
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$enr)){
                if(input$GSA_group_sel %in% DA_groups()){
                    if(length(scExp_cf()@metadata$enr$Underexpressed)>0){
                        DT::datatable(table_enriched_genes_scExp(
                            scExp_cf(), set = 'Underexpressed',  group = input$GSA_group_sel,
                            enr_class_sel = input$enr_class_sel) %>%
                                dplyr::filter(as.numeric(gsub(".*/","",Num_deregulated_genes)) >
                                                  input$pathway_size_GSA[1] &
                                                  as.numeric(gsub(".*/","",Num_deregulated_genes)) < 
                                                  input$pathway_size_GSA[2]),
                            options = list(pageLength = 10,  dom = 'tpi'), class = 'display', rownames = FALSE)
                    }
                }
                
            }
        }
    })
    
    output$download_enr_data <- downloadHandler(
        filename = paste(selected_filtered_dataset(), length(unique(scExp_cf()$cell_cluster)), input$qval.th, input$logFC.th, "enrichment_tables.zip", sep="_"),
        content = function(fname){
            fis <- c()
            for(i in 1:length(DA_groups())){
                if(!is.null(scExp_cf()@metadata$enr$Both[[i]])){
                    filename <- paste0(DA_groups()[i], "_significant_gene_sets.csv")
                    fis <- c(fis, filename)
                    write.table(scExp_cf()@metadata$enr$Both[[i]], file = filename, quote = FALSE, row.names = FALSE, sep=",")
                }
                if(!is.null(scExp_cf()@metadata$enr$Overexpressed[[i]])){
                    filename <- paste0(DA_groups()[i], "_enriched_gene_sets.csv")
                    fis <- c(fis, filename)
                    write.table(scExp_cf()@metadata$enr$Overexpressed[[i]], file = filename, quote = FALSE, row.names = FALSE, sep=",")
                }
                if(!is.null(scExp_cf()@metadata$enr$Underexpressed[[i]])){
                    filename <- paste0(DA_groups()[i], "_depleted_gene_sets.csv")
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
                s = selectizeInput(inputId = "gene_sel", label = "Select gene:", choices = NULL)
                most_diff = as.data.frame(SingleCellExperiment::rowData(scExp_cf())) %>% dplyr::select(ID, starts_with("qval."))
                most_diff[,"qval"] = Matrix::rowMeans(as.matrix(most_diff[,-1]))
                most_diff = dplyr::left_join(most_diff[order(most_diff$qval),], annotFeat_long(),by = c("ID"))
                most_diff = most_diff %>% dplyr::filter(!is.na(Gene)) %>%
                    dplyr::filter(distanceToTSS < 1000) 
                
                genes = base::intersect(most_diff$Gene,unique(GencodeGenes()))
                
                updateSelectizeInput(session = session, inputId = "gene_sel",
                                     label =  "Select gene:", choices = genes, server = TRUE)
                return(s)
            }
        }
    })
    
    observeEvent(input$save_plot_GSA,{
        req(gene_umap_p(), input$gene_sel, gene_violin_p())
        if(!dir.exists(plot_dir())) dir.create(plot_dir())
        pdf(file.path(plot_dir(), paste0("UMAP_", input$gene_sel,".pdf")))
        print(gene_umap_p())
        print(gene_violin_p())
        if(!is.null(pathways_mat())) print(pathways_umap_p())
        dev.off()
        shiny::showNotification(paste0("Plots saved in '",plot_dir(),"' !"),
                                duration = 15, closeButton = TRUE, type="message")
        updateActionButton(session = session, inputId = "save_plot_GSA", label = "Save HQ plots", icon = icon("check-circle"))
    })
    
    output$pathways_sel <- renderUI({ 
        shiny::selectizeInput("pathways_sel", label = "Select pathways:", choices = NULL)
    })
    
    pathways <- reactive({
        req(scExp_cf())
        if(!is.null(scExp_cf()@metadata$enr)){
            all_paths <- unique(as.character(unlist(sapply(1:length(DA_groups()), 
                                                           function(i) {
                                                               all = c()
                                                               if(length(scExp_cf()@metadata$enr$Overexpressed)>=i) 
                                                                   all = c(all, scExp_cf()@metadata$enr$Overexpressed[[i]]$Gene.Set)
                                                               if(length(scExp_cf()@metadata$enr$Underexpressed)>=i) 
                                                                   all = c(all, scExp_cf()@metadata$enr$Underexpressed[[i]]$Gene.Set)
                                                               if(length(scExp_cf()@metadata$enr$Both)>=i) 
                                                                   all = c(all, scExp_cf()@metadata$enr$Both[[i]]$Gene.Set)
                                                               return(all)
                                                           }))))
            shiny::updateSelectizeInput(session =  session, inputId = "pathways_sel",
                                        label = "Select pathways:", choices = all_paths, server = TRUE)
            
            all_paths
        }
    })
    
    pathways_mat <- reactiveVal(NULL)
    observeEvent(input$plot_pathways, {
        req(pathways(), annotation_id(), MSIG.classes())
        if(!is.null(scExp_cf()@metadata$enr)){
            progress <- shiny::Progress$new(session, min=0, max=1)
            on.exit(progress$close())
            progress$set(message='Loading Pathways for visualization...', value = 0.05)
            progress$set(detail='Initialization...', value = 0.15)
            normcounts = SingleCellExperiment::normcounts(scExp_cf())
            regions = SummarizedExperiment::rowRanges(scExp_cf())
            regions = regions[which(regions$distanceToTSS <= 1000)]
            database <- ChromSCape:::load_MSIGdb(annotation_id(), MSIG.classes())
            progress$set(message='Finished loading...', value = 0.65)
            progress$set(detail='Creating Pathways x Region mat...', value = 0.70)
            pathways_to_regions = sapply(seq_along(pathways()),
                                         function(i){
                                             genes = database$GeneSets[[pathways()[i]]]
                                             reg = unique(regions$ID[grep(paste(genes, collapse="|"),regions$Gene)])
                                         })
            names(pathways_to_regions) = pathways()
            pathways_mat. =  sapply(seq_along(pathways()),
                                    function(i){
                                        if(length(pathways_to_regions[[i]])>1) {
                                            r = colSums(normcounts[pathways_to_regions[[i]],])
                                        } else{
                                            r = rep(0, ncol(normcounts))
                                        }
                                        r = as.matrix(r)
                                        colnames(r) = pathways()[i]
                                        r
                                    })
            colnames(pathways_mat.) = pathways()
            rownames(pathways_mat.) = colnames(normcounts)
            progress$set(detail='Done !', value = 0.90)
            pathways_mat(pathways_mat.)
        }
    })
    
    pathways_umap_p <- reactive({
        req(pathways_mat())
        if(input$pathways_sel %in% colnames(pathways_mat())){
            cols = rev(viridis::magma(100))[round(changeRange(pathways_mat()[,input$pathways_sel],
                                                              newmin = 1, newmax = 100))]
            p <- ggplot(as.data.frame(SingleCellExperiment::reducedDim(scExp_cf(), "UMAP")),
                        aes(x = Component_1, y = Component_2)) +
                geom_point(alpha = 0.3,  color = cols, aes(shape = SummarizedExperiment::colData(scExp_cf())$cell_cluster)) +
                labs(color="Sum norm. count for region", title = "UMAP", shape="Cluster", x="Component 1", y="Component 2") +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour="black"),
                      panel.border = element_rect(colour="black", fill = NA))
        }
    })
    
    output$pathways_umap_plot <- plotly::renderPlotly({
        req(pathways_umap_p())
        plotly::ggplotly(pathways_umap_p(), tooltip="Sample", dynamicTicks = TRUE)
    })
    
    output$pathways_umap_UI <- renderUI({
        req(pathways_umap_p())
        plotly::plotlyOutput("pathways_umap_plot")
    })
    
    gene_umap_p <- reactive({
        req(input$gene_sel)
        p <- plot_reduced_dim_scExp(scExp_cf(),
                                    color_by = input$gene_sel,
                                    reduced_dim = "UMAP",
                                    transparency = input$options.dotplot_transparency,
                                    size = input$options.dotplot_size,
                                    max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                    downsample = input$options.dotplot_downsampling,
                                    min_quantile = input$options.dotplot_min_quantile,
                                    max_quantile = input$options.dotplot_max_quantile,
                                    annotate_clusters = input$label_cluster_umap_GSA
        )
        p
    })
    
    output$gene_umap_UI <- renderUI({
        req(input$gene_sel)
        
        output$gene_umap_plot <- renderPlot({
            gene_umap_p()
        })
        shiny::plotOutput("gene_umap_plot") %>% 
            shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)
        
    })
    
    gene_violin_p <- reactive({
        req(input$gene_sel, input$color_by_violin_GSA)
        p <- plot_percent_active_feature_scExp(scExp_cf(),
                                       gene = input$gene_sel,
                                       by = input$color_by_violin_GSA, 
                                       max_distanceToTSS = input$options.dotplot_max_distanceToTSS,
                                       downsample = input$options.dotplot_downsampling
        )
        p
    })
    output$gene_violin_UI <- renderUI({
        req(input$gene_sel)
        
        output$gene_violin_plot <- renderPlot({
            gene_violin_p()
        })
        shiny::plotOutput("gene_violin_plot") %>% 
            shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)
        
    })
    
    ###############################################################
    # 8. TF Enrichment using ChEA3
    ###############################################################
    
    output$use_peaks <- renderUI({
        if(canUsePeaks()){ checkboxInput("use_peaks", "use peak calling results to refine analysis", value = TRUE) }
    })
    
    output$TF_group_sel <- renderUI({ 
        if(!is.null(scExp_cf())){
            if(!is.null(DA_groups())){
                selectInput("TF_group_sel", "Select group:",
                            choices = DA_groups()) 
            }
        }
    })
    
    observeEvent(input$do_enrich_TF, {
        
        progress <- shiny::Progress$new(session, min=0, max=1)
        on.exit(progress$close())
        progress$set(message='Performing TF Analysis...', value = 0.05)
        progress$set(detail='Initialization...', value = 0.15)
        scExp_cf(enrich_TF_ChEA3_scExp(scExp_cf(), qval.th = input$qval.th,
                                       ref = annotation_id(), logFC.th = input$logFC.th,
                                       min.percent = input$min.percent,
                                       peak_distance = 1000, use_peaks = input$use_peaks,
                                       progress = progress))
        gc()
        progress$inc(detail='Finished TF Analysis - Saving...', amount = 0.0)
        data = list("scExp_cf" = getMainExperiment(scExp_cf()) )
        
        qs::qsave(data, file = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), "Diff_Analysis_Gene_Sets",
                                         paste0(input$selected_DA_GSA_dataset, ".qs")), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
        rm(data)
        gc()
        progress$inc(detail='Done !', amount = 0.1)
        
    })
    
    output$all_enrich_table_TF <- DT::renderDataTable({
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$TF_enrichment) && !is.null(input$TF_group_sel)){
                if(input$TF_group_sel %in% DA_groups()){
                    if(length(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]])>0){
                        if(!is.null(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]]$Differential)){
                            DT::datatable(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]]$Differential,
                                          options = list(pageLength = 10,  dom = 'tpi'), class = 'display', rownames = FALSE)
                        }
                    }
                }
            }
        }
    })
    
    output$over_enrich_table_TF <- DT::renderDataTable({
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$TF_enrichment) && !is.null(input$TF_group_sel)){
                if(input$TF_group_sel %in% DA_groups()){
                    if(length(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]])>0){
                        if(!is.null(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]]$Enriched)){
                            DT::datatable(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]]$Enriched,
                                          options = list(pageLength = 10,  dom = 'tpi'), class = 'display', rownames = FALSE)
                        }
                    }
                }
            }
        }
    })
    
    output$under_enrich_table_TF <- DT::renderDataTable({
        if(!is.null(scExp_cf())){
            if(!is.null(scExp_cf()@metadata$TF_enrichment) && !is.null(input$TF_group_sel)){
                if(input$TF_group_sel %in% DA_groups()){
                    if(length(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]])>0){
                        if(!is.null(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]]$Depleted)){
                            DT::datatable(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]]$Depleted,
                                          options = list(pageLength = 10,  dom = 'tpi'), class = 'display', rownames = FALSE)
                        }
                    }
                }
            }
        }
    })
    
    output$download_enr_data_TF <- downloadHandler(
        filename = paste(selected_filtered_dataset(), length(unique(scExp_cf()$cell_cluster)), input$qval.th, input$logFC.th, "TF_enrichment_tables.zip", sep="_"),
        content = function(fname){
            fis <- c()
            for(i in 1:length(DA_groups())){
                if(!is.null(scExp_cf()@metadata$TF_enrichment[[i]]$Differential)){
                    filename <- paste0(DA_groups()[i], "_significant_gene_sets.csv")
                    fis <- c(fis, filename)
                    write.table(scExp_cf()@metadata$TF_enrichment[[i]]$Differential, file = filename, quote = FALSE, row.names = FALSE, sep=",")
                }
                if(!is.null(scExp_cf()@metadata$TF_enrichment[[i]]$Enriched)){
                    filename <- paste0(DA_groups()[i], "_enriched_gene_sets.csv")
                    fis <- c(fis, filename)
                    write.table(scExp_cf()@metadata$TF_enrichment[[i]]$Enriched, file = filename, quote = FALSE, row.names = FALSE, sep=",")
                }
                if(!is.null(scExp_cf()@metadata$TF_enrichment[[i]]$Depleted)){
                    filename <- paste0(DA_groups()[i], "_depleted_gene_sets.csv")
                    fis <- c(fis, filename)
                    write.table(scExp_cf()@metadata$TF_enrichment[[i]]$Depleted, file = filename, quote = FALSE, row.names = FALSE, sep=",")
                }
            }
            zip(zipfile = fname, files = fis)},
        contentType = "application/zip"
    )
    
    TF_barplot_p <- reactive({
        req(scExp_cf(), input$TF_set, input$TF_group_sel, input$TF_set, input$TF_plot_y_axis,
            input$TF_plot_n_top)
        if(!is.null(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]][[input$TF_set]])){
            if(length(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]][[input$TF_set]]) > 1){
                if(nrow(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]][[input$TF_set]]) > 1){
                    p <- plot_top_TF_scExp(scExp_cf(),
                                           group = input$TF_group_sel,
                                           set = input$TF_set,
                                           type = input$TF_plot_y_axis,
                                           n_top = input$TF_plot_n_top)
                    p
                } else{
                    NULL
                }
            } else{
                NULL
            }
        } else{
            NULL
        }
        
        
    })
    
    output$TF_plot <- renderPlot({
        TF_barplot_p()
    })
    
    shiny::plotOutput("TF_plot") %>% 
        shinycssloaders::withSpinner(type=8,color="#434C5E",size = 0.75)
    
    output$TF_results_UI = renderUI({
        req(scExp_cf(), input$TF_group_sel, input$TF_group_sel, DA_groups())
        if(!is.null(scExp_cf()@metadata$TF_enrichment)){
            if(input$TF_group_sel %in% DA_groups()){
                if(length(scExp_cf()@metadata$TF_enrichment[[input$TF_group_sel]]) > 0){
                    shinydashboard::box(title="Results TF analysis per cluster", width = NULL, status="success", solidHeader = TRUE,
                                        column(4, align="left", selectInput("TF_set", "Select set:", choices = c("Differential", "Enriched", "Depleted"))),
                                        column(4, align="left", selectInput("TF_plot_y_axis", "Select Y axis:", choices = c("Score", "nTargets", "nTargets_over_TF", "nTargets_over_genes"))),
                                        column(4, align="left", sliderInput("TF_plot_n_top", "Select set:",min = 1, max = 200, value = 25), br()),
                                        column(12, align="left", plotOutput("TF_plot", height = "500px", width = "100%"))
                    )
                }
            }
        }
    })
    
    ###############################################################
    # 9. Close app
    ###############################################################
    
    output$analysis_deletion_info <- renderText({"The selected data set will be fully deleted from the computer, including all reduced data versions that have been produced so far for this set."})
    output$selected_delete_analysis <- renderUI({ selectInput("selected_delete_analysis", "Select data set :", choices = init$available_analyses) })
    
    observeEvent(input$delete_analysis, {  # delete selected dataset
        withProgress(message='Deleting data set', value = 0, {
            incProgress(amount=0.5, detail=paste("..."))
            unlink(file.path(init$data_folder, "ChromSCape_analyses", input$selected_delete_analysis), recursive=TRUE)
            init$available_analyses <- list.dirs(path=file.path(init$data_folder, "ChromSCape_analyses"), full.names=FALSE, recursive=FALSE)
            init$available_reduced_datasets <- get.available.reduced.datasets(analysis_name())
            incProgress(amount=0.5, detail=paste("... finished"))
        })
        showNotification("Data set successfully deleted.", duration=5, closeButton=TRUE, type="warning")
    })
    
    observeEvent(input$generate_report, {  # delete selected dataset
        withProgress(message='Creating HTML report', value = 0, {
            incProgress(amount=0.5, detail=paste("..."))
            print(file.path(init$data_folder, "ChromSCape_analyses", analysis_name()))
            generate_report(ChromSCape_directory = file.path(init$data_folder, "ChromSCape_analyses", analysis_name()))
            incProgress(amount=0.5, detail=paste("... finished"))
        })
        showNotification("Data set successfully deleted.", duration=5, closeButton=TRUE, type="warning")
    })
    
    observeEvent(input$close_and_save, {
        unlink(file.path("www", "images", "*"))
        unlink(file.path(".", "*.csv"))
        
        if(!is.null(scExp()) & !is.null(input$selected_reduced_dataset)){
            scExp = isolate(getMainExperiment(scExp()))
            qs::qsave(scExp, file = file.path(init$data_folder, "ChromSCape_analyses",
                                              analysis_name(), "Filtering_Normalize_Reduce",
                                              paste0(input$selected_reduced_dataset,".qs")), nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
            
        }
        if(!is.null(scExp_cf()) & !is.null(selected_filtered_dataset())){
            data = list()
            data$scExp_cf = isolate(getMainExperiment(scExp_cf()))
            
            if("consclust" %in% names(scExp_cf()@metadata)){
                if("diff" %in% names(scExp_cf()@metadata)){
                    dir =file.path(init$data_folder, "ChromSCape_analyses", analysis_name(), 
                                   "Diff_Analysis_Gene_Sets", paste0(
                                       selected_filtered_dataset(), "_",
                                       length(unique(scExp_cf()$cell_cluster)),
                                       "_", input$qval.th, "_", input$logFC.th, "_",
                                       input$de_type, ".qs"))
                } else{
                    dir = file.path(init$data_folder, "ChromSCape_analyses", analysis_name(),
                                    "correlation_clustering",
                                    paste0(selected_filtered_dataset(), 
                                           ".qs"))
                }
                scExp = getMainExperiment(scExp)
                qs::qsave(scExp, file = dir, nthreads = as.numeric(BiocParallel::bpworkers(CS_options.BPPARAM())))
            } 
            
        }
        
        lapply(names(resourcePaths()), removeResourcePath)
        
        print("Thank you & see you next time !")
        stopApp()
    })
    
})
