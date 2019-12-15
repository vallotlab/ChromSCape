
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=5000*1024^2) # allow upload of files with max 5GB
  
  ###############################################################
  # 0. Global variables and functions
  ###############################################################
  
  #Initializating user experience functions 
  js$init_directory() #Getting cookie for the directory
  
  tab_vector=c("pca_plots","cor_clustering","cons_clustering","peak_calling","diff_analysis","enrich_analysis") #list of all lockable tabs
  unlocked = reactiveValues(list=list(selected_reduced_dataset=FALSE,pca=FALSE,tsne=FALSE,cor_clust_plot=FALSE,filtered_datasets=FALSE,affectation=FALSE,diff_my_res=FALSE)) #list of all required items to unlock a tab
  for(tab in tab_vector){
    js$disableTab(tab); #Disabling all tabs but the first one
  }
  #Global reactives values
  cell_cov_df <- reactive ({data.frame(coverage=sort(unname(colSums(init$datamatrix)))) })  # used for plotting cell coverage on first page
  dataset_name <- reactive({ sub("_\\d+_\\d+(\\.\\d+)?_\\d+_[A-z]+$", "", input$selected_reduced_dataset) })
  annotation_id_norm <- reactive({ read.table(file.path(init$data_folder, 'datasets', input$selected_raw_dataset, 'annotation.txt'), header=FALSE, stringsAsFactors=FALSE)[[1]] })
  annotation_id <- reactive({ read.table(file.path(init$data_folder, 'datasets', dataset_name(), 'annotation.txt'), header=FALSE, stringsAsFactors=FALSE)[[1]] })
  
  #Global Functions
  init <- reactiveValues(data_folder="/var/lib/shiny-server/", datamatrix=data.frame(), annot_raw=data.frame(), available_raw_datasets = NULL,
                         available_reduced_datasets=NULL, available_filtered_datasets=NULL)
  reduced_datasets <- reactive({ if (is.null(init$available_reduced_datasets)) c() else gsub('.{6}$', '', basename(init$available_reduced_datasets)) })
  annotCol <- reactive({ c("sample_id","total_counts") })
  
  
  observeEvent(dataset_name(), { # application header (tells you which data set is selected)
    req(dataset_name())
    header <- paste0('Sample : ', dataset_name())
    shinyjs::html("pageHeader", header)
  })
  
  get.available.reduced.datasets <- function(){
    list.files(path=file.path(init$data_folder, "datasets"), full.names=FALSE, recursive=TRUE, pattern="[[:print:]]+_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?_[[:digit:]]+_(uncorrected).RData")
  }
  get.available.filtered.datasets <- function(name, preproc){
    list.files(path=file.path(init$data_folder, "datasets", name, "cor_filtered_data"), full.names=FALSE, recursive=FALSE, pattern=paste0(preproc, "_[[:digit:]]+_[[:digit:]]+(.[[:digit:]]+)?.RData"))
  }
  observe({
    init$available_raw_datasets <- list.dirs(path=file.path(init$data_folder, "datasets"), full.names=FALSE, recursive=FALSE)
    init$available_reduced_datasets <- get.available.reduced.datasets()
  })
  able_disable_tab <- function(variables_to_check, tab_id) {
      able_or_disable=c()
      for(var in variables_to_check){
        if (unlocked$list[[var]]==T) {
          able_or_disable=c(able_or_disable,TRUE)
        } else{
          able_or_disable=c(able_or_disable,FALSE)
        }}
        for (tab in tab_id) {
          if (FALSE %in% able_or_disable) {
            js$disableTab(tab_id)
          }
          else{
            js$enableTab(tab)
          }}
  }
  
  
  ###############################################################
  # 1. Select or upload dataset
  ###############################################################
  output$selected_raw_dataset <- renderUI({ selectInput("selected_raw_dataset", "Select data set :", choices=init$available_raw_datasets) })
  output$data_folder_info <- renderText({"Please select the directory which contains the 'dataset' folder that was built by the app. If this folder does not yet exist, the app will create it in the specified directory."})
  output$selected_reduced_dataset <- renderUI({ selectInput("selected_reduced_dataset", "Select filtered & normalized set :", choices=reduced_datasets()) })
  output$red_data_selection_info <- renderText({"The selected data set is automatically loaded and will be used for all subsequent analysis. If you just filtered a new data set, don't forget to select it here."})
  output$red_data_selection_format <- renderText({"The name of the preprocessed dataset is composed of the following information: data set name, min percentage of reads per cell, min percentage of cells to support a window, quantile of cell read counts to keep"})
  output$data_matrices_info <- renderText({"The filename of each selected matrix must be the sample id, e.g. T23_K4me3.txt (Must contain only alpha-numeric charachter and underscore !)"})
  output$data_deletion_info <- renderText({"The selected data set will be fully deleted from the computer, including all reduced data versions that have been produced so far for this set."})
  output$selected_delete_dataset <- renderUI({ selectInput("selected_delete_dataset", "Select data set :", choices=init$available_raw_datasets) })
  output$exclude_file <- renderUI({ if(input$exclude_regions){
    fileInput("exclude_file", ".bed file containing the regions to exclude from data set:", multiple=FALSE, accept=c(".bed"))
  }})
  
  
  #Look for existing cookie
  # observeEvent(
  #   ignoreNULL = TRUE,
  #   eventExpr = {
  #     input$path_cookie
  #   },
  #   handlerExpr = {
  #      if ( (input$path_cookie != "[null]") && !is.null(input$path_cookie) && !is.na(input$path_cookie)) {
  #         #Uploading the name displayed in Data Folder
  #         updateDirectoryInput(session, 'data_folder', value =  input$path_cookie)
  #         init$data_folder <- gsub(pattern = "\"|\\[|\\]|\\\\", "",as.character(input$path_cookie))
  # 
  #         init$available_raw_datasets <- list.dirs(path=file.path(init$data_folder, "datasets"), full.names=FALSE, recursive=FALSE)
  #         init$available_reduced_datasets <- get.available.reduced.datasets()
  # 
  #         unlink(file.path("www", "images", "*"))  # delete all images produced in the last run
  #         unlink(file.path(".", "*.csv"))
  #         file.copy(list.files(file.path(init$data_folder, "datasets"), ".pdf$", full.names=TRUE), file.path("www", "images"))  # copy saved images into app
  #       }
  #     })
  
  #Selecting a working directory using readDirectoryInput(input$data_folder) and saving cookie
  # observeEvent(
  #   ignoreNULL = TRUE,
  #   eventExpr = {
  #     input$data_folder
  #   },
  #   handlerExpr = {
  #     if (input$data_folder > 0) {
  #       # launch the directory selection dialog with initial path read from the widget
  #       print("OPENING DIR, LOOKING INTO SHINY SERVER")
  #       print(list.files("/home/"))
  #       folder = "/var/lib/shiny-server/" #default = readDirectoryInput(session, 'data_folder')
  #       print(list.files(folder))
  #       js$save_cookie(folder)
  #       if (!is.na(folder)){
  #         init$data_folder <- folder
  #         init$available_raw_datasets <- list.dirs(path=file.path(init$data_folder, "datasets"), full.names=FALSE, recursive=FALSE)
  #         init$available_reduced_datasets <- get.available.reduced.datasets()
  #         unlink(file.path("www", "images", "*"))  # delete all images produced in the last run
  #         unlink(file.path(".", "*.csv"))
  #         file.copy(list.files(file.path(init$data_folder, "datasets"), ".pdf$", full.names=TRUE), file.path("www", "images"))  # copy saved images into app
  # 
  #         updateDirectoryInput(session, 'data_folder', value = folder)
  #         js$save_cookie(folder)
  #       }
  #     }
  #   }
  # )
  
  observeEvent(input$compile_dataset, {  # save new dataset
    req(input$new_dataset_name, input$annotation, input$datafile_matrix)
    datamatrix <- NULL
    annot_raw <- NULL
    withProgress(message='Compiling new data set...',value = 0, {
      print(paste0("init$data_folder",init$data_folder))
      dir.create(file.path(init$data_folder, "datasets"))
      dir.create(file.path(init$data_folder, "datasets", input$new_dataset_name))
      dir.create(file.path(init$data_folder, "datasets", input$new_dataset_name, "reduced_data"))
      dir.create(file.path(init$data_folder, "datasets", input$new_dataset_name, "cor_filtered_data"))
      dir.create(file.path(init$data_folder, "datasets", input$new_dataset_name, "consclust"))
      dir.create(file.path(init$data_folder, "datasets", input$new_dataset_name, "supervised"))
      write.table(input$annotation, file.path(init$data_folder, 'datasets', input$new_dataset_name, 'annotation.txt'), row.names=FALSE, col.names=FALSE, quote=FALSE)
      incProgress(0.3, detail="reading data matrices")

      for(i in 1:dim(input$datafile_matrix)[1]){
        datamatrix_single <- read.table(input$datafile_matrix$datapath[i], header=TRUE, stringsAsFactors=FALSE)

        #perform some checks on data format
        # matchingRN <- grep("[[:alnum:]]+:[[:digit:]]+-[[:digit:]]+", rownames(datamatrix_single)) # check rowname format
        # if(length(matchingRN) < length(rownames(datamatrix_single))){
        #   showNotification(paste0(input$datafile_matrix$name, " contains ", (length(rownames(datamatrix_single))-length(matchingRN)), " rownames that do not conform to the required format. Please check your data matrix and try again."), duration=NULL, closeButton=TRUE, type="warning")
        #   if(length(matchingRN) < 5){ # almost all rownames are wrong
        #     showNotification("Maybe your rownames are contained in the first column instead? In this case, remove the header of this column so that they are interpreted as rownames.", duration=NULL, closeButton=TRUE, type="warning")
        #   }
        #   unlink(file.path(init$data_folder, "datasets", input$new_dataset_name), recursive=TRUE)
        #   return()
        # }
        # numericC <- sapply(datamatrix_single, is.numeric) # check if matrix is numeric
        # if(sum(numericC) < ncol(datamatrix_single)){
        #   showNotification(paste0(input$datafile_matrix$name, " contains non-numeric columns at the following indices: ", which(numericC==FALSE), ". Please check your data matrix and try again."), duration=NULL, closeButton=TRUE, type="warning")
        #   unlink(file.path(init$data_folder, "datasets", input$new_dataset_name), recursive=TRUE)
        #   return()
        # }
        #########################
        #TO REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # names = datamatrix_single$X0
        # datamatrix_single = datamatrix_single[,-1]
        # rownames(datamatrix_single) = names
        #TO REMOVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #########################
        datamatrix_single <- datamatrix_single[!duplicated(rownames(datamatrix_single)),] #put IN for new format
        
        
        rownames(datamatrix_single) <- gsub(":", "_", rownames(datamatrix_single))
        rownames(datamatrix_single) <- gsub("-", "_", rownames(datamatrix_single))
        

        total_cell <- length(datamatrix_single[1,])
        sample_name <- gsub('.{4}$', '', input$datafile_matrix$name[i])
        annot_single <- data.frame(barcode=colnames(datamatrix_single), cell_id=paste0(sample_name, "_c", 1:total_cell), sample_id=rep(sample_name, total_cell), batch_id=i)
        
        colnames(datamatrix_single) <- annot_single$cell_id
              
        if(is.null(datamatrix)){ datamatrix <- datamatrix_single
        }else{
          common_regions <- intersect(rownames(datamatrix), rownames(datamatrix_single))
          datamatrix <- cbind(datamatrix[common_regions,], datamatrix_single[common_regions,])
        }
        if(is.null(annot_raw)){ annot_raw <- annot_single} else{ annot_raw <- rbind(annot_raw, annot_single)}
      }
      incProgress(0.7, detail="saving raw data")
      
      #Removing weird chromosomes
      splitID <- sapply(rownames(datamatrix), function(x) strsplit(as.character(x), split="_"))
      normalChr <- which(sapply(splitID, length) <= 3) # weird chromosomes contain underscores in the name
      datamatrix <- datamatrix[normalChr,]
      
      #Remove chrM from mat if it is inside
      if(length(grep("chrM",rownames(datamatrix)))>0)  datamatrix <- datamatrix[-grep("chrM",rownames(datamatrix)),]
      

      save(datamatrix, annot_raw, file=file.path(init$data_folder, "datasets", input$new_dataset_name, "scChIP_raw.RData"))
      
      init$available_raw_datasets <- list.dirs(path=file.path(init$data_folder, "datasets"), full.names=FALSE, recursive=FALSE)
      init$datamatrix <- datamatrix
      init$annot_raw <- annot_raw
      
      updateActionButton(session, "compile_dataset", label="Added successfully", icon = icon("check-circle"))
    })
  })
  
  reactVal <- reactiveValues(annotColors=NULL, annotColors_filtered=NULL)
  
  observeEvent(input$selected_raw_dataset, {  # load precompiled dataset and update coverage plot
    if(!is.null(input$selected_raw_dataset) & input$selected_raw_dataset != ""){
      print("Selecting precompiled dataset and update coverage plot ")
      
      myData = new.env()
      load(file.path(init$data_folder,"datasets", input$selected_raw_dataset, "scChIP_raw.RData"), envir = myData)
      init$datamatrix <- myData$datamatrix
      init$annot_raw <- myData$annot_raw
          
      }
  })
  
  observeEvent(input$selected_raw_dataset,{
    req(input$selected_raw_dataset)
    
    if(!is.null(input$selected_reduced_dataset)){
      delay(1500, {
        if(gsub(pattern ="_\\d*_\\d*_\\d*_\\w*","",input$selected_reduced_dataset) != input$selected_raw_dataset){
          
          showNotification(paste0("Warning : Selected raw dataset '",input$selected_raw_dataset, "' is different from selected reduced dataset '", input$selected_reduced_dataset,"'"), duration=5, closeButton=TRUE, type="warning")
        }
      })
    }
  })
  observeEvent(input$dim_reduction, {  # perform QC filtering and dim. reduction
    annotationId <- annotation_id_norm()
    exclude_regions <- if(input$exclude_regions) setNames(read.table(input$exclude_file$datapath, header=FALSE, stringsAsFactors=FALSE), c("chr", "start", "stop")) else NULL
  
    callModule(moduleFiltering_and_Reduction, "Filtering_and_Reduction", reactive({input$selected_raw_dataset}), reactive({input$min_coverage_cell}),
               reactive({input$min_cells_window/100.0}), reactive({input$quant_removal}), reactive({init$datamatrix}), reactive({init$annot_raw}),
               reactive({init$data_folder}),reactive({annotationId}), reactive({exclude_regions}))
    init$available_reduced_datasets <- get.available.reduced.datasets()
    updateActionButton(session, "dim_reduction", label="Processed and saved successfully", icon = icon("check-circle"))
  })
  
  observeEvent(input$new_dataset_name, {  # reset label on actionButtion when new dataset is added
    updateActionButton(session, "compile_dataset", label="Compile dataset", icon=character(0))
  })
  
  observeEvent({input$selected_raw_dataset  # reset label on actionButtion when new filtering should be filtered
   input$min_coverage_cell
   input$quant_removal
   input$min_cells_window}, {
   updateActionButton(session, "dim_reduction", label="Apply and save data set", icon=character(0))
  })
  
  reduced_dataset <- eventReactive(input$selected_reduced_dataset, { # load reduced data set to work with on next pages
    req(input$selected_reduced_dataset)
    
    init$available_filtered_datasets <- get.available.filtered.datasets(dataset_name(), input$selected_reduced_dataset)
    file_index <- match(c(input$selected_reduced_dataset), reduced_datasets())
    filename_sel <- file.path(init$data_folder, "datasets", init$available_reduced_datasets[file_index])
    myData = new.env()
    load(filename_sel, envir = myData)
    myData
    
  })
  
  pca <- reactive({ reduced_dataset()$pca })
  annot <- reactive({ reduced_dataset()$annot })
  tsne <- reactive({ reduced_dataset()$tsne })
  
  observeEvent(input$selected_reduced_dataset, {  # load coloring used for PCA, tSNE etc.
    req(input$selected_reduced_dataset)
    
    updateSelectInput(session, "selected_raw_dataset", label = "Select data set :", choices = init$available_raw_datasets,selected = gsub(pattern = "_\\d*_\\d*_\\d*_\\w*", "",input$selected_reduced_dataset))
    
    if(file.exists(file.path(init$data_folder, "datasets", dataset_name(), "reduced_data", paste0(input$selected_reduced_dataset, "_annotColors.RData")))){
      myData = new.env()
      load(file.path(init$data_folder, "datasets", dataset_name(), "reduced_data", paste0(input$selected_reduced_dataset, "_annotColors.RData")), envir=myData)
      reactVal$annotColors <- myData$annotColors_save
    }else{
      reactVal$annotColors <- NULL
      tmp_meta <- data.frame(Sample=rownames(annot()), sample_id=annot()$sample_id) # modify if coloring should be possible for other columns
      reactVal$annotColors <- data.frame(sample_id=as.data.frame(anocol())$sample_id) %>% setNames(str_c(names(.), "_Color")) %>% rownames_to_column("Sample") %>% left_join(tmp_meta,. , by="Sample")
    }
  })
  
  cell_cov_plot <- reactive({ggplot(cell_cov_df(), aes(x=coverage)) + geom_histogram(color="black", fill="steelblue", bins=input$coverage_bins) + labs(x="read coverage per cell") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                                                                                                                                                                                            panel.background=element_blank(), axis.line=element_line(colour="black"),
                                                                                                                                                                                            panel.border=element_rect(colour="black", fill=NA))})
  output$cell_coverage <- renderPlotly( ggplotly(cell_cov_plot(), tooltip="Sample", dynamicTicks=T) )
  
  observeEvent(input$delete_dataset, {  # delete selected dataset
    withProgress(message='Deleting data set', value = 0, {
      incProgress(amount=0.5, detail=paste("..."))
      unlink(file.path(init$data_folder, "datasets", input$selected_delete_dataset), recursive=TRUE)
      file.remove(list.files(path=file.path(init$data_folder, "datasets"), pattern=paste0("consClust_", dataset_name(), "*"), full.names=TRUE))
      file.remove(list.files(path=file.path("www", "images"), pattern=paste0("consClust_", dataset_name(), "*"), full.names=TRUE))
      init$available_raw_datasets <- list.dirs(path=file.path(init$data_folder, "datasets"), full.names=FALSE, recursive=FALSE)
      init$available_reduced_datasets <- get.available.reduced.datasets()
      incProgress(amount=0.5, detail=paste("... finished"))
    })
    showNotification("Data set successfully deleted.", duration=5, closeButton=TRUE, type="warning")
  })
  
  output$num_cell_after_QC_filt <- function(){
    req(reduced_dataset())

    table <- as.data.frame(table(init$annot_raw$sample_id))
    table_filtered <- as.data.frame(table(reduced_dataset()$annot$sample_id))
    
    colnames(table) = c("Sample","#Cells Before Filtering")
    rownames(table) = NULL 
    colnames(table_filtered) = c("Sample","#Cells After Filtering")
    rownames(table_filtered) = NULL 

    table_both = left_join(table,table_filtered, by=c("Sample"))
    table_both[,1] = as.character(table_both[,1])
    table_both = table_both %>% 
      bind_rows(., tibble(Sample="",`#Cells Before Filtering`=sum(table_both[,2]),`#Cells After Filtering`=sum(table_both[,3]) ) )
    table_both %>% kable(escape=F, align="c") %>% kable_styling(c("striped", "condensed"), full_width = T) %>% group_rows("Total cell count", dim(table_both)[1], dim(table_both)[1])
  }
  
  output$table_QC_filt_box <- renderUI({
    if(!is.null(reduced_dataset())){
          column(12, align="left", tableOutput("num_cell_after_QC_filt"))
      
    }
  })
  
  observeEvent(input$selected_reduced_dataset,{if(nchar(input$selected_reduced_dataset)>0){unlocked$list$selected_reduced_dataset=T}else{for(i in names(unlocked$list)){unlocked$list[[i]]=F}}})
  observeEvent(unlocked$list,{able_disable_tab(("selected_reduced_dataset"),"pca_plots")}) # if conditions are met, unlock tab Dimensionality Reduction
  
  
  ###############################################################
  # 2. PCA
  ###############################################################
  
  output$pc_select_x <- renderUI({ selectInput("pc_select_x", "X",choices=paste0("PC", c(1:15)), selected="PC1") })
  output$pc_select_y <- renderUI({ selectInput("pc_select_y", "Y",choices=paste0("PC", c(1:15)), selected="PC2") })
  output$pca_color_2D <- renderUI({selectInput("pca_color_2D", "Color by", choices=c('sample_id', 'total_counts')) })
  output$pca_anno_2D <- renderUI({selectInput("pca_anno_2D", "Labels", choices=c('none', 'barcode', 'cell_id', 'sample_id', 'total_counts')) })
  output$tsne_color <- renderUI({ selectInput("tsne_color", "Color by", choices=c('sample_id', 'total_counts')) })
  output$tsne_anno <- renderUI({selectInput("tsne_anno", "Labels", choices=c('none', 'barcode', 'cell_id', 'sample_id', 'total_counts')) })
  

  
  pca_2D <- reactive({
    req(input$pca_anno_2D, input$pca_color_2D,levelsSelectedAnnot_pca())
    plot_df <- as.data.frame(cbind(pca(), annot()))
    p <- ggplot(plot_df, aes_string(x=input$pc_select_x, y=input$pc_select_y)) + geom_point(alpha=0.6, aes(color=annot()[, input$pca_color_2D])) +
      labs(color=input$pca_color_2D) +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.background=element_blank(), axis.line=element_line(colour="black"),
            panel.border=element_rect(colour="black", fill=NA))
    if(input$pca_anno_2D != 'none'){ p <- p + geom_text(aes(label=annot()[, input$pca_anno_2D])) }
    if(input$pca_color_2D == 'total_counts'){
      p <- p + scale_color_gradientn(colours = matlab.like(100))
    }else{
      cols <- paste0("c(", paste0("input$pca_col", sort(levelsSelectedAnnot_pca()), collapse = ", "), ")")
      cols <- eval(parse(text = cols))
      req(cols)
      lev <- sort(unique(levelsSelectedAnnot_pca()))
      names(cols) <- lev
      p <- p + scale_color_manual(values = cols)
    }
    levelsSelectedAnnot_pca = levelsSelectedAnnot_pca()
    save(levelsSelectedAnnot_pca,file=file.path(init$data_folder, "datasets", dataset_name(), "reduced_data", paste0(input$selected_reduced_dataset, "_levelsSelectedAnnot_pca.RData")) )
    unlocked$list$pca=T
    p
  })
  output$pca_plot_2D <- renderPlotly( ggplotly(pca_2D(), tooltip="Sample", dynamicTicks=T) )
  
  tsne_p <- reactive({
    req(input$tsne_anno, input$tsne_color,levelsSelectedAnnot_tSNE())
    p <- ggplot(as.data.frame(tsne()$Y), aes(x=V1, y=V2)) + geom_point(alpha=0.6, aes(color=annot()[, input$tsne_color])) +
      labs(color=input$tsne_color, x="t-SNE 1", y="t-SNE 2") +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.background=element_blank(), axis.line=element_line(colour="black"),
            panel.border=element_rect(colour="black", fill=NA))
    if(input$tsne_anno != 'none'){ p <- p + geom_text(aes(label=annot()[, input$tsne_anno])) }
    if(input$tsne_color == 'total_counts'){
      p <- p + scale_color_gradientn(colours = matlab.like(100))
    }else{
      cols <- paste0("c(", paste0("input$tsne_col", sort(levelsSelectedAnnot_tSNE()), collapse = ", "), ")")
      cols <- eval(parse(text = cols))
      req(cols)
      lev <- sort(unique(levelsSelectedAnnot_tSNE()))
      names(cols) <- lev
      p <- p + scale_color_manual(values = cols)
    }
    unlocked$list$tsne=T
    p
  })
  output$tsne_plot <- renderPlotly( ggplotly(tsne_p(), tooltip="Sample", dynamicTicks=T) )
  
  
  
  output$pca_color_box <- renderUI({
    req(input$pca_color_2D)
    if(input$pca_color_2D != 'total_counts'){
      box(title="PCA color settings", width = NULL, status = "success", solidHeader = T,
          column(6, htmlOutput("pca_colorPicker")),
          column(6 , br(), actionButton("PcaColReset", "Default colours", icon = icon("undo")),
                 br(), br(), actionButton("saveColors_pca", "Save colors & apply to all", icon = icon("save"))))
    }
  })
  
  output$tsne_color_box <- renderUI({
    req(input$tsne_color)
    if(input$tsne_color != 'total_counts'){
      box(title="tSNE color settings", width = NULL, status = "success", solidHeader = T,
          column(6, htmlOutput("tsne_colorPicker")),
          column(4 , br(), actionButton("tSNEColReset", "Default colours", icon = icon("undo")),
                 br(), br(), actionButton("saveColors_tsne", "Save colors & apply to all", icon = icon("save"))))
    }
  })
  
  observeEvent(input$saveColors_pca, {  # [pca] change modified colors for all plots where it applies and save RData locally
    cols <- paste0("c(", paste0("input$pca_col", sort(levelsSelectedAnnot_pca()), collapse = ", "), ")")
    cols <- eval(parse(text = cols))
    lev <- sort(unique(levelsSelectedAnnot_pca()))
    names(cols) <- lev
    annot_Color_Custom <- as.data.frame(cols) %>% rownames_to_column(input$pca_color_2D)%>% setNames(c(input$pca_color_2D, str_interp("${input$pca_color_2D}_Color") ))
    reactVal$annotColors <- reactVal$annotColors %>% select(-c(str_interp("${input$pca_color_2D}_Color"))) %>% right_join(annot_Color_Custom, by=input$pca_color_2D  )
    annotColors_save <- reactVal$annotColors
    save(annotColors_save, file=file.path(init$data_folder, "datasets", dataset_name(), "reduced_data", paste0(input$selected_reduced_dataset, "_annotColors.RData")))
  })
  
  observeEvent(input$saveColors_tsne, {  # [tsne] change modified colors for all plots where it applies and save RData locally
    cols <- paste0("c(", paste0("input$tsne_col", sort(levelsSelectedAnnot_tSNE()), collapse = ", "), ")")
    cols <- eval(parse(text = cols))
    lev <- sort(unique(levelsSelectedAnnot_tSNE()))
    names(cols) <- lev
    annot_Color_Custom <- as.data.frame(cols) %>% rownames_to_column(input$tsne_color)%>% setNames(c(input$tsne_color, str_interp("${input$tsne_color}_Color") ))
    reactVal$annotColors <- reactVal$annotColors %>% select(-c(str_interp("${input$tsne_color}_Color"))) %>% right_join(annot_Color_Custom, by=input$tsne_color  )
    annotColors_save <- reactVal$annotColors
    save(annotColors_save, file=file.path(init$data_folder, "datasets", dataset_name(), "reduced_data", paste0(input$selected_reduced_dataset, "_annotColors.RData")))
  })
  
  anocol <- reactive({
    req(annot(),dataset_name())
    ac <- geco.annotToCol2(annotS=annot()[, annotCol()], annotT=annot(), plotLegend=T, plotLegendFile=file.path(init$data_folder, "datasets", dataset_name(), "Annotation_legends.pdf"), categCol=NULL)
    if(!is.null(reactVal$annotColors)){
      lapply(colnames(ac), function(col){ if(col != "total_counts"){ ac[, col] <<- as.character(reactVal$annotColors[which(reactVal$annotColors$Sample %in% rownames(ac)), paste0(col, '_Color')]) } })
    }
    ac
  })
  
  levelsSelectedAnnot_pca <- reactive({
    annot() %>% pull(input$pca_color_2D) %>% unique() %>% as.vector()
  })
  
  levelsSelectedAnnot_tSNE <- reactive({
    annot() %>% pull(input$tsne_color) %>% unique() %>% as.vector()
  })
  
  observeEvent(input$PcaColReset, {
    lev <- sort(unique(levelsSelectedAnnot_pca()))
    cols <- gg_fill_hue(length(lev))
    lapply(seq_along(lev), function(i) {
      do.call(what="updateColourInput", args=list(session=session, inputId=paste0("pca_col", lev[i]), value=cols[i]))
    })
  })
  
  observeEvent(input$tSNEColReset, {
    lev <- sort(unique(levelsSelectedAnnot_tSNE()))
    cols <- gg_fill_hue(length(lev))
    lapply(seq_along(lev), function(i) {
      do.call(what="updateColourInput", args=list(session=session, inputId=paste0("tsne_col", lev[i]), value=cols[i]))
    })
  })
  
  output$pca_colorPicker <- renderUI({
    lev <- sort(unique(levelsSelectedAnnot_pca()))
    if(str_interp("${input$pca_color_2D}_Color") %in% colnames(reactVal$annotColors)) {
    colsModif <- reactVal$annotColors %>% dplyr::select(one_of( str_interp("${input$pca_color_2D}_Color")) ) %>% distinct() %>% pull( str_interp("${input$pca_color_2D}_Color"))
    
    lapply(seq_along(lev), function(i) {
      colourpicker::colourInput(inputId=paste0("pca_col", lev[i]), label=paste0("Choose colour for ", lev[i]), value=colsModif[i], returnName=TRUE) ## Add ", palette = "limited"" to get selectable color panel       
    })
    }
  })
  
  output$tsne_colorPicker <- renderUI({
    lev <- sort(unique(levelsSelectedAnnot_tSNE()))
    if(str_interp("${input$tsne_color}_Color") %in% colnames(reactVal$annotColors)) {
    colsModif <- reactVal$annotColors %>% dplyr::select(one_of( str_interp("${input$tsne_color}_Color")) ) %>% distinct() %>% pull( str_interp("${input$tsne_color}_Color"))
    
    lapply(seq_along(lev), function(i) {
      colourpicker::colourInput(inputId=paste0("tsne_col", lev[i]), label=paste0("Choose colour for ", lev[i]), value=colsModif[i], returnName=TRUE) ## Add ", palette = "limited"" to get selectable color panel       
    })
    }
  })
  observeEvent(unlocked$list,able_disable_tab(c("selected_reduced_dataset","pca","tsne"),"cor_clustering")) # if conditions are met, unlock tab Correlation Clustering
  
  
  ###############################################################
  # 3. Correlation clustering                                  ##
  ###############################################################
  
  corColors <- colorRampPalette(c("royalblue","white","indianred1"))(256)
  mati <- reactive({t(pca()[,1:50]) })
  
  hc_cor <- reactive({
    tmp = hclust(as.dist(1 - cor(mati())), method="ward.D") 
    tmp$labels = rep("",length(tmp$labels))
    tmp
    })
  
  mat.so.cor <- reactive({ mati()[,hc_cor()$order] })
  
  hc_pca_plot <- reactive({
    req(mat.so.cor(), hc_cor(), anocol())
    dev.set(4)  # this is a very weird fix, but it seems to work. If it causes problems, delete it and figure out another way.
    unlocked$list$cor_clust_plot=TRUE;
    geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor()),
                              hc=hc_cor(),
                              hmColors=corColors,
                              anocol=anocol()[hc_cor()$order,],
                              xpos=c(0.15, 0.9, 0.164, 0.885),
                              ypos=c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
                              dendro.cex=0.04,
                              xlab.cex=0.8,
                              hmRowNames=FALSE)
    
    })

    output$corr_clust_pca_plot <- renderPlot(hc_pca_plot())
  
  output$download_cor_clust_plot <- downloadHandler(
    filename=function(){ paste0("correlation_clustering_", input$selected_reduced_dataset, ".png")},
    content=function(file){
      png(file, width=1200, height=1400, res=300,pointsize = 8)
      geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor()),
                                                hc=hc_cor(),
                                                hmColors=corColors,
                                                anocol=anocol()[hc_cor()$order,],
                                                xpos=c(0.15, 0.9, 0.164, 0.885),
                                                ypos=c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
                                                dendro.cex=0.04,
                                                xlab.cex=0.5,
                                                hmRowNames=FALSE)
      dev.off()

  })
  
  correlation_values <- reactiveValues(limitC=vector(length=500))
  corChIP <- reactive({ cor(mati()) })
  z <- reactive({ matrix(sample(mati()), nrow=dim(mati())[1]) })
  thresh2 <- reactive({quantile(cor(z()), probs=seq(0,1,0.01))})
  limitC <- reactive({thresh2()[input$corr_thresh+1]})
  
  cell_cor_hist <- reactive({
    hist(corChIP(), prob=TRUE, col=alpha("steelblue", 0.8), breaks=50, ylim=c(0,4), main="Distribution of cell to cell correlation scores", xlab="Pearson Corr Scores")
    lines(density(corChIP()), col="blue", lwd=2)
    lines(density(cor(z())), col="black", lwd=2)
    abline(v=limitC(), lwd=2, col="red", lty=2)
    legend("topleft", legend=c("dataset", "randomized data", "correlation threshold"), col=c("blue", "black", "red"), lty=c(1, 1, 2), cex=0.8)
  })
  
  output$cell_cor_hist_plot <- renderPlot( cell_cor_hist() )
  
  output$download_cor_clust_hist_plot <- downloadHandler(
    filename=function(){ paste0("correlation_distribution_", input$selected_reduced_dataset, ".png")},
    content=function(file){
      png(file, width=2000, height=1400, res=300)
      hist(corChIP(), prob=TRUE, col=alpha("steelblue", 0.8), breaks=50, ylim=c(0,4),cex=0.4, main="Distribution of cell to cell correlation scores", xlab="Pearson Corr Scores")
      lines(density(corChIP()), col="blue", lwd=2)
      lines(density(cor(z())), col="black", lwd=2)
      abline(v=limitC(), lwd=2, col="red", lty=2)
      legend("topleft", legend=c("dataset", "randomized data", "correlation threshold"), col=c("blue", "black", "red"), lty=c(1, 1, 2), cex=0.8)
      dev.off()
  })
  
  cf <- reactiveValues(sel2=NULL, mati2=NULL, mat.so.cor2=NULL, hc_cor2=NULL, annot_sel=NULL, anocol_sel=NULL)
  
  observeEvent(input$filter_corr_cells, {  # retreiveing cells with low correlation score
    withProgress(message='Filtering correlated cells...', value = 0, {
      incProgress(amount=0.8, detail=paste("filtering"))
        for(i in 1:500){
        random_mat <-  matrix(sample(mati()), nrow=dim(mati())[1])
        thresh2 <- quantile(cor(random_mat), probs=seq(0,1,0.01))
        limitC <-  thresh2[input$corr_thresh+1]
        correlation_values$limitC[i] = limitC
      }
      
      correlation_values$limitC_mean = mean(correlation_values$limitC,na.rm = T)
      
      cf$sel2 <- (apply(corChIP(), 1, function(x) length(which(x>correlation_values$limitC_mean)))) > (input$percent_corr*0.01)*dim(corChIP())[1]
      cf$mati2 <- mati()[, cf$sel2]
      cf$hc_cor2 <- hclust(as.dist(1 - cor(cf$mati2)), method="ward.D")
      cf$hc_cor2$labels = rep("",length(cf$hc_cor2$labels))
      cf$mat.so.cor2 <- cf$mati2[, cf$hc_cor2$order]
      cf$annot_sel <- annot()[cf$sel2,]
      cf$anocol_sel <- geco.annotToCol2(annotS=cf$annot_sel[, annotCol()], annotT=cf$annot_sel, plotLegend=T, plotLegendFile=file.path(init$data_folder, "datasets", dataset_name(),"Annotation_legends_reclustering.pdf"), categCol=NULL)
      lapply(colnames(cf$anocol_sel), function(col){ if(col != "total_counts"){ cf$anocol_sel[, col] <<- as.character(reactVal$annotColors[which(reactVal$annotColors$Sample %in% rownames(cf$anocol_sel)), paste0(col, '_Color')]) } })
      
      mati2 <- cf$mati2
      annot_sel <- cf$annot_sel
      mat.so.cor2 <- cf$mat.so.cor2
      anocol_sel <- cf$anocol_sel

      
      limitC_mean= correlation_values$limitC_mean
      corChIP = corChIP()
      mati = mati()
      
      hc_cor2 <- cf$hc_cor2
      hc_cor = hc_cor()
      save(mati2, annot_sel, mat.so.cor2, anocol_sel,hc_cor,hc_cor2, file=file.path(init$data_folder, "datasets", dataset_name(), "cor_filtered_data", paste0(input$selected_reduced_dataset, "_", input$corr_thresh, "_", input$percent_corr, ".RData")))
      incProgress(amount=0.2, detail=paste("finished"))
      updateActionButton(session, "filter_corr_cells", label="Saved", icon = icon("check-circle"))
  
    })
  })
  
  observeEvent({input$selected_reduced_dataset  # reset label on actionButtion when new filtering should be filtered
    input$percent_corr}, {
      updateActionButton(session, "filter_corr_cells", label="Filter & save", icon=character(0))
    })
  
  output$corr_clust_filter_plot <- renderPlot(geco.hclustAnnotHeatmapPlot(x=cor(cf$mat.so.cor2),
                                                                          hc=cf$hc_cor2,
                                                                          hmColors=corColors,
                                                                          anocol=cf$anocol_sel[cf$hc_cor2$order,],
                                                                          xpos=c(0.15, 0.9, 0.164, 0.885),
                                                                          ypos=c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
                                                                          dendro.cex=0.04,
                                                                          xlab.cex=0.8,
                                                                          hmRowNames=FALSE))
  
  
  
  output$download_cor_clust_filtered_plot <- downloadHandler(
    filename=function(){ paste0("correlation_clustering_filtered_", input$selected_reduced_dataset, ".png")},
    content=function(file){
      png(file,  width=1200, height=1400, res=300,pointsize = 8)
      geco.hclustAnnotHeatmapPlot(x=cor(cf$mat.so.cor2),
                                                hc=cf$hc_cor2,
                                                hmColors=corColors,
                                                anocol=cf$anocol_sel[cf$hc_cor2$order,],
                                                xpos=c(0.15, 0.9, 0.164, 0.885),
                                                ypos=c(0.1, 0.5, 0.5, 0.6, 0.62, 0.95),
                                                dendro.cex=0.04,
                                                xlab.cex=0.8,
                                                hmRowNames=FALSE)
      dev.off()
  })
  
  output$corr_filtered_hist <- renderUI({
    if(!is.null(cf$anocol_sel)){
      column(12, align="left", plotOutput("corr_clust_filter_plot", height=500, width=500),
             downloadButton("download_cor_clust_filtered_plot", "Download image"), br(),
             tableOutput('num_cell_after_cor_filt')
             )
    }
  })
  
 
  
  output$num_cell_before_cor_filt <- function(){
    req(annot())
    table <- as.data.frame(table(annot()$sample_id))
    colnames(table) = c("Sample","#Cells")
    rownames(table) = NULL
    
    #Retrieve sample colors from user specified colors & add to table
    colors = unique(reactVal$annotColors[,c("sample_id","sample_id_Color")])
    colors = as.vector(as.character(left_join(table,colors,by=c("Sample"="sample_id"))[,"sample_id_Color"]))
    colors = c(col2hex(colors),"")
    

    table[,1] = as.character(table[,1])
    table = table %>% bind_rows(., tibble(Sample="",`#Cells`=sum(table[,-1])) )
    
    table %>% mutate(Sample=cell_spec(Sample, color="white", bold=T, background=colors)) %>% kable(escape=F, align="c") %>% kable_styling(c("striped", "condensed"), full_width = T) %>% group_rows("Total cell count", dim(table)[1], dim(table)[1])
  }
  
  output$num_cell_after_cor_filt <- function(){
    req(annot(),cf$annot_sel)
    table <- as.data.frame(table(annot()$sample_id))
    table_filtered <- as.data.frame(table(cf$annot_sel$sample_id))
    colnames(table) = c("Sample","#Cells Before Filtering")
    rownames(table) = NULL 
    colnames(table_filtered) = c("Sample","#Cells After Filtering")
    rownames(table_filtered) = NULL 
    
    #Retrieve sample colors from user specified colors & add to table
    colors = unique(reactVal$annotColors[,c("sample_id","sample_id_Color")])
    colors = as.vector(as.character(left_join(table,colors,by=c("Sample"="sample_id"))[,"sample_id_Color"]))
    colors = c(col2hex(colors),"")
    
    table_both = left_join(table,table_filtered, by=c("Sample"))
    table_both[,1] = as.character(table_both[,1])
    table_both = table_both %>% bind_rows(., tibble(Sample="",`#Cells Before Filtering`=sum(table_both[,2]),`#Cells After Filtering`=sum(table_both[,3]) ) )
    
    table_both %>% mutate(Sample=cell_spec(Sample, color="white", bold=T, background=colors)) %>% kable(escape=F, align="c") %>% kable_styling(c("striped", "condensed"), full_width = T) %>% group_rows("Total cell count", dim(table_both)[1], dim(table_both)[1])
  }
  
  
  observeEvent(init$available_filtered_datasets,ignoreNULL = T,{if(length(init$available_filtered_datasets)>0){unlocked$list$filtered_datasets =TRUE}else{unlocked$list$filtered_datasets =FALSE}})
  observeEvent(unlocked$list,able_disable_tab(c("selected_reduced_dataset","cor_clust_plot","filtered_datasets"),"cons_clustering")) # if conditions are met, unlock tab Consensus Clustering on Correlated cells
  
  ###############################################################
  # 4. Consensus clustering on correlated cells
  ###############################################################
  
  filtered_datasets <- reactive({ if (is.null(init$available_filtered_datasets)) c() else gsub('.{6}$', '', basename(init$available_filtered_datasets)) })
  output$selected_filtered_dataset <- renderUI({ selectInput("selected_filtered_dataset", "Select data filtered by correlated cells (saved on previous page):", choices=filtered_datasets()) })
  output$filtered_data_selection_format <- renderText({"The name of the filtered dataset is composed of the following information: data set name, min percentage of reads per cell, 
    min percentage of cells to support a window, quantile of cell read counts to keep, correlation threshold, percent of cell correlation. To work on a different dataset or different preprocessing state, select it on the first page."})
  
  filtered_dataset <- eventReactive(input$selected_filtered_dataset, {
    if(!is.null(input$selected_filtered_dataset) & nchar(input$selected_filtered_dataset) > 5){
      myData = new.env()
      load(file.path(init$data_folder, "datasets", dataset_name(), "cor_filtered_data", paste0(input$selected_filtered_dataset, ".RData")), envir = myData)
      clust$available_k=get_available_k()
      myData
    }
  })
  
  
  observeEvent(input$do_final_figs, {  # load coloring used for PCA, tSNE etc.
    req(input$selected_filtered_dataset)
    if(file.exists(file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_filtered_dataset,"_k_",input$nclust, "_annotColors_filtered.RData")))){
      myData = new.env()
      load(file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_filtered_dataset, "_k_",input$nclust,"_annotColors_filtered.RData")), envir=myData)
      reactVal$annotColors_filtered <- myData$annotColors_filtered_save
    }
    else{
      reactVal$annotColors_filtered <- NULL
    }
  })
  
  get_available_k <- function(){
    saved_clustFiles <- reactive({ list.files(path=file.path(init$data_folder, "datasets", dataset_name(), "consclust"), pattern=paste0(input$selected_filtered_dataset, '_affectation_k*'), full.names=FALSE) })
    k <- reactive({ gsub(".*_affectation_k(.+)\\.RData", "\\1", saved_clustFiles()) })
    isolate(k)
  }
  
  mati2 <- reactive({ filtered_dataset()$mati2 })
  annot_sel <- reactive({ filtered_dataset()$annot_sel })
  mat.so.cor2 <- reactive({ filtered_dataset()$mat.so.cor2 })
  anocol_sel <- reactive({ filtered_dataset()$anocol_sel })
  hc_cor2 <- reactive({ filtered_dataset()$hc_cor2 })
  clust <- reactiveValues(cc.col=NULL, consclust.mat=NULL, hc=NULL, tsne_corr=NULL, annot_sel2=NULL,
                          clust_pdf=NULL, available_k=get_available_k(), chi=NULL)
  
  observeEvent(input$selected_filtered_dataset,priority = 10, { clust$clust_pdf <- paste0("images/consClust_", input$selected_filtered_dataset, ".pdf")})
  
  observeEvent(input$do_cons_clust, {
    withProgress(message='Performing consensus clustering...', value = 0, {
      incProgress(amount=0.4, detail=paste("part one"))
      consclust <- ConsensusClusterPlus(mati2(), maxK=10, reps=1000, pItem=0.8, pFeature=1,
                                      title="Consensus_clustering_dir", clusterAlg="hc", distance="pearson",
                                      innerLinkage="ward.D", finalLinkage="ward.D", seed=3.14, plot="png")
      incProgress(amount=0.4, detail=paste("part two"))
      icl <- calcICL(consclust, plot="png", title="Consensus_clustering_dir")
      consPlots <- lapply(c('01','02','03','04','05','06','07','08','09','10','11','12'), function(i){ rasterGrob(readPNG(file.path("Consensus_clustering_dir", paste0("consensus0", i, ".png")), native=FALSE), interpolate=FALSE, width=c(1))})
      consPlots <- list.append(consPlots, rasterGrob(readPNG(file.path("Consensus_clustering_dir", "icl004.png"), native=FALSE), interpolate=FALSE, width=c(1)))
      layout <- matrix(rep(1, 13), ncol=1, nrow=13)
      pdf(paste0("www/images/consClust_", input$selected_filtered_dataset, ".pdf"))
      print(do.call("marrangeGrob", list(grobs=consPlots, ncol=1, nrow=13, layout_matrix=layout, top="")))
      dev.off()
      unlink("Consensus_clustering_dir", recursive=TRUE)
      file.copy(paste0("www/images/consClust_", input$selected_filtered_dataset, ".pdf"), file.path(init$data_folder, "datasets"))  # copy pdf into local dir on computer
      clust$clust_pdf <- NULL  # needed in order to update the pdf output
      clust$clust_pdf <- paste0("images/consClust_", input$selected_filtered_dataset, ".pdf")
      save(consclust, icl, file=file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_filtered_dataset, ".RData")))
      
      
      
      incProgress(amount=0.2, detail=paste("finished"))
    })
  })
  

  
  consData <- eventReactive(c(input$selected_filtered_dataset,  # load consclust and icl if it exists already for this filtered dataset
    input$do_cons_clust),{
      filename <- file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_filtered_dataset, ".RData"))
    if(file.exists(filename)){
      myData = new.env()
      load(filename, envir = myData)
      myData
    }else{
      NULL
    }
  })
  consclust <- reactive({ consData()$consclust })
  icl <- reactive({ consData()$icl })
  
  output$cons_clust_pdf <- renderUI({req(clust$clust_pdf);if(file.exists(file.path('www', clust$clust_pdf))){tags$iframe(style="height:500px; width:100%", src=clust$clust_pdf) }})
  output$nclust_selection_info <- renderText({"After performing the clustering and checking the results for different numbers of clusters, select here the preferred number of clusters to make additional annotated plots."})
  clustCol <- "ChromatinGroup"
  
  observeEvent(input$do_final_figs, {
    if(is.null(consclust())){
      showNotification("Could not find clustering results. Please make sure to perform the clustering before generating the final figures.", duration=7, closeButton=TRUE, type="warning")
    }else{
      withProgress(message='Preparing final figures...', value = 0, {
        incProgress(amount=0.2, detail=paste("annotating"))
        so <- colnames(mat.so.cor2())
        conscol <- c("#4285F4", "#F4B400","#DB4437", "#0F9D58", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#bd18ea", "#2ef4ca", "#f4cced", "#f4cc03", "#05188a", "#e5a25a", "#06f106", "#85848f", "#000000", "#076f25", "#93cd7f", "#4d0776", "#ffffff")
        cc <- consclust()[[as.integer(input$nclust)]]$consensusClass
        cc. <- lapply(unique(cc), function(z) names(which(cc==z)))
        mat.cc <- geco.groupMat(mati2(), margin=1, groups=cc., method="mean")
        hcc <- hclust(distPearson(t(mat.cc)), method="ward.D")
        clust$annot_sel2 <- isolate(annot_sel())
        clust$annot_sel2[, clustCol] <- paste("C", match(cc, hcc$order), sep="")
        clust$cc.col <- cbind(ChromatinGroup=conscol[match(clust$annot_sel2[so, clustCol], paste("C", 1:as.integer(input$nclust), sep=""))], anocol_sel()[so,])
        
        if( is.null(reactVal$annotColors_filtered) ){
          reactVal$annotColors_filtered <- NULL
          tmp_meta <- data.frame(Sample=rownames(clust$annot_sel2), sample_id=clust$annot_sel2$sample_id, ChromatinGroup=clust$annot_sel2[, clustCol]) # modify if coloring should be possible for other columns
          tmp_meta =data.frame(ChromatinGroup=as.data.frame(clust$cc.col)[,c("ChromatinGroup")]) %>% setNames(str_c(names(.), "_Color")) %>% rownames_to_column("Sample") %>% left_join(tmp_meta,. , by="Sample")
          reactVal$annotColors_filtered = reactVal$annotColors
          reactVal$annotColors_filtered <- inner_join(reactVal$annotColors_filtered,tmp_meta,by=c("Sample","sample_id"))
          
        } 
        
        # else{
        #   lapply(colnames(clust$cc.col), function(col){ if(col != "total_counts"){ clust$cc.col[, col] <<- as.character(reactVal$annotColors_filtered[which(reactVal$annotColors_filtered$Sample %in% rownames(clust$cc.col)), paste0(col, '_Color')]) } })
        #   }
        # 

        
        tmp <- lapply(paste("C", 1:as.integer(input$nclust), sep=""), function(k){
          colnames(mati2())[which(clust$annot_sel2[, clustCol]==k)]
        })
        affectation <- clust$annot_sel2[,c("barcode", "cell_id", "ChromatinGroup", "sample_id")]
        save(affectation, file=file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_filtered_dataset, '_affectation_k', input$nclust, ".RData")))
        clust$consclust.mat <- consclust()[[as.integer(input$nclust)]]$ml
        clust$hc <- hclust(as.dist(1 - clust$consclust.mat), method="ward.D")
        clust$consclust.mat <- clust$consclust.mat[clust$hc$order,]
        incProgress(amount=0.4, detail=paste("performing tSNE"))
        clust$tsne_corr <- Rtsne(t(mati2()), dims=2, pca=FALSE, theta=0.0, perplexity=choose_perplexity(t(mati2())), verbose=TRUE, max_iter = 1000)
        clust$available_k=get_available_k()
        incProgress(amount=0.4, detail=paste("finished"))
      })
    }
  })
  

    output$cons_clust_anno_plot <- renderPlot({
      colors <- cc.col()[names(consclust()[[as.integer(input$nclust)]]$consensusClass),"ChromatinGroup"]

      heatmap(clust$consclust.mat, Colv=as.dendrogram(clust$hc), Rowv=NA, symm=FALSE, scale="none",
              col=colorRampPalette(c("white", "blue"))(100), na.rm=TRUE, labRow=F, labCol=F, mar=c(5, 5),
              main=paste("consensus matrix k=", input$nclust, sep=""), ColSideCol=colors)
    })
  
  output$anno_cc_box <- renderUI({
    if(!is.null(clust$hc)){
      box(title="Annotated consensus clustering", width=NULL, status="success", solidHeader=T,
          column(12, align="left", plotOutput("cons_clust_anno_plot", height=500, width=500),
                 downloadButton("download_anno_cc_plot", "Download image")))
    }
  })
  
  output$download_anno_cc_plot <- downloadHandler(
    filename=function(){ paste0("consensus_clustering_k", input$nclust, "_", input$selected_filtered_dataset, ".png")},
    content=function(file){
      png(file, width=1200, height=800, res=300)
      colors <- cc.col()[names(consclust()[[as.integer(input$nclust)]]$consensusClass),"ChromatinGroup"]
      heatmap(clust$consclust.mat, Colv=as.dendrogram(clust$hc), Rowv=NA, symm=FALSE, scale="none",
              col=colorRampPalette(c("white", "blue"))(100), na.rm=TRUE, labRow=F, labCol=F, mar=c(5, 5),
              main = paste("consensus matrix k=", input$nclust, sep=""), ColSideCol=colors)
      
      dev.off()
  })
  
  output$cor_clust_anno_plot <- renderPlot(geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor2()), hc=hc_cor2(), hmColors=corColors,
                                                                       anocol=cc.col()[,ncol(cc.col()):1], xpos=c(0.15,0.9,0.164,0.885),
                                                                       ypos=c(0.1,0.5,0.5,0.6,0.62,0.95), dendro.cex=0.05,
                                                                       xlab.cex=0.5, hmRowNames=FALSE, hmRowNames.cex=0.01))
  
  output$anno_corc_box <- renderUI({
    if(!is.null(clust$hc)){
      box(title="Annotated correlation clustering", width=NULL, status="success", solidHeader=T,
          column(12, align="left", plotOutput("cor_clust_anno_plot", height=500, width=500)),
          column(12, align="left", tableOutput("hcor_kable")),
          column(5,offset = 2, align="left", br(), htmlOutput("chi_info"),br()),
          column(12, align="left", downloadButton("download_anno_hc_plot", "Download image")))
    }
  })
  
  output$download_anno_hc_plot <- downloadHandler(
    filename=function(){ paste0("hierarchical_clustering_k", input$nclust, "_", input$selected_filtered_dataset, ".png")},
    content=function(file){
      png(file,  width=1200, height=1400, res=300,pointsize = 8)
      geco.hclustAnnotHeatmapPlot(x=cor(mat.so.cor2()), hc=hc_cor2(), hmColors=corColors,
                                    anocol=cc.col()[,ncol(cc.col()):1], xpos=c(0.15,0.9,0.164,0.885),
                                  ypos=c(0.1,0.5,0.5,0.6,0.62,0.95), dendro.cex=0.05,
                                  xlab.cex=0.5, hmRowNames=FALSE, hmRowNames.cex=0.01)
      dev.off()
  })

    
  output$hcor_kable <- function(){
    req(cc.col(), anocol_sel())
    table_raw <- as.data.frame.matrix(table(cc.col()[,"ChromatinGroup"], cc.col()[,"sample_id"]))

     #Overall goodness of fit testing : how fairly are cells  allocated between the clusters ?
    clust$chi <- chisq.test(x=as.matrix(table_raw), correct=FALSE)
    
    # Cluster goodness of fit testing : how fairly are cells  allocated to a particular cluster ?
    chi_pvalues=c()
    for(i in 1:(dim(as.matrix(table_raw))[1])){
      contingency_tab = rbind(table_raw[i,], colSums(table_raw))
      chi <- chisq.test(x=contingency_tab, correct=FALSE)
      chi_pvalues[i]=chi$p.value
    }

    tab <- table_raw
    if(length(unique(annot_sel()$sample_id))==1){
      chi_pvalues=rep(1,dim(as.matrix(table_raw))[1])
      tab <- t(tab)
    }
    chi_pvalues= round(chi_pvalues, 5)
    chi_pvalues[which(chi_pvalues==0)] <- "<0.00001"
    chi_pvalues = c(chi_pvalues,"")

    #colnames(tab) <- sapply(colnames(table_raw), function(x){ annot_sel()[names(anocol_sel()[anocol_sel()[, "sample_id"]==x, "sample_id"][1]), "sample_id"] })
    
    colnames(tab) <- sapply(colnames(table_raw), function(x){  clust$annot_sel2[names(cc.col()[cc.col()[, "sample_id"]==x, "sample_id"][1]), "sample_id"] })
    colors = c(col2hex(rownames(tab)),"") # necessary as kable seems to not recognize all color names so it makes everything white
    rownames(tab) <- sapply(rownames(table_raw), function(x){  clust$annot_sel2[names(cc.col()[cc.col()[, "ChromatinGroup"]==x, "ChromatinGroup"][1]), "ChromatinGroup"] })
    
    
    tab <- rbind(tab, colSums(table_raw))
    tab <- cbind(tab, "#Cells"=c(rowSums(table_raw),sum(rowSums(table_raw))))
    tab <- cbind(tab, "p-value"=chi_pvalues)
    tab <- as.data.frame(cbind(Cluster=c(rownames(tab)[1:length(rownames(tab))-1],""), tab))
    rownames(tab) <- NULL
    if(length(unique(annot_sel()$sample_id))==1){
      tab[,-c(1,dim(tab)[2])] <- sapply(tab[,-c(1,dim(tab)[2])], function(x){ as.numeric(as.character(x)) })
    }else{
      tab[,-c(1,dim(tab)[2])] =apply(tab[,-c(1,dim(tab)[2])], 2, function(x){ as.numeric(x) })
    }
    print(tab)
    # ord = c(order(tab$Cluster[1:nrow(tab)-1]),nrow(tab))
    # tab = tab[ord,]
    tab= tab %>% mutate(Cluster=cell_spec(Cluster, color="white", bold=T, background=colors))
      
    tab %>% kable(escape=F, align="c") %>%
      kable_styling(c("striped", "condensed"), full_width = F) %>%
      group_rows("Total cell count", dim(tab)[1], dim(tab)[1])
  }
  
  output$chi_info <- renderUI({
    req(clust$chi)
    x2 <- round(clust$chi$statistic, 2)
    pval <- round(clust$chi$p.value, 5)
    pval <- if(pval==0) "< 0.00001" else paste0("= ", pval)
    HTML(paste0("Chi-squared test on<br/>the contingency table:<br/>X-squared = ", x2, "<br/>degrees of freedom = ",
           clust$chi$parameter, "<br/>p-value ", pval))
  })
  
  anno_tsne_p <- reactive({
    p <- ggplot(as.data.frame(clust$tsne_corr$Y), aes(x=V1, y=V2)) + geom_point(alpha=0.6, aes(color=clust$annot_sel2[, input$anno_tsne_color])) +
      labs(color=input$anno_tsne_color, x="t-SNE 1", y="t-SNE 2") +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.background=element_blank(), axis.line=element_line(colour="black"),
            panel.border=element_rect(colour="black", fill=NA))
    if(input$anno_tsne_anno != 'none'){ p <- p + geom_text(aes(label=annot()[, input$anno_tsne_anno])) }
    if(input$anno_tsne_color == 'total_counts'){
      p <- p + scale_color_gradientn(colours = matlab.like(100))
    }else{

      cols <- paste0("c(", paste0("input$anno_tsne_col", levelsSelectedAnnot_anno_tSNE(), collapse = ", "), ")")
      cols <- eval(parse(text = cols))
      req(cols)
      lev <- unique(levelsSelectedAnnot_anno_tSNE())
      names(cols) <- lev
      p <- p + scale_color_manual(values = cols)
    }
    p
  })
  output$anno_tsne_plot <- renderPlotly( ggplotly(anno_tsne_p(), tooltip="Sample", dynamicTicks=T) )
  
  levelsSelectedAnnot_anno_tSNE <- reactive({
    clust$annot_sel2 %>% pull(input$anno_tsne_color) %>% unique() %>% as.vector()
  })
  
  output$anno_tsne_box <- renderUI({
    if(!is.null(clust$tsne_corr)){
        box(title="Annotated tSNE", width=NULL, status="success", solidHeader=T,
            column(6, align="left", selectInput("anno_tsne_color", "Color by", choices=c('sample_id', 'total_counts', 'ChromatinGroup'))),
            column(6, align="left", selectInput("anno_tsne_anno", "Labels", choices=c('none', 'barcode', 'cell_id', 'sample_id', 'total_counts'))),
            column(12, align="left", plotlyOutput("anno_tsne_plot")))
      
    }
  })
  
  output$anno_tsne_color_box <- renderUI({
    req(input$anno_tsne_color, clust$tsne_corr)
    if(input$anno_tsne_color != 'total_counts'){
      box(title="tSNE color settings", width = NULL, status = "success", solidHeader = T,
          column(6, htmlOutput("anno_tsne_colorPicker")),
          column(4 , br(), actionButton("anno_tSNEColReset", "Default colours", icon = icon("undo")),
                 br(), br(), actionButton("saveColors_anno_tsne", "Save colors & apply to all", icon = icon("save"))))
    }
  })
  
  cc.col <- reactive({
    if(!is.null(clust$cc.col) && !is.null(reactVal_tsneAnno$annotColors_filtered()) ){
      
      tmp_col = reactVal_tsneAnno$annotColors_filtered() %>% select(Sample,sample_id_Color,ChromatinGroup_Color)
      tmp = as.data.frame(clust$cc.col)
      names <-rownames(tmp)
      tmp[,"Sample"] = names
      tmp = left_join(tmp,tmp_col,by=c("Sample")) %>% select(ChromatinGroup_Color,sample_id_Color,total_counts)
      colnames(tmp) =c("ChromatinGroup","sample_id","total_counts")
      tmp = as.matrix(tmp)
      rownames(tmp)= names
      cc.col=tmp
      save(cc.col, file=file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_reduced_dataset, "cc.col.RData")))
      tmp
    }
  })

  observeEvent(input$saveColors_anno_tsne, {  # [anno_tsne] change modified colors for all plots where it applies and save RData locally

      cols <- paste0("c(", paste0("input$anno_tsne_col", levelsSelectedAnnot_anno_tSNE(), collapse = ", "), ")")
      cols <- eval(parse(text = cols))
      lev <- unique(levelsSelectedAnnot_anno_tSNE())
      names(cols) <- lev
      annot_Color_Custom <- as.data.frame(cols) %>% rownames_to_column(input$anno_tsne_color) %>% setNames(c(input$anno_tsne_color, str_interp("${input$anno_tsne_color}_Color") ))
      reactVal$annotColors_filtered <- isolate(reactVal_tsneAnno$annotColors_filtered()) %>% select(-c(str_interp("${input$anno_tsne_color}_Color"))) %>% right_join(annot_Color_Custom, by=input$anno_tsne_color  )
      annotColors_filtered_save <- reactVal$annotColors_filtered
      save(annotColors_filtered_save, file=file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_filtered_dataset,"_k_",input$nclust, "_annotColors_filtered.RData")))

  })
  
  observeEvent(input$anno_tSNEColReset, {
    lev <- unique(levelsSelectedAnnot_anno_tSNE())
    cols <- gg_fill_hue(length(lev))
    lapply(seq_along(lev), function(i) {
      do.call(what="updateColourInput", args=list(session=session, inputId=paste0("anno_tsne_col", lev[i]), value=cols[i]))
    })
  })
  
  output$anno_tsne_colorPicker <- renderUI({
   if(str_interp("${input$anno_tsne_color}_Color") %in% colnames(isolate(reactVal_tsneAnno$annotColors_filtered()))) {
    lev <- unique(levelsSelectedAnnot_anno_tSNE())

      colsModif <- isolate(reactVal_tsneAnno$annotColors_filtered()) %>% dplyr::select(one_of( str_interp("${input$anno_tsne_color}_Color")) ) %>% distinct() %>% pull( str_interp("${input$anno_tsne_color}_Color"))
        
      lapply(seq_along(lev), function(i) {
        colourpicker::colourInput(inputId=paste0("anno_tsne_col", lev[i]), label=paste0("Choose colour for ", lev[i]), value=colsModif[i], returnName=TRUE) ## Add ", palette = "limited"" to get selectable color panel       
      })
    }
  })
  
  anocol_sel_tsne <- reactive({ geco.annotToCol2(clust$annot_sel2[, c(annotCol(), "ChromatinGroup")], annotT=clust$annot_sel2, plotLegend=T, plotLegendFile=file.path(init$data_folder, "datasets", dataset_name(),"Annotation_legends_tsne_anno.pdf"), categCol=NULL) })
  
  annot_tsneAnno_Color <- reactive({
    df <- reactVal$annotColors_filtered[reactVal$annotColors_filtered$Sample %in% rownames(clust$annot_sel2), ]
    df
  })
  
  
  reactVal_tsneAnno <- reactiveValues(annotColors_filtered=reactive({annot_tsneAnno_Color()}) )#,annotColors_filtered_Chrom=reactive({annot_tsneAnno_Color_Chrom_Orig()}))
  
  observeEvent(input$selected_filtered_dataset ,{if(length(input$selected_filtered_dataset)>0 & !is.null(clust$available_k))
                   {unlocked$list$affectation =T}else{unlocked$list$affectation =F}})

  observeEvent(unlocked$list,{able_disable_tab(c("selected_reduced_dataset","affectation"),c("peak_calling","diff_analysis"))}) # if conditions are met, unlock tab Dimensionality Reduction
  
  
  ###############################################################
  # 5. Peak calling [optional]
  ###############################################################
  
  output$peak_calling_info <- renderText({"This module is optional, but recommended in order to obtain more precise results for enrichment analysis. 
    Based on the BAM files for the samples in your project, peaks will be called so that the counts can be mapped to the gene TSS more specifically."})
  output$pc_k_selection <- renderUI({
    selectInput("pc_k_selection", "Select number of clusters:", choices=clust$available_k())
  })
  output$bam_upload <- renderUI({
    sample_ids <- unique(annot_sel()$sample_id)
    lapply(sample_ids, function(x){
      textInput(paste0('bam_', x), paste0('Full path to .bam file for sample ', x, ':'), value="/home/example/sample.bam")
    })
  })
  
  observeEvent(input$do_pc, {
    if(is.null(input$pc_k_selection) | input$pc_k_selection == ''){
      showNotification("No cluster number selected. It seems that you didn't select a clustering on the previous page, please go back and check.", duration=7, closeButton=TRUE, type="warning")
    }else{
      myData = new.env()
      load(file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_filtered_dataset, '_affectation_k', input$pc_k_selection, ".RData")), envir=myData)
      load(file.path(init$data_folder, "datasets", dataset_name(), "reduced_data", paste0(input$selected_reduced_dataset, "_annotFeat.RData")), envir=myData)
      dir.create(file.path(init$data_folder, "datasets", dataset_name(), "peaks"), showWarnings=FALSE)
      dir.create(file.path(init$data_folder, "datasets", dataset_name(), "peaks", paste0(input$selected_filtered_dataset, "_k", input$pc_k_selection)), showWarnings=FALSE)
      odir <- file.path(init$data_folder, "datasets", dataset_name(), "peaks", paste0(input$selected_filtered_dataset, "_k", input$pc_k_selection))
      sample_ids <- unique(annot_sel()$sample_id)
      inputBams <- sapply(sample_ids, function(x){ input[[paste0('bam_', x)]] })
      checkBams <- sapply(inputBams, function(x){ if(file.exists(x)){ 0 } else{ showNotification(paste0("Could not find file ", x, ". Please make sure to give a full path including the file name."), duration=7, closeButton=TRUE, type="warning"); 1}  })
      stat.value <- if(input$pc_stat=="p.value") paste("-p", input$pc_stat_value) else paste("-q", input$pc_stat_value)
      if(sum(checkBams)==0){
        subsetBamAndCallPeaks(myData$affectation, myData$annotFeat, odir, inputBams, stat.value, annotation_id(),input$peak_distance_to_merge)
        pc$new <- Sys.time()
        updateActionButton(session, "do_pc", label="Finished successfully", icon = icon("check-circle"))
      }
    }
  })
  
  pc <- reactiveValues(new="")
  
  observeEvent({input$selected_filtered_dataset  # reset label on actionButtion when new peak calling should be performed
    input$pc_k_selection
    input$pc_stat
    input$pc_stat_value}, {
      pc$available_pc <- has_available_pc()
      updateActionButton(session, "do_pc", label="Start", icon=character(0))
  })
  
  has_available_pc <- reactive({
    print(paste0("new peaks done: ", pc$new))
    file.exists(file.path(init$data_folder, "datasets", dataset_name(), "peaks", paste0(input$selected_filtered_dataset, "_k", input$pc_k_selection), "pm.annot.gene.bed"))
      })
  available_pc_plots <- reactive({
    print(paste0("new peaks done: ", pc$new))
    fe <- sapply(c(1:input$pc_k_selection), function(i){file.exists(file.path(init$data_folder, "datasets", dataset_name(), "peaks", paste0(input$selected_filtered_dataset, "_k", input$pc_k_selection), paste0("C", i, "_model.r")))})
    which(fe==TRUE)
  })

  output$pc_plot_box <- renderUI({
    if(has_available_pc()){
      box(title="Peak calling visualization", width=NULL, status="success", solidHeader=T,
          column(8, align="left", selectInput("pc_cluster","Select cluster (only those shown for which plots are available):", choices=paste0("C", available_pc_plots()))),
          column(12, align="left",
                 plotOutput("peak_model_plot", height=500, width=500),
                 plotOutput("cross_corr_plot", height=500, width=500)))
    }
  })
  
  pc_data_p <- reactive({ req(has_available_pc()); as.numeric(unlist(read.table(file.path(init$data_folder, "datasets", dataset_name(), "peaks", paste0(input$selected_filtered_dataset, "_k", input$pc_k_selection), paste0(input$pc_cluster, "_data_p.txt"))))) })
  pc_data_m <- reactive({ req(has_available_pc()); as.numeric(unlist(read.table(file.path(init$data_folder, "datasets", dataset_name(), "peaks", paste0(input$selected_filtered_dataset, "_k", input$pc_k_selection), paste0(input$pc_cluster, "_data_m.txt"))))) })
  pc_data_xcorr <- reactive({ req(has_available_pc()); as.numeric(unlist(read.table(file.path(init$data_folder, "datasets", dataset_name(), "peaks", paste0(input$selected_filtered_dataset, "_k", input$pc_k_selection), paste0(input$pc_cluster, "_data_xcorr.txt"))))) })
  pc_data_ycorr <- reactive({ req(has_available_pc()); as.numeric(unlist(read.table(file.path(init$data_folder, "datasets", dataset_name(), "peaks", paste0(input$selected_filtered_dataset, "_k", input$pc_k_selection), paste0(input$pc_cluster, "_data_ycorr.txt"))))) })
  peak_model_p <- reactive({
    req(has_available_pc())
    x <- seq.int((length(pc_data_p())-1)/2*-1,(length(pc_data_p())-1)/2)
    plot(x,pc_data_p(),type='l',col=c('red'),main='Peak Model',xlab='Distance to the middle',ylab='Percentage')
    lines(x,pc_data_m(),col=c('blue'))
    legend('topleft',c('forward tags','reverse tags'),lty=c(1,1,1),col=c('red','blue'))
  })
  output$peak_model_plot <- renderPlot(peak_model_p())
  cross_corr_p <- reactive({
    req(has_available_pc())
    altd  <- c(297)
    plot(pc_data_xcorr(),pc_data_ycorr(),type='l',col=c('black'),main='Cross-Correlation',xlab='Lag between + and - tags',ylab='Correlation')
    abline(v=altd,lty=2,col=c('red'))
    legend('topleft','alternative lag(s)',lty=2,col='red')
    legend('topright','alt lag(s) : 297',bty='n')
  })
  output$cross_corr_plot <- renderPlot(cross_corr_p())
  
  
  ###############################################################
  # 6. Differential analysis
  ###############################################################
  
  output$diff_analysis_info <- renderText({"Differential analysis will be performed using the cluster assignment obtained on the Consensus clustering page. To use a different number of clusters, 
    go to this page and first perform the clustering, then select the preferred number of clusters in the box on the right in order to display and save the data. 
    It will then appear here for selection."})
  output$selected_k <- renderUI({ selectInput("selected_k", "Select number of clusters:", choices=clust$available_k()) })
  output$contrib_hist <- renderUI({ if(input$only_contrib_cells){ plotOutput("contrib_hist_p", height=250, width=500) }})
  output$contrib_hist_p <- renderPlot( contrib_hist_plot() )
  contrib_hist_plot <- reactive({
    maxCons <- tapply(icl()$itemConsensus$itemConsensus[icl()$itemConsensus$k==input$selected_k], icl()$itemConsensus$item[icl()$itemConsensus$k==input$selected_k], max)
    hist(maxCons, col="steelblue", breaks=80, main="Max cluster contribution per cell", xlab="", ylab="number of cells")
    abline(v=input$contrib_thresh, lwd=2, col="red", lty=2)
    legend("topleft", legend=c("cluster contribution threshold"), col=c("red"), lty=c(2), cex=0.8)
  })
  output$contrib_thresh <- renderUI({ if(input$only_contrib_cells){ sliderInput("contrib_thresh", "Select minimum cluster contribution for cells:", min=0.6, max=1, value=0.9, step=0.01) }})
  output$contrib_info <- renderUI({ if(input$only_contrib_cells){ textOutput("contrib_info_text") }})
  output$contrib_info_text <- renderText({
    total_cells <- length(unique(icl()$itemConsensus$item[icl()$itemConsensus$k==input$selected_k]))
    sel_cells <- length(unique(icl()$itemConsensus$item[icl()$itemConsensus$k==input$selected_k & icl()$itemConsensus$itemConsensus >= input$contrib_thresh]))
    paste("Selected top", sel_cells, "cells out of", total_cells)
  })
  
  affectation <- reactive({
    req(input$selected_filtered_dataset)
    myData = new.env()
    load(file.path(init$data_folder, "datasets", dataset_name(), "consclust", paste0(input$selected_filtered_dataset, '_affectation_k', input$selected_k, ".RData")), envir=myData)
    myData$affectation
  })
  affectation_filtered <- reactive({
    if(input$only_contrib_cells){
      cells <- unique(icl()$itemConsensus$item[icl()$itemConsensus$k==input$selected_k & icl()$itemConsensus$itemConsensus >= input$contrib_thresh])
      affectation()[cells,]
    }else{
      affectation()
    }
  })
  norm_mat <- reactive({
    myData = new.env()
    load(file.path(init$data_folder, "datasets", dataset_name(), "reduced_data", paste0(input$selected_reduced_dataset, "_normMat.RData")), envir=myData)
    myData$norm_mat
  })
  diff <- reactiveValues(my.res=NULL, summary=NULL, groups=NULL, refs=NULL)
  
  observeEvent({input$selected_filtered_dataset  # load saved diff analysis results if they exist for these parameters
    input$selected_k
    input$qval.th
    input$cdiff.th
    input$de_type}, {
      req(input$new_dataset_name, input$qval.th, input$cdiff.th, init$data_folder, dataset_name(), input$selected_filtered_dataset)
      filename <- file.path(init$data_folder, "datasets", dataset_name(), "supervised", paste0(input$selected_filtered_dataset, "_", input$selected_k, "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".RData"))
      if(file.exists(filename)){
        myData = new.env()
        load(filename, envir=myData)
        diff$my.res <- myData$my.res_save
        diff$summary <- myData$summary_save
        diff$groups <- myData$groups_save
        diff$refs <- myData$refs_save
      }
  })
  
  observeEvent(input$do_wilcox, {  # perform differential analysis based on wilcoxon test
    withProgress(message='Performing differential analysis...', value = 0, {
      incProgress(amount=0.1, detail=paste("initializing"))
     
      Counts <- norm_mat()[, rownames(affectation_filtered())]
      feature <- data.frame(ID=annotFeat()$ID, chr=annotFeat()$chr, start=annotFeat()$start, end=annotFeat()$end)
      
      if(input$de_type=="one_vs_rest"){ # compare each cluster to all the rest
        mygps <- lapply(1:input$selected_k, function(i){ affectation_filtered()[which(affectation_filtered()$ChromatinGroup==paste0("C", i)), "cell_id"]})
        names(mygps) <- paste0('C', 1:input$selected_k)
        groups <- names(mygps)
        myrefs <- lapply(1:input$selected_k, function(i){ affectation_filtered()[which(affectation_filtered()$ChromatinGroup!=paste0("C", i)), "cell_id"]})
        names(myrefs) <- paste0('notC', 1:input$selected_k)
        refs <- names(myrefs)
        incProgress(amount=0.5, detail=paste("performing wilcoxon rank tests"))
        diff$my.res <- geco.CompareWilcox(dataMat=Counts, annot=affectation_filtered(), ref=myrefs, groups=mygps, featureTab=feature)
        
      }else{ # pairwise one-vs-one testing for each cluster
        incProgress(amount=0.5, detail=paste("performing wilcoxon rank tests"))
        diff$my.res <- feature
        count_save <- data.frame(ID=feature$ID)
        single_results <- list()
        pairs <- setNames(data.frame(matrix(ncol=2, nrow=0)), c("group", "ref"))
        for(i in 1:(as.integer(input$selected_k)-1)){
          mygps <- list(affectation_filtered()[which(affectation_filtered()$ChromatinGroup==paste0("C", i)), "cell_id"])
          names(mygps) <- paste0('C', i)
          groups <- names(mygps)
          for(j in (i+1):as.integer(input$selected_k)){
            myrefs <- list(affectation_filtered()[which(affectation_filtered()$ChromatinGroup==paste0("C", j)), "cell_id"])
            names(myrefs) <- paste0('C', j)
            refs <- names(myrefs)
            tmp_result <- geco.CompareWilcox(dataMat=Counts, annot=affectation_filtered(), ref=myrefs, groups=mygps, featureTab=feature)
            tmp_result <- tmp_result[tmp_result$ID==count_save$ID,]
            rownames(tmp_result) <- tmp_result$ID
            tmp_result[5] <- NULL #remove rank because it will be added later
            colnames(tmp_result)[5:8] <- c("count", "cdiff", "p.val", "adj.p.val")
            count_save[, paste0('C', i)] <- tmp_result$count
            single_results <- list.append(single_results, tmp_result)
            pairs[nrow(pairs)+1,] <- list(groups[1], refs[1])
            tmp_mirror <- tmp_result
            tmp_mirror$cdiff <- tmp_mirror$cdiff*(-1.0)
            tmp_mirror$count <- 0 #not correct, but doesn't matter because it won't be used
            single_results <- list.append(single_results, tmp_mirror)
            pairs[nrow(pairs)+1,] <- list(refs[1], groups[1])
          }
        }
        #get count for last group as the loop doesn't cover it
        tmp.gp <- list(affectation_filtered()[which(affectation_filtered()$ChromatinGroup==paste0("C", input$selected_k)), "cell_id"])
        count_save[, paste0('C', input$selected_k)] <- apply(Counts, 1, function(x) mean(x[as.character(tmp.gp[[1]])]))
        combinedTests <- combineMarkers(de.lists=single_results, pairs=pairs, pval.field="p.val",
                                        effect.field="cdiff", pval.type="any", log.p.in=FALSE,
                                        log.p.out=FALSE, output.field="stats", full.stats=TRUE)
        for(i in 1:as.integer(input$selected_k)){
          cdiffs <- sapply(1:(as.integer(input$selected_k)-1), function(k){ combinedTests[[paste0("C", i)]][feature$ID, k+3]$cdiff })
          diff$my.res[, paste0("Rank.C", i)] <- combinedTests[[paste0("C", i)]][feature$ID, "Top"]
          diff$my.res[, paste0("Count.C", i)] <- as.numeric(count_save[, paste0("C", i)])
          diff$my.res[, paste0("cdiff.C", i)] <- rowMeans(cdiffs)
          diff$my.res[, paste0("pval.C", i)] <- combinedTests[[paste0("C", i)]][feature$ID, "p.value"]
          diff$my.res[, paste0("qval.C", i)] <- combinedTests[[paste0("C", i)]][feature$ID, "FDR"]
        }
        groups <- paste0('C', 1:input$selected_k) #needed for following code
        refs <- paste0('pairedTest', 1:input$selected_k)
      }
      
      incProgress(amount=0.3, detail=paste("building summary"))
      diff$summary <- matrix(nrow=3, ncol=length(groups), dimnames=list(c("differential", "over", "under"), groups))
      for(k in 1:length(groups)){
        gpsamp <- groups[k]
        
        #For log2(x1/x2) > 1 || log2(x1/x2) > -1
        diff$summary["differential", gpsamp] <- sum(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & abs(diff$my.res[, paste("cdiff", gpsamp, sep=".")]) > input$cdiff.th, na.rm=T)
        diff$summary["over", gpsamp] <- sum(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gpsamp, sep=".")] > input$cdiff.th, na.rm=T)
        diff$summary["under", gpsamp] <- sum(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gpsamp, sep=".")] < -input$cdiff.th, na.rm=T)


        #For mean(x1) - 2 mean(x2) > 0 || mean(x1) - 0.5 mean(x2) < 0
        # diff$summary["over", gpsamp] <- sum(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff1", gpsamp, sep=".")] > input$cdiff.th, na.rm=T)
        # diff$summary["under", gpsamp] <- sum(diff$my.res[, paste("qval", gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff2", gpsamp, sep=".")] < -input$cdiff.th, na.rm=T)
        # diff$summary["differential", gpsamp] <- diff$summary["over", gpsamp] + diff$summary["under", gpsamp]

        }
      diff$groups <- groups
      diff$refs <- refs
      my.res_save <- diff$my.res
      summary_save <- diff$summary
      groups_save <- diff$groups
      refs_save <- diff$refs
      save(my.res_save, summary_save, groups_save, refs_save, file=file.path(init$data_folder, "datasets", dataset_name(), "supervised", paste0(input$selected_filtered_dataset, "_", input$selected_k, "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".RData")))
      incProgress(amount=0.1, detail=paste("finished"))
    })
  })
  
  output$da_summary_box <- renderUI({
    if(!is.null(diff$groups)){
      box(title="Number of differentially bound regions", width=NULL, status="success", solidHeader=T,
          column(5, align="left", br(), tableOutput("da_summary_kable")),
          column(7, align="left", plotOutput("da_barplot", height=270, width=250)),
          column(12, align="left", div(style = 'overflow-x: scroll', DT::dataTableOutput('da_table')), br()),
          column(4, align="left", downloadButton("download_da_barplot", "Download barplot")),
          column(4, align="left", downloadButton("download_da_table", "Download table")))
    }
  })
  
  output$da_summary_kable <- function(){
    req(diff$summary)
    diff$summary %>%
      kable(escape=F, align="c") %>%
      kable_styling(c("striped", "condensed"), full_width = F)
  }
  
  output$da_barplot <- renderPlot({
    myylim <- range(c(diff$summary["over",], -diff$summary["under", ]))
    barplot(diff$summary["over",], col="red", las=1, ylim=myylim, main="Differentially bound regions",
            ylab="Number of regions", axes=F)
    barplot(-diff$summary["under", ], col="forestgreen", ylim=myylim, add=T, axes=F, names.arg="")
    z <- axis(2, pos=-10)
    axis(2, at=z, labels=abs(z), las=1)
  })
  
  output$download_da_barplot <- downloadHandler(
    filename=function(){ paste0("diffAnalysis_numRegions_barplot_", input$selected_filtered_dataset, "_", input$selected_k, "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".png")},
    content=function(file){
      png(file, width=500, height=400, res=300)
      myylim <- range(c(diff$summary["over",], -diff$summary["under", ]))
      barplot(diff$summary["over",], col="red", las=1, ylim=myylim, main="Differentially bound regions",
              ylab="Number of regions", axes=F)
      barplot(-diff$summary["under", ], col="forestgreen", ylim=myylim, add=T, axes=F, names.arg="")
      z <- axis(2, pos=-10)
      axis(2, at=z, labels=abs(z), las=1)
      dev.off()
    })
  
  output$da_table <- DT::renderDataTable({
    table <- diff$my.res[, -c(1)]
    rownames(table) <- NULL
    for(i in seq(from=0, to=(dim(table)[2]-8)/5, by=1)){
      table[, (5*i+5)] <- round(table[, (5*i+5)], 3) #counts
      table[, (5*i+6)] <- round(table[, (5*i+6)], 3) #cdiff
        #table[, (5*i+7)] <- round(table[, (5*i+7)], 8) #pval
        #table[, (5*i+8)] <- round(table[, (5*i+8)], 8) #qval
    }
    table <- table[order(table$Rank.C1),]
    DT::datatable(table, options=list(dom='tpi'))
  })
  
  output$download_da_table <- downloadHandler(
    filename=function(){ paste0("diffAnalysis_data_", input$selected_filtered_dataset, "_", input$selected_k, "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, ".csv")},
    content=function(file){
      write.table(diff$my.res, file, row.names=F, quote=F, sep=",")
    })
  
  output$da_visu_box <- renderUI({
    if(!is.null(diff$groups)){
      box(title="Visualization", width=NULL, status="success", solidHeader=T,
          column(4, align="left", selectInput("gpsamp", "Select cluster:", choices=diff$groups)),
          column(12, align="left", plotOutput("h1_prop", height=300, width=500),
                 plotOutput("da_volcano", height=500, width=500)),
          column(4, align="left", downloadButton("download_h1_plot", "Download histogram")),
          column(4, align="left", downloadButton("download_da_volcano", "Download volcano plot")))
    }
  })
  
  output$h1_prop <- renderPlot({
    tmp <- geco.H1proportion(diff$my.res[, paste("pval", input$gpsamp, sep=".")])
    hist(diff$my.res[, paste("pval", input$gpsamp, sep=".")], breaks=seq(0, 1, by=0.05), xlab="P-value",
         ylab="Frequency", main=paste(input$gpsamp, "vs the rest", "\n", "H1 proportion:", round(tmp, 3)))
  })
  
  output$download_h1_plot <- downloadHandler(
    filename=function(){ paste0("diffAnalysis_numRegions_barplot_", input$selected_filtered_dataset, "_", input$selected_k, "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, "_", input$gpsamp, ".png")},
    content=function(file){
      png(file, width=1000, height=600, res=300)
      tmp <- geco.H1proportion(diff$my.res[, paste("pval", input$gpsamp, sep=".")])
      hist(diff$my.res[, paste("pval", input$gpsamp, sep=".")], breaks=seq(0, 1, by=0.05), xlab="P-value",
           ylab="Frequency", main=paste(input$gpsamp, "vs the rest", "\n", "H1 proportion:", round(tmp, 3)))
      dev.off()
  })
  
  output$da_volcano <- renderPlot({
    mycol <- rep("black", nrow(diff$my.res))
    mycol[which(diff$my.res[, paste("qval", input$gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", input$gpsamp, sep=".")] > input$cdiff.th)] <- "red"
    mycol[which(diff$my.res[, paste("qval", input$gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", input$gpsamp, sep=".")] < -input$cdiff.th)] <- "forestgreen"
    plot(diff$my.res[, paste("cdiff", input$gpsamp, sep=".")], -log10(diff$my.res[, paste("qval", input$gpsamp, sep=".")]),
         col=mycol, cex=0.7, pch=16, xlab="count difference", ylab="-log10(adjusted p-value)", las=1,
         main=paste(input$gpsamp, "vs the rest","\n", diff$summary["over",input$gpsamp], "enriched,", diff$summary["under", input$gpsamp], "depleted"))
    abline(v=input$cdiff.th, lty=2)
    abline(h=-log10(input$qval.th), lty=2)
    abline(v=-input$cdiff.th, lty=2)
  })
  
  output$download_da_volcano <- downloadHandler(
    filename=function(){ paste0("diffAnalysis_numRegions_barplot_", input$selected_filtered_dataset, "_", input$selected_k, "_", input$qval.th, "_", input$cdiff.th, "_", input$de_type, "_", input$gpsamp, ".png")},
        content=function(file){
        png(file, width=900, height=900, res=300)
        mycol <- rep("black", nrow(diff$my.res))
        mycol[which(diff$my.res[, paste("qval", input$gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", input$gpsamp, sep=".")] > input$cdiff.th)] <- "red"
        mycol[which(diff$my.res[, paste("qval", input$gpsamp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", input$gpsamp, sep=".")] < -input$cdiff.th)] <- "forestgreen"
        plot(diff$my.res[, paste("cdiff", input$gpsamp, sep=".")], -log10(diff$my.res[, paste("qval", input$gpsamp, sep=".")]),
             col=mycol, cex=0.7, pch=16, xlab="count difference", ylab="-log10(adjusted p-value)", las=1,
             main=paste(input$gpsamp, "vs the rest","\n", diff$summary["over",input$gpsamp], "enriched,", diff$summary["under", input$gpsamp], "depleted"))
        abline(v=input$cdiff.th, lty=2)
        abline(h=-log10(input$qval.th), lty=2)
        abline(v=-input$cdiff.th, lty=2)
        dev.off()
    })
    observeEvent(unlocked$list,{able_disable_tab(c("selected_reduced_dataset","diff_my_res"),"enrich_analysis")}) # if conditions are met, unlock tab Dimensionality Reduction
    observeEvent(diff$my.res,{if(dim(diff$my.res)[1]>0){unlocked$list$diff_my_res=TRUE}else{unlocked$list$diff_my_res=FALSE}}) # if conditions are met, unlock tab Dimensionality Reduction
    
  
  ###############################################################
  # 7. Enrichment analysis
  ###############################################################
  
  msig_db <- reactive({
    myData = new.env()
    load(file.path('annotation', annotation_id(), 'MSigDB.RData'), envir=myData)
    myData
  })
  MSIG.ls <- reactive({ msig_db()$MSIG.ls })
  MSIG.gs <- reactive({ msig_db()$MSIG.gs })
  MSIG.classes <- reactive({unique(MSIG.gs()$Class)})
  
  annotFeat <- reactive({
    myData = new.env()
    load(file.path(init$data_folder, 'datasets', dataset_name(), "reduced_data", paste0(input$selected_reduced_dataset, "_annotFeat.RData")), envir=myData)
    myData$annotFeat
  })
  annotFeat_long <- reactive({ as.data.frame(cSplit(annotFeat(), splitCols="Gene", sep=", ", direction="long")) })
  
  output$enr_info <- renderText({"Enrichment will be performed based on the significant genes per cluster that were computed on the previous page. 
    Please make sure that you have run the differential analysis on the clustering that you prefer before running the enrichment analysis."})
  canUsePeaks <- reactive({
    print("File peak calling exist ?")
    print(paste(init$data_folder, 'datasets', dataset_name(), 'peaks', paste0(input$selected_filtered_dataset, '_k', input$selected_k), 'pm.annot.gene.bed',sep="/"))
    print(file.exists(file.path(init$data_folder, 'datasets', dataset_name(), 'peaks', paste0(input$selected_filtered_dataset, '_k', input$selected_k), 'pm.annot.gene.bed')))
    file.exists(file.path(init$data_folder, 'datasets', dataset_name(), 'peaks', paste0(input$selected_filtered_dataset, '_k', input$selected_k), 'pm.annot.gene.bed')) })
  
  output$use_peaks <- renderUI({
    if(canUsePeaks()){ checkboxInput("use_peaks", "use peak calling results to refine analysis", value=TRUE) }
  })
  enr <- reactiveValues(Both=NULL, Overexpressed=NULL, Underexpressed=NULL)
  
  peak_window <- reactive({
    req(input$use_peaks==TRUE)
    pm.annot.window.file <- file.path(init$data_folder, 'datasets', dataset_name(), 'peaks', paste0(input$selected_filtered_dataset, '_k', input$selected_k), 'pm.annot.window.bed')
    #pm.annot.window.file <- "/home/depic2/Downloads/BC976_mm10_peak_merge.annot.window.bed"
    peak_window <- read.table(pm.annot.window.file, sep="\t", header=F)
    colnames(peak_window) <- c("chr", "start", "end", "w_chr", "w_start", "w_end")
    peak_window$peak_ID <- paste(peak_window$chr, peak_window$start, peak_window$end, sep="_")
    peak_window$window_ID <- paste(peak_window$w_chr, peak_window$w_start, peak_window$w_end, sep="_")
    peak_window
  })
  peak_gene <- reactive({
    req(input$use_peaks==TRUE)
    pm.annot.gene.file <- file.path(init$data_folder, 'datasets', dataset_name(), 'peaks', paste0(input$selected_filtered_dataset, '_k', input$selected_k), 'pm.annot.gene.bed')
    #pm.annot.gene.file <- "/home/depic2/Downloads/BC976_mm10_peak_merge.annot.gene.bed"
    peak_gene <- read.table(pm.annot.gene.file, sep="\t", header=F)
    peak_gene <- peak_gene[, c(1:3, 7:8)]
    colnames(peak_gene) <- c("chr","start", "end", "Gene_name", "dist_TSS")
    peak_gene$peak_ID <- paste(peak_gene$chr, peak_gene$start, peak_gene$end, sep="_")
    peak_gene
  })
  
  GencodeGenes <- reactive({
    Gencode <- read.table(file.path('annotation', annotation_id(), 'Gencode_TSS_pc_lincRNA_antisense.bed'))
    unique(Gencode$V4)
  })
  
  observeEvent(input$do_enrich, {
    withProgress(message='Performing enrichment analysis...', value = 0, {

      incProgress(amount=0.1, detail=paste("initializing"))
      unlink(file.path(".", "*.csv"))
      Both <- list()
      Overexpressed <- list()
      Underexpressed <- list()
      for(i in 1:length(diff$groups)){
        gp <- diff$groups[i]
        incProgress(amount=(0.8/length(diff$groups)), detail=paste("working on cluster", gp))
        ref <- diff$refs[i]
        signific <- diff$my.res$ID[which(diff$my.res[, paste("qval", gp, sep=".")] <= input$qval.th & abs(diff$my.res[, paste("cdiff", gp, sep=".")]) > input$cdiff.th)]
        significG <- unique(annotFeat_long()$Gene[annotFeat_long()$distance < 1000 & annotFeat_long()$ID %in% signific])
        over <- diff$my.res$ID[which(diff$my.res[, paste("qval", gp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gp, sep=".")] > input$cdiff.th)]
        overG <- unique(annotFeat_long()$Gene[annotFeat_long()$distance < 1000 & annotFeat_long()$ID %in% over])
        under <- diff$my.res$ID[which(diff$my.res[, paste("qval", gp, sep=".")] <= input$qval.th & diff$my.res[, paste("cdiff", gp, sep=".")] < -input$cdiff.th)]
        underG <- unique(annotFeat_long()$Gene[annotFeat_long()$distance < 1000 & annotFeat_long()$ID %in% under])

        if( !is.null(input$use_peaks)){
          if(input$use_peaks ==TRUE){
          signific_associated_peak <- peak_window()[peak_window()$window_ID %in% signific, 8]
          over_associated_peak <- peak_window()[peak_window()$window_ID %in% over, 8]
          under_associated_peak <- peak_window()[peak_window()$window_ID %in% under, 8]
          signific_associated_gene <- peak_gene()[peak_gene()$peak_ID %in% signific_associated_peak & peak_gene()$dist_TSS<1000, ]
          over_associated_gene <- peak_gene()[peak_gene()$peak_ID %in% over_associated_peak & peak_gene()$dist_TSS<1000, ] 
          under_associated_gene <- peak_gene()[peak_gene()$peak_ID %in% under_associated_peak & peak_gene()$dist_TSS<1000, ]
          significG <- unique(signific_associated_gene$Gene_name)
          overG <- unique(over_associated_gene$Gene_name)
          underG <- unique(under_associated_gene$Gene_name)
          }
        }
        
        if(length(significG)){
          enrich.test <- geco.enrichmentTest(gene.sets=MSIG.ls(), mylist=significG, possibleIds=GencodeGenes())
          enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
          enrich.test <- merge(subset(MSIG.gs(), select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
          enrich.test <- enrich.test[order(enrich.test$`p-value`),]
          ind <- which(enrich.test$`q-value`<= 0.1)
          if(!length(ind)){ind <- 1:10}
          Both[[i]] <- enrich.test[ind,]
        }
        if(length(overG)){
          enrich.test <- geco.enrichmentTest(gene.sets=MSIG.ls(), mylist=overG, possibleIds=GencodeGenes())
          enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
          enrich.test <- merge( subset(MSIG.gs(), select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
          enrich.test <- enrich.test[order(enrich.test$`p-value`),]
          ind <- which(enrich.test$`q-value`<= 0.1)
          if(!length(ind)){ind <- 1:10}
          Overexpressed[[i]]  <- enrich.test[ind,]
        }
        if(length(underG)){
          enrich.test <- geco.enrichmentTest(gene.sets=MSIG.ls(), mylist=underG, possibleIds=GencodeGenes())
          enrich.test <- data.frame(Gene_set_name=rownames(enrich.test), enrich.test, check.names=FALSE)
          enrich.test <- merge( subset(MSIG.gs(), select=-Genes), enrich.test, by.x="Gene.Set", by.y="Gene_set_name", all.y=TRUE, sort=FALSE ) ## Get class of gene set
          enrich.test <- enrich.test[order(enrich.test$`p-value`),]
          ind <- which(enrich.test$`q-value`<= 0.1)
          if(!length(ind)){ind <- 1:10}
          Underexpressed[[i]] <- enrich.test[ind,]
        }
      }
      enr$Both <- Both
      enr$Overexpressed <- Overexpressed
      enr$Underexpressed <- Underexpressed
      incProgress(amount=0.1, detail=paste("finished"))
    })
  })
  
  output$enr_clust_sel <- renderUI({ selectInput("enr_clust_sel", "Select cluster:", choices=diff$groups) })
  output$enr_class_sel <- renderUI({shiny::checkboxGroupInput(inputId = "enr_class_sel",inline = T,label =  "Select classes to display:",selected = MSIG.classes() ,choiceNames=MSIG.classes(), choiceValues=MSIG.classes())})
  
  output$all_enrich_table <- DT::renderDataTable({
    req(enr$Both,input$enr_class_sel)
    table <- enr$Both[[match(input$enr_clust_sel, diff$groups)]]
    table <- table[which(table[,"Class"] %in% input$enr_class_sel),]
    if(is.null(table)){ return(setNames(data.frame(matrix(ncol=6, nrow=0)), c("Gene_set","Class", "Num_deregulated_genes", "p.value", "q.value", "Deregulated_genes"))) }
    table <- unite(table, "dereg_genes", c("Nb_of_deregulated_genes", "Nb_of_genes"), sep="/")
    colnames(table) <- c("Gene_set","Class", "Num_deregulated_genes", "p.value", "adj.p.value", "Deregulated_genes")
    table[, 4] <- round(table[, 4], 9)
    table[, 5] <- round(table[, 5], 9)
    table <- table[order(table$adj.p.value, table$p.value),]
    DT::datatable(table, options=list(dom='tpi'))
  })
  
  output$over_enrich_table <- DT::renderDataTable({
    req(enr$Overexpressed,input$enr_class_sel)
    table <- enr$Overexpressed[[match(input$enr_clust_sel, diff$groups)]]
    table <- table[which(table[,"Class"] %in% input$enr_class_sel),]
    if(is.null(table)){ return(setNames(data.frame(matrix(ncol=6, nrow=0)), c("Gene_set","Class", "Num_deregulated_genes", "p.value", "q.value", "Deregulated_genes"))) }
    table <- unite(table, "dereg_genes", c("Nb_of_deregulated_genes", "Nb_of_genes"), sep="/")
    colnames(table) <- c("Gene_set","Class", "Num_deregulated_genes", "p.value", "adj.p.value", "Deregulated_genes")
    table[, 4] <- round(table[, 4], 9)
    table[, 5] <- round(table[, 5], 9)
    table <- table[order(table$adj.p.value, table$p.value),]
    DT::datatable(table, options=list(dom='tpi'))
  })
  
  output$under_enrich_table <- DT::renderDataTable({
    req(enr$Underexpressed,input$enr_class_sel)
    table <- enr$Underexpressed[[match(input$enr_clust_sel, diff$groups)]]
    table <- table[which(table[,"Class"] %in% input$enr_class_sel),]
    if(is.null(table)){ return(setNames(data.frame(matrix(ncol=6, nrow=0)), c("Gene_set","Class", "Num_deregulated_genes", "p.value", "q.value", "Deregulated_genes"))) }
    table <- unite(table, "dereg_genes", c("Nb_of_deregulated_genes", "Nb_of_genes"), sep="/")
    colnames(table) <- c("Gene_set","Class", "Num_deregulated_genes", "p.value", "adj.p.value", "Deregulated_genes")
    table[, 4] <- round(table[, 4], 9)
    table[, 5] <- round(table[, 5], 9)
    table <- table[order(table$adj.p.value, table$p.value),]
    DT::datatable(table, options=list(dom='tpi'))
  })
  
  output$download_enr_data <- downloadHandler(
    filename = paste(input$selected_filtered_dataset, input$selected_k, input$qval.th, input$cdiff.th, "enrichment_tables.zip", sep="_"),
    content = function(fname){
      fs <- c()
      for(i in 1:length(diff$groups)){
        if(!is.null(enr$Both[[i]])){
          filename <- paste0(diff$groups[i], "_significant_gene_sets.csv")
          fs <- c(fs, filename)
          write.table(enr$Both[[i]], file=filename, quote=FALSE, row.names=FALSE, sep=",")
        }
        if(!is.null(enr$Overexpressed[[i]])){
          filename <- paste0(diff$groups[i], "_enriched_gene_sets.csv")
          fs <- c(fs, filename)
          write.table(enr$Overexpressed[[i]], file=filename, quote=FALSE, row.names=FALSE, sep=",")
        }
        if(!is.null(enr$Underexpressed[[i]])){
          filename <- paste0(diff$groups[i], "_depleted_gene_sets.csv")
          fs <- c(fs, filename)
          write.table(enr$Underexpressed[[i]], file=filename, quote=FALSE, row.names=FALSE, sep=",")
        }
      }
      zip(zipfile=fname, files=fs)},
    contentType = "application/zip"
  )
  
  output$gene_sel <- renderUI({
    req(enr$Both, enr$Overexpressed, enr$Underexpressed)
    
    most_diff = diff$my.res %>% select(ID,starts_with("qval."))
    most_diff[,"qval"]= rowMeans2(as.matrix(most_diff[,-1]))
    most_diff = left_join(most_diff[order(most_diff$qval),],annotFeat_long(),by = c("ID") )
    most_diff = most_diff %>% filter(!is.na(Gene)) 
    genes <- intersect(most_diff$Gene,unique(GencodeGenes()))
    
    # genes = genes[which(genes %in% peak_gene()$Gene_name)]
    selectizeInput(inputId = "gene_sel", "Select gene:",options= list(maxOptions = 250),genes)
  })
  
  output$region_sel <- renderUI({
    req(input$gene_sel)

    subset <- annotFeat_long()[annotFeat_long()$Gene==input$gene_sel, ]
    subset <- subset[order(subset$distance),]
    regions <- paste0(subset$ID, " (distance to gene TSS: ", subset$distance, ")")
    selectInput("region_sel", "Select associated genomic region:", choices=regions)
  })
  
  enr_tsne_res <- reactive({
    tsne_res <- NULL
    withProgress(message='Performing tSNE...', value = 0, {
      incProgress(amount=0.5, detail=paste("running"))
      tsne_res <- Rtsne(t(mati2()[, rownames(affectation_filtered())]), dims=2, pca=FALSE, theta=0.0, perplexity=choose_perplexity(t(mati2()[, rownames(affectation_filtered())])), verbose=TRUE, max_iter=600)
      incProgress(amount=0.5, detail=paste("finished"))
    })
    tsne_res
  })
  
  gene_tsne_p <- reactive({
    req(input$gene_sel, input$region_sel)
    region <- strsplit(input$region_sel, " ")[[1]][1]
    if(region %in% rownames(norm_mat())){
      p <- ggplot(as.data.frame(enr_tsne_res()$Y), aes(x=V1, y=V2)) +
        geom_point(alpha=0.5, aes(color=norm_mat()[region, rownames(affectation_filtered())], shape=affectation_filtered()$ChromatinGroup)) +
        labs(color="norm. count for region", shape="Cluster", x="t-SNE 1", y="t-SNE 2") +
        theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              panel.background=element_blank(), axis.line=element_line(colour="black"),
              panel.border=element_rect(colour="black", fill=NA)) +
        scale_color_viridis(direction=-1)
    }
  })
  output$gene_tsne_plot <- renderPlotly({
      req(input$gene_sel, input$region_sel)
      region <- strsplit(input$region_sel, " ")[[1]][1]
      if(region %in% rownames(norm_mat())){
      ggplotly(gene_tsne_p(), tooltip="Sample", dynamicTicks=T)
      }
  })
  
  
  ###############################################################
  # 8. Close app
  ###############################################################
  
  observeEvent(input$close, {
    unlink(file.path("www", "images", "*"))
    unlink(file.path(".", "*.csv"))
    js$closeWindow()
    stopApp()
  })
  
}