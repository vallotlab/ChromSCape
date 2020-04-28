// initialize an introjs instance          
var intro = introJs();

$(document).ready(function() {
  tabs = $('a[data-toggle="tab"]')
    intro.setOptions({
        steps: [
          {
            intro: 'Welcome to ChromSCape ! This application allows you to analyze and vizualize single cell epigenomic data (scChIP-seq, scATAC-seq, …). Follow the steps to learn how to use ChromSCape.'
            },
          {
            element: '#data_folder',
          intro: 'The first thing you need to do is to define a location where your analysis will be saved. This also allows to re-load an analysis that you might want to vizualize or try some other parameters.'
            
          },
          {
            element: '#new_dataset_name',
            intro:'Give a name to your analysis & select a reference genome (hg38 or mm10).',
            position:'auto'},
          {
            element:'#datafile_matrix_progress',
            intro:'Browse your computer to upload one or multiple matrice(s). The matrices should be tab-separated and in .tsv or .txt format. The first column must be cell-barcode names, the first row must be  genomic location (chr:start-end) or gene name. The matrices must all be placed in the same folder and uploaded simultaneously.',position:'auto'},{step:5,element:'#compile_dataset',
            intro:'Click here to load & compile your dataset.',
            position:'auto'},
          {
            element:'#selected_raw_dataset',
            intro:'All compiled raw dataset will appear here.',
            position:'auto'},
          {
            element:'#table_QC_filt_box',
            intro:'Adjust filtering parameters. If you try multiple combinations of parameters, all the results will be saved to allow you to select the optimal combination of parameter for your datasets.  You can first try with default parameters.',
            position:'auto'},
          {
            element:'#exclude_regions',
            intro:'If your samples contain known amplified / deleted regions, you might want to exclude some of the regions by uploading a BED file of those regions.',
            position:'auto'},
          {
            element:'#do_batch_corr',
            intro:'If you observe strong differences between known batch, you might want to specify to which batches each sample belongs and perform batch correction (using fastMNN algorithm from scran).',
            position:'auto'},
          {
            element:'#dim_reduction',
            intro:'Click here to perform QC filtering, and you are good to go !',
            position:'auto'
            
          },
          {
            element:'#selected_reduced_dataset',
            intro:'All the datasets you already filtered will be available here. Select the filtered dataset of your choice and continue to all the steps one by one.',
            position:'auto'
            
          },
          {
            element: tabs[1],
            intro: "<b>Dimensionality Reduction:</b> <br> On this tab, you will be able to vizualize your cells in different feature spaces: PCA, T-SNE or UMAP space.",
            position: 'right'
          },
          {
            element: tabs[2],
            intro: "<b>Correlation Clustering:</b> <br> Some of the cells are not correlated with enough other cells and will disturb clustering, filter them out in this tab.",
            position: 'right'
          },
          {
            element: tabs[3],
            intro: "<b>Consensus Clustering:</b> <br> Cluster your cells in an unsupervised manner to find the optimal number of epigenomic subpopulations.",
            position: 'right'
          },
          {
            element: tabs[4],
            intro: "<b>Peak Calling:</b> <br> This tab is reserved for advanced users on a Unix machine that have installed MACS2 & samtools. Using BAM sequencing files as input for peak calling of each cell cluster, this allows to refine the annotation of genes to genomic bins.",
            position: 'right'
          },
          {
            element: tabs[5],
            intro: "<b>Differential Analysis: </b> <br> Identify the most differential features between the different epigenomic cell populations.",
            position: 'right'
          },
          {
            element: tabs[6],
            intro: "<b>Enrichment Analysis: </b> <br> Find significant gene sets associated to differential features and vizualize signal of features in cells.",
            position: 'right'
          },
          {
            intro: "You are good to go ! Take a look at https://github.com/vallotlab/ChromSCape for more details or questions.",
            position: 'right'
          }
          ]
});

});

window.addEventListener('load', (event) => {
    console.log('The page has fully loaded');
});

window.addEventListener('load', function () {
    var doneTour = localStorage.getItem('EventTour') === 'Completed';
    console.log("doneTour")
    console.log(doneTour)
    console.log(intro)
    if (doneTour) {
        return;
    }
    else {  
      console.log("intro.start()")
        intro.start();
        intro.oncomplete(function () {
             console.log("Completed - complete");
            localStorage.setItem('EventTour', 'Completed');
        });

        intro.onexit(function () {
            console.log("Completed -exit");
            localStorage.setItem('EventTour', 'Completed');
        });
    }
});

// handler 1
Shiny.addCustomMessageHandler("setHelpContent",

  // callback function. 
  // note: message is passed by shiny and contains the tour data
  function(message){
    console.log("setHelpContent 2 ");
    // load data 
    intro.setOptions({steps: message.steps });
    
  }
  
);

// handler 2
Shiny.addCustomMessageHandler("startHelp",
  
  // callback function
  function(message) {
    console.log("startHelp 2 ");
    // start intro.js
    // note: we don't need information from shiny, just start introJS
    intro.start();
  }
  
);