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
            element: '#new_analysis_name',
            intro:'Give a name to your analysis & select a reference genome (hg38 or mm10).',
            position:'auto'},
          {
            element:document.getElementsByClassName("input-group-btn")[0],
            intro:'Browse your computer to upload your data. ChromSCape accepts a wide variety of format, either processed (count matrices / Sparse matrices) or raw (fragment file, single-cell BED or BAM files). For the latter 3, place the files of each condition in a separate folder, and places all those folders in the same directory.',position:'auto'},
          {
            element:'#data_choice_box',
            intro:'Upload Dense count matrix, Sparse matrix (10X format), Fragment File raw data (10X format, requires "Seurat" package installed), Single-cell BAM or BED raw data (one file per cell). The latter 3 options will open a menu to enable counting on feature of your choice, e.g. genomic bins, peaks or genes.',
            position:'auto'},
            {
            element:'#create_analysis',
            intro:'Click here to create analysis & upload the selected matrices.',
            position:'auto'},
            {
            element:'#add_to_current_analysis',
            intro:'Add multiple layers of features in your object, e.g. genomic bins, peaks or genes on the same cells by ticking this box. Switch between features by clicking on top of the application ("Features" dropdown menu) that will appear once you loaded multi-feature analysis.',
            position:'auto'},
          {
            element:'#selected_analysis',
            intro:'All the analysis you created in the output directory will be available here.',
            position:'auto'},
          {
            element: tabs[1],
            intro: "<b>Filter & Normalize:</b> <br> On this tab, you will apply various filters to remove lowly covered cells, also normalize and reduce dimensions of your dataset.",
            position: 'right'
          },
          {
            element: tabs[2],
            intro: "<b>Visualize Cells:</b> In this tab, you can visualize cells dispersion in reduced feature space by methods such as PCA, UMAP or TSNE.<br> ",
            position: 'right'
          },
          {
            element: tabs[3],
            intro: "<b>Cluster Cells:</b> <br> Cluster your cells in an unsupervised manner to find the optimal number of epigenomic subpopulations. Optionally filter lowly correlated cells before clustering.",
            position: 'right'
          },
          {
            element: tabs[4],
            intro: "<b>Coverage:</b> <br> Generate coverage tracks for each one of your cluster and browse for genes / loci of intereset in this tab. Coverage tracks must be generated from raw data (scBED).",
            position: 'right'
          },
          {
            element: tabs[5],
            intro: "<b>Peak Call:</b> <br> The peak calling & coverage tabs will refine the signal based on cell clusters defined in previous step. It will also creates pseudo-bulk coverage of each cluster for you to explore the intra-cluster heterogeneity anywhere in the genome !.",
            position: 'right'
          },
          {
            element: tabs[6],
            intro: "<b>Differential Analysis: </b> <br> Identify the most differential features between the different epigenomic cell populations.",
            position: 'right'
          },
          {
            element: tabs[7],
            intro: "<b>Enrichment Analysis: </b> <br> Find significant gene sets associated to differential features and vizualize signal of features in cells.",
            position: 'right'
          },
          {
            intro: "You are good to go ! Take a look at the <a href='https://vallotlab.github.io/ChromSCape/ChromSCape_guide.html' target='_blank'>User Guide</a> for more details, and please report issues at <a href='https://github.com/vallotlab/ChromSCape'  target='_blank'>ChromSCape Github</a>.",
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