set_ChromSCape_options <- function() {
    if(is.null(getOption("ChromSCape_options"))){
        
    ChromSCape_options = list(
        heatmap_downsample = 1000,
        dotplot_downsample = 5000,
        dotplot_transparency = 0.6,
        dotplot_size = 1,
        dotplot_max_distanceToTSS = 1000,
        dotplot_min_quantile = 0.01,
        dotplot_max_quantile = 0.99
    )
    options("ChromSCape_options" = ChromSCape_options)
    } else{
        ChromSCape_options. = list(
            heatmap_downsample = 1000,
            dotplot_downsample = 5000,
            dotplot_transparency = 0.6,
            dotplot_size = 1,
            dotplot_max_distanceToTSS = 1000,
            dotplot_min_quantile = 0.01,
            dotplot_max_quantile = 0.99
        )
        
        ChromSCape_options = getOption("ChromSCape_options")
        for(item in names(ChromSCape_options.)){
            if(!item %in% names(ChromSCape_options)){
                ChromSCape_options[[item]] = ChromSCape_options.[[item]]
            }
        }
        options("ChromSCape_options" = ChromSCape_options)
    }
}