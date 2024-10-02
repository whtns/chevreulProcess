#' Unite metadata
#'
#' @param object An object
#' @param group_bys A feature or variable to combine
#'
#' @return an object with Idents formed from concatenation of group_bys
#' @export
#'
#' @examples
#' 
#' 
#' data(small_example_dataset)
#' unite_metadata(small_example_dataset, "Mutation_Status")
#'
unite_metadata <-
    function(object, group_bys) {
        newcolname <- paste(group_bys, collapse = "_by_")
        newdata <- colData(object)[group_bys] |>
            as.data.frame() |>
            unite(!!newcolname, any_of(group_bys)) |>
            deframe()

        return(object)
    }


#' Plotly settings
#'
#' Change settings of a plotly plot
#'
#' @param plotly_plot  A plotly plot
#' @param width Default set to '600'
#' @param height Default set to '700'
#'
#' @noRd
plotly_settings <- function(plotly_plot, width = 600, height = 700) {
    plotly_plot |>
        layout(dragmode = "lasso") |>
        config(toImageButtonOptions = list(format = "svg", 
                                           filename = "myplot", 
                                           width = width, height = height)) |>
        identity()
}




#' Annotate Cell Cycle
#'
#' Annotate Cell Cycle for Gene and Transcript SingleCellExperiment Objects
#'
#' @param object A SingleCellExperiment object
#'
#' @return a SingleCellExperiment object
annotate_cell_cycle <- function(object) {
    
    data_env <- new.env(parent = emptyenv())
    data("cc.genes.cyclone", envir = data_env, package = "chevreul")
    cc.genes.cyclone <- data_env[["cc.genes.cyclone"]]
    
   assignments <- cyclone(object, cc.genes.cyclone, 
                          gene.names = rownames(object))
    colData(object)[colnames(assignments$scores)] <- assignments$scores
    colData(object)["Phase"] <- assignments$phases
    return(object)
}

