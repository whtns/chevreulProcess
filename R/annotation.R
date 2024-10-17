#' Gene Symbols to Ensembl Transcript Ids
#'
#' convert hgnc gene symbols to ensembl transcript ids
#'
#' @param symbols character vector of gene symbols
#' @param organism mouse or human
#'
#' @return a vector of transcripts
#' @export
#' @importFrom GenomicFeatures transcripts
#'
#' @examples
#'
#' genes_to_transcripts("NRL")
genes_to_transcripts <- function(symbols, organism = "human") {

ensdb <- switch(organism, human = EnsDb.Hsapiens.v86, mouse = EnsDb.Mmusculus.v79)

    feature_table <- transcripts(ensdb,
            columns = c("gene_name", "gene_biotype", "gene_id"),
            return.type = "DataFrame"
        )

    feature_table[(feature_table[["gene_name"]] %in% symbols), "tx_id"]
}

#' Ensembl Transcript Ids to Gene Symbols
#'
#' Convert ensembl transcript ids to hgnc gene symbols
#'
#' @param transcripts transcripts
#' @param organism human or mouse
#'
#' @return a vector of gene symbols
#' @export
#'
#' @examples
#'
#' NRL_transcripts_hs <-
#' c("ENST00000359842", "ENST00000470566", "ENST00000465764")
#'
#'  data("grch38_tx2gene")
#'  data("grch38")
#' transcripts_to_genes(transcripts = NRL_transcripts_hs)
#'
transcripts_to_genes <- function(transcripts, organism = "human") {
    
    data_env <- new.env(parent = emptyenv())
    data("grch38", envir = data_env, package = "chevreul")
    data("grch38_tx2gene", envir = data_env, package = "chevreul")
    grch38 <- data_env[["grch38"]]
    grch38_tx2gene <- data_env[["grch38_tx2gene"]]

    gene_table <- switch(organism, human = grch38, mouse = grcm38)

    transcript_table = switch(organism, 
    human = grch38_tx2gene, 
    mouse = grcm38_tx2gene)

    tibble(enstxp = transcripts) |>
        left_join(transcript_table, by = "enstxp") |>
        left_join(gene_table, by = "ensgene") |>
        pull("symbol") |>
        identity()
}

#' Annotate percent mitochondrial reads per cell
#'
#'  Add a Percentage of Mitochondrial Read Count Categorical Variable to the
#'  Object (based on nCount_RNA)
#'
#' @param object A object
#' @param experiment gene
#'
#' @return a single cell object with
#' cell metadata column containing mitochondrial percentage
#' @export
#' @examples
#' 
#' data(small_example_dataset)
#' add_percent_mito(small_example_dataset)
#'
add_percent_mito <- function(object, experiment = "gene") {
    is.mito <- grepl("^MT-*", rownames(object))
    object <- addPerCellQCMetrics(object, subsets = list(Mito = is.mito))
    return(object)
}
