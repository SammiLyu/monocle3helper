#' Propagate cell annotations from a subset cds to a full cds
#'
#' Copies cell-level annotations from a subset \code{cell_data_set}
#' back into a full \code{cell_data_set}, matching cells by column names.
#'
#'
#' @param cds_full Full Monocle3 \code{cell_data_set}.
#' @param cds_sub Subset \code{cell_data_set} containing refined annotations.
#' @param from_col Column name in \code{cds_sub} to copy from.
#' @param to_col Column name in \code{cds_full} to write to.
#' @param overwrite Logical; overwrite existing values in \code{to_col}.
#'
#' @return Updated \code{cell_data_set} with propagated annotations.
#'
#'#' @importFrom SummarizedExperiment colData
#' @export
#'
#' @examples
#' \dontrun{
#' cds <- propagate_annotations(
#'   cds_full = cds,
#'   cds_sub  = cds_epi,
#'   from_col = "epi_annotation"
#' )
#' }
propagate_annotations <- function(
    cds_full,
    cds_sub,
    from_col,
    to_col = from_col,
    overwrite = TRUE
) {
  
  common_cells <- intersect(
    colnames(cds_full),
    colnames(cds_sub)
  )
  
  if (length(common_cells) == 0) {
    stop("No overlapping cell IDs between cds_full and cds_sub")
  }
  
  if (!to_col %in% colnames(colData(cds_full))) {
    colData(cds_full)[[to_col]] <- NA
  }
  
  colData(cds_full)[[to_col]] <- as.character(colData(cds_full)[[to_col]])
  sub_vals <- as.character(colData(cds_sub)[[from_col]])
  names(sub_vals) <- colnames(cds_sub)
  
  if (overwrite) {
    colData(cds_full)[common_cells, to_col] <- sub_vals[common_cells]
  } else {
    idx <- common_cells[is.na(colData(cds_full)[common_cells, to_col])]
    colData(cds_full)[idx, to_col] <- sub_vals[idx]
  }
  
  return(cds_full)
}