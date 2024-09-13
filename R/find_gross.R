#' @title Find the gross outliers.
#' 
#' @inheritParams find_elbow
#' @inheritParams outcast_gmm
#' @param k_neighbours Number of neighbours for dbscan::glosh
#' 
#' @return List:
#' * $choice: a numeric value indicating the elbow's location.
#' * $bool: a logical vector identifying the gross outliers.
#' * $plot
#'
#' @export
find_gross <- function(
    x, max_out, search_centre, k_neighbours = 10, choice = NULL
  ) {
  
  x_glosh <- dbscan::glosh(x, k_neighbours)
  
  elbow <- find_elbow(-sort(-x_glosh)[seq_len(max_out)], search_centre, FALSE)
  
  if (is.null(choice)) {
    choice <- elbow$choice
  }
  
  bool <- rank(-x_glosh) <= choice
  
  return(list(choice = choice, bool = bool, plot = elbow$plot))
}

#' @export
plot_gross <- function(x, max_out, k_neighbours = 10) {
  
  x_glosh <- dbscan::glosh(x, k_neighbours)
  
  gg <- cbind(
    outliers = seq_len(max_out), glosh = -sort(-x_glosh)[seq_len(max_out)]
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = outliers, y = glosh)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = paste0(
        "k = ", k_neighbours
      )
    )
  
  return(gg)
}
  