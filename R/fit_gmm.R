#' @title Fit a Gaussian mixture model given the number of outliers and a ranking of 
#' the candidate outliers.
#'
#' @inheritParams outcast_gmm
#' @inheritParams outback_gmm
#' @param outlier_num Number of outliers.
#' @param init_z Initial assignment probability matrix.
#'
#' @return List:
#' * $labels
#' * $mix
#' 
#' @export
fit_gmm <- function(
    x,
    outlier_rank, 
    outlier_num, 
    comp_num, 
    mnames = "VVV", 
    seed = 123,
    init_z = NULL
  ) {
  
  x0 <- as.matrix(x)
  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0
  
  labels <- rep(0, nrow(x0))
  
  if (is.null(init_z)) {
    init_z <- init_kmpp(x0[!outlier_bool], comp_num, seed)
  }
  
  mix <- mixture::gpcm(
    x0[!outlier_bool],
    G = comp_num,
    mnames = mnames,
    start = init_z,
    seed = seed
  )
  
  labels[!outlier_bool] <- mix$map
  
  return(list(labels = labels, mix = mix))
}