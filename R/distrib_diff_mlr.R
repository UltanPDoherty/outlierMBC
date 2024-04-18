#' distrib_diff_mlr
#'
#' @param student_resids Data.
#' @param df Component assignment probability matrix.
#'
#' @return List of
#' * distrib_diff
#' * sq_st_res
#' * choice_id
distrib_diff_mlr <- function(student_resids, df) {

  sq_st_res <- student_resids^2 / df

  eps <- 0.001
  check_seq <- seq(eps, 1 - eps, eps)
  checkpoints <- stats::qbeta(check_seq, 1 / 2, (df - 1) / 2)

  # sq_st_res_ecdf_func <- stats::ecdf(sq_st_res)
  sq_st_res_ecdf_func <- spatstat.geom::ewcdf(sq_st_res, 1)

  sq_st_res_ecdf <- sq_st_res_ecdf_func(checkpoints)
  beta_cdf <- stats::pbeta(checkpoints, 1 / 2, (df - 1) / 2)
  distrib_diff <- mean(abs(sq_st_res_ecdf - beta_cdf))

  choice_id <- which.max(sq_st_res)

  return(list(distrib_diff = distrib_diff,
              sq_st_res = sq_st_res,
              choice_id = choice_id))
}
