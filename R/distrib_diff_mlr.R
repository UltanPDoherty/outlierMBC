#' distrib_diff_mlr
#'
#' @param lm_mod Output from `lm` function.
#'
#' @return List of
#' * distrib_diff
#' * scsqst_res
#' * choice_id
distrib_diff_mlr <- function(lm_mod) {

  student_resids <- stats::rstandard(lm_mod)

  scsqst_res <- student_resids^2 / lm_mod$df.residual

  eps <- 0.001
  check_seq <- seq(eps, 1 - eps, eps)
  checkpoints <- stats::qbeta(check_seq, 1 / 2, (lm_mod$df.residual - 1) / 2)

  scsqst_res_ecdf_func <- spatstat.geom::ewcdf(scsqst_res, 1)

  scsqst_res_ecdf <- scsqst_res_ecdf_func(checkpoints)
  beta_cdf <- stats::pbeta(checkpoints, 1 / 2, (lm_mod$df.residual - 1) / 2)
  distrib_diff <- mean(abs(scsqst_res_ecdf - beta_cdf))

  choice_id <- which.min(dnorm(lm_mod$residuals, sd = stats::sigma(lm_mod)))

  return(list(distrib_diff = distrib_diff,
              scsqst_res = scsqst_res,
              choice_id = choice_id))
}
