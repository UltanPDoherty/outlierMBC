#' ombc_lcwm
#'
#' @description
#' Iterative Detection & Identification of Outliers for a Linear
#' Cluster-Weighted Model
#'
#'
#' @param xy `data.frame` containing covariates and response.
#' @param x Covariate data only.
#' @param formula Regression formula.
#' @param comp_num Number of components.
#' @param max_out Maximum number of outliers.
#' @param mnames Model names for flexCWM::cwm.
#' @param seed Seed.
#' @param print_interval How frequently the iteration count is printed.
#' @param alpha Factor to control relative importance of covariate and response.
#'
#' @return List of
#' * distrib_diffs
#' * outlier_bool
#' * outlier_num
#' * outlier_rank
#' * gmm_labels
#' @export
#'
#' @examples
#' ombc_p1 <- simulate_ombc(
#' n = c(1000, 1000),
#' mu = list(c(-1), c(+1)),
#' sigma = list(as.matrix(0.2), as.matrix(0.2)),
#' beta = list(c(1, 0), c(1, 3)),
#' error_sd = c(1, 1),
#' outlier_num = c(25, 25),
#' outlier_type = "x_and_y",
#' seed = 123,
#' crit_val = 0.9999,
#' range_multipliers = c(1.5, 2)
#' )
#'
#' ombc_p1_lcwm <- ombc_lcwm(
#' xy = ombc_p1[, -3],
#' x = ombc_p1$X1,
#' formulaY = Y ~ X1,
#' comp_num = 2,
#' max_out = 200,
#' mnames = "V",
#' seed = 123,
#' print_interval = 100,
#' alpha = 0.5
#' )
#'
#' plot(
#' ombc_p1[, c("X1", "Y")], pch = ombc_p1$G + 1, col = ombc_p1_lcwm$gmm_labels
#' )
ombc_lcwm <- function(
  xy,
  x,
  y_formula,
  comp_num,
  max_out,
  mnames = "VVV",
  seed = 123,
  print_interval = Inf,
  alpha = 0.5
) {
  xy0 <- xy
  x <- as.matrix(x)
  x0 <- x

  var_num <- ncol(x)

  z0 <- NULL
  cwm_init <- "kmeans"

  distrib_diffs <- c()
  outlier_rank <- rep(0, nrow(x))
  for (i in seq_len(max_out + 1)) {
    if (i %% print_interval == 0) cat("i = ", i, "\n")

    set.seed(seed)
    invisible(utils::capture.output(lcwm <- flexCWM::cwm(
      formulaY = y_formula,
      familyY = stats::gaussian(link = "identity"),
      data = xy,
      Xnorm = x,
      modelXnorm = mnames,
      k = comp_num,
      initialization = cwm_init,
      start.z = z0,
      seed = seed
    )))

    prop <- colSums(lcwm$models[[1]]$posterior)

    if (any(prop < (var_num + 1))) {
      warning(paste0(
        "One of the components became too small after removing ",
        i - 1, " outliers."
      ))
      break()
    }

    mod_list <- lapply(lcwm$models[[1]]$GLModel, function(x) x$model)
    names(mod_list) <- paste0("comp.", seq_along(mod_list))
    y_sigma <- vapply(lcwm$models[[1]]$GLModel, function(x) x$sigma, double(1L))

    out <- distrib_diff_lcwm(
      x,
      lcwm$models[[1]]$posterior,
      prop,
      lcwm$models[[1]]$concomitant$normal.mu,
      lcwm$models[[1]]$concomitant$normal.Sigma,
      mod_list,
      y_sigma,
      alpha
    )

    distrib_diffs[i] <- out$distrib_diff

    outlier_rank[!outlier_rank][out$choice_id] <- i
    x <- x[-out$choice_id, , drop = FALSE]
    xy <- xy[-out$choice_id, , drop = FALSE]
    z0 <- lcwm$models[[1]]$posterior[-out$choice_id, , drop = FALSE]
    cwm_init <- "manual"
  }

  outlier_num <- which.min(distrib_diffs) - 1

  outlier_bool <- outlier_rank <= outlier_num & outlier_rank != 0

  set.seed(seed)
  lcwm <- flexCWM::cwm(
    formulaY = y_formula,
    familyY = stats::gaussian(link = "identity"),
    data = xy0[!outlier_bool, ],
    Xnorm = x0[!outlier_bool, ],
    modelXnorm = mnames,
    k = comp_num,
    initialization = "kmeans",
    seed = seed
  )

  gmm_labels <- rep(1, nrow(x0))
  gmm_labels[!outlier_bool] <- 1 + lcwm$models[[1]]$cluster

  return(list(
    distrib_diffs = distrib_diffs,
    outlier_bool = outlier_bool,
    outlier_num = outlier_num,
    outlier_rank = outlier_rank,
    gmm_labels = gmm_labels
  ))
}
