#' Modified, minimal version of the `cwm` function from the flexCWM package
#'
#' @description
#' Maximum likelihood fitting of the cluster-weighted model by the EM algorithm.
#'
#' @param formula_y
#'  an optional object of class "\code{\link[stats:formula]{formula}}" (or
#'  one that can be coerced to that class): a symbolic description of the
#'  model to be fitted.
#' @param data
#'  an optional \code{\link[base:data.frame]{data.frame}},
#'  \code{\link[base:list]{list}}, or
#'  \code{\link[base:environment]{environment}} with the variables needed to
#'  use \code{formula_y}.
#' @param x_norm
#'  an optional matrix containing variables to be used for marginalization
#'  having normal distributions.
#' @param model_x_norm
#'  an optional vector of character strings indicating the parsimonious
#'  models to be fitted for variables in \code{x_norm}. The default is
#'  \code{c("E", "V")} for a single continuous covariate, and
#'  \code{c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE",
#'  "EEV", "VVE", "VEV", "EVV", "VVV")} for multivariate continuous
#'  covariates (see \code{\link[mixture:gpcm]{mixture:gpcm} for details}).
#' @param k
#'  an optional vector containing the numbers of mixture components to be tried.
#'  Default value is \code{1:3}.
#' @param start_z
#'  matrix of soft or hard classification. Only models with the same number of
#'  mixture components as the columns of this matrix will be fit.
#' @param seed
#'  an optional scalar. It sets the seed for the random number generator, when
#'  random initializations are used; if \code{NULL}, current seed is not
#'  changed. Default value is \code{NULL}.
#' @param iter_max
#'  an optional scalar. It sets the maximum number of iterations in the
#'  EM-algorithm. Default value is 200.
#' @param threshold
#'  an optional scalar. It sets the threshold for the Aitken acceleration
#'  procedure. Default value is 1.0e-04.
#' @param eps
#'  an optional scalar. It sets the smallest value for eigenvalues of covariance
#'  matrices for \code{x_norm}. Default value is 1e-100.
#' @param pwarning When \code{TRUE}, warnings are printed.
#'
#' @details
#' When several models have been estimated, methods \code{summary} and
#' \code{print} consider the best model according to the information criterion
#' in \code{criterion}, among the estimated models having a number of components
#' among those in \code{k} an error distribution among those in \code{family_y}
#' and a parsimonious model among those in \code{model_x_norm}.
#'
#' @returns
#'  This function returns a class \code{cwm} object, which is a list of values
#'  related to the model selected. It contains:
#'  \item{call}{an object of class \code{call}.}
#'  \item{formula_y}{
#'    an object of class \code{\link[stats:formula]{formula}} containing a
#'    symbolic description of the model fitted.}
#'  \item{data}{
#'    a \code{\link[base:data.frame]{data.frame}} with the variables needed to
#'    use \code{formula_y}.}
#'  \item{concomitant}{
#'    a list containing \code{x_norm}
#'  }
#'  \item{models}{
#'    a list; each element is related to one of the models fitted. Each element
#'    is a list and contains:
#'  }
#'    \itemize{
#'      \item{\code{posterior} posterior probabilities}
#'      \item{\code{iter} number of iterations performed in EM algorithm}
#'      \item{\code{k} number of (fitted) mixture components.}
#'      \item{\code{size} estimated size of the groups.}
#'      \item{\code{cluster} classification vector}
#'      \item{\code{loglik} final log-likelihood value}
#'      \item{\code{df} overall number of estimated parameters}
#'      \item{\code{prior} weights for the mixture components}
#'      \item{\code{ic} list containing values of the information criteria }
#'      \item{\code{converged} logical; \code{TRUE} if EM algorithm converged}
#'      \item{\code{gl_model}
#'        a list;  each element is related to a mixture component and contains:
#'      }
#'      \itemize{
#'        \item{\code{model} a "\code{\link[stats:glm]{glm}}" class object.}
#'        \item{
#'          \code{concomitant} a list with estimated concomitant variables
#'          parameters for each mixture component
#'        }
#'        \item{\code{normal.d} marginal distribution of concomitant variables}
#'        \item{\code{normal.mu} mixture component means for \code{x_norm}}
#'        \item{
#'          \code{normal.sigma} mixture component covariance matrices for
#'          \code{x_norm}
#'        }
#'        \item{\code{normal.model} models fitted for \code{x_norm}}
#'      }
#'    }
#'
#' @references {
#' Mazza, A., Ingrassia, S., and Punzo, A. (2018).
#' {flexCWM}: A Flexible Framework for Cluster-Weighted Models.
#' \emph{Journal of Statistical Software}, \bold{86}(2), 1-30.
#'
#' Ingrassia, S., Minotti, S. C., and Vittadini, G. (2012).
#' Local Statistical Modeling via the Cluster-Weighted Approach with Elliptical
#' Distributions.
#' \emph{Journal of Classification}, \bold{29}(3), 363-401.
#'
#' Ingrassia, S., Minotti, S. C., and Punzo, A. (2014).
#' Model-based clustering via linear cluster-weighted models.
#' \emph{Computational Statistics and Data Analysis}, \bold{71}, 159-182.
#'
#' Ingrassia, S., Punzo, A., and Vittadini, G. (2015).
#' The Generalized Linear Mixed Cluster-Weighted Model.
#' \emph{Journal of Classification}, \bold{32}(forthcoming)
#'
#' McCullagh, P. and Nelder, J. (1989).
#' Generalized Linear Models.
#' Chapman & Hall, Boca Raton, 2nd edition
#'
#' Punzo, A. (2014).
#' Flexible Mixture Modeling with the Polynomial Gaussian Cluster-Weighted
#' Model.
#' \emph{Statistical Modelling}, \bold{14}(3), 257-291.
#' }
#'
#' @author Mazza A., Punzo A., Ingrassia S.
fit_cwm <- function(
  formula_y = NULL,
  data,
  x_norm = NULL,
  model_x_norm = NULL,
  k,
  start_z = NULL,
  seed = NULL,
  iter_max = 1000,
  threshold = 1.0e-04,
  eps = 1.0e-100,
  pwarning = FALSE
) {
  if (is.null(formula_y) && is.null(x_norm)) {
    stop("Either formula_y or x_norm must be provided.")
  }

  if (!is.null(formula_y)) {
    if (is.null(data)) {
      data <- stats::get_all_vars(stats::as.formula(formula_y))
    } else {
      data <- stats::get_all_vars(stats::as.formula(formula_y), data)
    }
    y <- as.matrix(stats::model.response(stats::model.frame(formula_y, data)))
  } else {
    y <- data <- NULL
  }

  col_x_n <- 0
  if (!is.null(x_norm)) {
    x_norm <- as.matrix(x_norm)
    stopifnot("Data length mismatch" = nrow(x_norm) == nrow(y))
    col_x_n <- ncol(x_norm)
    if (is.null(colnames(x_norm))) {
      colnames(x_norm) <- paste0("x_norm", seq_len(col_x_n))
    }
  }

  if (col_x_n > 0) {
    if (col_x_n == 1) {
      model_x_names <- c("E", "V")
      eqmod <- data.frame(name = model_x_names, number = c(1, 1))
    } else {
      model_x_names <- c(
        "EII", "VII", "EEI", "VEI", "EVI", "VVI",
        "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV"
      )
      eqmod <- data.frame(name = model_x_names, number = rep(1:3, c(2, 4, 8)))
    }
    if (is.null(model_x_norm)) {
      model_x_norm <- model_x_names
    } else {
      model_x_norm <- match.arg(model_x_norm, model_x_names, several.ok = TRUE)
    }
  } else {
    model_x_norm <- eqmod <- NULL
  }

  mm <- expand.grid(list("k" = k, "model_x_norm" = model_x_norm))
  if (1 %in% as.integer(ceiling(k))) {
    mm <- merge(mm, eqmod, by.x = "model_x_norm", by.y = "name")
    mm1 <- mm[mm$k == 1, ]
    mm2 <- mm1[!duplicated(subset(mm1, select = -model_x_norm, drop = FALSE)), ]

    if (nrow(mm1) > nrow(mm2)) {
      cat(paste0(
        "\nWith k = 1, some models in model_x_norm are equivalent, ",
        "so only the first is estimated.\n"
      ))
    }
    mm <- rbind(mm2, mm[mm$k != 1, ])
  }
  mm <- mm[order(mm$k), , drop = FALSE]
  job <- function(i) {
    cat(paste0(
      "\nEstimating model with k = ", mm$k[i],
      ", x_norm = ", mm$model_x_norm[i], " "
    ))
    cwm2(
      formula_y = formula_y, data = data, y = y,
      x_norm = x_norm, n = nrow(y), col_x_n = col_x_n,
      k = mm$k[i], model_x_norm = mm$model_x_norm[i],
      start_z = start_z, iter_max = iter_max, threshold = threshold,
      eps = eps, pwarning = pwarning
    )
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  models <- lapply(seq_len(nrow(mm)), \(i) job(i))

  concomitant <- data.frame(x_norm)
  for (i in seq_len(ncol(x_norm))) {
    attr(concomitant[[i]], "type") <- "x_norm"
  }

  cat("\n")
  res <- structure(
    list(
      "call" = match.call(),
      "formula_y" = formula_y,
      "data" = data,
      "concomitant" = concomitant,
      "models" = models
    ),
    class = "cwm"
  )
  print(res)
  invisible(res)
}

cwm2 <- function(
  formula_y, data, y, x_norm, col_x_n, n, k, model_x_norm,
  start_z, iter_max, threshold, eps, pwarning
) {
  prior <- numeric(length = k)

  stats_gaussian <- stats::gaussian()

  if (!is.null(data)) {
    m_y <- rep(1, nrow(y))
  }
  a <- try(
    {
      z <- (start_z + 1e-4) / rowSums(start_z + 1e-4)
      exit_while <- FALSE
      iter <- 1
      loglik <- NULL
      aloglik <- c(0, 0)
      a <- c(0, 0)

      while (!exit_while) {
        prior <- colMeans(z)
        norm <- .p_x_norm(
          col_x_n = col_x_n, x_norm = x_norm, model_x_norm = model_x_norm,
          z = z, k = k, n = n, eps = eps
        )

        if (!is.null(data)) {
          l <- do.call(
            ".family_y_gaussian",
            list(
              "stats_gaussian" = stats_gaussian,
              "k" = k, "y" = y, "m_y" = m_y, "n" = n,
              "data" = data, "z" = z, "formula_y" = formula_y
            )
          )
        } else {
          l <- list(p_y = 1)
        }

        dens <- sweep(norm$p_x * l$p_y, MARGIN = 2, STATS = prior, FUN = "*")
        loglik[iter] <- sum(log(rowSums(dens)))

        if (iter > 2 & k > 1) {
          loglik_diff12 <- loglik[iter - 1] - loglik[iter - 2]
          if (!is.na(loglik_diff12) && !is.na(loglik[iter])) {
            if (abs(loglik_diff12) > 0) {
              loglik_diff01 <- loglik[iter] - loglik[iter - 1]
              a[iter - 1] <- loglik_diff01 / loglik_diff12
              a_ratio <- 1 / (1 - a[iter - 1])
              aloglik[iter] <- loglik[iter - 1] + a_ratio * loglik_diff01

              exit_while <- abs(aloglik[iter] - loglik[iter]) < threshold
            }
          } else {
            exit_while <- TRUE
          }
        }

        exit_while <- exit_while || iter == iter_max || k == 1

        cat("*")
        iter <- iter + 1

        dens[dens == 0] <- 1e-320
        z <- dens / rowSums(dens)

        if (min(colSums(z)) <= col_x_n) {
          z <- (z + 1e-7) / rowSums(z + 1e-7)
        }
      }

      finalloglik <- loglik[iter - 1]

      df_model <- length(stats::coef(l$lmodel_y[[1]])) * k
      if (!is.null(model_x_norm)) {
        n_cov_par <- ncovpar(
          modelname = strtrim(paste0(model_x_norm, "II"), 3), p = col_x_n, g = k
        )
        df_model <- df_model + n_cov_par
      }
      df <- df_model + (k - 1) + col_x_n * k
      if (df >= n) {
        warning(paste0(
          df, "parameters with only ", n,
          " observations; AICc and AICu are set to -Inf."
        ), call. = FALSE)
      }

      gl_model <- NULL
      if (!is.null(data)) {
        gl_model <- list()
        for (i in seq_len(k)) {
          gl_model[[i]] <- c(
            list("model" = l$lmodel_y[[i]]),
            list("sigma" = sqrt(l$disp_y[i]))
          )
        }
        names(gl_model) <- paste0("comp.", 1:k)
      }

      if (is.null(x_norm)) {
        concomitant <- NULL
      } else {
        concomitant <- c("normal" = list(
          "d" = norm$p_x,
          "mu" = norm$mu,
          "sigma" = norm$sigma,
          "model" = model_x_norm
        ))
      }

      cluster <- apply(z, 1, base::which.max)

      result <- list(
        "posterior" = z,
        "iter" = iter,
        "k" = k,
        "size" = table(cluster),
        "cluster" = cluster,
        "logLik" = finalloglik,
        "df" = df,
        "prior" = prior,
        "ic" = compute_ic(finalloglik, n, df, z),
        "converged" = TRUE,
        "gl_model" = gl_model,
        "concomitant" = concomitant
      )
    },
    silent = !pwarning
  )
  if (!length(a) > 1) {
    result <- list(iter = iter, k = k, ic = list(BIC = NA))
    cat(" Not estimated")
  }
  class(result) <- "cwm"

  result
}

ncovpar <- function(modelname = NULL, p = NULL, g = NULL) {
  if (is.null(p)) {
    stop("p is null")
  } else if (is.null(g)) {
    stop("g is null")
  } else if (is.null(modelname)) {
    stop("modelname is null")
  } else if (modelname == "EII") {
    npar <- 1
  } else if (modelname == "VII") {
    npar <- g
  } else if (modelname == "EEI") {
    npar <- p
  } else if (modelname == "VEI") {
    npar <- p + g - 1
  } else if (modelname == "EVI") {
    npar <- p * g - g + 1
  } else if (modelname == "VVI") {
    npar <- p * g
  } else if (modelname == "EEE") {
    npar <- p * (p + 1) / 2
  } else if (modelname == "EEV") {
    npar <- g * p * (p + 1) / 2 - (g - 1) * p
  } else if (modelname == "VEV") {
    npar <- g * p * (p + 1) / 2 - (g - 1) * (p - 1)
  } else if (modelname == "VVV") {
    npar <- g * p * (p + 1) / 2
  } else if (modelname == "EVE") {
    npar <- p * (p + 1) / 2 + (g - 1) * (p - 1)
  } else if (modelname == "VVE") {
    npar <- p * (p + 1) / 2 + (g - 1) * p
  } else if (modelname == "VEE") {
    npar <- p * (p + 1) / 2 + (g - 1)
  } else if (modelname == "EVV") {
    npar <- g * p * (p + 1) / 2 - (g - 1)
  } else {
    stop("modelname is not correctly defined")
  }

  npar
}

#' @importFrom mclust mstepE mstepV mstepEEE mstepVVV
#' @importFrom mclust mstepEII mstepVII mstepEEI mstepVEI mstepEVI mstepVVI
#' @importFrom mclust mstepEEV mstepEVE mstepVEE mstepVEV mstepEVV mstepVVE
.p_x_norm <- function(col_x_n, x_norm, model_x_norm, z, k, n, eps) {
  if (col_x_n > 0) {
    comp_names <- paste0("comp.", seq_len(k))

    p_x_norm <- array(1, c(n, k), dimnames = list(seq_len(n), comp_names))
    mu_x <- array(
      0, c(col_x_n, k),
      dimnames = list(dimnames(x_norm)[[2]], comp_names)
    )
    sigma <- array(
      0, c(col_x_n, col_x_n, k),
      dimnames = list(dimnames(x_norm)[[2]], dimnames(x_norm)[[2]], comp_names)
    )
    if (col_x_n > 1) {
      fit_m <- ContaminatedMixt::m.step(
        x = x_norm, modelname = model_x_norm, z = z, mtol = 1e-25, mmax = 10
      )
      mu_x <- fit_m$mu
      sigma <- fit_m$sigma
      p_x_norm <- fit_m$px
    } else {
      fit_m <- mclust::mstep(
        modelName = substr(model_x_norm, 1, 1), data = x_norm, z = z
      )
      mu_x[1, ] <- fit_m$parameters$mean
      if (substr(model_x_norm, 1, 1) == "E") {
        sigma[1, 1, ] <- rep(fit_m$parameters$variance$sigmasq, k)
      }
      if (substr(model_x_norm, 1, 1) == "V") {
        sigma[1, 1, ] <- fit_m$parameters$variance$sigmasq
      }
      for (j in 1:k) {
        p_x_norm[, j] <- stats::dnorm(x_norm, mu_x[1, j], sqrt(sigma[1, 1, j]))
      }
    }
  } else {
    p_x_norm <- 1
    mu_x <- sigma <- NULL
  }
  list("p_x" = p_x_norm, "mu" = mu_x, "sigma" = sigma)
}

.family_y_gaussian <- function(
  stats_gaussian, k, data, y, n, z, formula_y, ...
) {
  lmodel_y <- list()
  mu_y <- p_y <- matrix(NA, nrow = n, ncol = k)
  disp_y <- c()
  for (h in seq_len(k)) {
    lmodel_y[[h]] <- do.call(
      stats::glm,
      list(
        "formula" = formula_y,
        "family" = stats_gaussian,
        "data" = data,
        "weights" = z[, h],
        "model" = FALSE,
        "x" = FALSE,
        "y" = FALSE
      )
    )
    mu_y[, h] <- stats::fitted(lmodel_y[[h]])
    disp_y[h] <- as.numeric(
      crossprod(sqrt(z[, h] / sum(z[, h])) * (y - mu_y[, h]))
    )
    p_y[, h] <- stats::dnorm(y, mean = mu_y[, h], sd = sqrt(disp_y[h]))
  }
  colnames(mu_y) <- colnames(p_y) <- paste("comp.", 1:k, sep = "")

  list("lmodel_y" = lmodel_y, "mu_y" = mu_y, "disp_y" = disp_y, "p_y" = p_y)
}

compute_ic <- function(ll, n, df, z) {
  ic <- list()
  ic$AIC <- 2 * ll - df * 2
  ic$BIC <- 2 * ll - df * log(n)
  ic$AIC3 <- 2 * ll - df * 3
  ic$CAIC <- 2 * ll - df * (1 + log(n))
  ic$AWE <- 2 * ll - 2 * df * (3 / 2 + log(n))
  if (n - df - 1 > 0) {
    ic$AICc <- ic$AIC - (2 * df * (df + 1)) / (n - df - 1)
    ic$AICu <- ic$AICc - n * log(n / (n - df - 1))
  } else {
    ic$AICc <- ic$AICu <- -Inf
  }

  z_const <- ifelse(z < 1e-322, 1e-322, z)
  hard_z <- t(1 * apply(z, 1, \(u) u == max(u)))
  ecm <- sum(hard_z * log(z_const))
  ic$ICL <- ic$BIC + ecm

  ic
}
