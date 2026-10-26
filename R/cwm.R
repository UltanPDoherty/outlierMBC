#' `cwm` function from Angelo Mazza's flexCWM package
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
#' @param initialization
#'  an optional character string. It sets the initialization strategy for the
#'  EM-algorithm. It can be:
#'  \itemize{
#'    \item \code{"random.soft"}
#'    \item \code{"random.hard"}
#'    \item \code{"kmeans"}
#'    \item \code{"mclust"}
#'    \item \code{"manual"}
#'  }
#'  Default value is \code{"random.soft"}.
#' @param start_z
#'  matrix of soft or hard classification: it is used only if
#'  \code{initialization = "manual"}. Only models with the same number of
#'  mixture components as the columns of this matrix will be fit.
#' @param seed
#'  an optional scalar. It sets the seed for the random number generator, when
#'  random initializations are used; if \code{NULL}, current seed is not
#'  changed. Default value is \code{NULL}.
#' @param max_init number of initializations to be tried. Default value is 1.
#' @param iter_max
#'  an optional scalar. It sets the maximum number of iterations in the
#'  EM-algorithm. Default value is 200.
#' @param threshold
#'  an optional scalar. It sets the threshold for the Aitken acceleration
#'  procedure. Default value is 1.0e-04.
#' @param eps
#'  an optional scalar. It sets the smallest value for eigenvalues of covariance
#'  matrices for \code{x_norm}. Default value is 1e-100.
#' @param parallel
#'  When \code{TRUE}, the package
#'  \code{\link[parallel:parallel-package]{parallel}} is used for parallel
#'  computation. When several models are estimated, computational time is
#'  reduced. The number of cores to use may be set with the global option
#'  \code{cl.cores}; default value is detected using
#'  \code{\link[parallel:detectCores]{detectCores()}}.
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
cwm <- function(
  formula_y = NULL,
  data = NULL,
  x_norm = NULL,
  model_x_norm = NULL,
  k = 1:3,
  initialization = c(
    "random.soft", "random.hard", "kmeans", "mclust", "manual"
  ),
  start_z = NULL,
  seed = NULL,
  max_init = 1,
  iter_max = 1000,
  threshold = 1.0e-04,
  eps = 1.0e-100,
  parallel = FALSE,
  pwarning = FALSE
) {
  # Preliminary checks and init ----------------------------------------------
  if (is.null(formula_y) && is.null(x_norm)) {
    stop("No data were entered.")
  }

  initialization <- match.arg(initialization)
  lm <- list(k = k)
  # Get regression terms ----------------------------------------
  if (!is.null(formula_y)) {
    family_y <- list(family_y = stats::gaussian())
    if (is.null(data)) {
      data <- stats::get_all_vars(stats::as.formula(formula_y))
    } else {
      data <- stats::get_all_vars(stats::as.formula(formula_y), data)
    }
    y <- as.matrix(stats::model.response(stats::model.frame(formula_y, data)))
    lm$family_y <- seq_along(family_y)
  } else {
    y <- data <- family_y <- NULL
  }

  n <- unlist(sapply(list(y, x_norm), function(x) nrow(cbind(x))))
  if (max(n) != min(n)) stop("Data length mismatch")
  n <- n[1]

  # Fix variables to use for marginal --------------------------
  col_x_n <- 0
  x_mod <- list()
  if (!is.null(x_norm)) {
    x_norm <- as.matrix(x_norm)
    col_x_n <- ncol(x_norm)
    if (is.null(colnames(x_norm))) {
      colnames(x_norm) <- paste0("x_norm", 1:col_x_n)
    }
  }

  # More checks ------------------------------------------
  if (
    (initialization == "mclust" || initialization == "kmeans") &&
      is.null(x_norm) && is.null(data)
  ) {
    stop(paste0(
      "when initialization is '", initialization, "', numeric variables are
      needed."
    ))
  }

  if (col_x_n > 0) {
    if (col_x_n == 1) {
      model_x_norm_names <- c("E", "V")
      eqmod <- data.frame(name = model_x_norm_names, number = c(1, 1))
    } else {
      model_x_norm_names <- c(
        "EII", "VII", "EEI", "VEI", "EVI", "VVI",
        "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV"
      )
      eqmod <- data.frame(
        name = model_x_norm_names,
        number = c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3)
      )
    }
    if (is.null(model_x_norm)) {
      model_x_norm <- model_x_norm_names
    } else {
      model_x_norm <- match.arg(
        model_x_norm, model_x_norm_names,
        several.ok = TRUE
      )
    }
    lm$model_x_norm <- model_x_norm
  } else {
    model_x_norm <- eqmod <- NULL
  }

  k <- as.integer(ceiling(k))
  par <- list()
  # Compute ----------------------------------------------
  mm <- expand.grid(lm)
  if (!is.null(model_x_norm) && 1 %in% k) {
    mm <- merge(mm, eqmod, by.x = "model_x_norm", by.y = "name")
    mm1 <- mm[mm$k == 1, ]
    mm2 <- mm1[!duplicated(subset(mm1, select = -model_x_norm, drop = FALSE)), ]

    if (nrow(mm1) > nrow(mm2)) {
      cat(
        "\nWith k=1, some models in model_x_norm are equivalent, so only the
        first is estimated.\n"
      )
    }
    mm <- rbind(mm2, mm[mm$k != 1, ])
  }
  mm <- mm[order(mm$k), , drop = FALSE]
  job <- function(i) {
    cat("\nEstimating model")
    cat(paste0(" with k=", mm$k[i]))
    if (!is.null(mm$model_x_norm[i])) {
      cat(paste0(", x_norm=", mm$model_x_norm[i]))
    }
    if (!is.null(mm$family_y[i])) {
      cat(paste0(", family_y=", family_y[[mm$family_y[i]]][1]))
    }
    cat(" ")
    cwm2(
      formula_y = formula_y, data = data, y = y,
      x_norm = x_norm, n = n, col_x_n = col_x_n, x_mod = x_mod,
      k = mm$k[i], model_x_norm = mm$model_x_norm[i],
      family_y = family_y[[mm$family_y[i]]], method = "Nelder-Mead",
      initialization = initialization,
      start_z = start_z, iter_max = iter_max, threshold = threshold,
      max_init = max_init, eps = eps, pwarning = pwarning
    )
  }
  start_z <- start_z
  if (parallel) {
    cores <- getOption("cl.cores", parallel::detectCores())
    cat(paste("Using", cores, "cores\n"))
    cl <- parallel::makeCluster(
      cores,
      outfile = file.path(tempdir(), "temp.out")
    )
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl = cl, iseed = seed)
    par <- parallel::parLapply(cl = cl, seq_len(nrow(mm)), function(i) {
      foo <- job(i)
      gc()
      foo
    })
    parallel::stopCluster(cl)
  } else {
    if (!is.null(seed)) set.seed(seed)
    par <- lapply(seq_len(nrow(mm)), function(i) job(i))
  }

  cat("\n")
  res <- structure(list(
    call = match.call(),
    formula_y = formula_y,
    data = data,
    concomitant = list(x_norm = as_df(x_norm, "x_norm")),
    models = par
  ), class = "cwm")
  print(res)
  invisible(res)
}

as_df <- function(ma, type) {
  if (!is.null(ma)) {
    df <- data.frame(ma)
    for (i in seq_len(ncol(ma))) attr(df[[i]], "type") <- type
    df
  } else {
    NULL
  }
}

cwm2 <- function(
  formula_y, data, y, x_norm, col_x_n,
  n, x_mod, k, model_x_norm, family_y,
  method, initialization, start_z, iter_max, threshold, loglikplot,
  max_init, eps, pwarning
) {
  if (initialization == "manual" || initialization == "mclust") {
    max_init <- 1
  }

  # Parameters definition ------------------------------------------------------
  prior <- numeric(k) # weights

  if (!is.null(data)) {
    m_y <- rep(1, nrow(y))
  }
  result <- list()
  # Do Estimation --------------------------------------------------------------
  for (r in 1:max_init) {
    a <- try(
      {
        if (max_init > 1) cat(paste("\nrandom init.", r))
        z <- .post_init(
          initialization = initialization, y = y,
          x_norm = x_norm, n = n, k = k,
          start_z = start_z
        )

        # EM algorithm --------------------------------------------------
        # Preliminary definition of convergence criterions
        check <- 0
        iteration <- 1
        loglik <- NULL
        aloglik <- c(0, 0)
        a <- c(0, 0)

        while (check < 1) {
          # M-step --------------------------------------------------
          # Weights #
          prior <- colMeans(z)
          # Posterior probabilities
          norm <- .p_x_norm(
            col_x_n = col_x_n, x_norm = x_norm, model_x_norm = model_x_norm,
            z = z, k = k, n = n, eps = eps
          )

          # Y|x pag 45 of McCullag & Nelder (1989) -----------------------------
          if (!is.null(data)) {
            l <- do.call(
              ".family_y_gaussian",
              list(
                family_y = family_y, k = k, y = y, m_y = m_y, n = n,
                data = data, z = z, method = method, formula_y = formula_y
              )
            )
          } else {
            l <- list(p_y = 1)
          }
          # Global - Observed-data log-likelihood ##
          dens <- matrix(rep(prior, n), n, k, byrow = TRUE) * norm$p_x * l$p_y
          loglik[iteration] <- sum(log(rowSums(dens)))

          # Aitkane's Acceleration-Based Stopping Criterion # #
          if (iteration > 2 & k > 1) {
            if (
              !is.na(loglik[iteration - 1] - loglik[iteration - 2]) &&
                !is.na(loglik[iteration])
            ) {
              if (abs(loglik[iteration - 1] - loglik[iteration - 2]) > 0) {
                loglik_diff01 <- loglik[iteration] - loglik[iteration - 1]
                loglik_diff12 <- loglik[iteration - 1] - loglik[iteration - 2]
                a[iteration - 1] <- loglik_diff01 / loglik_diff12
                aloglik[iteration] <- loglik[iteration - 1] +
                  (1 / (1 - a[iteration - 1]) * loglik_diff01)
                if (abs(aloglik[iteration] - loglik[iteration]) < threshold) {
                  check <- 1
                }
              } else {
                check <- 1
              }
            }
          }
          if (iteration == iter_max | k == 1) check <- 1

          cat("*")
          iteration <- iteration + 1

          # E-Step -------------------------------------------------------------
          dens[dens == 0] <- 1E-320
          z <- dens / matrix(rep(rowSums(dens), k), ncol = k) # (n x k)

          # z adjustment
          if (min(colSums(z)) <= col_x_n) {
            z <- (z + 0.0000001) / rowSums(z + 0.0000001)
          }
        }

        finalloglik <- loglik[iteration - 1]

        # EM-algorithm is finished ---------------------------------------------

        # Check the number of parameters -------------------------------
        df_model <- length(stats::coef(l$lmodel_y[[1]])) * k
        if (!is.null(model_x_norm)) {
          df_model <-
            df_model + ncovpar(
              modelname = strtrim(paste0(model_x_norm, "II"), 3),
              p = col_x_n, g = k
            )
        }
        df <- df_model + (k - 1) + col_x_n * k
        if (df >= n) {
          warning(paste0(
            df, "parameters with only ", n,
            " observations; AICc and AICu are set to -Inf."
          ), call. = FALSE)
        }

        # Classification Matrix ------------------------------------------------
        cluster <- apply(z, 1, base::which.max)

        # Information criteria -------------------------------------------------
        ic <- list()
        ic$AIC <- 2 * finalloglik - df * 2
        ic$BIC <- 2 * finalloglik - df * log(n)
        ic$AIC3 <- 2 * finalloglik - df * 3
        ic$CAIC <- 2 * finalloglik - df * (1 + log(n))
        ic$AWE <- 2 * finalloglik - 2 * df * (3 / 2 + log(n))
        if (n - df - 1 > 0) {
          ic$AICc <- ic$AIC - (2 * df * (df + 1)) / (n - df - 1)
          ic$AICu <- ic$AICc - n * log(n / (n - df - 1))
        } else {
          ic$AICc <- ic$AICu <- -Inf
        }

        z_const <- (z < 10^(-322)) * 10^(-322) + (z > 10^(-322)) * z
        # vincolo per evitare i NaN nel calcolo di tau*log(tau)
        hard_z <-
          (matrix(rep(apply(z, 1, max), k), n, k, byrow = FALSE) == z) * 1
        ecm <- sum(hard_z * log(z_const))
        ic$ICL <- ic$BIC + ecm
        ## GLM parameters
        gl_model <- NULL
        if (!is.null(data)) {
          gl_model <- list()
          for (i in seq_len(k)) {
            parameters <- list(model = l$lmodel_y[[i]])
            parameters <- c(parameters, list(sigma = sqrt(l$var_y[1, i])))
            gl_model[[i]] <- parameters
          }
          names(gl_model) <- paste0("comp.", 1:k)
        }
        ## Concomitant
        concomitant <- NULL
        if (!is.null(x_norm)) {
          concomitant <- list()
          normal <- list(
            d = norm$p_x,
            mu = norm$mu,
            sigma = norm$sigma,
            model = model_x_norm
          )
          concomitant <- c(concomitant, normal = normal)
        }
        result[[r]] <- list(
          posterior = z,
          iter = iteration,
          k = k,
          size = table(cluster),
          cluster = cluster,
          logLik = finalloglik,
          df = df,
          prior = prior,
          ic = ic,
          converged = TRUE, # Logical, TRUE if EM algorithm converged.
          gl_model = gl_model,
          concomitant = concomitant
        )
      },
      silent = !pwarning
    )
    if (!length(a) > 1) {
      result[[r]] <- list(
        iter = iteration,
        k = k,
        ic = list(BIC = NA)
      )
      cat(" Not estimated")
    }
  } # max_init
  m <- which.max(sapply(1:max_init, function(r) result[[r]]$ic$BIC))
  if (length(m) == 0) m <- 1
  result <- result[[m]]
  class(result) <- "cwm"

  result
}

ncovpar <- function(modelname = NULL, p = NULL, g = NULL) {
  if (is.null(p)) {
    stop("p is null")
  }
  if (is.null(g)) {
    stop("g is null")
  }
  if (is.null(modelname)) {
    stop("modelname is null")
  }
  if (modelname == "EII") {
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

.post_init <- function(
  initialization, y, x_norm, n, k, start_z
) {
  if (k == 1) {
    z <- matrix(rep(1, n), ncol = 1)
    return(z)
  }
  if (initialization == "mclust") {
    data <- cbind(y[, 1], x_norm)
    z <- mclust::Mclust(
      data = data, modelNames = ifelse(ncol(data) == 1, "V", "VVV"), G = k
    )$z
  }
  if (initialization == "kmeans") {
    clusters <- stats::kmeans(
      x = cbind(y, x_norm), centers = k
    )
    # clusters on D
    z <- mclust::unmap(clusters$cluster, 1:k)
  }
  if (initialization == "random.soft") {
    z <- array(stats::runif(n * k), c(n, k))
    # soft posterior probabilities (no-normalized) (n x k)
    z <- z / rowSums(z)
    # soft posterior probabilities (n x k)
  }
  if (initialization == "random.hard") {
    z <- t(stats::rmultinom(n, size = 1, prob = rep(1 / k, k)))
    # hard posterior probabilities (n x k)
  }
  if (initialization == "manual") {
    # z.start can be both soft and hard initialization
    z <- start_z # posterior probabilities (n x k) no-normalized
  }
  z <- (z + 0.0001) / rowSums(z + 0.0001)
  z
}

#' @importFrom mclust mstepE
#' @importFrom mclust mstepV
#' @importFrom mclust mstepEII
#' @importFrom mclust mstepVII
#' @importFrom mclust mstepEEI
#' @importFrom mclust mstepVEI
#' @importFrom mclust mstepEVI
#' @importFrom mclust mstepVVI
#' @importFrom mclust mstepEEE
#' @importFrom mclust mstepEVE
#' @importFrom mclust mstepVEV
#' @importFrom mclust mstepEVV
#' @importFrom mclust mstepVVV
.p_x_norm <- function(col_x_n, x_norm, model_x_norm, z, k, n, eps) {
  if (col_x_n > 0) {
    p_x_norm <- array(
      1, c(n, k),
      dimnames = list(1:n, paste("comp.", 1:k, sep = ""))
    )
    mu_x <- array(
      0, c(col_x_n, k),
      dimnames = list(dimnames(x_norm)[[2]], paste("comp.", 1:k, sep = ""))
    )
    sigma <- array(
      0, c(col_x_n, col_x_n, k),
      dimnames = list(
        dimnames(x_norm)[[2]],
        dimnames(x_norm)[[2]],
        paste("comp.", 1:k, sep = "")
      )
    )
    if (col_x_n > 1) { #  More than one normal concomitant
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
  list(p_x = p_x_norm, mu = mu_x, sigma = sigma)
}

.family_y_gaussian <- function(family_y, k, data, y, n, z, formula_y, ...) {
  beta <- var_fun_y <- var_y <- mu_y <- p_y <- disp_y <- NULL
  lmodel_y <- list()
  for (h in 1:k) {
    model_y <- do.call(
      stats::glm,
      list(
        formula = formula_y, data = data, family = family_y, weights = z[, h],
        y = FALSE, model = FALSE, x = FALSE
      )
    )
    sigma_hat <- sqrt(
      crossprod(sqrt(z[, h] / sum(z[, h])) * (y - stats::fitted(model_y)))
    )
    beta <- rbind(beta, stats::coef(model_y))
    mu_y <- cbind(mu_y, stats::fitted(model_y))
    disp_y <- rbind(disp_y, sigma_hat^2)
    var_fun_y <- cbind(var_fun_y, rep(1, n))
    var_y <- cbind(var_y, as.vector(sigma_hat^2) * rep(1, n))
    p_y <- cbind(
      p_y, stats::dnorm(y, mean = stats::fitted(model_y), sd = sigma_hat)
    )
    lmodel_y[[h]] <- model_y
  }
  dimnames(var_fun_y) <- dimnames(var_y) <- dimnames(mu_y) <- dimnames(p_y) <-
    list(1:n, paste("comp.", 1:k, sep = ""))

  list(
    lmodel_y = lmodel_y, mu_y = mu_y, disp_y = as.vector(disp_y),
    var_fun_y = var_fun_y, var_y = var_y, p_y = p_y
  )
}
