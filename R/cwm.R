#' `cwm` function from Angelo Mazza's flexCWM package
#'
#' @description
#' Maximum likelihood fitting of the cluster-weighted model by the EM algorithm.
#'
#' @param formulaY
#'  an optional object of class "\code{\link[stats:formula]{formula}}" (or
#'  one that can be coerced to that class): a symbolic description of the
#'  model to be fitted.
#' @param familyY
#'  a description of the error distribution and link function to be used for the
#'  conditional distribution of \eqn{Y} in each mixture component. This can be a
#'  character string naming a \code{\link[stats:family]{family function}}, a
#'  family function or the result of a call to a family function.
#'  The following family functions are supported:
#'  \itemize{
#'    \item \code{binomial(link = "logit")}
#'    \item \code{gaussian(link = "identity")}
#'    \item \code{Gamma(link = "log")}
#'    \item \code{inverse.gaussian(link = "1/mu^2")}
#'    \item \code{poisson(link = "log")}
#'    \item \code{student.t(link = "identity")}
#'  }
#'  Default value is \code{gaussian(link = "identity")}.
#' @param data
#'  an optional \code{\link[base:data.frame]{data.frame}},
#'  \code{\link[base:list]{list}}, or
#'  \code{\link[base:environment]{environment}} with the variables needed to
#'  use \code{formulaY}.
#' @param Xnorm
#'  an optional matrix containing variables to be used for marginalization
#'  having normal, binomial, Poisson and multinomial distributions.
#' @param Xbin
#'  an optional matrix containing variables to be used for marginalization
#'  having normal, binomial, Poisson and multinomial distributions.
#' @param Xpois
#'  an optional matrix containing variables to be used for marginalization
#'  having normal, binomial, Poisson and multinomial distributions.
#' @param Xmult
#'  an optional matrix containing variables to be used for marginalization
#'  having normal, binomial, Poisson and multinomial distributions.
#' @param modelXnorm
#'  an optional vector of character strings indicating the parsimonious
#'  models to be fitted for variables in \code{Xnorm}. The default is
#'  \code{c("E", "V")} for a single continuous covariate, and
#'  \code{c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE",
#'  "EEV", "VVE", "VEV", "EVV", "VVV")} for multivariate continuous
#'  covariates (see \code{\link[mixture:gpcm]{mixture:gpcm} for details}).
#' @param Xbtrials
#'  an optional vector containing the number of trials for each column in
#'  \code{Xbin}. If omitted, the maximum of each column in \code{Xbin} is
#'  used.
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
#' @param start.z
#'  matrix of soft or hard classification: it is used only if
#'  \code{initialization = "manual"}. Only models with the same number of
#'  mixture components as the columns of this matrix will be fit.
#' @param seed
#'  an optional scalar. It sets the seed for the random number generator, when
#'  random initializations are used; if \code{NULL}, current seed is not
#'  changed. Default value is \code{NULL}.
#' @param maxR number of initializations to be tried. Default value is 1.
#' @param iter.max
#'  an optional scalar. It sets the maximum number of iterations in the
#'  EM-algorithm. Default value is 200.
#' @param threshold
#'  an optional scalar. It sets the threshold for the Aitken acceleration
#'  procedure. Default value is 1.0e-04.
#' @param eps
#'  an optional scalar. It sets the smallest value for eigenvalues of covariance
#'  matrices for \code{Xnorm}. Default value is 1e-100.
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
#' When \code{familyY = binomial}, the response variable must be a matrix with
#' two columns, where the first column is the number of "successes" and the
#' second column is the number of "failures". When several models have been
#' estimated, methods \code{summary} and \code{print} consider the best model
#' according to the information criterion in \code{criterion}, among the
#' estimated models having a number of components among those in \code{k} an
#' error distribution among those in \code{familyY} and a parsimonious model
#' among those in \code{modelXnorm}.
#'
#' @returns
#'  This function returns a class \code{cwm} object, which is a list of values
#'  related to the model selected. It contains:
#'  \item{call}{an object of class \code{call}.}
#'  \item{formulaY}{
#'    an object of class \code{\link[stats:formula]{formula}} containing a
#'    symbolic description of the model fitted.}
#'  \item{familyY}{
#'    the distribution used for the conditional distribution of \eqn{Y} in each
#'    mixture component.}
#'  \item{data}{
#'    a \code{\link[base:data.frame]{data.frame}} with the variables needed to
#'    use \code{formulaY}.}
#'  \item{concomitant}{
#'    a list containing \code{Xnorm}, \code{Xbin}, \code{Xpois}, \code{Xmult}.
#'  }
#'  \item{\code{Xbtrials}}{number of trials used for \code{Xbin}.}
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
#'      \item{\code{IC} list containing values of the information criteria }
#'      \item{\code{converged} logical; \code{TRUE} if EM algorithm converged}
#'    }
#'  \item{\code{GLModels}}{
#'    a list;  each element is related to a mixture component and contains:
#'  }
#'    \itemize{
#'      \item{\code{model} a "\code{\link[stats:glm]{glm}}" class object.}
#'      \item{
#'        \code{sigma} estimated local scale parameters of the conditional
#'        distribution of \eqn{Y}, when \code{familyY} is \code{gaussian} or
#'        \code{student.t}
#'      }
#'      \item{
#'        \code{t_df} estimated degrees of freedom of the t distribution, when
#'        \code{familyY} is \code{student.t}
#'      }
#'      \item{
#'        \code{nuY} estimated shape parameter, when \code{familyY} is
#'        \code{Gamma}. The gamma distribution is parameterized according to
#'        McCullagh & Nelder (1989, p. 30)
#'      }
#'      \item{
#'        \code{concomitant} a list with estimated concomitant variables
#'        parameters for each mixture component
#'      }
#'      \item{
#'        \code{normal.d, multinomial.d, poisson.d, binomial.d} marginal
#'        distribution of concomitant variables
#'      }
#'      \item{\code{normal.mu} mixture component means for \code{Xnorm}}
#'      \item{
#'        \code{normal.Sigma} mixture component covariance matrices for
#'        \code{Xnorm}
#'      }
#'      \item{\code{normal.model} models fitted for \code{Xnorm} }
#'      \item{
#'        \code{multinomial.probs} multinomial distribution probabilities for
#'        \code{Xmult}
#'      }
#'      \item{\code{poisson.lambda} lambda parameters for \code{Xpois}}
#'      \item{\code{binomial.p} binomial probabilities for \code{Xbin}}
#'    }
#'
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
  formulaY = NULL,
  familyY = stats::gaussian,
  data = NULL,
  Xnorm = NULL,
  Xbin = NULL,
  Xpois = NULL,
  Xmult = NULL,
  modelXnorm = NULL,
  Xbtrials = NULL,
  k = 1:3,
  initialization = c(
    "random.soft", "random.hard", "kmeans", "mclust", "manual"
  ),
  start.z = NULL,
  seed = NULL,
  maxR = 1,
  iter.max = 1000,
  threshold = 1.0e-04,
  eps = 1.0e-100,
  parallel = FALSE,
  pwarning = FALSE
) {
  # Preliminary checks and init ----------------------------------------------
  if (
    is.null(formulaY) && is.null(Xnorm) && is.null(Xbin) && is.null(Xpois) &&
      is.null(Xmult)
  ) {
    stop("No data were entered.")
  }

  initialization <- match.arg(initialization)
  lm <- list(k = k)
  # Get regression terms ----------------------------------------
  if (!is.null(formulaY)) {
    familyY <- .familyY(familyY)
    if (is.null(data)) {
      data <- stats::get_all_vars(stats::as.formula(formulaY))
    } else {
      data <- stats::get_all_vars(stats::as.formula(formulaY), data)
    }
    Y <- as.matrix(stats::model.response(stats::model.frame(formulaY, data)))
    lm$familyY <- seq_along(familyY)
  } else {
    Y <- data <- familyY <- NULL
  }

  n <- unlist(
    sapply(
      list(Y, Xnorm, Xpois, Xbin, Xmult),
      function(x) nrow(cbind(x))
    )
  )
  if (max(n) != min(n)) stop("Data length mismatch")
  n <- n[1]


  # Fix variables to use for marginal --------------------------
  colXn <- colXm <- m <- colXb <- colXp <- 0
  Xmod <- list()
  if (!is.null(Xmult)) {
    Xmult <- as.data.frame(Xmult)
    colXm <- ncol(Xmult)
    for (i in 1:colXm) {
      Xmult[[i]] <- as.factor(Xmult[[i]])
      m[i] <- length(levels(Xmult[, i]))
      Xmod[[i]] <- mclust::unmap(Xmult[, i], levels(Xmult)[i])
    }
    names(Xmod) <- colnames(Xmult)
  }
  if (!is.null(Xnorm)) {
    Xnorm <- as.matrix(Xnorm)
    colXn <- ncol(Xnorm)
    if (is.null(colnames(Xnorm))) colnames(Xnorm) <- paste0("Xnorm", 1:colXn)
  }
  if (!missing(Xpois) && !is.null(Xpois)) {
    Xpois <- as.matrix(Xpois)
    colXp <- ncol(Xpois)
    if (is.null(colnames(Xpois))) colnames(Xpois) <- paste0("Xpois", 1:colXp)
  }
  if (!is.null(Xbin)) {
    Xbin <- as.matrix(Xbin)
    colXb <- ncol(Xbin)
    if (is.null(colnames(Xbin))) colnames(Xbin) <- paste0("Xbin", 1:colXb)
    if (!is.null(Xbtrials)) {
      Xbtrials <- as.vector(Xbtrials)
      if (length(Xbtrials) != colXb) {
        stop(paste0(
          "If provided, Xbtrials has to be a vector of length equal to the
          number of columns of Xbin."
        ))
      }
    } else {
      Xbtrials <- sapply(1:colXb, function(i) max(Xbin[, i]))
    }
  }
  # More checks ------------------------------------------
  if (
    (initialization == "mclust" || initialization == "kmeans") &&
      is.null(Xnorm) && is.null(Xpois) && is.null(Xbin) && is.null(data)
  ) {
    stop(paste0(
      "when initialization is '", initialization, "', numeric variables are
      needed."
    ))
  }

  if (colXn > 0) {
    if (colXn == 1) {
      modelXnormNames <- c("E", "V")
      eqmod <- data.frame(name = modelXnormNames, number = c(1, 1))
    } else {
      modelXnormNames <- c(
        "EII", "VII", "EEI", "VEI", "EVI", "VVI",
        "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV"
      )
      eqmod <- data.frame(
        name = modelXnormNames,
        number = c(1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3)
      )
    }
    if (is.null(modelXnorm)) {
      modelXnorm <- modelXnormNames
    } else {
      modelXnorm <- match.arg(modelXnorm, modelXnormNames, several.ok = TRUE)
    }
    lm$modelXnorm <- modelXnorm
  } else {
    modelXnorm <- eqmod <- NULL
  }

  k <- as.integer(ceiling(k))
  par <- list()
  # Compute ----------------------------------------------
  mm <- expand.grid(lm)
  if (!is.null(modelXnorm) && 1 %in% k) {
    mm <- merge(mm, eqmod, by.x = "modelXnorm", by.y = "name")
    mm1 <- mm[mm$k == 1, ]
    mm2 <- mm1[!duplicated(subset(mm1, select = -modelXnorm, drop = FALSE)), ]

    if (nrow(mm1) > nrow(mm2)) {
      cat(
        "\nWith k=1, some models in modelXnorm are equivalent, so only the
        first is estimated.\n"
      )
    }
    mm <- rbind(mm2, mm[mm$k != 1, ])
  }
  mm <- mm[order(mm$k), , drop = FALSE]
  job <- function(i) {
    cat("\nEstimating model")
    cat(paste0(" with k=", mm$k[i]))
    if (!is.null(mm$modelXnorm[i])) {
      cat(paste0(", Xnorm=", mm$modelXnorm[i]))
    }
    if (!is.null(mm$familyY[i])) {
      cat(paste0(", familyY=", familyY[[mm$familyY[i]]][1]))
    }
    cat(" ")

    cwm2(
      formulaY = formulaY, data = data, Y = Y,
      Xnorm = Xnorm, Xmult = Xmult, Xbin = Xbin, Xpois = Xpois,
      Xbtrials = Xbtrials, n = n, m = m,
      colXn = colXn, colXp = colXp, colXb = colXb, colXm = colXm,
      Xmod = Xmod,
      k = mm$k[i], modelXnorm = mm$modelXnorm[i],
      familyY = familyY[[mm$familyY[i]]], method = "Nelder-Mead",
      initialization = initialization,
      start.z = start.z, iter.max = iter.max, threshold = threshold,
      maxR = maxR, eps = eps, pwarning = pwarning
    )
  }
  start.z <- start.z
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
    formulaY = formulaY,
    data = data,
    concomitant = list(
      Xnorm = as.df(Xnorm, "Xnorm"), Xbin = as.df(Xbin, "Xbin"),
      Xpois = as.df(Xpois, "Xpois"), Xmult = as.df(Xmult, "Xmult")
    ),
    Xbtrials = Xbtrials,
    models = par
  ), class = "cwm")
  print(res)
  invisible(res)
}

as.df <- function(ma, type) {
  if (!is.null(ma)) {
    df <- data.frame(ma)
    for (i in seq_len(ncol(ma))) attr(df[[i]], "type") <- type
    df
  } else {
    NULL
  }
}

student.t <- function(link = "identity") {
  x <- stats::gaussian(link)
  x$family <- "student.t"
  x
}

cwm2 <- function(
  formulaY, data, Y,
  Xnorm, Xmult, Xpois, Xbin,
  colXm, colXn, colXb, colXp,
  Xbtrials, n, m, Xmod, k, modelXnorm, familyY,
  method, initialization, start.z, iter.max, threshold, loglikplot,
  maxR, eps, pwarning
) {
  # colXm number of factors in data
  # colXn number of Xnorm variables in data
  # m vector with number of levels of factors
  t_df <- vY <- NULL # only for the t distribution
  familyYname <- familyY$family

  if (initialization == "manual" || initialization == "mclust") maxR <- 1


  # Xnorm Parameters definition ------------------------------------------------
  if (colXn > 0) {}
  # Parameters definition ------------------------------------------------------
  prior <- numeric(k) # weights


  if (!is.null(data)) {
    if (familyYname == "student.t") {
      vY <- array(
        1, c(n, k),
        dimnames = list(
          1:n,
          paste("comp.", 1:k, sep = "")
        )
      )
      t_df <- rep(30, k) # only for the t distribution
    }
    if (familyYname == "binomial" && !is.factor(Y)) {
      mY <- rowSums(Y)
    } else {
      mY <- rep(1, nrow(Y))
    }
  }
  result <- list()
  # Do Estimation --------------------------------------------------------------
  for (r in 1:maxR) {
    a <- try(
      {
        if (maxR > 1) cat(paste("\nrandom init.", r))
        z <- .postInit(
          initialization = initialization, Y = Y,
          Xnorm = Xnorm, Xpois = Xpois, Xbin = Xbin, n = n, k = k,
          start.z = start.z
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
          nj <- colSums(z) # for multinomial variables
          # Posterior probabilities
          norm <- .PX_norm(
            colXn = colXn, Xnorm = Xnorm, modelXnorm = modelXnorm,
            z = z, k = k, n = n, eps = eps
          )
          bin <- .PX_bin(
            k = k, X = Xbin, weights = z, Xbtrials = Xbtrials, n = n
          )
          pois <- .PX_Poisson(k = k, X = Xpois, weights = z, n = n)
          multi <- .PX_multi(
            colXm = colXm, Xmod = Xmod, z = z, nj = nj, m = m, k = k, n = n
          )

          # Y|x pag 45 of McCullag & Nelder (1989) -----------------------------
          if (!is.null(data)) {
            l <- do.call(
              paste0(".familyY.", familyYname),
              list(
                familyY = familyY, k = k, Y = Y, mY = mY, n = n, data = data,
                z = z, method = method, vY = vY, t_df = t_df,
                formulaY = formulaY
              )
            )
          } else {
            l <- list(PY = 1)
          }
          # Global - Observed-data log-likelihood ##
          dens <- matrix(rep(prior, n), n, k, byrow = TRUE) *
            norm$PX * multi$PX * pois$PX * bin$PX * l$PY
          loglik[iteration] <- sum(log(rowSums(dens)))

          # Aitkane's Acceleration-Based Stopping Criterion # #

          if (iteration > 2 & k > 1) {
            if (
              !is.na(loglik[iteration - 1] - loglik[iteration - 2]) &&
                !is.na(loglik[iteration])
            ) {
              if (abs(loglik[iteration - 1] - loglik[iteration - 2]) > 0) {
                a[iteration - 1] <-
                  (loglik[iteration] - loglik[iteration - 1]) /
                    (loglik[iteration - 1] - loglik[iteration - 2])
                aloglik[iteration] <- loglik[iteration - 1] +
                  (1 / (1 - a[iteration - 1]) *
                    (loglik[iteration] - loglik[iteration - 1]))
                if (abs(aloglik[iteration] - loglik[iteration]) < threshold) {
                  check <- 1
                }
              } else {
                check <- 1
              }
            }
          }
          if (iteration == iter.max | k == 1) check <- 1

          cat("*")
          iteration <- iteration + 1

          # E-Step -------------------------------------------------------------
          dens[dens == 0] <- 1E-320
          z <- dens / matrix(rep(rowSums(dens), k), ncol = k) # (n x k)

          # z adjustment
          if (min(colSums(z)) <= colXn) {
            z <- (z + 0.0000001) / rowSums(z + 0.0000001)
          }

          if (!is.null(data)) {
            if (familyYname == "student.t") {
              t_df <- l$t_df
              for (j in 1:k) {
                vY[, j] <-
                  (t_df[j] + 1) /
                    (t_df[j] + (Y - stats::fitted(l$lmodelY[[j]]))^2 / l$sig[j])
              }
            }
          }
        }

        finalloglik <- loglik[iteration - 1]

        # EM-algorithm is finished ---------------------------------------------

        # Check the number of parameters -------------------------------
        df_model <- length(stats::coef(l$lmodelY[[1]])) * k
        if (!is.null(data)) {
          if (familyYname == "gaussian" | familyYname == "Gamma") {
            df_model <- df_model + k
          }
          if (familyYname == "student.t") df_model <- df_model + 2 * k
        }
        if (!is.null(modelXnorm)) {
          df_model <-
            df_model + ncovpar(
              modelname = strtrim(paste0(modelXnorm, "II"), 3), p = colXn, G = k
            )
        }
        df <- df_model + (k - 1) + colXb * k + colXp * k + colXn * k +
          k * (sum(m) - colXm)
        if (df >= n) {
          warning(paste0(
            df, "parameters with only ", n,
            " observations; AICc and AICu are set to -Inf."
          ), call. = FALSE)
        }


        # Classification Matrix ------------------------------------------------
        cluster <- apply(z, 1, which.max)

        # Information criteria -------------------------------------------------
        IC <- list()
        IC$AIC <- 2 * finalloglik - df * 2
        IC$BIC <- 2 * finalloglik - df * log(n)
        IC$AIC3 <- 2 * finalloglik - df * 3
        IC$CAIC <- 2 * finalloglik - df * (1 + log(n))
        IC$AWE <- 2 * finalloglik - 2 * df * (3 / 2 + log(n))
        if (n - df - 1 > 0) {
          IC$AICc <- IC$AIC - (2 * df * (df + 1)) / (n - df - 1)
          IC$AICu <- IC$AICc - n * log(n / (n - df - 1))
        } else {
          IC$AICc <- IC$AICu <- -Inf
        }

        z.const <- (z < 10^(-322)) * 10^(-322) + (z > 10^(-322)) * z
        # vincolo per evitare i NaN nel calcolo di tau*log(tau)
        hard.z <-
          (matrix(rep(apply(z, 1, max), k), n, k, byrow = FALSE) == z) * 1
        ECM <- sum(hard.z * log(z.const))
        IC$ICL <- IC$BIC + ECM
        ## GLM parameters
        GLModel <- NULL
        if (!is.null(data)) {
          GLModel <- list()
          for (i in seq_len(k)) {
            parameters <- list(model = l$lmodelY[[i]])
            if (familyYname == "gaussian" | familyYname == "student.t") {
              parameters <- c(parameters, list(sigma = sqrt(l$VarY[1, i])))
            }
            if (familyYname == "student.t") {
              parameters <- c(parameters, list(t_df = t_df[i]))
            }
            if (familyYname == "Gamma") {
              parameters <- c(parameters, list(nuY = l$nuY[i]))
            }
            GLModel[[i]] <- parameters
          }
          names(GLModel) <- paste0("comp.", 1:k)
        }
        ## Concomitant
        concomitant <- NULL
        if (!is.null(Xnorm)) {
          concomitant <- list()
          normal <- list(
            d = norm$PX,
            mu = norm$mu,
            Sigma = norm$Sigma,
            model = modelXnorm
          )
          concomitant <- c(concomitant, normal = normal)
        }
        if (!is.null(Xmult)) {
          multinomial <- list(
            d         = multi$PX,
            probs     = multi$alpha
          )
          concomitant <- c(concomitant, multinomial = multinomial)
        }
        if (!is.null(Xpois)) {
          poisson <- list(
            d         = pois$PX,
            lambda    = pois$lambda
          )
          concomitant <- c(concomitant, poisson = poisson)
        }
        if (!is.null(Xbin)) {
          binomial <- list(
            d    = bin$PX,
            p    = bin$p
          )
          concomitant <- c(concomitant, binomial = binomial)
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
          IC = IC,
          converged = TRUE, # Logical, TRUE if EM algorithm converged.
          GLModel = GLModel,
          concomitant = concomitant
        )
      },
      silent = !pwarning
    )
    if (!length(a) > 1) {
      result[[r]] <- list(
        iter = iteration,
        k = k,
        IC = list(BIC = NA)
      )
      cat(" Not estimated")
    }
  } # maxr
  m <- which.max(sapply(1:maxR, function(r) result[[r]]$IC$BIC))
  if (length(m) == 0) m <- 1
  result <- result[[m]]
  class(result) <- "cwm"

  result
}

ncovpar <- function(modelname = NULL, p = NULL, G = NULL) {
  if (is.null(p)) {
    stop("p is null")
  }
  if (is.null(G)) {
    stop("G is null")
  }
  if (is.null(modelname)) {
    stop("modelname is null")
  }
  if (modelname == "EII") {
    npar <- 1
  } else if (modelname == "VII") {
    npar <- G
  } else if (modelname == "EEI") {
    npar <- p
  } else if (modelname == "VEI") {
    npar <- p + G - 1
  } else if (modelname == "EVI") {
    npar <- p * G - G + 1
  } else if (modelname == "VVI") {
    npar <- p * G
  } else if (modelname == "EEE") {
    npar <- p * (p + 1) / 2
  } else if (modelname == "EEV") {
    npar <- G * p * (p + 1) / 2 - (G - 1) * p
  } else if (modelname == "VEV") {
    npar <- G * p * (p + 1) / 2 - (G - 1) * (p - 1)
  } else if (modelname == "VVV") {
    npar <- G * p * (p + 1) / 2
  } else if (modelname == "EVE") {
    npar <- p * (p + 1) / 2 + (G - 1) * (p - 1)
  } else if (modelname == "VVE") {
    npar <- p * (p + 1) / 2 + (G - 1) * p
  } else if (modelname == "VEE") {
    npar <- p * (p + 1) / 2 + (G - 1)
  } else if (modelname == "EVV") {
    npar <- G * p * (p + 1) / 2 - (G - 1)
  } else {
    stop("modelname is not correctly defined")
  }

  npar
}

.familyY <- function(familyY) {
  if (!is.list(familyY) || !is.null(familyY$family)) {
    familyY <- list(familyY = familyY)
  }

  for (i in seq_len(familyY)) {
    if (is.character(familyY[[i]])) {
      familyYname <- match.arg(
        familyY[[i]],
        c(
          "gaussian", "poisson", "binomial", "Gamma", "student.t",
          "inverse.gaussian"
        )
      )
      familyY[[i]] <- do.call(familyYname, list())
    } else if (is.function(familyY[[i]])) {
      familyY[[i]] <- familyY[[i]]()
    }
  }

  familyY
}

.postInit <- function(initialization, Y, Xnorm, Xpois, Xbin, n, k, start.z) {
  if (k == 1) {
    z <- matrix(rep(1, n), ncol = 1)
    return(z)
  }
  if (initialization == "mclust") {
    data <- cbind(Y[, 1], Xnorm, Xpois, Xbin)
    z <- mclust::Mclust(
      data = data, modelNames = ifelse(ncol(data) == 1, "V", "VVV"), G = k
    )$z
  }
  if (initialization == "kmeans") {
    clusters <- stats::kmeans(x = cbind(Y, Xnorm, Xpois, Xbin), centers = k)
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
    z <- start.z # posterior probabilities (n x k) no-normalized
  }
  z <- (z + 0.0001) / rowSums(z + 0.0001)
  z
}

.PX_multi <- function(colXm, Xmod, z, nj, m, k, n) {
  alpha <- NULL
  if (colXm > 0) {
    PXmult <- array(
      1, c(n, k),
      dimnames = list(1:n, paste("comp.", 1:k, sep = ""))
    )
    for (u in 1:colXm) {
      alpha[[u]] <- array(
        (t(Xmod[[u]]) %*% z) / matrix(rep(nj, m[u]), nrow = m[u], byrow = TRUE),
        c(m[u], k),
        dimnames = list(
          paste("level.", 1:m[u], sep = ""), paste("comp.", 1:k, sep = "")
        )
      )
    }
    names(alpha) <- names(Xmod)
    for (j in 1:k) {
      for (i in 1:n) {
        PXmult[i, j] <- 1
        for (u in 1:colXm) {
          PXmult[i, j] <- PXmult[i, j] *
            stats::dmultinom(Xmod[[u]][i, ], size = 1, prob = alpha[[u]][, j])
        }
      }
    }
  } else {
    PXmult <- 1
  }
  list(PX = PXmult, alpha = alpha)
}
.PX_norm <- function(colXn, Xnorm, modelXnorm, z, k, n, eps) {
  if (colXn > 0) {
    PXnorm <- array(
      1, c(n, k),
      dimnames = list(1:n, paste("comp.", 1:k, sep = ""))
    )
    muX <- array(
      0, c(colXn, k),
      dimnames = list(dimnames(Xnorm)[[2]], paste("comp.", 1:k, sep = ""))
    )
    Sigma <- array(
      0, c(colXn, colXn, k),
      dimnames = list(
        dimnames(Xnorm)[[2]],
        dimnames(Xnorm)[[2]],
        paste("comp.", 1:k, sep = "")
      )
    )
    if (colXn > 1) { #  More than one normal concomitant
      fitM <- ContaminatedMixt::m.step(
        X = Xnorm, modelname = modelXnorm, z = z, mtol = 1e-25, mmax = 10
      )
      muX <- fitM$mu
      Sigma <- fitM$Sigma
      PXnorm <- fitM$px
    } else {
      fitM <- mclust::mstep(
        modelName = substr(modelXnorm, 1, 1), data = Xnorm, z = z
      )
      muX[1, ] <- fitM$parameters$mean
      if (substr(modelXnorm, 1, 1) == "E") {
        Sigma[1, 1, ] <- rep(fitM$parameters$variance$sigmasq, k)
      }
      if (substr(modelXnorm, 1, 1) == "V") {
        Sigma[1, 1, ] <- fitM$parameters$variance$sigmasq
      }
      for (j in 1:k) {
        PXnorm[, j] <- stats::dnorm(Xnorm, muX[1, j], sqrt(Sigma[1, 1, j]))
      }
    }
  } else {
    PXnorm <- 1
    muX <- Sigma <- NULL
  }
  list(PX = PXnorm, mu = muX, Sigma = Sigma)
}

.fixSigma <- function(sigma, eps) {
  es <- eigen(sigma)
  es$values[es$values < eps] <- eps
  es$vectors %*% diag(es$values) %*% t(es$vectors)
}

.PX_Poisson <- function(k, X, weights, n) {
  if (!is.null(X)) {
    lambdaX <- array(
      0, c(ncol(X), k),
      dimnames = list(dimnames(X)[[2]], paste("comp.", 1:k, sep = ""))
    )
    PX <- array(1, c(n, k), dimnames = list(1:n, paste("comp.", 1:k, sep = "")))
    for (i in seq_len(ncol(X))) {
      for (h in 1:k) {
        lambdaX[i, h] <- .PoissonX(X[, i], weights = weights[, h])$lambdaX
        PX[, h] <- stats::dpois(X[, i], lambda = lambdaX[i, h])
      }
    }
    list(PX = PX, lambda = lambdaX)
  } else {
    list(PX = 1, lambda = NULL)
  }
}

.PX_bin <- function(k, X, weights, Xbtrials, n) {
  if (!is.null(X)) {
    p <- array(
      0, c(ncol(X), k),
      dimnames = list(dimnames(X)[[2]], paste("comp.", 1:k, sep = ""))
    )
    PX <- array(1, c(n, k), dimnames = list(1:n, paste("comp.", 1:k, sep = "")))
    for (i in seq_len(ncol(X))) {
      for (h in 1:k) {
        p[i, h] <- .BinomialX(X[, i], weights = weights[, h], m = Xbtrials[i])$p
        # stima p
        PX[, h] <- stats::dbinom(X[, i], size = Xbtrials[i], prob = p[i, h])
      }
    }
    list(PX = PX, p = p)
  } else {
    list(PX = 1, p = NULL)
  }
}

.BinomialX <- function(X, weights, m) {
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the
  # weights to use for elements of X.

  f <- function(par, X, weights, m) {
    -sum(weights * stats::dbinom(x = X, size = m, prob = par, log = TRUE))
  }

  res <- stats::optimize(
    f = f, interval = c(0, 1), X = X, weights = weights, m = m
  )

  list(loglik = -res$objective, p = res$minimum)
}

.PoissonX <- function(X, weights) {
  # X: a numerical vector for the unique covariate
  # weights: a numerical vector of weights the same length as X giving the
  # weights to use for elements of X.
  lambda.obs <- stats::weighted.mean(x = X, w = weights)

  f <- function(par, X, weights) {
    par
  }

  res <- stats::optimize(
    f = f, interval = c(0, 2 * lambda.obs), X = X, weights = weights
  )

  loglik <- -res$objective
  lambda.hat <- res$minimum

  list(loglik = loglik, lambdaX = lambda.hat)
}
