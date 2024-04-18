idioClust
================
Ultán P. Doherty
2024-04-18

## Iterative Detection & Identification of Outliers while Clustering

- `idio_gmm` - Identify multivariate outliers while clustering the data
  with a Gaussian mixture model.
- `idio_mlr` - Identify response variable outliers while fitting a
  multiple linear regression model to the data.
- `simulate_noisy_gmm` - Simulate data from a Gaussian mixture model
  with multivariate outliers.
- `simulate_noisy_mlr` - Simulate data from a multiple linear regression
  model with response variable outliers.

``` r
devtools::load_all()
```

    ## ℹ Loading idioClust

### `simulate_noisy_gmm`

``` r
n_vec <- c(2000, 1000, 1000)
mu_list <- list(c(-1, 0), c(+1, -1), c(+1, +1))
sigma_list <- list(diag(c(0.2, 4 * 0.2)),
                  diag(c(0.2, 0.2)),
                  diag(c(0.2, 0.2)))
noisy_gmm_p2g3 <- simulate_noisy_gmm(
 n_vec, mu_list, sigma_list,
 outlier_num = 40, seed = 123, crit_val = 0.9999,
 unif_range_multiplier = 1.5
)
par(mfrow = c(1, 1))
plot(noisy_gmm_p2g3[, 1:2],
     col = 1 + noisy_gmm_p2g3[, 3], pch = 1 + noisy_gmm_p2g3[, 3])
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### `idio_gmm`

``` r
idio_gmm_p2g3 <- idio_gmm(noisy_gmm_p2g3[, 1:2], comp_num = 3, max_out = 100,
                         print_interval = Inf)
par(mfrow = c(1, 2))
plot(0:100, idio_gmm_p2g3$distrib_diffs, type = "l")
abline(v = idio_gmm_p2g3$outlier_num)
plot(noisy_gmm_p2g3[, 1:2], col = idio_gmm_p2g3$gmm_labels,
   pch = 1 + noisy_gmm_p2g3[, 3])
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### `simulate_noisy_mlr`

``` r
n_vec <- c(1000)
mu_list <- list(+1)
sigma_list <- list(as.matrix(0.1))
beta_list <- list(c(1, 1))
error_sd_vec <- c(0.1)
noisy_mlr_p1 <- simulate_noisy_mlr(n_vec, mu_list, sigma_list, beta_list,
                                   error_sd_vec,
                                   outlier_num = 20, seed = 123,
                                   crit_val = 0.9999)
par(mfrow = c(1, 1))
plot(x = noisy_mlr_p1$covariates[, 1], y = noisy_mlr_p1$responses,
     col = 1 + noisy_mlr_p1$labels, pch = 1 + noisy_mlr_p1$labels)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### `idio_mlr`

``` r
idio_mlr_p1 <- idio_mlr(noisy_mlr_p1$covariates, noisy_mlr_p1$responses,
                        max_out = 50, print_interval = Inf)
par(mfrow = c(1, 2))
plot(0:50, idio_mlr_p1$distrib_diffs, type = "l")
abline(v = idio_mlr_p1$outlier_num)
plot(x = noisy_mlr_p1$covariates[, 1], y = noisy_mlr_p1$responses,
     pch = 1 + noisy_mlr_p1$labels, col = 1 + idio_mlr_p1$outlier_bool)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
