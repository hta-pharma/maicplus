#' Generate time-to-event data
#'
#' Function generates time-to-event data. We follow the data-generating mechanism from
#' Bender et al. and Remiro-Azocar et al. 2021 to simulate Weibull-distribued survival times
#' under a proportional hazards parameterization. As per Remiro-Azocar et al, we set the
#' default scale parameter (lambda) to 8.5 and shape parameter (nu) to 1.3 to reflect frequently
#' observed mortality trends in metastatic cancer patients. Censoring times are generated from
#' exponential distribution where the rate parameter is 0.96 is selected to achieve a censoring rate
#' of 35% under the active treatment at baseline (with the values of the covariates set to zero).
#'
#' @param N sample size
#' @param ncov number of covariates to generate
#' @param lambda scale parameter in baseline hazard function
#' @param nu shape parameter in baseline hazard function
#' @param censor_rate rate parameter of the exponential distribution of censoring times
#' @param cor correlation between covariates
#' @param marginal_mean marginal mean of covariates
#' @param marginal_sd marginal standard deviation of covariates
#' @param beta vector of coefficients for prognostic factors. Defaults is `-log(0.67)` for all covariates.
#' @param gamma vector of coefficients for effect modifiers. Default is `0` for half of the covariates and
#'  `-log(0.67)` for the others.
#' @param d0 coefficient for the treatment effect for patients with `treat = 0`
#' @param seed seed used for simulation reproducibility
#'
#' @return data.frame with survival time, censor (event) status, treatment indicator, and covariates
#'
#' @references Remiro-Azocar et al. 2021
#'
#' @export

generate_survival_data <- function(N = 150, ncov = 4, lambda = 8.5, nu = 1.3, censor_rate = 0.96,
                                   cor = 0.35, marginal_mean = 0.6, marginal_sd = 0.2,
                                   beta = rep(-log(0.67), ncov),
                                   gamma = sort(rep(c(0, -log(0.67)), length = ncov)),
                                   d0 = log(0.25), seed = 1) {
  set.seed(seed)

  # if beta and gamma is not specified, give it default values
  if (is.null(beta)) {
    beta <- rep(-log(0.67), ncov)
  }

  if (is.null(gamma)) {
    gamma <- rep(-log(0.67), ncov)
    gamma[1:round(ncov / 2)] <- 0
  }

  # generate ncov covariates (X) using a multivariate normal
  sigma <- matrix(cor, nrow = ncov, ncol = ncov)
  diag(sigma) <- 1
  sigma <- sigma * marginal_sd^2

  covariates <- MASS::mvrnorm(N, mu = rep(marginal_mean, ncov), Sigma = sigma)

  treat <- sample(x = c(0, 1), size = N, replace = TRUE, prob = c(0.5, 0.5))

  # linear predictor
  lp <- covariates %*% beta + covariates %*% gamma * treat + d0 * treat

  # Weibull latent event times
  u <- runif(n = N)
  latent_time <- (-log(u) / (lambda * exp(lp)))^(1 / nu)

  # exponential censoring times
  censor_time <- rexp(n = N, rate = censor_rate)

  # follow-up times and event indicators
  time <- pmin(latent_time, censor_time)
  event <- as.numeric(latent_time <= censor_time)

  colnames(covariates) <- paste0("x", seq_len(ncov))
  survival_data <- data.frame(USUBJID = seq_len(N), time = time, event = event, treat = treat, x = covariates)
  return(survival_data)
}
