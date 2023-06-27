
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
#' @param lambda scale parameter in baseline hazard function
#' @param nu shape parameter in baseline hazard function
#' @param rateC rate parameter of the exponential distribution of censoring times
#' @param cor_X correlation between covariates
#' @param marginal_X_mean marginal mean of covariates
#' @param marginal_X_sd marginal standard deviation of covariates
#' @param beta vector of coefficients for prognostic factors
#' @param gamma vector of coefficients for effect modifiers
#' @param d0 coefficient for the treatment effect when X = 0
#' @param seed seed used for simulation reproducibility
#'
#' @return data.frame with survival time, censor (event) status, treatment indicator, and covariates
#'
#' @references Remiro-Azocar et al. 2021
#'
#' @export


generate_survival_data <- function(N = 150, ncov = 4, lambda = 8.5, nu = 1.3, rateC = 0.96, 
                                   cor_X = 0.35, marginal_X_mean = 0.6, marginal_X_sd = 0.2,
                                   beta = NULL, gamma = NULL, d0 = log(0.25), seed = 1){
  
  set.seed(seed)
  
  #if beta and gamma is not specified, give it default values
  if(is.null(beta)){
    beta <- rep(-log(0.67), ncov)
  }
  
  if(is.null(gamma)){
    gamma <- rep(-log(0.67), ncov)
    gamma[1:round(ncov/2)] <- 0
  }
  
  # generate ncov covariates (X) using a multivariate normal
  Sigma <- matrix(cor_X, nrow = ncov, ncol = ncov)
  diag(Sigma) <- 1
  Sigma <- Sigma * marginal_X_sd^2 
  
  x <- MASS::mvrnorm(N, mu = rep(marginal_X_mean, ncov), Sigma = Sigma)
  
  treat <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))
  
  # linear predictor
  lp <- x %*% beta + x %*% gamma * treat + d0 * treat
  
  # Weibull latent event times
  u <- runif(n=N)
  Tlat <- (-log(u) / (lambda * exp(lp)))^(1 / nu)
  
  # exponential censoring times
  C <- rexp(n=N, rate=rateC)
  
  # follow-up times and event indicators
  Time <- pmin(Tlat, C)
  Event <- as.numeric(Tlat <= C)
  
  survival_data <- data.frame(USUBJID=1:N, Time=Time, Event=Event, treat = treat, x = x)
  colnames(survival_data)[(ncol(survival_data)-ncov+1):ncol(survival_data)] <- paste0("x", 1:ncov)  
  return(survival_data)
}