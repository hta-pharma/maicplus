
#' Preprocess data
#'
#' Preprocess data before estimating weights for matching-adjusted indirect comparison (MAIC).
#'
#' @param intervention_input A data frame containing individual patient data from
#'   the intervention study.
#' @param target_pop A data frame containing aggregate dataset for the target population. 
#'   Variables are followed by one of the following suffixes to denote the type of summary: 
#'   varname_mean, varname_sd, varname_median, varname_prop. 
#'   After preprocessing these summary suffixes, intervention_input is 
#'   centered using the aggregate data averages. 
#' @return A list containing 2 objects. First, a data frame named analysis_data
#'   containing intervention_data with additional columns named wt (weights) and
#'   wt_rs (rescaled weights). Second, a vector called matching_vars of the
#'   names of the centered matching variables used.
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#'
#' @export

preprocess_data <- function(intervention_input, target_pop){
  
  # Check intervention_data and target_pop are data frame and match_cov is a character vector
  if(!is.data.frame(intervention_input)){stop("intervention_input is expected to be a data frame")}
  if(!is.data.frame(target_pop)){stop("target_pop is expected to be a data frame")}
  
  # Check if target_pop is 1 row of aggregate data
  if(nrow(target_pop)!=1){stop("target_pop should have exactly 1 row")}
  
  # Strip off naming convention in the aggregate data
  varnames <- gsub("_([^_]+)$","", names(target_pop))
  vartype <- gsub("^.*_","", names(target_pop))
  
  #Preprocess standard deviation
  for(i in 1:dim(target_pop)[2]){
    if(vartype[i] == "sd"){
      
      # retrieve which variable sd was specified
      varwithsd <- varnames[i]
      
      if(!paste0(varwithsd, "_mean") %in% names(target_pop)){
        stop(paste0("Also need to provide mean for ", varwithsd, " when specifying sd"))
      }
      
      # derive squared mean term
      target_pop[,paste0(varwithsd, "_squared_mean")] <- target_pop[,paste0(varwithsd, "_mean")]^2 + target_pop[,paste0(varwithsd, "_sd")]^2
      
      # remove standard deviation from the data frame
      target_pop <- target_pop[,-which(colnames(target_pop) == paste0(varwithsd, "_sd"))]
    }
  }
  
  # Preprocess median
  for(i in 1:dim(target_pop)[2]){
    if(vartype[i] == "median"){
      
      # retrieve which variable median was specified
      varwithmedian <- varnames[i]
      
      # make median into binary category
      intervention_input[,varwithmedian] <- ifelse(intervention_input[,varwithmedian] > target_pop[,paste0(varwithmedian, "_median")], 1, 0)
      target_pop[,paste0(varwithmedian, "_median")] <- 0.5
    }
  }
  
  # Remove everything that is not mean, median, or prop from target_pop
  if(!is.null(target_pop$N)){
    N <- target_pop$N  
  }
  
  vartype <- gsub("^.*_","", names(target_pop))
  target_pop <- target_pop[,which(vartype %in% c("mean", "median", "prop"))]
  
  varnames <- gsub("_([^_]+)$","", names(target_pop))
  if(any(duplicated(varnames))){stop("Cannot have more than 1 summary stat for each variable")}
  names(target_pop) <- varnames
  
  # intervention_input is centered using the aggregate data averages.
  intervention_data <- intervention_input
  for(i in varnames){
    intervention_data[,paste0(i, "_centered")] <- intervention_input[,i] - target_pop[,i]
  }
  
  # Add back in N
  if(!is.null(N)){
    target_pop$N <- N 
  }
  
  return(list(intervention_data = intervention_data, target_pop = target_pop))
}



#' Estimate MAIC propensity weights
#'
#' Estimate propensity weights for matching-adjusted indirect comparison (MAIC).
#'
#' @param intervention_data A data frame containing individual patient data from
#'   the intervention study. Intervention_data is assumed to have been preprocessed using 
#'   preprocess_data (i.e. centered using aggregate data means)
#' @param match_cov A character vector giving the names of the covariates to
#'   use in matching. These names must match the column names in intervention_data.
#' @param method The method used for optimisation - The default is method =
#'   "BFGS". Refer to \code{\link[stats]{optim}} for options.
#' @param startVal a scalar, the starting value for all coefficients of the propensity score 
#'   regression
#' @param ... Additional arguments to be passed to optimisation functions such
#'   as the method for maximum likelihood optimisation. Refer to \code{\link[stats]{optim}} 
#'   for options.
#' @return a list with 4 elements,
#' \describe{
#'   \item wt - a numeric vector of unscaled individual weights.
#'   \item wt.rs - a numerical vector of rescaled individual weights, with summation equaling to sample size (# rows of input \code{EM})
#'   \item ess - effective sample size, square of sum divided by sum of squares
#'   \item opt - R object returned by \code{base::optim()}, for assess convergence and other details
#' }
#' @references NICE DSU TECHNICAL SUPPORT DOCUMENT 18: METHODS FOR
#'   POPULATION-ADJUSTED INDIRECT COMPARISONS IN SUBMSISSIONS TO NICE, REPORT BY
#'   THE DECISION SUPPORT UNIT, December 2016
#' @seealso \code{\link{optim}}
#' @export

estimate_weights <- function(intervention_data, match_cov, startVal = 0, 
                             method = "BFGS", ...){
  
  # Check intervention_data is a data frame and match_cov is a character vector
  if(!is.data.frame(intervention_data)){stop("intervention_data is expected to be a data frame")}
  if(!is.character(match_cov)){stop("match_cov is expected to be a character vector")}
  
  # Check if match_cov name is included in the IPD data
  
  
  # Check if there is any missingness in intervention_data
  missing <- apply(intervention_data[,match_cov], 1, function(x) any(is.na(x)))
  if(any(missing)){
    stop(paste0("Following rows have missing values: ", paste(which(missing), collapse = ",")))
  } 
  
  for(i in match_cov){
    # Check that match_vars is in one of the columns of intervention_data
    if(!paste0(i, "_centered") %in% colnames(intervention_data)){
      stop(paste0("Variable ", i, " is not one of intervention_data column names"))
    }
    
    # Check whether intervention_data has not been centered by the aggregate data means
    # by looking at whether binary variables have only values of 0 and 1 
    if(all(unique(intervention_data[,i]) == 2 & unique(intervention_data[,i]) %in% c(0,1))){
      stop("intervention_data does not seem to be centered by the aggregate data means")
    }
  }
  
  # Objective function
  objfn <- function(a1, X){
    sum(exp(X %*% a1))
  }
  
  # Gradient function
  gradfn <- function(a1, X){
    colSums(sweep(X, 1, exp(X %*% a1), "*"))
  }
  
  # Optimise Q(b) using Newton-Raphson techniques
  opt1 <- stats::optim(par = rep(startVal,dim(intervention_data[,match_cov])[2]),
                       fn = objfn,
                       gr = gradfn,
                       X = as.matrix(intervention_data[,paste0(match_cov,"_centered")]),
                       method = method,
                       control = list(maxit = 300, trace = 2),
                       ...)
  
  alpha <- opt1$par
  wt <- as.vector(exp(as.matrix(intervention_data[,paste0(match_cov,"_centered")]) %*% alpha))
  wt_rs <- (wt / sum(wt)) * nrow(intervention_data)
  
  output <- list(
    wt = wt,
    wt_rs = wt_rs,
    ess = sum(wt)^2 / sum(wt^2),
    opt = opt1
  )
  return(output)
}

summarize_wts <- function(weights){
  
  with(weights, {
    summary <- data.frame(
      type = c("Weights", "Rescaled weights"),
      mean = c(mean(wt), mean(wt_rs)),
      sd = c(stats::sd(wt), stats::sd(wt_rs)),
      median = c(stats::median(wt), stats::median(wt_rs)),
      min = c(min(wt), min(wt_rs)),
      max = c(max(wt), max(wt_rs))
    )
    return(summary)
  })
}


check_weights <- function(intervention_data, target_pop, weights, match_cov){
  
  if(is.null(target_pop$N)){
    stop("Assumes target_pop stores reported sample size (N)")
  }
  
  ARM <- c("Intervention", "Intervention_weighted", "Comparator")
  ESS <- round(c(nrow(intervention_data), weights$ess,
                 target_pop$N))
  
  unweighted_cov <- intervention_data %>% summarise_at(match_cov, list(~ mean(.)))
  weighted_cov <- intervention_data %>% summarise_at(match_cov, list(~ weighted.mean(., weights$wt)))
  comparator_cov <- select(target_pop, all_of(match_cov))
  
  cov <- rbind(unweighted_cov, weighted_cov, comparator_cov)
  baseline_summary <- cbind(ARM, ESS, cov)
  
  return(baseline_summary)
}


#' Bootstrapping for MAIC weighted hazard ratios
#'
#' @param intervention_data A data frame containing individual patient data
#'   from the intervention study.
#' @param match_cov A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param model A model formula in the form 'Surv(Time, Event==1) ~ ARM'.
#'   Variable names need to match the corresponding columns in intervention_data.
#' @param comparator_input A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   The outcome variables names must match intervention_data.
#' @param min_weight A numeric value that defines the minimum weight allowed. 
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted hazard ratio (HR) from a
#'   Cox proportional hazards model. This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The HR as a numeric value.
#' @export

bootstrap_HR <- function(intervention_data, i, match_cov, comparator_input, model, min_weight = 0.0001){
  
  # Samples the data
  bootstrap_data <- intervention_data[i,]
  
  # Estimates weights
  weights <- estimate_weights(intervention_data=bootstrap_data, match_cov=match_cov)
  
  bootstrap_data$wt <- weights$wt
  bootstrap_data$ARM <- "Intervention"
  bootstrap_data <- bootstrap_data[,c("Time", "Event", "wt", "ARM")]
  
  # Give comparator data weights of 1
  comparator_input$wt <- 1
  comparator_input$ARM <- "Comparator"
  comparator_input <- comparator_input[,c("Time", "Event", "wt", "ARM")] 
  
  # Add the comparator data
  combined_data <- rbind(bootstrap_data, comparator_input)
  combined_data$ARM <- stats::relevel(as.factor(combined_data$ARM), ref="Comparator")
  
  # set weights that are below min_weight to min_weight to avoid issues with 0 values
  combined_data$wt <- ifelse(combined_data$wt < min_weight, min_weight, combined_data$wt)
  
  # survival data stat
  cox_model <- survival::coxph(model, data = combined_data, weights = wt)
  HR <- exp(cox_model$coefficients)
}




#' Bootstrapping for MAIC weighted odds ratios
#'
#' @param intervention_data  A data frame containing individual patient data
#'   from the intervention study.
#' @param i Index used to select a sample within \code{\link{boot}}.
#' @param match_cov A character vector giving the names of the covariates to use
#'   in matching. These names must match the column names in intervention_data.
#' @param comparator_data A data frame containing pseudo individual patient data
#'   from the comparator study needed to derive the relative treatment effect.
#'   The outcome variables names must match intervention_data.
#' @param model A model formula in the form 'endpoint ~ treatment_var'. Variable
#'   names need to match the corresponding columns in intervention_data.
#' @param min_weight A numeric value that defines the minimum weight allowed. 
#'   This value (default 0.0001) will replace weights estimated at 0 in a sample.
#'
#' @details This function is intended to be used in conjunction with the
#'   \code{\link{boot}} function to return the statistic to be
#'   bootstrapped. In this case by performing MAIC weighting using
#'   {\link{estimate_weights}} and returning a weighted odds ratio (OR) from a
#'   logistic regression model. This is used as the 'statistic' argument in
#'   the boot function.
#'
#' @return The OR as a numeric value.
#'
#' @seealso \code{\link{estimate_weights}}, \code{\link{boot}}
#'
#' @export

bootstrap_OR <- function(intervention_data, i, match_cov, comparator_input, model, min_weight = 0.0001){
  
  # Samples the data
  bootstrap_data <- intervention_data[i,]
  
  # Estimates weights
  weights <- estimate_weights(intervention_data=bootstrap_data, match_cov=match_cov)
  
  bootstrap_data$wt <- weights$wt
  bootstrap_data$ARM <- "Intervention"
  bootstrap_data <- bootstrap_data[,c("response", "wt", "ARM")]
  
  # Give comparator data weights of 1
  comparator_input$wt <- 1
  comparator_input$ARM <- "Comparator"
  comparator_input <- comparator_input[,c("response", "wt", "ARM")]
  
  # Add the comparator data
  combined_data <- rbind(bootstrap_data, comparator_input)
  combined_data$ARM <- stats::relevel(as.factor(combined_data$ARM), ref="Comparator")
  
  # set weights that are below min_weight to min_weight to avoid issues with 0 values
  combined_data$wt <- ifelse(combined_data$wt < min_weight, min_weight, combined_data$wt)
  
  # Perform logistic regression and extract the OR estimate
  logistic.regr <- suppressWarnings(stats::glm(formula = model, family=stats::binomial(link="logit"), data = combined_data, weight = wt))
  OR <- exp(as.numeric(stats::coef(logistic.regr)[2]))
}
