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
  if(!is.null(target_pop$N)){
    target_pop$N <- N 
  }
  
  return(list(intervention_data = intervention_data, target_pop = target_pop))
}


estimate_weights <- function(intervention_data, match_cov, startVal = 0, 
                             method = "BFGS", ...){
  
  # Check intervention_data is a data frame and match_cov is a character vector
  if(!is.data.frame(intervention_data)){stop("intervention_data is expected to be a data frame")}
  if(!is.character(match_cov)){stop("match_cov is expected to be a character vector")}
  # Check if there is any missingness in intervention_data
  missing <- apply(intervention_data[,match_cov], 1, function(x) any(is.na(x)))
  if(any(missing)){
    stop(paste0("Following rows have missing values: ", paste(which(missing), collapse = ",")))
  } 
  for(i in match_cov){
    # Check that match_vars is in one of the columns of intervention_data
    if(!i %in% colnames(intervention_data)){
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


check_weights <- function(intervention_data, weights, match_cov, target_pop_N = NULL){
  
  ARM <- c("Intervention", "Intervention_weighted", "Comparator")
  ESS <- round(c(nrow(intervention_data), weights$ess,
                 target_pop_N))
  
  weighted_cov <- intervention_data %>% summarise_at(match_cov, list(~ weighted.mean(., weights$wt)))
  unweighted_cov <- intervention_data %>% summarise_at(match_cov, list(~ mean(.)))
  comparator_cov <- select(target_pop, all_of(match_cov))
  
  cov <- rbind(unweighted_cov, weighted_cov, comparator_cov)
  baseline_summary <- cbind(ARM, ESS, cov)
  
  return(baseline_summary)
}