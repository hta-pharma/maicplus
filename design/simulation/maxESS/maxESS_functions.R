# Define MAIC function.
# This function computes conventional MAIC weights.
# It also computes our new "optimal" weights.

jackson_MAIC<-function (the_data, id, match_mean, agg_means, agg_sds, accuracy=sqrt(.Machine$double.eps)) 
{
  ### MAIC is a function of 5 arguments
  ### the_data is a dataframe containing the data
  ### id is the name of the patient identifier in the dataset
  ### match_mean is a vector containing the names of the covariates used in the matching.
  ### match_mean must contain column headings of the_data (argument 1)
  ### agg_means is a vector of the aggregate data's mean values corresponding to match_mean
  ### agg_means are therefore the average covariate values used in the matching
  ### agg_sds is a corresponding vector of standard deviations: NA means match mean only
  
  the_data_o=the_data
  ### Create a copy of the original data that we can manipulate
  ### Without changing the user's arguments
  
  match_sds=which(1*(is.na(agg_sds))==0)
  also_match_var=match_mean[match_sds]
  ### We have now identified the covariates for which we also match sd/variance
  
  vars_keep=c(id, match_mean)
  the_data_m=select(the_data, all_of(vars_keep))
  complete_case=complete.cases(the_data_m)
  the_data_m=na.omit(the_data_m)
  ### the_data_m is a modified version of the_data
  ### the_data_m contains only the covariates that we match on (plus patient identifier)
  ### the_data_m only contains complete cases for these covariates
  ### this is because patients with missing covariates needed in the matching will be given weight=0
  
  if(length(also_match_var)>0)
  {
    the_data_v=select(the_data_m, also_match_var)
    the_data_v=the_data_v^2
    sq_list <- paste(also_match_var,"sq",sep="_")
    colnames(the_data_v) = sq_list
    the_data_m=cbind(the_data_m, the_data_v)
  }
  ### the_data_m now also contains extra columns with squared values of covariates that we also match sd/variance. 
  ### the column headings of the squared values have the appropriate variable name followed by _sq
  
  M1=matrix(agg_means, nrow=nrow(the_data_m), ncol=length(match_mean), byrow=T)
  ### Matrix M1 currently contains columns of data containing the average covariate values used in the matching 
  
  if(length(also_match_var)>0)
  {
    mom2=agg_means^2+agg_sds^2; mom2=mom2[!is.na(mom2)]
    M2=matrix(mom2, nrow=nrow(the_data_m), ncol=length(also_match_var), byrow=T)
    M1=cbind(M1, M2)
  }
  ### Now M1 also contains columns of containing average values of second moments
  
  the_data_m[,2:ncol(the_data_m)]=the_data_m[,2:ncol(the_data_m)] - M1
  X=as.matrix(the_data_m[,2:ncol(the_data_m)])
  
  ### Recall that the first column of the_data_m is patient identifiers and should be ignored above
  ### X is now a matrix containing CENTRED matching covariates (including squared terms).
  
  start=rep(0, ncol(X)); Answers=weight(start, X, accuracy)
  ### Answers contains the weights, calculated using the weights function 
  
  the_data_o$weight.maic=0
  the_data_o$weight.maic[which(1*(complete_case)==1)]=Answers$w
  the_data_o$weight.maic=the_data_o$weight.maic/sum(the_data_o$weight.maic)
  ### Weights for incomplete cases are 0.
  ### Weights are now calculated. 
  
  ### Now calculate optimum weights (constrained to be positive)
  
  n_complete_obs=nrow(X)
  A=diag(n_complete_obs)
  B=rep(0, n_complete_obs)
  E=rbind(rep(1,n_complete_obs),t(X))
  F=c(1, rep(0, ncol(X)))
  G=A
  H=rep(0, n_complete_obs)
  
  weights_optimal=lsei(A=A, B=B, E=E, F=F, G=G, H=H)$X
  the_data_o$weight.opt=0
  the_data_o$weight.opt[which(1*(complete_case)==1)]=weights_optimal
  return(list(the_data=the_data_o))
}

# Also define a function weight that is needed by the MAIC function to compute conventional MAIC weights.

weight<-function(start, X, accuracy) {  
  objective <- function(x, X) return( sum(exp(X %*% x)) ) 
  gradient <- function(x, X) return(t((exp(X %*% x))) %*% X)
  max <- optim(par=start, fn=objective, gr=gradient, method="BFGS", X=X, control=list(maxit=100000, reltol=accuracy))
  weights <- exp(X %*% max$par)
  L <- list(weights, max$par, max$convergence)
  names(L) <- c("w", "alpha", "converge")
  return(L)
}
