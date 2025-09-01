# Load R packages
rm(list=ls())
library(dplyr)
library(sandwich)
library(limSolve)

# Define MAIC function.
# This function computes conventional MAIC weights.
# It also computes our new "optimal" weights.

MAIC<-function (the_data, id, match_mean, agg_means, agg_sds, accuracy=sqrt(.Machine$double.eps)) 
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
  the_data_m=select(the_data, vars_keep)
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

# Perform numerical example in section 4.1. 

set.seed(1)
X<-rnorm(1000)

my.data<-data.frame(id=1:1000, X=X)

# Match to zero.
Zero<-MAIC(my.data, "id", "X", 0, NA, 0.000000000001)
# Check the matching.
sum(Zero$the_data$X*Zero$the_data$weight.maic)/sum(Zero$the_data$weight.maic)
sum(Zero$the_data$X*Zero$the_data$weight.opt)/sum(Zero$the_data$weight.opt)

# Now match to 0.5, 1, 1.5 and 2. 


point5<-MAIC(my.data, "id", "X", 0.5, NA, 0.000000000001)

sum(point5$the_data$X*point5$the_data$weight.maic)/sum(point5$the_data$weight.maic)
sum(point5$the_data$X*point5$the_data$weight.opt)/sum(point5$the_data$weight.opt)


One<-MAIC(my.data, "id", "X", 1, NA, 0.000000000001)

sum(One$the_data$X*One$the_data$weight.maic)/sum(One$the_data$weight.maic)
sum(One$the_data$X*One$the_data$weight.opt)/sum(One$the_data$weight.opt)

Onep5<-MAIC(my.data, "id", "X", 1.5, NA, 0.000000000001)

sum(Onep5$the_data$X*Onep5$the_data$weight.maic)/sum(Onep5$the_data$weight.maic)
sum(Onep5$the_data$X*Onep5$the_data$weight.opt)/sum(Onep5$the_data$weight.opt)



Two<-MAIC(my.data, "id", "X", 2, NA, 0.000000000001)

sum(Two$the_data$X*Two$the_data$weight.maic)/sum(Two$the_data$weight.maic)
sum(Two$the_data$X*Two$the_data$weight.opt)/sum(Two$the_data$weight.opt)



# Provide function for computing ESS for a set of weights. 

ESS<-function(weights)
{
  ((sum(weights))^2)/sum(weights^2)  
}  

# Produce Table of results. 
results<-matrix(nrow=4, ncol=5)
results[1,1]<-ESS(Zero$the_data$weight.maic)
results[1,2]<-ESS(point5$the_data$weight.maic)
results[1,3]<-ESS(One$the_data$weight.maic)
results[1,4]<-ESS(Onep5$the_data$weight.maic)
results[1,5]<-ESS(Two$the_data$weight.maic)

results[2,1]<-ESS(Zero$the_data$weight.opt)
results[2,2]<-ESS(point5$the_data$weight.opt)
results[2,3]<-ESS(One$the_data$weight.opt)
results[2,4]<-ESS(Onep5$the_data$weight.opt)
results[2,5]<-ESS(Two$the_data$weight.opt)

results[3,1]<-max(Zero$the_data$weight.maic)
results[3,2]<-max(point5$the_data$weight.maic)
results[3,3]<-max(One$the_data$weight.maic)
results[3,4]<-max(Onep5$the_data$weight.maic)
results[3,5]<-max(Two$the_data$weight.maic)

results[4,1]<-max(Zero$the_data$weight.opt)
results[4,2]<-max(point5$the_data$weight.opt)
results[4,3]<-max(One$the_data$weight.opt)
results[4,4]<-max(Onep5$the_data$weight.opt)
results[4,5]<-max(Two$the_data$weight.opt)


rownames(results)<-c("ESS(MAIC)", "ESS(Optimum)", "Largest w(MAIC)", "Largest w(Optimum)")
colnames(results)<-c("0", "0.5", "1", "1.5", "2")

results[1,]<-round(results[1,], 0)
results[2,]<-round(results[2,], 0)
results[3,]<-round(results[3,], 3)
results[4,]<-round(results[4,], 3)

results

# Observe the linear trends mentioned in the paper
# That is, the weights are linear in X but with some weights set to zero.
plot(X, Zero$the_data$weight.opt, xlim=c(-4,4))
plot(X, point5$the_data$weight.opt, xlim=c(-4,4))
plot(X, One$the_data$weight.opt, xlim=c(-4,4))
plot(X, Onep5$the_data$weight.opt, xlim=c(-4,4))
plot(X, Two$the_data$weight.opt, xlim=c(-4,4))


plot(Two$the_data$weight.maic, Two$the_data$weight.opt, xlim=c(0, 0.20), ylim=c(0, 0.04),
     xlab="Conventional MAIC", ylab="Alternative method", main="Alternative versus conventional weights when matching the numerical example 
     in Section 4.1 to an average covariate value of two")
lines(c(0,1), c(0,1))


sum(Two$the_data$weight.opt==0)
sum(Two$the_data$weight.maic<0.0001)

# Do "real example in section 4.2"
# Read simulated data in (this is "pre-simulated")
AB.IPD<-read.csv("AB.csv")
AC.IPD<-read.csv("AC.csv")

# Extract aggregrate level data from AC trial
AC.AgD <-
  cbind(
    # Trial level stats: mean and sd of age, number and proportion of males
    summarise(AC.IPD, age.mean = mean(age), age.sd = sd(age),
              N.male = sum(gender=="Male"), prop.male = mean(gender=="Male")),
    # Summary outcomes for A arm
    filter(AC.IPD, trt == "A") %>%
      summarise(y.A.sum = sum(y), y.A.bar = mean(y), N.A = n()),
    # Summary outcomes for C arm
    filter(AC.IPD, trt == "C") %>%
      summarise(y.C.sum = sum(y), y.C.bar = mean(y), N.C = n())
    
  )


AC.AgD

MAIC_res<-MAIC(AB.IPD, "ID", "age", as.numeric(AC.AgD[1]), as.numeric(AC.AgD[2]), 0.000000000001)

the_data<-MAIC_res$the_data

ESS(the_data$weight.maic)

ESS(the_data$weight.opt)

# Check with weights of zero removed. 
Pos<-subset(the_data, weight.opt>0)
ESS(Pos$weight.opt)
length(Pos$weight.opt)

# Produce aggregate level outcome data for competitor trial.

d.AC <- with(AC.AgD, log(y.C.sum * (N.A - y.A.sum) / (y.A.sum * (N.C - y.C.sum))))
var.d.AC <- with(AC.AgD, 1/y.A.sum + 1/(N.A - y.A.sum) + 1/y.C.sum + 1/(N.C - y.C.sum))



# Check the matching (age, mean and standard deviation)
as.numeric(AC.AgD[1])
sum(the_data$age*the_data$weight.maic)/sum(the_data$weight.maic)
sum(the_data$age*the_data$weight.opt)/sum(the_data$weight.opt)
#Mean of age is matched.
as.numeric(AC.AgD[2])
sqrt(sum(the_data$age^2*the_data$weight.maic)/sum(the_data$weight.maic) - as.numeric(AC.AgD[1])^2)

sqrt(sum(the_data$age^2*the_data$weight.opt)/sum(the_data$weight.opt) - as.numeric(AC.AgD[1])^2)

#Standard deviation of age is also matched.

# Calculate metrics
X_age<-the_data$age - as.numeric(AC.AgD[1])

X_age_sq<-the_data$age^2 - (as.numeric(AC.AgD[1])^2 + as.numeric(AC.AgD[2])^2)

metric_age<-nrow(AB.IPD)*(1-((mean(X_age))^2)/(mean(X_age^2)))

metric_age_sq<-nrow(AB.IPD)*(1-((mean(X_age_sq))^2)/(mean(X_age_sq^2)))

metric_age

metric_age_sq

# Do unadjusted analysis
the_data$y0<-1-the_data$y
fit_un<-glm(cbind(y,y0) ~ trt, family = binomial, data=the_data)
V_un <- vcov(fit_un)
print(d.AB_un <- coef(fit_un)["trtB"])
print(var.d.AB_un <- V_un["trtB","trtB"])
print(d.BC_un <- d.AC - d.AB_un)
print(var.d.BC_un <- var.d.AC + var.d.AB_un)
print(sqrt(var.d.BC_un))
c(d.BC_un-1.96*sqrt(var.d.BC_un), d.BC_un+1.96*sqrt(var.d.BC_un)) 

# Do conventional analysis
fit_maic<-glm(cbind(y,y0) ~ trt, family = binomial, weights = weight.maic, data=the_data)
V_maic <- vcovHC(fit_maic)
print(d.AB.MAIC <- coef(fit_maic)["trtB"])
print(var.d.AB.MAIC <- V_maic["trtB","trtB"])
print(d.BC.MAIC <- d.AC - d.AB.MAIC)
print(var.d.BC.MAIC <- var.d.AC + var.d.AB.MAIC)
print(sqrt(var.d.BC.MAIC))

c(d.BC.MAIC-1.96*sqrt(var.d.BC.MAIC), d.BC.MAIC+1.96*sqrt(var.d.BC.MAIC)) 

# Do novel analysis
fit_opt<-glm(cbind(y,y0) ~ trt, family = binomial, weights = weight.opt, data=subset(the_data, weight.opt>0))
V_opt <- vcovHC(fit_opt)
print(d.AB_opt <- coef(fit_opt)["trtB"])
print(var.d.AB_opt <- V_opt["trtB","trtB"])
print(d.BC_opt <- d.AC - d.AB_opt)
print(var.d.BC_opt <- var.d.AC + var.d.AB_opt)
print(sqrt(var.d.BC_opt))
c(d.BC_opt-1.96*sqrt(var.d.BC_opt), d.BC_opt+1.96*sqrt(var.d.BC_opt))

plot(the_data$weight.maic, the_data$weight.opt, xlab="Conventional MAIC", ylab="Alternative method", xlim=c(0, 0.007), ylim=c(0, 0.007), 
     main="Alternative versus conventional weights for the numerical example in Section 4.2")
lines(c(0,0.007), c(0, 0.007), lty=2)