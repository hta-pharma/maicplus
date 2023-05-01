#' Derive treatment effect in the analysis step of MAIC
#'
#' @param Ntrt 
#' @param ntrt 
#' @param Ncom 
#' @param ncom 
#' @param link 
#'
#' @return
#' @export
#'
#' @examples
dummy_mod = function(Ntrt,ntrt,Ncom,ncom,link=poisson(link="identity")){
  tmpdat = data.frame(res=unlist(mapply(rep,times=c(ntrt,Ntrt-ntrt,ncom,Ncom-ncom), x=c(1,0,1,0))),
                      trt=unlist(mapply(rep,times=c(Ntrt,Ncom),x=c("Treatment","Placebo")))
  )
  tmpdat$trt = factor(tmpdat$trt,levels=c("Placebo","Treatment"))
  tmpdat$id = paste0("ID",1:nrow(tmpdat))
  mod = glm(res~trt,tmpdat,family=link)
  vmod = clubSandwich::vcovCR(mod,cluster=tmpdat$id,type="CR2")
  coef_res = clubSandwich::conf_int(mod,vmod,coef=2)
  list(mod=summary(mod),
       est=coef_res$beta*100,
       ci_l=round(coef_res$CI_L*100,2),
       ci_u=round(coef_res$CI_U*100,2),
       se = coef_res$SE*100)
}
ipd_mod = function(mod,tmpdat){
  vmod = clubSandwich::vcovCR(mod,cluster=tmpdat$USUBJID,type="CR2")
  coef_res = clubSandwich::conf_int(mod,vmod,coef=2)
  list(mod=summary(mod),
       est=coef_res$beta*100,
       ci_l=round(coef_res$CI_L*100,2),
       ci_u=round(coef_res$CI_U*100,2),
       se = coef_res$SE*100)
}

