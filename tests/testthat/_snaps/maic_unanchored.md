# test binary case

    Code
      testout$descriptive$summary
    Output
        trt_ind treatment            type   n   events events_pct
      1       B         B Before matching 400 280.0000   70.00000
      2       A         A Before matching 500 390.0000   78.00000
      3       B         B  After matching 400 280.0000   70.00000
      4       A         A  After matching 500 142.8968   28.57935

---

    Code
      testout$inferential$summary
    Output
               case       RR       LCL      UCL      pval
      1          AB 1.114286 0.9557015 1.299185 0.1671206
      2 adjusted_AB 1.023627 0.9123647 1.148457 0.6908096

---

    Code
      testout$inferential$fit
    Output
      $model_before
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = dat)
      
      Coefficients:
      (Intercept)         ARMA  
          -0.3567       0.1082  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    395.5 
      Residual Deviance: 393.5 	AIC: 1738
      
      $model_after
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = dat, 
          weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
         -0.35667      0.02335  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    295.1 
      Residual Deviance: 295 	AIC: 1145
      
      $res_AB
      $res_AB$est
      [1] 1.023627
      
      $res_AB$se
      [1] 0.06025155
      
      $res_AB$ci_l
      [1] 0.9123647
      
      $res_AB$ci_u
      [1] 1.148457
      
      $res_AB$pval
      [1] 0.6908096
      
      
      $res_AB_unadj
      $res_AB_unadj$est
      [1] 1.114286
      
      $res_AB_unadj$se
      [1] 0.08768422
      
      $res_AB_unadj$ci_l
      [1] 0.9557015
      
      $res_AB_unadj$ci_u
      [1] 1.299185
      
      $res_AB_unadj$pval
      [1] 0.1671206
      
      
      $boot_res
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      testout2$descriptive
    Output
      $summary
        trt_ind treatment            type   n   events events_pct
      1       B         B Before matching 400 280.0000   70.00000
      2       A         A Before matching 500 390.0000   78.00000
      3       B         B  After matching 400 280.0000   70.00000
      4       A         A  After matching 500 142.8968   28.57935
      

---

    Code
      testout2$inferential$fit
    Output
      $model_before
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = dat)
      
      Coefficients:
      (Intercept)         ARMA  
          -0.3567       0.1082  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    395.5 
      Residual Deviance: 393.5 	AIC: 1738
      
      $model_after
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = dat, 
          weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
         -0.35667      0.02335  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    295.1 
      Residual Deviance: 295 	AIC: 1145
      
      $res_AB
      $res_AB$est
      [1] 1.023627
      
      $res_AB$se
      [1] 0.06025155
      
      $res_AB$ci_l
      [1] 0.9123647
      
      $res_AB$ci_u
      [1] 1.148457
      
      $res_AB$pval
      [1] 0.6908096
      
      
      $res_AB_unadj
      $res_AB_unadj$est
      [1] 1.114286
      
      $res_AB_unadj$se
      [1] 0.08768422
      
      $res_AB_unadj$ci_l
      [1] 0.9557015
      
      $res_AB_unadj$ci_u
      [1] 1.299185
      
      $res_AB_unadj$pval
      [1] 0.1671206
      
      
      $boot_res
      
      STRATIFIED BOOTSTRAP
      
      
      Call:
      boot(data = boot_ipd, statistic = stat_fun, R = R, strata = weights_object$boot_strata, 
          w_obj = weights_object, pseudo_ipd = pseudo_ipd, normalize = normalize_weights)
      
      
      Bootstrap Statistics :
            original       bias     std. error
      t1* 0.02335185 1.366708e-02 0.0538052765
      t2* 0.01056949 8.109601e-05 0.0005263733
      
      $boot_res_AB
      $boot_res_AB$est
      [1] 1.023627
      
      $boot_res_AB$se
      [1] NA
      
      $boot_res_AB$ci_l
      [1] 0.9086715
      
      $boot_res_AB$ci_u
      [1] 1.122032
      
      $boot_res_AB$pval
      [1] NA
      
      

---

    Code
      testout2$inferential$summary
    Output
               case       RR       LCL      UCL      pval
      1          AB 1.114286 0.9557015 1.299185 0.1671206
      2 adjusted_AB 1.023627 0.9123647 1.148457 0.6908096

# test time to event case

    Code
      testout$descriptive$summary
    Output
        trt_ind treatment            type records    n.max  n.start    events
      1       B         B Before matching     300 300.0000 300.0000 178.00000
      2       A         A Before matching     500 500.0000 500.0000 190.00000
      3       B         B  After matching     300 300.0000 300.0000 178.00000
      4       A         A  After matching     500 199.4265 199.4265  65.68878
            rmean se(rmean)    median  0.95LCL   0.95UCL
      1  4.303551 0.3367260  2.746131 2.261125  3.320857
      2  8.709690 0.3551477  7.587627 6.278691 10.288538
      3  4.303551 0.3367260  2.746131 2.261125  3.320857
      4 10.166029 0.5499915 11.900015 7.815275 14.873786

---

    Code
      testout$inferential$summary
    Output
               case        HR       LCL       UCL         pval
      1          AB 0.3748981 0.3039010 0.4624815 5.245204e-20
      2 adjusted_AB 0.2834780 0.2074664 0.3873387 2.473442e-15

---

    Code
      testout$inferential$fit
    Output
      $km_before
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = dat, conf.type = km_conf_type)
      
              n events median 0.95LCL 0.95UCL
      ARM=B 300    178   83.6    68.8     101
      ARM=A 500    190  230.9   191.1     313
      
      $km_after
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = dat, weights = dat$weights, 
          conf.type = km_conf_type)
      
            records   n events median 0.95LCL 0.95UCL
      ARM=B     300 300  178.0   83.6    68.8     101
      ARM=A     500 199   65.7  362.2   237.9     453
      
      $model_before
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = dat)
      
              coef exp(coef) se(coef)      z      p
      ARMA -0.9811    0.3749   0.1071 -9.159 <2e-16
      
      Likelihood ratio test=80.62  on 1 df, p=< 2.2e-16
      n= 800, number of events= 368 
      
      $model_after
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = dat, weights = weights, 
          robust = TRUE)
      
              coef exp(coef) se(coef) robust se      z        p
      ARMA -1.2606    0.2835   0.1504    0.1593 -7.915 2.47e-15
      
      Likelihood ratio test=80.4  on 1 df, p=< 2.2e-16
      n= 800, number of events= 368 
      
      $res_AB
      $res_AB$est
      [1] 0.283478
      
      $res_AB$se
      [1] 0.04601759
      
      $res_AB$ci_l
      [1] 0.2074664
      
      $res_AB$ci_u
      [1] 0.3873387
      
      $res_AB$pval
      [1] 2.473442e-15
      
      
      $res_AB_unadj
      $res_AB_unadj$est
      [1] 0.3748981
      
      $res_AB_unadj$se
      [1] 0.0405065
      
      $res_AB_unadj$ci_l
      [1] 0.303901
      
      $res_AB_unadj$ci_u
      [1] 0.4624815
      
      $res_AB_unadj$pval
      [1] 5.245204e-20
      
      
      $boot_res
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      testout2$descriptive$summary
    Output
        trt_ind treatment            type records    n.max  n.start    events
      1       B         B Before matching     300 300.0000 300.0000 178.00000
      2       A         A Before matching     500 500.0000 500.0000 190.00000
      3       B         B  After matching     300 300.0000 300.0000 178.00000
      4       A         A  After matching     500 199.4265 199.4265  65.68878
            rmean se(rmean)    median  0.95LCL   0.95UCL
      1  4.303551 0.3367260  2.746131 2.261125  3.320857
      2  8.709690 0.3551477  7.587627 6.278691 10.288538
      3  4.303551 0.3367260  2.746131 2.261125  3.320857
      4 10.166029 0.5499915 11.900015 7.815275 14.873786

---

    Code
      testout2$inferential$summary
    Output
               case        HR       LCL       UCL         pval
      1          AB 0.3748981 0.3039010 0.4624815 5.245204e-20
      2 adjusted_AB 0.2834780 0.2074664 0.3873387 2.473442e-15

---

    Code
      testout2$inferential$fit
    Output
      $km_before
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = dat, conf.type = km_conf_type)
      
              n events median 0.95LCL 0.95UCL
      ARM=B 300    178   83.6    68.8     101
      ARM=A 500    190  230.9   191.1     313
      
      $km_after
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = dat, weights = dat$weights, 
          conf.type = km_conf_type)
      
            records   n events median 0.95LCL 0.95UCL
      ARM=B     300 300  178.0   83.6    68.8     101
      ARM=A     500 199   65.7  362.2   237.9     453
      
      $model_before
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = dat)
      
              coef exp(coef) se(coef)      z      p
      ARMA -0.9811    0.3749   0.1071 -9.159 <2e-16
      
      Likelihood ratio test=80.62  on 1 df, p=< 2.2e-16
      n= 800, number of events= 368 
      
      $model_after
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = dat, weights = weights, 
          robust = TRUE)
      
              coef exp(coef) se(coef) robust se      z        p
      ARMA -1.2606    0.2835   0.1504    0.1593 -7.915 2.47e-15
      
      Likelihood ratio test=80.4  on 1 df, p=< 2.2e-16
      n= 800, number of events= 368 
      
      $res_AB
      $res_AB$est
      [1] 0.283478
      
      $res_AB$se
      [1] 0.04601759
      
      $res_AB$ci_l
      [1] 0.2074664
      
      $res_AB$ci_u
      [1] 0.3873387
      
      $res_AB$pval
      [1] 2.473442e-15
      
      
      $res_AB_unadj
      $res_AB_unadj$est
      [1] 0.3748981
      
      $res_AB_unadj$se
      [1] 0.0405065
      
      $res_AB_unadj$ci_l
      [1] 0.303901
      
      $res_AB_unadj$ci_u
      [1] 0.4624815
      
      $res_AB_unadj$pval
      [1] 5.245204e-20
      
      
      $boot_res
      
      STRATIFIED BOOTSTRAP
      
      
      Call:
      boot(data = boot_ipd, statistic = stat_fun, R = R, strata = weights_object$boot_strata, 
          w_obj = weights_object, pseudo_ipd = pseudo_ipd, normalize = normalize_weights)
      
      
      Bootstrap Statistics :
             original       bias    std. error
      t1* -1.26062079 0.0024513461 0.131388233
      t2*  0.02536718 0.0005819358 0.002670424
      
      $boot_res_AB
      $boot_res_AB$est
      [1] 0.283478
      
      $boot_res_AB$se
      [1] NA
      
      $boot_res_AB$ci_l
      [1] 0.2185832
      
      $boot_res_AB$ci_u
      [1] 0.3658412
      
      $boot_res_AB$pval
      [1] NA
      
      

