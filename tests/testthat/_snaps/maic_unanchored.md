# test binary case

    Code
      testout_RR$descriptive$summary
    Output
        trt_ind treatment            type        n   events events_pct
      1       B         B Before matching 400.0000 280.0000   70.00000
      2       A         A Before matching 500.0000 390.0000   78.00000
      3       B         B  After matching 400.0000 280.0000   70.00000
      4       A         A  After matching 199.4265 142.8968   71.65386

---

    Code
      testout_RR$inferential$summary
    Output
               case       RR       LCL      UCL        pval
      1          AB 1.114286 1.0293724 1.206204 0.007455267
      2 adjusted_AB 1.023627 0.9123647 1.148457 0.690809607

---

    Code
      testout_RR$inferential$fit
    Output
      $model_before
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = dat)
      
      Coefficients:
      (Intercept)         ARMA  
          -0.3567       0.1082  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    1023 
      Residual Deviance: 1016 	AIC: 1020
      
      $model_after
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = dat, weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
         -0.35667      0.02335  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    726.7 
      Residual Deviance: 726.5 	AIC: 712.5
      
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
      [1] 0.04511891
      
      $res_AB_unadj$ci_l
      [1] 1.029372
      
      $res_AB_unadj$ci_u
      [1] 1.206204
      
      $res_AB_unadj$pval
      [1] 0.007455267
      
      
      $boot_res
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      testout_RD$descriptive$summary
    Output
        trt_ind treatment            type        n   events events_pct
      1       B         B Before matching 400.0000 280.0000   70.00000
      2       A         A Before matching 500.0000 390.0000   78.00000
      3       B         B  After matching 400.0000 280.0000   70.00000
      4       A         A  After matching 199.4265 142.8968   71.65386

---

    Code
      testout_RD$inferential$summary
    Output
               case       RD       LCL       UCL        pval
      1          AB 8.000000  2.224920 13.775080 0.006626293
      2 adjusted_AB 1.653865 -6.532173  9.839902 0.692119096

---

    Code
      testout_RD$inferential$fit
    Output
      $model_before
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = dat)
      
      Coefficients:
      (Intercept)         ARMA  
             0.70         0.08  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    1023 
      Residual Deviance: 1016 	AIC: 1020
      
      $model_after
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = dat, weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
          0.70000      0.01654  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    726.7 
      Residual Deviance: 726.5 	AIC: 712.5
      
      $res_AB
      $res_AB$est
      [1] 1.653865
      
      $res_AB$se
      [1] 4.176627
      
      $res_AB$ci_l
      [1] -6.532173
      
      $res_AB$ci_u
      [1] 9.839902
      
      $res_AB$pval
      [1] 0.6921191
      
      
      $res_AB_unadj
      $res_AB_unadj$est
      [1] 8
      
      $res_AB_unadj$se
      [1] 2.946523
      
      $res_AB_unadj$ci_l
      [1] 2.22492
      
      $res_AB_unadj$ci_u
      [1] 13.77508
      
      $res_AB_unadj$pval
      [1] 0.006626293
      
      
      $boot_res
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      testout_OR$descriptive$summary
    Output
        trt_ind treatment            type        n   events events_pct
      1       B         B Before matching 400.0000 280.0000   70.00000
      2       A         A Before matching 500.0000 390.0000   78.00000
      3       B         B  After matching 400.0000 280.0000   70.00000
      4       A         A  After matching 199.4265 142.8968   71.65386

---

    Code
      testout_OR$inferential$summary
    Output
               case       OR       LCL      UCL        pval
      1          AB 1.519481 1.1247154 2.052805 0.006417064
      2 adjusted_AB 1.083350 0.7268601 1.614683 0.694183560

---

    Code
      testout_OR$inferential$fit
    Output
      $model_before
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = dat)
      
      Coefficients:
      (Intercept)         ARMA  
           0.8473       0.4184  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    1023 
      Residual Deviance: 1016 	AIC: 1020
      
      $model_after
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = dat, weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
          0.84730      0.08006  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    726.7 
      Residual Deviance: 726.5 	AIC: 712.5
      
      $res_AB
      $res_AB$est
      [1] 1.08335
      
      $res_AB$se
      [1] 0.2275624
      
      $res_AB$ci_l
      [1] 0.7268601
      
      $res_AB$ci_u
      [1] 1.614683
      
      $res_AB$pval
      [1] 0.6941836
      
      
      $res_AB_unadj
      $res_AB_unadj$est
      [1] 1.519481
      
      $res_AB_unadj$se
      [1] 0.2373883
      
      $res_AB_unadj$ci_l
      [1] 1.124715
      
      $res_AB_unadj$ci_u
      [1] 2.052805
      
      $res_AB_unadj$pval
      [1] 0.006417064
      
      
      $boot_res
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      print(testout_boot_RR$descriptive$summary, digits = 5)
    Output
        trt_ind treatment            type      n events events_pct
      1       B         B Before matching 400.00  280.0     70.000
      2       A         A Before matching 500.00  390.0     78.000
      3       B         B  After matching 400.00  280.0     70.000
      4       A         A  After matching 199.43  142.9     71.654

---

    Code
      print(testout_boot_RR$inferential$summary, digits = 5)
    Output
               case     RR     LCL    UCL      pval
      1          AB 1.1143 1.02937 1.2062 0.0074553
      2 adjusted_AB 1.0236 0.91236 1.1485 0.6908096

---

    Code
      print(testout_boot_RR$inferential$fit, digits = 5)
    Output
      $model_before
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = dat)
      
      Coefficients:
      (Intercept)         ARMA  
         -0.35667      0.10821  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    1023 
      Residual Deviance: 1015.6 	AIC: 1019.6
      
      $model_after
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = dat, weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
        -0.356675     0.023352  
      
      Degrees of Freedom: 899 Total (i.e. Null);  898 Residual
      Null Deviance:	    726.66 
      Residual Deviance: 726.48 	AIC: 712.47
      
      $res_AB
      $res_AB$est
      [1] 1.0236
      
      $res_AB$se
      [1] 0.060252
      
      $res_AB$ci_l
      [1] 0.91236
      
      $res_AB$ci_u
      [1] 1.1485
      
      $res_AB$pval
      [1] 0.69081
      
      
      $res_AB_unadj
      $res_AB_unadj$est
      [1] 1.1143
      
      $res_AB_unadj$se
      [1] 0.045119
      
      $res_AB_unadj$ci_l
      [1] 1.0294
      
      $res_AB_unadj$ci_u
      [1] 1.2062
      
      $res_AB_unadj$pval
      [1] 0.0074553
      
      
      $boot_res
      
      STRATIFIED BOOTSTRAP
      
      
      Call:
      boot(data = boot_ipd, statistic = stat_fun, R = R, strata = weights_object$boot_strata, 
          w_obj = weights_object, pseudo_ipd = pseudo_ipd, normalize = normalize_weights)
      
      
      Bootstrap Statistics :
           original      bias    std. error
      t1* 0.0233518  0.01366708  0.05380528
      t2* 0.0030551 -0.00004315  0.00038235
      
      $boot_res_AB
      $boot_res_AB$est
      [1] 1.0236
      
      $boot_res_AB$se
      [1] NA
      
      $boot_res_AB$ci_l
      [1] 0.90867
      
      $boot_res_AB$ci_u
      [1] 1.122
      
      $boot_res_AB$pval
      [1] NA
      
      

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
      print(testout2$descriptive$summary, digits = 5)
    Output
        trt_ind treatment            type records  n.max n.start  events   rmean
      1       B         B Before matching     300 300.00  300.00 178.000  4.3036
      2       A         A Before matching     500 500.00  500.00 190.000  8.7097
      3       B         B  After matching     300 300.00  300.00 178.000  4.3036
      4       A         A  After matching     500 199.43  199.43  65.689 10.1660
        se(rmean)  median 0.95LCL 0.95UCL
      1   0.33673  2.7461  2.2611  3.3209
      2   0.35515  7.5876  6.2787 10.2885
      3   0.33673  2.7461  2.2611  3.3209
      4   0.54999 11.9000  7.8153 14.8738

---

    Code
      print(testout2$inferential$summary, digits = 5)
    Output
               case      HR     LCL     UCL       pval
      1          AB 0.37490 0.30390 0.46248 5.2452e-20
      2 adjusted_AB 0.28348 0.20747 0.38734 2.4734e-15

---

    Code
      print(testout2$inferential$fit, digits = 5)
    Output
      $km_before
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = dat, conf.type = km_conf_type)
      
              n events  median 0.95LCL 0.95UCL
      ARM=B 300    178  83.585  68.823  101.08
      ARM=A 500    190 230.948 191.108  313.16
      
      $km_after
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = dat, weights = dat$weights, 
          conf.type = km_conf_type)
      
            records      n  events  median 0.95LCL 0.95UCL
      ARM=B     300 300.00 178.000  83.585  68.823  101.08
      ARM=A     500 199.43  65.689 362.207 237.877  452.72
      
      $model_before
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = dat)
      
               coef exp(coef) se(coef)       z         p
      ARMA -0.98110   0.37490  0.10712 -9.1589 < 2.2e-16
      
      Likelihood ratio test=80.62  on 1 df, p=< 2.22e-16
      n= 800, number of events= 368 
      
      $model_after
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = dat, weights = weights, 
          robust = TRUE)
      
               coef exp(coef) se(coef) robust se      z         p
      ARMA -1.26062   0.28348  0.15035   0.15927 -7.915 2.473e-15
      
      Likelihood ratio test=80.4  on 1 df, p=< 2.22e-16
      n= 800, number of events= 368 
      
      $res_AB
      $res_AB$est
      [1] 0.28348
      
      $res_AB$se
      [1] 0.046018
      
      $res_AB$ci_l
      [1] 0.20747
      
      $res_AB$ci_u
      [1] 0.38734
      
      $res_AB$pval
      [1] 2.4734e-15
      
      
      $res_AB_unadj
      $res_AB_unadj$est
      [1] 0.3749
      
      $res_AB_unadj$se
      [1] 0.040506
      
      $res_AB_unadj$ci_l
      [1] 0.3039
      
      $res_AB_unadj$ci_u
      [1] 0.46248
      
      $res_AB_unadj$pval
      [1] 5.2452e-20
      
      
      $boot_res
      
      STRATIFIED BOOTSTRAP
      
      
      Call:
      boot(data = boot_ipd, statistic = stat_fun, R = R, strata = weights_object$boot_strata, 
          w_obj = weights_object, pseudo_ipd = pseudo_ipd, normalize = normalize_weights)
      
      
      Bootstrap Statistics :
           original     bias    std. error
      t1* -1.260621 0.00245135   0.1313882
      t2*  0.025367 0.00058194   0.0026704
      
      $boot_res_AB
      $boot_res_AB$est
      [1] 0.28348
      
      $boot_res_AB$se
      [1] NA
      
      $boot_res_AB$ci_l
      [1] 0.21858
      
      $boot_res_AB$ci_u
      [1] 0.36584
      
      $boot_res_AB$pval
      [1] NA
      
      

