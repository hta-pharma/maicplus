# maic_anchored works for TTE

    Code
      testout$descriptive$summary
    Output
        trt_ind treatment                 type records    n.max  n.start    events
      1       C         C IPD, before matching     500 500.0000 500.0000 500.00000
      2       A         A IPD, before matching     500 500.0000 500.0000 190.00000
      3       C         C  IPD, after matching     500 173.3137 173.3137 173.31374
      4       A         A  IPD, after matching     500 173.3137 173.3137  55.37398
      5       C         C        AgD, external     500 500.0000 500.0000 500.00000
      6       B         B        AgD, external     300 300.0000 300.0000 178.00000
          events%     rmean  se(rmean)    median   0.95LCL   0.95UCL
      1 100.00000  2.564797 0.11366994  1.836467  1.644765  2.045808
      2  38.00000  8.709690 0.35514766  7.587627  6.278691 10.288538
      3 100.00000  2.363509 0.18354814  1.641770  1.052627  2.093468
      4  31.95014 10.584609 0.57397937 12.166430 10.244293        NA
      5 100.00000  2.455272 0.09848888  1.851987  1.670540  2.009650
      6  59.33333  4.303551 0.33672602  2.746131  2.261125  3.320857

---

    Code
      print(testout$inferential$summary, digits = 5)
    Output
               case      HR     LCL     UCL       pval
      1          AC 0.22166 0.18672 0.26314 2.1366e-66
      2 adjusted_AC 0.13675 0.09450 0.19789 4.9838e-26
      3          BC 0.57180 0.48120 0.67946 2.1437e-10
      4          AB 0.38765 0.30393 0.49443 2.2704e-14
      5 adjusted_AB 0.23916 0.15906 0.35959 6.1920e-12

---

    Code
      testout$inferential$fit
    Output
      $km_before_ipd
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = ipd, conf.type = km_conf_type)
      
              n events median 0.95LCL 0.95UCL
      ARM=C 500    500   55.9    50.1    62.3
      ARM=A 500    190  230.9   191.1   313.2
      
      $km_after_ipd
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = ipd, weights = ipd$weights, 
          conf.type = km_conf_type)
      
            records   n events median 0.95LCL 0.95UCL
      ARM=C     500 173  173.3     50      32    63.7
      ARM=A     500 173   55.4    370     312      NA
      
      $km_agd
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = pseudo_ipd, 
          conf.type = km_conf_type)
      
              n events median 0.95LCL 0.95UCL
      ARM=C 500    500   56.4    50.8    61.2
      ARM=B 300    178   83.6    68.8   101.1
      
      $model_before_ipd
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = ipd)
      
               coef exp(coef) se(coef)      z      p
      ARMA -1.50662   0.22166  0.08753 -17.21 <2e-16
      
      Likelihood ratio test=341.2  on 1 df, p=< 2.2e-16
      n= 1000, number of events= 690 
      
      $model_after_ipd
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = ipd, weights = weights, 
          robust = TRUE)
      
              coef exp(coef) se(coef) robust se      z      p
      ARMA -1.9896    0.1368   0.1691    0.1886 -10.55 <2e-16
      
      Likelihood ratio test=175.6  on 1 df, p=< 2.2e-16
      n= 1000, number of events= 690 
      
      $model_agd
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = pseudo_ipd)
      
               coef exp(coef) se(coef)      z        p
      ARMB -0.55897   0.57180  0.08802 -6.351 2.14e-10
      
      Likelihood ratio test=43.66  on 1 df, p=3.91e-11
      n= 800, number of events= 678 
      
      $res_AC
      $res_AC$est
      [1] 0.1367514
      
      $res_AC$se
      [1] 0.1885548
      
      $res_AC$ci_l
      [1] 0.09450039
      
      $res_AC$ci_u
      [1] 0.1978927
      
      $res_AC$pval
      [1] 4.98377e-26
      
      
      $res_AC_unadj
      $res_AC_unadj$est
      [1] 0.2216588
      
      $res_AC_unadj$se
      [1] 0.08752989
      
      $res_AC_unadj$ci_l
      [1] 0.1867151
      
      $res_AC_unadj$ci_u
      [1] 0.2631423
      
      $res_AC_unadj$pval
      [1] 2.13665e-66
      
      
      $res_BC
      $res_BC$est
      [1] 0.5718004
      
      $res_BC$se
      [1] 0.0880166
      
      $res_BC$ci_l
      [1] 0.4811989
      
      $res_BC$ci_u
      [1] 0.6794607
      
      $res_BC$pval
      [1] 2.14366e-10
      
      
      $res_AB
                  result             pvalue 
      "0.24[0.16; 0.36]"           "<0.001" 
      
      $res_AB_unadj
                  result             pvalue 
      "0.39[0.30; 0.49]"           "<0.001" 
      
      $boot_res
      NULL
      
      $boot_res_AC
      NULL
      
      $boot_res_AB_mc
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      print(testout2$descriptive$summary, digits = 5)
    Output
        trt_ind treatment                 type records  n.max n.start  events events%
      1       C         C IPD, before matching     500 500.00  500.00 500.000 100.000
      2       A         A IPD, before matching     500 500.00  500.00 190.000  38.000
      3       C         C  IPD, after matching     500 173.31  173.31 173.314 100.000
      4       A         A  IPD, after matching     500 173.31  173.31  55.374  31.950
      5       C         C        AgD, external     500 500.00  500.00 500.000 100.000
      6       B         B        AgD, external     300 300.00  300.00 178.000  59.333
          rmean se(rmean)  median 0.95LCL 0.95UCL
      1  2.5648  0.113670  1.8365  1.6448  2.0458
      2  8.7097  0.355148  7.5876  6.2787 10.2885
      3  2.3635  0.183548  1.6418  1.0526  2.0935
      4 10.5846  0.573979 12.1664 10.2443      NA
      5  2.4553  0.098489  1.8520  1.6705  2.0096
      6  4.3036  0.336726  2.7461  2.2611  3.3209

---

    Code
      print(testout2$inferential$summary, digits = 5)
    Output
               case      HR     LCL     UCL       pval
      1          AC 0.22166 0.18672 0.26314 2.1366e-66
      2 adjusted_AC 0.13675 0.09450 0.19789 4.9838e-26
      3          BC 0.57180 0.48120 0.67946 2.1437e-10
      4          AB 0.38765 0.30393 0.49443 2.2704e-14
      5 adjusted_AB 0.23916 0.15906 0.35959 6.1920e-12

---

    Code
      print(testout2$inferential$fit, digits = 5)
    Output
      $km_before_ipd
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = ipd, conf.type = km_conf_type)
      
              n events  median 0.95LCL 0.95UCL
      ARM=C 500    500  55.897  50.063  62.269
      ARM=A 500    190 230.948 191.108 313.157
      
      $km_after_ipd
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = ipd, weights = ipd$weights, 
          conf.type = km_conf_type)
      
            records      n  events  median 0.95LCL 0.95UCL
      ARM=C     500 173.31 173.314  49.971  32.039   63.72
      ARM=A     500 173.31  55.374 370.316 311.811      NA
      
      $km_agd
      Call: survfit(formula = Surv(TIME, EVENT) ~ ARM, data = pseudo_ipd, 
          conf.type = km_conf_type)
      
              n events median 0.95LCL 0.95UCL
      ARM=C 500    500 56.370  50.847  61.169
      ARM=B 300    178 83.585  68.823 101.079
      
      $model_before_ipd
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = ipd)
      
               coef exp(coef) se(coef)       z         p
      ARMA -1.50662   0.22166  0.08753 -17.213 < 2.2e-16
      
      Likelihood ratio test=341.18  on 1 df, p=< 2.22e-16
      n= 1000, number of events= 690 
      
      $model_after_ipd
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = ipd, weights = weights, 
          robust = TRUE)
      
               coef exp(coef) se(coef) robust se       z         p
      ARMA -1.98959   0.13675  0.16905   0.18855 -10.552 < 2.2e-16
      
      Likelihood ratio test=175.56  on 1 df, p=< 2.22e-16
      n= 1000, number of events= 690 
      
      $model_agd
      Call:
      coxph(formula = Surv(TIME, EVENT) ~ ARM, data = pseudo_ipd)
      
                coef exp(coef)  se(coef)       z         p
      ARMB -0.558965  0.571800  0.088017 -6.3507 2.144e-10
      
      Likelihood ratio test=43.66  on 1 df, p=3.9096e-11
      n= 800, number of events= 678 
      
      $res_AC
      $res_AC$est
      [1] 0.13675
      
      $res_AC$se
      [1] 0.18855
      
      $res_AC$ci_l
      [1] 0.0945
      
      $res_AC$ci_u
      [1] 0.19789
      
      $res_AC$pval
      [1] 4.9838e-26
      
      
      $res_AC_unadj
      $res_AC_unadj$est
      [1] 0.22166
      
      $res_AC_unadj$se
      [1] 0.08753
      
      $res_AC_unadj$ci_l
      [1] 0.18672
      
      $res_AC_unadj$ci_u
      [1] 0.26314
      
      $res_AC_unadj$pval
      [1] 2.1366e-66
      
      
      $res_BC
      $res_BC$est
      [1] 0.5718
      
      $res_BC$se
      [1] 0.088017
      
      $res_BC$ci_l
      [1] 0.4812
      
      $res_BC$ci_u
      [1] 0.67946
      
      $res_BC$pval
      [1] 2.1437e-10
      
      
      $res_AB
                  result             pvalue 
      "0.24[0.16; 0.36]"           "<0.001" 
      
      $res_AB_unadj
                  result             pvalue 
      "0.39[0.30; 0.49]"           "<0.001" 
      
      $boot_res
      
      STRATIFIED BOOTSTRAP
      
      
      Call:
      boot(data = boot_ipd, statistic = stat_fun, R = R, strata = weights_object$boot_strata, 
          w_obj = weights_object, normalize = normalize_weights)
      
      
      Bootstrap Statistics :
           original    bias    std. error
      t1* -1.422586 0.0150047    0.287870
      t2*  0.043300 0.0044132    0.010425
      t3*  0.208086 0.0092974    0.023911
      t4* -1.989591 0.0098093    0.250584
      t5*  0.188555 0.0099839    0.026183
      t6*  0.035553 0.0044132    0.010425
      
      $boot_res_AC
      $boot_res_AC$est
      [1] 0.13675
      
      $boot_res_AC$se
      [1] NA
      
      $boot_res_AC$ci_l
      [1] 0.082866
      
      $boot_res_AC$ci_u
      [1] 0.22129
      
      $boot_res_AC$pval
      [1] NA
      
      
      $boot_res_AB_mc
      $boot_res_AB_mc$est
      [1] 0.24109
      
      $boot_res_AB_mc$se
      [1] NA
      
      $boot_res_AB_mc$ci_l
      [1] 0.13509
      
      $boot_res_AB_mc$ci_u
      [1] 0.41754
      
      $boot_res_AB_mc$pval
      [1] NA
      
      
      $boot_res_AB
      $boot_res_AB$est
      [1] 0.23916
      
      $boot_res_AB$se
      [1] NA
      
      $boot_res_AB$ci_l
      [1] 0.14211
      
      $boot_res_AB$ci_u
      [1] 0.40249
      
      $boot_res_AB$pval
      [1] NA
      
      

# maic_anchored for binary case gives the expected result

    Code
      testout_OR$descriptive$summary
    Output
        trt_ind treatment                 type        n   events events_pct
      1       C         C IPD, before matching 500.0000 338.0000   67.60000
      2       A         A IPD, before matching 500.0000 390.0000   78.00000
      3       C         C  IPD, after matching 199.4265 134.0364   67.21094
      4       A         A  IPD, after matching 199.4265 142.8968   71.65386
      5       C         C        AgD, external 320.0000 120.0000   37.50000
      6       B         B        AgD, external 480.0000 280.0000   58.33333

---

    Code
      testout_OR$inferential$summary
    Output
               case        OR       LCL       UCL         pval
      1          AC 1.6993007 1.2809976 2.2541985 2.354448e-04
      2 adjusted_AC 1.2332036 0.7710134 1.9724576 3.817109e-01
      3          BC 2.3333333 1.7458092 3.1185794 1.035032e-08
      4          AB 0.7282717 0.4857575 1.0918611 1.248769e-01
      5 adjusted_AB 0.5285158 0.3043103 0.9179084 2.356848e-02

---

    Code
      testout_OR$inferential$fit
    Output
      $model_before_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = ipd)
      
      Coefficients:
      (Intercept)         ARMA  
           0.7354       0.5302  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    1170 
      Residual Deviance: 1157 	AIC: 1161
      
      $model_after_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = ipd, weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
           0.7177       0.2096  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    491.1 
      Residual Deviance: 490.1 	AIC: 458.4
      
      $model_agd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = pseudo_ipd)
      
      Coefficients:
      (Intercept)         ARMB  
          -0.5108       0.8473  
      
      Degrees of Freedom: 799 Total (i.e. Null);  798 Residual
      Null Deviance:	    1109 
      Residual Deviance: 1075 	AIC: 1079
      
      $res_AC
      $res_AC$est
      [1] 1.233204
      
      $res_AC$se
      [1] 0.3085377
      
      $res_AC$ci_l
      [1] 0.7710134
      
      $res_AC$ci_u
      [1] 1.972458
      
      $res_AC$pval
      [1] 0.3817109
      
      
      $res_AC_unadj
      $res_AC_unadj$est
      [1] 1.699301
      
      $res_AC_unadj$se
      [1] 0.2488482
      
      $res_AC_unadj$ci_l
      [1] 1.280998
      
      $res_AC_unadj$ci_u
      [1] 2.254199
      
      $res_AC_unadj$pval
      [1] 0.0002354448
      
      
      $res_BC
      $res_BC$est
      [1] 2.333333
      
      $res_BC$se
      [1] 0.3510631
      
      $res_BC$ci_l
      [1] 1.745809
      
      $res_BC$ci_u
      [1] 3.118579
      
      $res_BC$pval
      [1] 1.035032e-08
      
      
      $res_AB
                  result             pvalue 
      "0.53[0.30; 0.92]"            "0.024" 
      
      $res_AB_unadj
                  result             pvalue 
      "0.73[0.49; 1.09]"            "0.125" 
      
      $boot_res
      NULL
      
      $boot_res_AC
      NULL
      
      $boot_res_AB_mc
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      testout_RR$descriptive$summary
    Output
        trt_ind treatment                 type        n   events events_pct
      1       C         C IPD, before matching 500.0000 338.0000   67.60000
      2       A         A IPD, before matching 500.0000 390.0000   78.00000
      3       C         C  IPD, after matching 199.4265 134.0364   67.21094
      4       A         A  IPD, after matching 199.4265 142.8968   71.65386
      5       C         C        AgD, external 320.0000 120.0000   37.50000
      6       B         B        AgD, external 480.0000 280.0000   58.33333

---

    Code
      testout_RR$inferential$summary
    Output
               case        RR       LCL       UCL         pval
      1          AC 1.1538462 1.0688892 1.2455556 2.451956e-04
      2 adjusted_AC 1.0661042 0.9234409 1.2308077 3.824937e-01
      3          BC 1.5555555 1.3250564 1.8261510 6.678781e-08
      4          AB 0.7417582 0.6210074 0.8859883 9.832938e-04
      5 adjusted_AB 0.6853527 0.5525932 0.8500075 5.832632e-04

---

    Code
      testout_RR$inferential$fit
    Output
      $model_before_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = ipd)
      
      Coefficients:
      (Intercept)         ARMA  
          -0.3916       0.1431  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    1170 
      Residual Deviance: 1157 	AIC: 1161
      
      $model_after_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = ipd, weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
         -0.39733      0.06401  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    491.1 
      Residual Deviance: 490.1 	AIC: 458.4
      
      $model_agd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = pseudo_ipd)
      
      Coefficients:
      (Intercept)         ARMB  
          -0.9808       0.4418  
      
      Degrees of Freedom: 799 Total (i.e. Null);  798 Residual
      Null Deviance:	    1109 
      Residual Deviance: 1075 	AIC: 1079
      
      $res_AC
      $res_AC$est
      [1] 1.066104
      
      $res_AC$se
      [1] 0.07845782
      
      $res_AC$ci_l
      [1] 0.9234409
      
      $res_AC$ci_u
      [1] 1.230808
      
      $res_AC$pval
      [1] 0.3824937
      
      
      $res_AC_unadj
      $res_AC_unadj$est
      [1] 1.153846
      
      $res_AC_unadj$se
      [1] 0.04507631
      
      $res_AC_unadj$ci_l
      [1] 1.068889
      
      $res_AC_unadj$ci_u
      [1] 1.245556
      
      $res_AC_unadj$pval
      [1] 0.0002451956
      
      
      $res_BC
      $res_BC$est
      [1] 1.555556
      
      $res_BC$se
      [1] 0.1279272
      
      $res_BC$ci_l
      [1] 1.325056
      
      $res_BC$ci_u
      [1] 1.826151
      
      $res_BC$pval
      [1] 6.678781e-08
      
      
      $res_AB
                  result             pvalue 
      "0.69[0.55; 0.85]"            "0.001" 
      
      $res_AB_unadj
                  result             pvalue 
      "0.74[0.62; 0.89]"            "0.001" 
      
      $boot_res
      NULL
      
      $boot_res_AC
      NULL
      
      $boot_res_AB_mc
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      testout_RD$descriptive$summary
    Output
        trt_ind treatment                 type        n   events events_pct
      1       C         C IPD, before matching 500.0000 338.0000   67.60000
      2       A         A IPD, before matching 500.0000 390.0000   78.00000
      3       C         C  IPD, after matching 199.4265 134.0364   67.21094
      4       A         A  IPD, after matching 199.4265 142.8968   71.65386
      5       C         C        AgD, external 320.0000 120.0000   37.50000
      6       B         B        AgD, external 480.0000 280.0000   58.33333

---

    Code
      testout_RD$inferential$summary
    Output
               case         RD        LCL       UCL         pval
      1          AC  10.400000   4.921741 15.878259 1.985755e-04
      2 adjusted_AC   4.442927  -5.499175 14.385028 3.811014e-01
      3          BC  20.833333  13.934963 27.731704 3.235832e-09
      4          AB -10.433333 -19.242354 -1.624313 2.026711e-02
      5 adjusted_AB -16.390407 -28.491353 -4.289460 7.937461e-03

---

    Code
      testout_RD$inferential$fit
    Output
      $model_before_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = ipd)
      
      Coefficients:
      (Intercept)         ARMA  
            0.676        0.104  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    1170 
      Residual Deviance: 1157 	AIC: 1161
      
      $model_after_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = ipd, weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
          0.67211      0.04443  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    491.1 
      Residual Deviance: 490.1 	AIC: 458.4
      
      $model_agd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = pseudo_ipd)
      
      Coefficients:
      (Intercept)         ARMB  
           0.3750       0.2083  
      
      Degrees of Freedom: 799 Total (i.e. Null);  798 Residual
      Null Deviance:	    1109 
      Residual Deviance: 1075 	AIC: 1079
      
      $res_AC
      $res_AC$est
      [1] 4.442927
      
      $res_AC$se
      [1] 5.072594
      
      $res_AC$ci_l
      [1] -5.499175
      
      $res_AC$ci_u
      [1] 14.38503
      
      $res_AC$pval
      [1] 0.3811014
      
      
      $res_AC_unadj
      $res_AC_unadj$est
      [1] 10.4
      
      $res_AC_unadj$se
      [1] 2.795081
      
      $res_AC_unadj$ci_l
      [1] 4.921741
      
      $res_AC_unadj$ci_u
      [1] 15.87826
      
      $res_AC_unadj$pval
      [1] 0.0001985755
      
      
      $res_BC
      $res_BC$est
      [1] 20.83333
      
      $res_BC$se
      [1] 3.519641
      
      $res_BC$ci_l
      [1] 13.93496
      
      $res_BC$ci_u
      [1] 27.7317
      
      $res_BC$pval
      [1] 3.235832e-09
      
      
      $res_AB
                       result                  pvalue 
      "-16.39[-28.49; -4.29]"                 "0.008" 
      
      $res_AB_unadj
                       result                  pvalue 
      "-10.43[-19.24; -1.62]"                 "0.020" 
      
      $boot_res
      NULL
      
      $boot_res_AC
      NULL
      
      $boot_res_AB_mc
      NULL
      
      $boot_res_AB
      NULL
      

---

    Code
      print(testout_boot_OR$descriptive$summary, digits = 5)
    Output
        trt_ind treatment                 type      n events events_pct
      1       C         C IPD, before matching 500.00 338.00     67.600
      2       A         A IPD, before matching 500.00 390.00     78.000
      3       C         C  IPD, after matching 199.43 134.04     67.211
      4       A         A  IPD, after matching 199.43 142.90     71.654
      5       C         C        AgD, external 320.00 120.00     37.500
      6       B         B        AgD, external 480.00 280.00     58.333

---

    Code
      print(testout_boot_OR$inferential$summary, digits = 5)
    Output
               case      OR     LCL     UCL       pval
      1          AC 1.69930 1.28100 2.25420 2.3544e-04
      2 adjusted_AC 1.23320 0.77101 1.97246 3.8171e-01
      3          BC 2.33333 1.74581 3.11858 1.0350e-08
      4          AB 0.72827 0.48576 1.09186 1.2488e-01
      5 adjusted_AB 0.52852 0.30431 0.91791 2.3568e-02

---

    Code
      print(testout_boot_OR$inferential$fit, digits = 5)
    Output
      $model_before_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = ipd)
      
      Coefficients:
      (Intercept)         ARMA  
          0.73545      0.53022  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    1170.5 
      Residual Deviance: 1156.8 	AIC: 1160.8
      
      $model_after_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = ipd, weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
          0.71774      0.20962  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    491.07 
      Residual Deviance: 490.14 	AIC: 458.35
      
      $model_agd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = binomial(link = glm_link), 
          data = pseudo_ipd)
      
      Coefficients:
      (Intercept)         ARMB  
         -0.51083      0.84730  
      
      Degrees of Freedom: 799 Total (i.e. Null);  798 Residual
      Null Deviance:	    1109 
      Residual Deviance: 1075.4 	AIC: 1079.4
      
      $res_AC
      $res_AC$est
      [1] 1.2332
      
      $res_AC$se
      [1] 0.30854
      
      $res_AC$ci_l
      [1] 0.77101
      
      $res_AC$ci_u
      [1] 1.9725
      
      $res_AC$pval
      [1] 0.38171
      
      
      $res_AC_unadj
      $res_AC_unadj$est
      [1] 1.6993
      
      $res_AC_unadj$se
      [1] 0.24885
      
      $res_AC_unadj$ci_l
      [1] 1.281
      
      $res_AC_unadj$ci_u
      [1] 2.2542
      
      $res_AC_unadj$pval
      [1] 0.00023544
      
      
      $res_BC
      $res_BC$est
      [1] 2.3333
      
      $res_BC$se
      [1] 0.35106
      
      $res_BC$ci_l
      [1] 1.7458
      
      $res_BC$ci_u
      [1] 3.1186
      
      $res_BC$pval
      [1] 1.035e-08
      
      
      $res_AB
                  result             pvalue 
      "0.53[0.30; 0.92]"            "0.024" 
      
      $res_AB_unadj
                  result             pvalue 
      "0.73[0.49; 1.09]"            "0.125" 
      
      $boot_res
      
      STRATIFIED BOOTSTRAP
      
      
      Call:
      boot(data = boot_ipd, statistic = stat_fun, R = R, strata = weights_object$boot_strata, 
          w_obj = weights_object, eff_measure = eff_measure, normalize = normalize_weights)
      
      
      Bootstrap Statistics :
           original   bias    std. error
      t1* -0.754756 0.661173   0.2591600
      t2*  0.069346 0.010683   0.0040404
      t3*  0.263336 0.019472   0.0071395
      t4*  0.209615 0.515383   0.2915831
      t5*  0.217810 0.023141   0.0083805
      t6*  0.047441 0.010683   0.0040404
      
      $boot_res_AC
      $boot_res_AC$est
      [1] 1.2332
      
      $boot_res_AC$se
      [1] NA
      
      $boot_res_AC$ci_l
      [1] 0.41592
      
      $boot_res_AC$ci_u
      [1] 1.3044
      
      $boot_res_AC$pval
      [1] NA
      
      
      $boot_res_AB_mc
      $boot_res_AB_mc$est
      [1] 0.47013
      
      $boot_res_AB_mc$se
      [1] NA
      
      $boot_res_AB_mc$ci_l
      [1] 0.14604
      
      $boot_res_AB_mc$ci_u
      [1] 0.40334
      
      $boot_res_AB_mc$pval
      [1] NA
      
      
      $boot_res_AB
      $boot_res_AB$est
      [1] 0.52852
      
      $boot_res_AB$se
      [1] NA
      
      $boot_res_AB$ci_l
      [1] 0.27843
      
      $boot_res_AB$ci_u
      [1] 1.0032
      
      $boot_res_AB$pval
      [1] NA
      
      

