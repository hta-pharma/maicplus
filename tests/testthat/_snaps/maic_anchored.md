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
      3 100.00000  2.679473 0.20670827  1.815795  1.457222  2.292484
      4  31.95014 10.584609 0.57397937 12.166430 10.244293        NA
      5 100.00000  2.455272 0.09848888  1.851987  1.670540  2.009650
      6  59.33333  4.303551 0.33672602  2.746131  2.261125  3.320857

---

    Code
      testout$inferential$summary
    Output
               case        HR       LCL       UCL         pval
      1          AC 0.2216588 0.1867151 0.2631423 2.136650e-66
      2 adjusted_AC 0.1631852 0.1113815 0.2390829 1.361531e-20
      3          BC 0.5718004 0.4811989 0.6794607 2.143660e-10
      4          AB 0.3876507 0.3039348 0.4944253 2.270430e-14
      5 adjusted_AB 0.2853885 0.1876867 0.4339497 4.509575e-09

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
      ARM=C     500 173  173.3   55.3    44.4    69.8
      ARM=A     500 173   55.4  370.3   311.8      NA
      
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
      ARMA -1.8129    0.1632   0.1602    0.1949 -9.303 <2e-16
      
      Likelihood ratio test=155.2  on 1 df, p=< 2.2e-16
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
      [1] 0.1631852
      
      $res_AC$se
      [1] 0.194863
      
      $res_AC$ci_l
      [1] 0.1113815
      
      $res_AC$ci_u
      [1] 0.2390829
      
      $res_AC$pval
      [1] 1.361531e-20
      
      
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
      "0.29[0.19; 0.43]"           "<0.001" 
      
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
      testout2$descriptive$summary
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
      3 100.00000  2.679473 0.20670827  1.815795  1.457222  2.292484
      4  31.95014 10.584609 0.57397937 12.166430 10.244293        NA
      5 100.00000  2.455272 0.09848888  1.851987  1.670540  2.009650
      6  59.33333  4.303551 0.33672602  2.746131  2.261125  3.320857

---

    Code
      testout2$inferential$summary
    Output
               case        HR       LCL       UCL         pval
      1          AC 0.2216588 0.1867151 0.2631423 2.136650e-66
      2 adjusted_AC 0.1631852 0.1113815 0.2390829 1.361531e-20
      3          BC 0.5718004 0.4811989 0.6794607 2.143660e-10
      4          AB 0.3876507 0.3039348 0.4944253 2.270430e-14
      5 adjusted_AB 0.2853885 0.1876867 0.4339497 4.509575e-09

---

    Code
      testout2$inferential$fit
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
      ARM=C     500 173  173.3   55.3    44.4    69.8
      ARM=A     500 173   55.4  370.3   311.8      NA
      
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
      ARMA -1.8129    0.1632   0.1602    0.1949 -9.303 <2e-16
      
      Likelihood ratio test=155.2  on 1 df, p=< 2.2e-16
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
      [1] 0.1631852
      
      $res_AC$se
      [1] 0.194863
      
      $res_AC$ci_l
      [1] 0.1113815
      
      $res_AC$ci_u
      [1] 0.2390829
      
      $res_AC$pval
      [1] 1.361531e-20
      
      
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
      "0.29[0.19; 0.43]"           "<0.001" 
      
      $res_AB_unadj
                  result             pvalue 
      "0.39[0.30; 0.49]"           "<0.001" 
      
      $boot_res
      
      STRATIFIED BOOTSTRAP
      
      
      Call:
      boot(data = boot_ipd, statistic = stat_fun, R = R, strata = weights_object$boot_strata, 
          w_obj = weights_object, normalize = normalize_weights)
      
      
      Bootstrap Statistics :
             original       bias    std. error
      t1* -1.24586465 -0.161716794  0.28787014
      t2*  0.04571853  0.001994459  0.01042480
      t3*  0.21381891  0.003564636  0.02391077
      t4* -1.81286926 -0.166912210  0.25058383
      t5*  0.19486304  0.003675676  0.02618306
      t6*  0.03797160  0.001994459  0.01042480
      
      $boot_res_AC
      $boot_res_AC$est
      [1] 0.1631852
      
      $boot_res_AC$se
      [1] NA
      
      $boot_res_AC$ci_l
      [1] 0.1179974
      
      $boot_res_AC$ci_u
      [1] 0.3151137
      
      $boot_res_AC$pval
      [1] NA
      
      
      $boot_res_AB_mc
      $boot_res_AB_mc$est
      [1] 0.287692
      
      $boot_res_AB_mc$se
      [1] NA
      
      $boot_res_AB_mc$ci_l
      [1] 0.1923646
      
      $boot_res_AB_mc$ci_u
      [1] 0.5945605
      
      $boot_res_AB_mc$pval
      [1] NA
      
      
      $boot_res_AB
      $boot_res_AB$est
      [1] 0.2853885
      
      $boot_res_AB$se
      [1] NA
      
      $boot_res_AB$ci_l
      [1] 0.1695758
      
      $boot_res_AB$ci_u
      [1] 0.4802958
      
      $boot_res_AB$pval
      [1] NA
      
      

# maic_anchored for binary case gives the expected result

    Code
      testout$descriptive$summary
    Output
        trt_ind treatment                 type   n   events events_pct
      1       C         C IPD, before matching 500 338.0000   67.60000
      2       A         A IPD, before matching 500 390.0000   78.00000
      3       C         C  IPD, after matching 500 131.2892   26.25784
      4       A         A  IPD, after matching 500 142.8968   28.57935
      5       C         C        AgD, external 320 120.0000   37.50000
      6       B         B        AgD, external 480 280.0000   58.33333

---

    Code
      testout$inferential$summary
    Output
               case        OR       LCL       UCL         pval
      1          AC 1.3119021 0.8210000 2.0963303 2.562849e-01
      2 adjusted_AC 1.6993007 1.2809976 2.2541985 2.354448e-04
      3          BC 2.3333333 1.7458092 3.1185794 1.035032e-08
      4          AB 0.5622438 0.3239933 0.9756933 4.061296e-02
      5 adjusted_AB 0.7282717 0.4857575 1.0918611 1.248769e-01

---

    Code
      testout$inferential$fit
    Output
      $model_before_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = ipd)
      
      Coefficients:
      (Intercept)         ARMA  
           0.7354       0.5302  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    1170 
      Residual Deviance: 1157 	AIC: 1161
      
      $model_after_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = ipd, 
          weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
           0.6559       0.2715  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    495.5 
      Residual Deviance: 493.9 	AIC: 454.5
      
      $model_agd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = pseudo_ipd)
      
      Coefficients:
      (Intercept)         ARMB  
          -0.5108       0.8473  
      
      Degrees of Freedom: 799 Total (i.e. Null);  798 Residual
      Null Deviance:	    1109 
      Residual Deviance: 1075 	AIC: 1079
      
      $res_AC
      $res_AC$est
      [1] 1.311902
      
      $res_AC$se
      [1] 0.3275028
      
      $res_AC$ci_l
      [1] 0.821
      
      $res_AC$ci_u
      [1] 2.09633
      
      $res_AC$pval
      [1] 0.2562849
      
      
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
      "0.56[0.32; 0.98]"            "0.041" 
      
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
      testout2$descriptive$summary
    Output
        trt_ind treatment                 type   n   events events_pct
      1       C         C IPD, before matching 500 338.0000   67.60000
      2       A         A IPD, before matching 500 390.0000   78.00000
      3       C         C  IPD, after matching 500 131.2892   26.25784
      4       A         A  IPD, after matching 500 142.8968   28.57935
      5       C         C        AgD, external 320 120.0000   37.50000
      6       B         B        AgD, external 480 280.0000   58.33333

---

    Code
      testout2$inferential$summary
    Output
               case        OR       LCL       UCL         pval
      1          AC 1.3119021 0.8210000 2.0963303 2.562849e-01
      2 adjusted_AC 1.6993007 1.2809976 2.2541985 2.354448e-04
      3          BC 2.3333333 1.7458092 3.1185794 1.035032e-08
      4          AB 0.5622438 0.3239933 0.9756933 4.061296e-02
      5 adjusted_AB 0.7282717 0.4857575 1.0918611 1.248769e-01

---

    Code
      testout2$inferential$fit
    Output
      $model_before_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = ipd)
      
      Coefficients:
      (Intercept)         ARMA  
           0.7354       0.5302  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    1170 
      Residual Deviance: 1157 	AIC: 1161
      
      $model_after_ipd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = ipd, 
          weights = weights)
      
      Coefficients:
      (Intercept)         ARMA  
           0.6559       0.2715  
      
      Degrees of Freedom: 999 Total (i.e. Null);  998 Residual
      Null Deviance:	    495.5 
      Residual Deviance: 493.9 	AIC: 454.5
      
      $model_agd
      
      Call:  glm(formula = RESPONSE ~ ARM, family = glm_link, data = pseudo_ipd)
      
      Coefficients:
      (Intercept)         ARMB  
          -0.5108       0.8473  
      
      Degrees of Freedom: 799 Total (i.e. Null);  798 Residual
      Null Deviance:	    1109 
      Residual Deviance: 1075 	AIC: 1079
      
      $res_AC
      $res_AC$est
      [1] 1.311902
      
      $res_AC$se
      [1] 0.3275028
      
      $res_AC$ci_l
      [1] 0.821
      
      $res_AC$ci_u
      [1] 2.09633
      
      $res_AC$pval
      [1] 0.2562849
      
      
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
      "0.56[0.32; 0.98]"            "0.041" 
      
      $res_AB_unadj
                  result             pvalue 
      "0.73[0.49; 1.09]"            "0.125" 
      
      $boot_res
      
      STRATIFIED BOOTSTRAP
      
      
      Call:
      boot(data = boot_ipd, statistic = stat_fun, R = R, strata = weights_object$boot_strata, 
          w_obj = weights_object, eff_measure = eff_measure, normalize = normalize_weights)
      
      
      Bootstrap Statistics :
             original     bias    std. error
      t1* -0.69289322 0.59931037 0.259160031
      t2*  0.06888563 0.01114340 0.004040355
      t3*  0.26246071 0.02034770 0.007139486
      t4*  0.27147808 0.45352064 0.291583054
      t5*  0.21675070 0.02420063 0.008380500
      t6*  0.04698086 0.01114340 0.004040355
      
      $boot_res_AC
      $boot_res_AC$est
      [1] 1.311902
      
      $boot_res_AC$se
      [1] NA
      
      $boot_res_AC$ci_l
      [1] 0.4706998
      
      $boot_res_AC$ci_u
      [1] 1.476168
      
      $boot_res_AC$pval
      [1] NA
      
      
      $boot_res_AB_mc
      $boot_res_AB_mc$est
      [1] 0.500127
      
      $boot_res_AB_mc$se
      [1] NA
      
      $boot_res_AB_mc$ci_l
      [1] 0.1652744
      
      $boot_res_AB_mc$ci_u
      [1] 0.4564577
      
      $boot_res_AB_mc$pval
      [1] NA
      
      
      $boot_res_AB
      $boot_res_AB$est
      [1] 0.5622438
      
      $boot_res_AB$se
      [1] NA
      
      $boot_res_AB$ci_l
      [1] 0.2962009
      
      $boot_res_AB$ci_u
      [1] 1.067242
      
      $boot_res_AB$pval
      [1] NA
      
      

