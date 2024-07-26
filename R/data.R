# unanchored datasets ------

#' Patient data from single arm study
#' @format a data frame with 500 rows and 8 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm.}
#'     \item{AGE}{Age in years at baseline.}
#'     \item{SEX}{Sex of patient recorded as character `"Male"`/`"Female"`.}
#'     \item{SMOKE}{Smoking status at baseline as integer `1`/`0`.}
#'     \item{ECOG0}{Indicator of ECOG score = 0 at baseline as integer `1`/`0`.}
#'     \item{N_PR_THER}{Number of prior therapies received as integer `1, 2, 3, 4`.}
#'     \item{SEX_MALE}{Indicator of `SEX == "Male"` as numeric `1`/`0`.}
#'
#'   }
#' @keywords dataset
#' @family unanchored datasets
"adsl_sat"

#' Survival data from single arm trial
#' @format A data frame with 500 rows and 10 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm, `"A"`.}
#'     \item{AVAL}{Analysis value which in this dataset overall survival time in days.}
#'     \item{AVALU}{Unit of `AVAL`.}
#'     \item{PARAMCD}{Paramater code of `AVAL`, `"OS"`.}
#'     \item{PARAM}{Parameter name of `AVAL`, `"Overall Survival`.}
#'     \item{CNSR}{Censoring indicator `0`/`1`.}
#'     \item{TIME}{Survival time in days.}
#'     \item{EVENT}{Event indicator `0`/`1`.}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"adtte_sat"


#' Pseudo individual patient survival data from published study
#' @format A data frame with 300 rows and 3 columns:
#'   \describe{
#'     \item{TIME}{Survival time in days.}
#'     \item{EVENT}{Event indicator `0`/`1`.}
#'     \item{ARM}{Assigned treatment arm, `"B"`.}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"pseudo_ipd_sat"


#' Centered patient data from single arm trial
#' @format A data frame with 500 rows and 14 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm.}
#'     \item{AGE}{Age in years at baseline.}
#'     \item{SEX}{Sex of patient recorded as character `"Male"`/`"Female"`.}
#'     \item{SMOKE}{Smoking status at baseline as integer `1`/`0`.}
#'     \item{ECOG0}{Indicator of ECOG score = 0 at baseline as integer `1`/`0`.}
#'     \item{N_PR_THER}{Number of prior therapies received as integer `1, 2, 3, 4`.}
#'     \item{SEX_MALE}{Indicator of `SEX == "Male"` as numeric `1`/`0`.}
#'     \item{AGE_CENTERED}{Age in years at baseline relative to average in aggregate data [agd].}
#'     \item{AGE_MEDIAN_CENTERED}{`AGE` greater/less than `MEDIAN_AGE` in [agd] coded as `1`/`0` and then centered at
#'      0.5.}
#'     \item{AGE_SQUARED_CENTERED}{`AGE` squared and centered with respect to the `AGE` in [agd]. The squared age in the
#'       aggregate data is derived from the \eqn{E(X^2)} term in the variance formula.}
#'     \item{SEX_MALE_CENTERED}{`SEX_MALE` centered by the proportion of male patients in [agd]}
#'     \item{ECOG0_CENTERED}{`ECOG0` centered by the proportion of `ECOG0` in [agd]}
#'     \item{SMOKE_CENTERED}{`SMOKE` centered by the proportion of `SMOKE` in [agd]}
#'     \item{N_PR_THER_MEDIAN_CENTERED}{`N_PR_THER` centered by the median in [agd].}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"centered_ipd_sat"

#' Binary outcome data from single arm trial
#' @format A data frame with 500 rows and 5 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm.}
#'     \item{AVAL}{Analysis value, in this dataset an indicator of response.}
#'     \item{PARAM}{Parameter type of `AVAL`.}
#'     \item{RESPONSE}{Indicator of response.}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"adrs_sat"

#' Weighted object for single arm trial data
#' @format A `maicplus_estimate_weights` object created by [estimate_weights()] containing
#'   \describe{
#'     \item{data}{patient level data with weights}
#'     \item{centered_colnames}{Columns used in MAIC}
#'     \item{nr_missing}{Number of observations with missing data}
#'     \item{ess}{Expected sample size}
#'     \item{opt}{Information from `optim` from weight calculation}
#'     \item{boot}{Parameters and bootstrap sample weights, `NULL` in this object}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"weighted_sat"

# aggregate data ------

#' Aggregate effect modifier data from published study
#'
#' This data is formatted to be used in [center_ipd()].
#'
#' @format A data frame with 3 rows and 9 columns:
#'   \describe{
#'     \item{STUDY}{The study name, Study_XXXX}
#'     \item{ARM}{Study arm name or total}
#'     \item{N}{Number of observations in study arm}
#'     \item{AGE_MEAN}{Mean age in study arm}
#'     \item{AGE_MEDIAN}{Median age in study arm}
#'     \item{AGE_SD}{Standard deviation of age in study arm}
#'     \item{SEX_MALE_COUNT}{Number of male patients}
#'     \item{ECOG0_COUNT}{Number of patients with ECOG score = 0}
#'     \item{SMOKE_COUNT}{Number of smokers}
#'     \item{N_PR_THER_MEDIAN}{Median number of prior therapies}
#'   }
#' @family unanchored datasets
#' @family anchored datasets
#' @keywords dataset
"agd"


# anchored datasets -------

#' Patient data from two arm trial
#' @format A data frame with 1000 rows and 8 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm.}
#'     \item{AGE}{Age in years at baseline.}
#'     \item{SEX}{Sex of patient recorded as character "Male"/"Female"}
#'     \item{SMOKE}{Smoking status at baseline as integer `1`/`0`.}
#'     \item{ECOG0}{Indicator of ECOG score = 0 at baseline as integer `1`/`0`.}
#'     \item{N_PR_THER}{Number of prior therapies received as integer `1, 2, 3, 4`.}
#'     \item{SEX_MALE}{Indicator of SEX == "Male" as numeric 1/0}
#'   }
#' @family anchored datasets
#' @keywords dataset
"adsl_twt"


#' Survival data from two arm trial
#' @format A data frame with 1000 rows and 10 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm, `"A"`, `"C"`.}
#'     \item{AVAL}{Analysis value which in this dataset overall survival time in days.}
#'     \item{AVALU}{Unit of `AVAL`.}
#'     \item{PARAMCD}{Parameter code of `AVAL`, `"OS"`.}
#'     \item{PARAM}{Parameter name of `AVAL`, `"Overall Survival`.}
#'     \item{CNSR}{Censoring indicator `0`/`1`.}
#'     \item{TIME}{Survival time in days.}
#'     \item{EVENT}{Event indicator `0`/`1`.}
#'   }
#' @family anchored datasets
#' @keywords dataset
"adtte_twt"

#' Binary outcome data from two arm trial
#' @format A data frame with 1000 rows and 5 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm, `"A"`, `"C"`.}
#'     \item{AVAL}{Analysis value, in this dataset an indicator of response.}
#'     \item{PARAM}{Parameter type of `AVAL`.}
#'     \item{RESPONSE}{Indicator of response.}
#'   }
#' @family anchored datasets
#' @keywords dataset
"adrs_twt"

#' Pseudo individual patient survival data from published two arm study
#' @format A data frame with 800 rows and 3 columns:
#'   \describe{
#'     \item{TIME}{Survival time in days.}
#'     \item{EVENT}{Event indicator `0`/`1`.}
#'     \item{ARM}{Assigned treatment arm, `"B"`, `"C"`.}
#'   }
#' @family anchored datasets
#' @keywords dataset
"pseudo_ipd_twt"


#' Centered patient data from two arm trial
#' @format A data frame with 1000 rows and 14 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm.}
#'     \item{AGE}{Age in years at baseline.}
#'     \item{SEX}{Sex of patient recorded as character `"Male"`/`"Female"`.}
#'     \item{SMOKE}{Smoking status at baseline as integer `1`/`0`.}
#'     \item{ECOG0}{Indicator of ECOG score = 0 at baseline as integer `1`/`0`.}
#'     \item{N_PR_THER}{Number of prior therapies received as integer `1, 2, 3, 4`.}
#'     \item{SEX_MALE}{Indicator of `SEX == "Male"` as numeric `1`/`0`.}
#'     \item{AGE_CENTERED}{Age in years at baseline relative to average in aggregate data [agd].}
#'     \item{AGE_MEDIAN_CENTERED}{`AGE` greater/less than `MEDIAN_AGE` in [agd] coded as `1`/`0` and then centered at
#'      0.5.}
#'     \item{AGE_SQUARED_CENTERED}{`AGE` squared and centered with respect to the `AGE` in [agd]. The squared age in the
#'       aggregate data is derived from the \eqn{E(X^2)} term in the variance formula.}
#'     \item{SEX_MALE_CENTERED}{`SEX_MALE` centered by the proportion of male patients in [agd]}
#'     \item{ECOG0_CENTERED}{`ECOG0` centered by the proportion of `ECOG0` in [agd]}
#'     \item{SMOKE_CENTERED}{`SMOKE` centered by the proportion of `SMOKE` in [agd]}
#'     \item{N_PR_THER_MEDIAN_CENTERED}{`N_PR_THER` centered by the median in [agd].}
#'   }
#' @keywords dataset
#' @family anchored datasets
"centered_ipd_twt"

#' Weighted object for two arm trial data
#' @format A `maicplus_estimate_weights` object created by [estimate_weights()] containing
#'   \describe{
#'     \item{data}{patient level data with weights}
#'     \item{centered_colnames}{Columns used in MAIC}
#'     \item{nr_missing}{Number of observations with missing data}
#'     \item{ess}{Expected sample size}
#'     \item{opt}{Information from `optim` from weight calculation}
#'     \item{boot}{Parameters and bootstrap sample weights for the 100 samples}
#' }
#' @description
#' The weighted patient data for a two arm trial generated from the centered patient data ([centered_ipd_twt]).
#' It has weights calculated for 100 bootstrap samples.
#'
#' The object is generated using the following code:
#' ```
#' estimate_weights(
#'   data = centered_ipd_twt,
#'   centered_colnames = c(
#'     "AGE_CENTERED",
#'     "AGE_MEDIAN_CENTERED",
#'     "AGE_SQUARED_CENTERED",
#'     "SEX_MALE_CENTERED",
#'     "ECOG0_CENTERED",
#'     "SMOKE_CENTERED"
#'     ),
#'   n_boot_iteration = 100
#'  )
#' ```
#'
#' @keywords dataset
#' @family anchored datasets
"weighted_twt"


if (FALSE) {
  make_roxygen_data <- function(df) {
    cn <- colnames(df)
    cat("#' @format A data frame with", nrow(df), "rows and", ncol(df), "columns:\n")
    cat("#'   \\describe{\n")
    for (i in cn) cat("#'     \\item{", i, "}{}\n", sep = "")
    cat("#'   }")
  }
}
