# unanchored datasets ------

#' Patient data from single arm study
#' @format a data frame with 500 rows and 8 columns:
#'   \describe{
#'     \item{USUBJID}{Unique subject identifiers for patients.}
#'     \item{ARM}{Assigned treatment arm.}
#'     \item{AGE}{Age in years at baseline.}
#'     \item{SEX}{Sex of patient recorded as character `"Male"`/`"Female"`.}
#'     \item{SMOKE}{Smoking status at baseline as integer `0`/`1`.}
#'     \item{ECOG0}{ECOG score at baseline as integer `0`/`1`.}
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
#'     \item{EVNT}{Event indicator `0`/`1`.}
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
#'     \item{Time}{Survival time in days.}
#'     \item{Event}{Event indicator `0`/`1`.}
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
#'     \item{SMOKE}{Smoking status at baseline as integer `0`/`1`.}
#'     \item{ECOG0}{ECOG score at baseline as integer `0`/`1`.}
#'     \item{N_PR_THER}{Number of prior therapies received as integer `1, 2, 3, 4`.}
#'     \item{SEX_MALE}{Indicator of `SEX == "Male"` as numeric `1`/`0`.}
#'     \item{AGE_CENTERED}{Age in years at baseline relative to average in aggregate data [agd].}
#'     \item{AGE_MEDIAN_CENTERED}{Age greater/less than median age in [agd] coded as  TODO .}
#'     \item{AGE_SQUARED_CENTERED}{TODO}
#'     \item{SEX_MALE_CENTERED}{TODO}
#'     \item{ECOG0_CENTERED}{TODO}
#'     \item{SMOKE_CENTERED}{TODO}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"centered_ipd_sat"

#' Binary outcome data from single arm trial
#' @format A data frame with 500 rows and 5 columns:
#'   \describe{
#'     \item{USUBJID}{}
#'     \item{ARM}{}
#'     \item{AVAL}{}
#'     \item{PARAM}{}
#'     \item{RESPONSE}{}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"adrs_sat"



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
#'     \item{AGE_SD}{TODO}
#'     \item{SEX_MALE_COUNT}{Number of male patients}
#'     \item{ECOG0_COUNT}{Number of patients with ECOG score TODO}
#'     \item{SMOKE_COUNT}{Number of smokers}
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
#'     \item{SMOKE}{Smoking status at baseline as integer 0/1}
#'     \item{ECOG0}{ECOG score at baseline as integer 0/1}
#'     \item{N_PR_THER}{Number of prior therapies received as integer 1, 2, 3, 4}
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
#'     \item{EVNT}{Event indicator `0`/`1`.}
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


#' Pseudo individual patient survival data from published two arm study
#' @format A data frame with 800 rows and 3 columns:
#'   \describe{
#'     \item{Time}{Survival time in days.}
#'     \item{Event}{Event indicator `0`/`1`.}
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
#'     \item{SMOKE}{Smoking status at baseline as integer `0`/`1`.}
#'     \item{ECOG0}{ECOG score at baseline as integer `0`/`1`.}
#'     \item{N_PR_THER}{Number of prior therapies received as integer `1, 2, 3, 4`.}
#'     \item{SEX_MALE}{Indicator of `SEX == "Male"` as numeric `1`/`0`.}
#'     \item{AGE_CENTERED}{Age in years at baseline relative to average in aggregate data [agd].}
#'     \item{AGE_MEDIAN_CENTERED}{Age greater/less than median age in [agd] coded as  TODO .}
#'     \item{AGE_SQUARED_CENTERED}{TODO}
#'     \item{SEX_MALE_CENTERED}{TODO}
#'     \item{ECOG0_CENTERED}{TODO}
#'     \item{SMOKE_CENTERED}{TODO}
#'   }
#' @keywords dataset
#' @family anchored datasets
"centered_ipd_twt"


if (FALSE) {
  make_roxygen_data <- function(df) {
    cn <- colnames(df)
    cat("#' @format A data frame with", nrow(df), "rows and", ncol(df), "columns:\n")
    cat("#'   \\describe{\n")
    for (i in cn) cat("#'     \\item{", i, "}{}\n", sep = "")
    cat("#'   }")
  }
}
