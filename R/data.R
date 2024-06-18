# unanchored datasets ------

#' Patient data from single arm study
#' @format a data frame with 500 rows and 8 columns:
#'   \describe{
#'     \item{USUBJID}{}
#'     \item{ARM}{}
#'     \item{AGE}{}
#'     \item{SEX}{}
#'     \item{SMOKE}{}
#'     \item{ECOG0}{}
#'     \item{N_PR_THER}{}
#'     \item{SEX_MALE}{}
#'   }
#' @keywords dataset
#' @family unanchored datasets
"adsl_sat"

#' Survival data from single arm trial
#' @format A data frame with 500 rows and 10 columns:
#'   \describe{
#'     \item{USUBJID}{}
#'     \item{ARM}{}
#'     \item{EVNT}{}
#'     \item{AVAL}{}
#'     \item{AVALU}{}
#'     \item{PARAMCD}{}
#'     \item{PARAM}{}
#'     \item{CNSR}{}
#'     \item{TIME}{}
#'     \item{EVENT}{}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"adtte_sat"


#' Pseudo individual patient survival data from published study
#' @format A data frame with 300 rows and 3 columns:
#'   \describe{
#'     \item{Time}{}
#'     \item{Event}{}
#'     \item{ARM}{}
#'   }
#' @family unanchored datasets
#' @keywords dataset
"pseudo_ipd_sat"


#' Centered patient data from single arm trial
#' @format A data frame with 500 rows and 14 columns:
#'   \describe{
#'     \item{USUBJID}{}
#'     \item{ARM}{}
#'     \item{AGE}{}
#'     \item{SEX}{}
#'     \item{SMOKE}{}
#'     \item{ECOG0}{}
#'     \item{N_PR_THER}{}
#'     \item{SEX_MALE}{}
#'     \item{AGE_CENTERED}{}
#'     \item{AGE_MEDIAN_CENTERED}{}
#'     \item{AGE_SQUARED_CENTERED}{}
#'     \item{SEX_MALE_CENTERED}{}
#'     \item{ECOG0_CENTERED}{}
#'     \item{SMOKE_CENTERED}{}
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
#'     \item{STUDY}{}
#'     \item{ARM}{}
#'     \item{N}{}
#'     \item{AGE_MEAN}{}
#'     \item{AGE_MEDIAN}{}
#'     \item{AGE_SD}{}
#'     \item{SEX_MALE_COUNT}{}
#'     \item{ECOG0_COUNT}{}
#'     \item{SMOKE_COUNT}{}
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
#'     \item{USUBJID}{}
#'     \item{ARM}{}
#'     \item{EVNT}{}
#'     \item{AVAL}{}
#'     \item{AVALU}{}
#'     \item{PARAMCD}{}
#'     \item{PARAM}{}
#'     \item{CNSR}{}
#'     \item{TIME}{}
#'     \item{EVENT}{}
#'   }
#' @family anchored datasets
#' @keywords dataset
"adtte_twt"


#' Pseudo individual patient survival data from published two arm study
#' @format A data frame with 800 rows and 3 columns:
#'   \describe{
#'     \item{Time}{}
#'     \item{Event}{}
#'     \item{ARM}{}
#'   }
#' @family anchored datasets
#' @keywords dataset
"pseudo_ipd_twt"


#' Centered patient data from two arm trial
#' @format A data frame with 1000 rows and 14 columns:
#'   \describe{
#'     \item{USUBJID}{}
#'     \item{ARM}{}
#'     \item{AGE}{}
#'     \item{SEX}{}
#'     \item{SMOKE}{}
#'     \item{ECOG0}{}
#'     \item{N_PR_THER}{}
#'     \item{SEX_MALE}{}
#'     \item{AGE_CENTERED}{}
#'     \item{AGE_MEDIAN_CENTERED}{}
#'     \item{AGE_SQUARED_CENTERED}{}
#'     \item{SEX_MALE_CENTERED}{}
#'     \item{ECOG0_CENTERED}{}
#'     \item{SMOKE_CENTERED}{}
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
