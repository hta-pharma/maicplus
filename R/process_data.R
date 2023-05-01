# Functions for pre-processing data before conduct MAIC

# functions to be exported ---------------------------------------

#' Process and Check Aggregated data format
#'
#' @param raw.agd
#'
#' @return
#' @export
process.agd <- function(raw.agd) {

  raw.agd <- as.data.frame(raw.agd)
  # make all column names to be capital, avoid different style
  names(raw.agd) <- toupper(names(raw.agd))

  # regulaized column name patterns
  must_exist <- c("STUDY","ARM", "N")
  legal_suffix <- c("MEAN","MEDIAN","SD","COUNT")

  # swap "TREATMENT" column to "ARM", if applicable
  if("TREATMENT"%in%names(raw.agd) & (!"ARM"%in%names(raw.agd))){
    raw.agd$ARM <- raw.agd$TREATMENT
    raw.agd <- raw.agd[,names(raw.agd)!="TREATMENT"]
    warning("'TREATMENT' column is renamed as 'ARM'")
  }

  # check: must exist
  if(!all(must_exist%in%names(raw.agd))){
    stop("At least 1 of the must-exists columns (STUDY, ARM, N) cannot be found in raw.agd!")
  }

  # check: legal suffix
  other_colnames <- setdiff(names(raw.agd),must_exist)
  ind1 <- grepl("_",other_colnames,fixed=TRUE)
  ind2 <- sapply(other_colnames,function(xx){
      tmp <- unlist(strsplit(xx,split="_"))
      tmp[length(tmp)] # this deployment is robust to the cases that there are multiple _ in the column name
  })
  ind2 <- (ind2%in%legal_suffix)

  use_cols <- other_colnames[ind1&ind2]
  use_agd <- raw.agd[,c(must_exist,use_cols),drop=FALSE]
  if(!all(other_colnames%in%use_cols)){
    warning(paste0("following columns are ignored since it does not following the naming convention:",
                   paste(setdiff(other_colnames,use_cols),collapse = ",")))
  }

  # complete statistics for ARM=="Total"
  if(!"total"%in%tolower(use_agd$ARM)){
     use_agd <- complete_agd(use_agd)
  }

  # calculate percentage columns
  ind <- grepl("_COUNT$",names(use_agd))
  if(any(ind)){
    for(i in which(ind)){
      use_agd[[i]] <- use_agd[[i]]/use_agd$N
    }
    names(use_agd) <- gsub("_COUNT$","_PROP",names(use_agd))
  }

  # output
  with(use_agd,{ use_agd[tolower(ARM)=="total",,drop=FALSE] })
}

#' Process Individual Patient data to dummize categorical effective modifiers
#'
#' @param raw.ipd
#' @param dummize.cols
#' @param dummize.ref.level
#'
#' @return
#' @export

process.ipd <- function(raw.ipd, dummize.cols, dummize.ref.level){
   for(i in 1:length(dummize.cols)){
     yy <- raw.ipd[[ dummize.cols[i] ]]
     yy_levels <- na.omit(unique(yy))
     yy <- factor( as.character(yy), levels=c(dummize.ref.level[i], setdiff(yy_levels,dummize.ref.level[i])) )
     new_yy <- sapply(levels(yy)[-1],function(j){
                   as.numeric(yy == j)
               })
     new_yy <- as.data.frame(new_yy)
     names(new_yy) <- toupper(paste(dummize.cols[i],levels(yy)[-1], sep="_"))
     raw.ipd <- cbind(raw.ipd,new_yy)
   }
   raw.ipd
}


#' Center effect modifiers
#'
#' @param ipd
#' @param agd
#'
#' @return
#' @export
center.ipd <- function(ipd,agd){
  # regulaized column name patterns
  must_exist <- c("STUDY","ARM", "N")
  legal_suffix <- c("MEAN","MEDIAN","SD","PROP")
  suffix_pat <- paste(paste0("_",legal_suffix,"$"),collapse = "|")

  for(i in 1:nrow(agd)){# study i
     study_id <- agd$STUDY[i]
     use_agd <- agd[i,!names(agd)%in%must_exist,drop=FALSE]
     param_id <- gsub(suffix_pat,"",names(use_agd))

     for(j in 1:ncol(use_agd)){# effect modifier j
       if(is.na(use_agd[[j]])) next

       ipd_param <- param_id[j]

       if(grepl("_MEAN$|_PROP$",names(use_agd)[j])){

            ipd[[paste0(ipd_param,"_","CENTERED")]] <- ipd[[ipd_param]] - use_agd[[j]]

       }else if(grepl("_MEDIAN$",names(use_agd)[j])){

            ipd[[paste0(ipd_param,"_MED_","CENTERED")]] <- ipd[[ipd_param]] > use_agd[[j]]
            ipd[[paste0(ipd_param,"_MED_","CENTERED")]] <- ipd[[paste0(ipd_param,"_MED_","CENTERED")]] - 0.5

       }else if(grepl("_SD$",names(use_agd)[j])){

            ipd[[paste0(ipd_param,"_SD_","CENTERED")]] <- ipd[[ipd_param]]^2
            tmp_aim <- use_agd[[j]]^2 + (use_agd[[paste0(ipd_param,"_MEAN")]]^2)
            ipd[[paste0(ipd_param,"_SD_","CENTERED")]] <- ipd[[paste0(ipd_param,"_SD_","CENTERED")]] - tmp_aim

       }
     } # end of j
  } # end of i

  # output
  ipd
}










#' helper function: transform TTE ADaM data to suitable input for survival R pkg
#'
#' @param dd data frame, ADTTE read via haven::read_sas
#' @param time_scale a character string, 'year', 'month', 'week' or 'day', time unit of median survival time
#'
#' @return a data frame can be used as input to survival::Surv
ext_tte_transfer <- function(dd, time_scale = "month", trt = NULL) {
  timeUnit <- list("year" = 365.24, "month" = 30.4367, "week" = 7, "day" = 1)

  if (!time_scale %in% names(timeUnit)) stop("time_scale has to be 'year', 'month', 'week' or 'day'")

  if ("CENSOR" %in% names(dd)) {
    dd <- dd[!is.na(dd$CENSOR), ]
    dd$status <- as.numeric(dd$CENSOR == 0)
  }
  if ("EVENT" %in% names(dd)) {
    dd$status <- as.numeric(as.character(dd$EVENT))
  }
  if ("TIME" %in% names(dd)) {
    dd$AVAL <- as.numeric(as.character(dd$TIME))
  }

  dd$time <- dd$AVAL * timeUnit[[time_scale]]
  if (!is.null(trt)) dd$treatment <- trt
  as.data.frame(dd)
}



# functions NOT to be exported ---------------------------------------

#' Calculate pooled arm statistics in AgD based on arm-specific statistics
#'
#' @param use_agd
#'
#' @return
complete.agd <- function(use.agd) {
  use.agd <- as.data.frame(use.agd)
  use.agd <- with(use.agd, {use.agd[tolower(ARM)!="total",,drop=FALSE]})

  if(nrow(use.agd)<2) stop("error in call complete_agd: need to have at least 2 rows that ARM!='total' ")

  rowId <- nrow(use.agd)+1
  use.agd[rowId, ] <- NA
  use.agd$STUDY[rowId] <-   use.agd$STUDY[1]
  use.agd$ARM[rowId] <- "total"

  # complete N and count
  NN <- use.agd$N[rowId] <- sum(use.agd$N, na.rm=TRUE)
  nn <- use.agd$N[-rowId]
  for(i in grep("_COUNT$",names(use.agd),value=TRUE)){
     use.agd[[i]][rowId] <- sum(use.agd[[i]], na.rm=TRUE)
  }

  # complete MEAN
  for(i in grep("_MEAN$",names(use.agd),value=TRUE)){
    use.agd[[i]][rowId] <- sum(use.agd[[i]]*nn)/NN
  }

  # complete SD
  for(i in grep("_SD$",names(use.agd),value=TRUE)){
    use.agd[[i]][rowId] <- sqrt( sum(use.agd[[i]]^2*(nn-1))/(N-1) )
  }

  # complete MEDIAN, approximately!!
  for(i in grep("_MEDIAN$",names(use.agd),value=TRUE)){
    use.agd[[i]][rowId] <- mean(use.agd[[i]][-rowId])
  }

  # output
  use.agd
}


