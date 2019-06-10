####################
##
##
##
###################


take_only_WTs <- function(all_data_dt) 
  #for given data dt only WT individuals, assumes the column "Genotype"
{
  res <- try(if (!any(class(all_data_dt) == "data.table") |
                 (!"Genotype" %in% colnames(all_data_dt))) {
    stop("Input should be a data.table with a column 'Genotype'", call. = FALSE)
  } else{
    all_data_dt[Genotype %in% c("+/+", "+/Y", "WT")]
  })
  return(res)
}

get_percentiles_range_variable <- function(var, percentile) 
  #returns symmetrical percentile tails from the distribution var
{
  a <- quantile(var, c(percentile, 1 - percentile), na.rm = TRUE)
  names(a) <- c("lower", "upper")
  return(a)
}

add_zygosity <- function(dt,
                         Genotype_col = "Genotype",
                         new_format = FALSE) 
  #Adds zygosity info
  #new_format: Genotype:Zygosity, WT
  #old_format: Genotype/Genotype, Genotype/+, Genotype/Y, +/+
{
  res <- try(if (!any(class(dt) == "data.table") |
                 (!Genotype_col %in% colnames(dt))) {
    stop("Input should be a data.table with a column 'Genotype'", call. = FALSE)
  } else{
    if (new_format) {
      dt[, zygosity := sub(".*:", "", get(Genotype_col))]
      dt[zygosity == "Hemi", zygosity := "Hom"]
    }
    else{
      dt[, zygosity := case_when(
        get(Genotype_col) %in% c("+/+", "+/Y") ~ "WT",
        grepl("/\\+", get(Genotype_col)) ~ "Het",
        grepl("\\+/", get(Genotype_col)) ~ "Het",
        is.character(get(Genotype_col)) ~ "Hom"
      )]
    }
  })
  return(dt)
}



find_RR_days <- function(minInd,
                         samples_on_days,
                         days_you_want_RR_for = NULL) 
  # given number of WT samples available for each day  and required N for RR  ('minInd') returns start and end dates for computing RR;
  # days_you_want_RR_for allow to compute RR also for days not included in samples_on_days
  # if a choice of dates not provided, takes all days within 'samples_on_days'
  # samples_on_days is a dt with columns 'Assay_Date' and 'avail_WTs' (number of WT samples available on the given date)
{ 
  if (is.null(days_you_want_RR_for)) {
    days_you_want_RR_for <-
      data.frame(Assay_Date = samples_on_days[, Assay_Date])
  }
  
  # print(days_you_want_RR_for)
  #  print(sum(samples_on_days$avail_WTs) )
  if (sum(samples_on_days$avail_WTs, na.rm = TRUE) < minInd) {
    dates_for_RR <- NULL
    print("Not enough individuals: total number of animals lower than selected minInd")
    return(dates_for_RR)
  } else {
    dates_for_RR <-
      sapply(days_you_want_RR_for$Assay_Date, function(x) {
        #samples_available before the day x
        for_cumsum_before <-
          samples_on_days[Assay_Date < x, ] 
        #samples_available after the day x
        for_cumsum_after <-
          samples_on_days[Assay_Date > x, ] 
        
        #start collecting samples: samples from the day x
        if (length(samples_on_days[Assay_Date == x, avail_WTs]) > 0) {
          on_day <- samples_on_days[Assay_Date == x, avail_WTs]
        } else{
          on_day <- 0
        }
        
        #how many samples (smallest possible number) to get from before and from after the day x
        min_before <- ceiling((minInd - on_day) / 2)
        min_after <- minInd - on_day - min_before
        
        #identify the start of the RR e
        if (nrow(for_cumsum_before) == 0) {
          start_before <- "Not_enough_mice_before"
        } else{
          for_cumsum_before[, csb := rev(cumsum(rev(avail_WTs)))]
          #Here- what to do if there is not enough individuals before
          if (for_cumsum_before[1, csb] < min_before) {
            start_before <- "Not_enough_mice_before"
          } else{
            start_before <- max(for_cumsum_before[csb >= min_before, Assay_Date])
          }
        }
        
        #identify the end of the RR 
        if (nrow(for_cumsum_after) == 0) {
          start_after <- "Not_enough_mice_after"
        } else{
          for_cumsum_after[, csa := cumsum(avail_WTs)]
          #print(for_cumsum_after)
          if (for_cumsum_after[nrow(for_cumsum_after), csa] < min_after) {
            start_after <- "Not_enough_mice_after"
          } else{
            start_after <- min(for_cumsum_after[csa >= min_after, Assay_Date])
          }
        }
        
        return(
          c(
            Assay_Date = as.character(x),
            start_before = as.character(start_before),
            start_after = as.character(start_after)
          )
        )
      }) %>% t()
  #  print(dates_for_RR)
    return(dates_for_RR)
    
  }
  
  # Identify days where additional animals from another end of RR are needed 
  dates_to_change_before <-
    which(dates_for_RR[, "start_before"] == "Not_enough_mice_before")
  first_good_RR <-
    which(dates_for_RR[, "start_before"] != "Not_enough_mice_before")[1]
  
  dates_to_change_after <-
    which(as.character(dates_for_RR[, "start_after"]) == "Not_enough_mice_after")
  last_good_RR <-
    which(dates_for_RR[, "start_after"] == "Not_enough_mice_after")[1] - 1
  

  if (first_good_RR > last_good_RR) {
    dates_for_RR[, "start_after"] <- dates_for_RR[last_good_RR, 3]
    dates_for_RR[, "start_before"] <- dates_for_RR[first_good_RR, 2]
  } else{
    dates_for_RR[dates_to_change_before, 2] <- dates_for_RR[first_good_RR, 2]
    dates_for_RR[dates_to_change_before, 3] <- dates_for_RR[first_good_RR, 3]
    dates_for_RR[dates_to_change_after, 2] <- dates_for_RR[last_good_RR, 2]
    dates_for_RR[dates_to_change_after, 3] <- dates_for_RR[last_good_RR, 3]
  }
  return(dates_for_RR)
}


compute_RR_from_date <- function(WT_data,
                                 days_for_RR,
                                 percentile_cutoff,
                                 par,
                                 minInd) 
  #helper function for compute_RR_befaft
{
  RR <- apply(days_for_RR, 1, function(x) {
    data_to_make_RR <-
      WT_data[Assay_Date >= x["start_before"] &
                Assay_Date <= x["start_after"]]
    
    if (data_to_make_RR[, sum(!is.na(get(par)))] < minInd) {
  #    print(paste0(par, ": not enough individuals to prepare RR"))
      RR <-
        data.frame(
          Assay_Date = x["Assay_Date"],
          upper = NA,
          lower = NA,
          ind = NA,
          startRR = NA,
          endRR = NA
        )
    } else{
      a <-
        get_percentiles_range_variable(data_to_make_RR[, par, with = F], perc = percentile_cutoff)
      RR <-
        data.frame(
          Assay_Date = NA,
          upper = NA,
          lower = NA,
          ind = NA,
          startRR = NA,
          endRR = NA
        )
      RR["Assay_Date"] <- x["Assay_Date"]
      RR["upper"] <- a["upper"]
      RR["lower"] <- a["lower"]
      RR["ind"] <- sum(!is.na(data_to_make_RR[, par, with = F]))
      RR["startRR"] <- x["start_before"]
      RR["endRR"] <- x["start_after"]
    }
    return(RR)
  }) %>% bind_rows()
  return(RR)
}

tableh <- function(x) 
  #summary: number of samples above/below RR 
  {
  return(list(
    above = sum(x == 1, na.rm = T),
    below = sum(x == -1, na.rm = T),
    samples = sum(!is.na(x))
  ))
}


compute_RR_simple <-
  function(WT_data,
           days_without_WT = NULL,
           split_by_sex = TRUE,
           par,
           percentile = 0.025) 
    #computes RR for the whole time range, for one parameter
    {
    
    WT_dt <-
      data.table(WT_data)[, c("Assay_Date", "Sex", par), with = F][order(Assay_Date)][!is.na((par))]
     
    if (split_by_sex) {
      WT_male <- WT_dt[Sex == "Male", ]
      WT_female <- WT_dt[Sex == "Female", ]
      miss_males <-
        setdiff(WT_female$Assay_Date, WT_male$Assay_Date)
      if (length(miss_males) > 0) {
        class(miss_males) <- "Date"
        days_without_WT <- c(days_without_WT, miss_males)
      }
      miss_females <-
        setdiff(WT_male$Assay_Date, WT_female$Assay_Date)
      if (length(miss_females) > 0) {
        class(miss_females) <- "Date"
        days_without_WT <- c(days_without_WT, miss_females)
      }
      
      RR_male  <-
        compute_RR_simple(
          WT_male,
          perc = percentile,
          par = par,
          split_by_sex = FALSE,
          days_without_WT = days_without_WT
        )
      if (!is.null(RR_male)) {
        RR_male <- data.frame(RR_male, Sex = "Male")
      }
      
      RR_female <-
        compute_RR_simple(
          WT_female,
          perc = percentile,
          par = par,
          split_by_sex = FALSE,
          days_without_WT = days_without_WT
        )
      
      if (!is.null(RR_female)) {
        RR_female <- data.frame(RR_female, Sex = "Female")
      }

      if (!(is.null(RR_female) | is.null(RR_male))) {
        RRres <- rbind(RR_male, RR_female)
        RRres <- RRres[order(RRres$Assay_Date), ]
      } else{
        RRres <- NA
      }
      return(RRres)
    }
    
    
    perc_values <-
      get_percentiles_range_variable(WT_dt[, par, with = F], perc = percentile)
   
    rr_dates <- sort(unique(c(WT_dt$Assay_Date, days_without_WT)))
    class(rr_dates) <- "Date"
    RR <- data.frame(
      Assay_Date = rr_dates,
      upper = perc_values["upper"],
      lower  = perc_values["lower"],
      ind = sum(!is.na(WT_dt[, par, with = F])),
      startRR = min(rr_dates),
      endRR = max(rr_dates)
    )
    return(RR)
  }



prepare_RR_simple_dt <- function(all_data_dt,
                                 cols_for_RR,
                                 split_by_sex = TRUE) 
  #Prepare naive RR for every day, sex if needed, and parameter (wrapper for compute_RR_simple)
  { 
  all_params <-  lapply(cols_for_RR, function(x){
    a <- data.table(
      compute_RR_simple(
        take_only_WTs(all_data_dt),
        split_by_sex = split_by_sex,
        days_without_WT = as.Date(setdiff(
          as.character(all_data_dt$Assay_Date),
          as.character(take_only_WTs(all_data_dt)$Assay_Date)
        )),
        par = x
      ),
      par = x
    )
    return(a)
  }) %>%bind_rows()
  return(all_params)
}




compute_RR_befaft <- function(WT_data,
                              days_without_WT=NULL,
                              split_by_sex=TRUE,
                              minInd,
                              par,
                              percentile=0.025,
                              start=NULL,
                              end=NULL)
  #Compute time-window RR for the parameter par
  {
  WT_dt <- data.table(WT_data)[,c("Assay_Date", "Sex", par), with=F][order(Assay_Date)][!is.na((par))]
  
  if(!is.null(start)){
    WT_dt <- WT_dt[Assay_Date >=start] 
  }
  
  if(!is.null(end)){
    WT_dt <- WT_dt[Assay_Date <= end] 
  } 
  #add checks: if Assay_Date is class Date, Sex has aproprpaite name 
  
  if(split_by_sex){
    WT_male <- WT_dt[Sex=="Male",]
    WT_female <- WT_dt[Sex=="Female",]
    #identify dats missing one sex
    #days with female dates only
    miss_males <- setdiff(WT_female$Assay_Date , WT_male$Assay_Date )
    if(length(miss_males)>0){
      class(miss_males) <-"Date"
      days_without_WT <- c(days_without_WT, miss_males )
    }
    miss_females <- setdiff(WT_male$Assay_Date , WT_female$Assay_Date )
    if(length(miss_females)>0){
      class(miss_females) <-"Date"
      #days with male dates only
      days_without_WT <- c(days_without_WT,miss_females)
    }
    
    RR_male  <- compute_RR_befaft(WT_male,
                                  minInd=minInd,
                                  start=start,
                                  end=end,
                                  perc=percentile,
                                  par=par,
                                  split_by_sex=FALSE,
                                  days_without_WT=days_without_WT)
    if(!is.null(RR_male)){RR_male <- data.frame( RR_male, Sex="Male")}
    #  print(WT_female[,get(par)])
    
    RR_female <- compute_RR_befaft(WT_female,
                                   minInd=minInd,
                                   start=start,
                                   end=end,
                                   perc=percentile,
                                   par=par,
                                   split_by_sex=FALSE,
                                   days_without_WT=days_without_WT)
    if(!is.null(RR_female)){RR_female<-data.frame(RR_female, Sex="Female")}
    # print(RR_female)
    # print(RR_male)
    if(!(is.null(RR_female)|is.null(RR_male))){
      RRres <-rbind(RR_male, RR_female)
      RRres <-RRres[order(RRres$Assay_Date),]
    }else{
      RRres <- NA
    }  
    return(RRres)
  }
  
  WT_dt_sum <- WT_dt[,.(avail_WTs=sum(!is.na(get(par)))), by=Assay_Date]
  # print(WT_dt_sum)
  if(!is.null(days_without_WT)){print("adding days")
    days_without_WT <- setdiff(days_without_WT, WT_dt_sum$Assay_Date)
    class(days_without_WT) <-"Date"
    WT_dt_sum <- rbind(WT_dt_sum, data.table(Assay_Date=days_without_WT, avail_WTs=0) )[order(Assay_Date)]
  }
  #print(WT_dt_sum)
  days_for_RR <- find_RR_days(WT_dt_sum, minInd=minInd)
  #  print(days_for_RR )
  #  
  #change into dates again
  if(!is.null(days_for_RR )){
    days_for_RR <- data.table(data.frame(sapply(data.frame(days_for_RR),function(x) as.Date(x, format="%Y-%m-%d" ), simplify=F)))
    #  print(days_for_RR )
    RR <- compute_RR_from_date(WT_data=WT_data, days_for_RR=days_for_RR, percentile=percentile, par=par, minInd=minInd )
  }else{RR=NULL}
  
  return(RR)
}



prepare_RR <- function(cols_for_RR,
                       WT_dt=NULL,
                       all_data_dt,
                       minInd = 70) 
  #Prepare time-window RR for every day, sex if needed, and parameter (wrapper for compute_RR_befaft)
  {
  if (is.null(WT_dt)) {
    WT_dt <- take_only_WTs(all_data_dt)
  }
  RR_full <- lapply(cols_for_RR, function(par) {
    print(par)
    days_without_WT = setdiff(all_data_dt[, Assay_Date], WT_dt[, Assay_Date])
    class(days_without_WT) <- "Date"
    compute_RR_befaft(
      WT_dt,
      split_by_sex = TRUE,
      minInd = minInd,
      par = par,
      percentile = 0.025,
      days_without_WT = days_without_WT,
      start = NULL,
      end = NULL
    )
  })
  names(RR_full) <- cols_for_RR
  
  RR_full <- sapply(names(RR_full), function(x) {
    if (!all(is.na(RR_full[[x]]))) {
      if (nrow(RR_full[[x]]) > 0) {
        a <- RR_full[[x]]
        a$par <- x
      }
    } else{
      a = NULL
    }
    return(a)
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  RR_all_dt <- RR_full %>%
    bind_rows()%>%
    data.table()%>%
    unique()%>%
    .[order(Assay_Date)]
  
  return(RR_all_dt)
}


check_if_mouse_outside_dt <-
  function(all_data_dt, RR_all_dt, pars = NULL) 
    #Compare mouse values with RR and return info whether it is below/above/or within RR
    {
    if (is.null(pars)) {
      pars <- unique(grep("%", colnames(all_data_dt), val = T))
    }
 #   print(pars)
    per_mouse <- apply(all_data_dt[, c("Assay_Date", "Sex", pars), with = F], 1,
                       function(mouse) {
      datem <- mouse["Assay_Date"]
      sex <- mouse["Sex"]
      a <- RR_all_dt[Assay_Date == datem & Sex == sex]
      setkey(a, "par")
      upper_limits <- unlist(a[pars, upper])
      lower_limits <- unlist(a[pars, lower])
      above <- as.numeric(mouse[pars]) > upper_limits
      below <-  as.numeric(mouse[pars]) < lower_limits
      res <- rep(0, length(pars))
      res[(above)] <- 1
      res[(below)] <- -1
      res[is.na(above)] <- NA
      res[is.na(below)] <- NA
      if (length(pars) != length(res)) {
        print("HELP")
      }
      return(res)
    })%>%
      t()
#    print(dim(per_mouse))
 #   print(length(pars))
    colnames(per_mouse) <- pars
    rownames(per_mouse) <- all_data_dt$Mouse
    return(per_mouse)
  }

combine_outside <-
  function(manual_fract,
           outside_auto,
           outside_post_manual,
           outside_pre_manual,
           isMatrix = TRUE,
           colnames_for_RR) 
    #combines info from auto/premanual/postmanual 'outside info' in this order. Only for mice present in full manual (manual_fract, all mice done). Returns a dt with info about being outside of RR and another one with info about source of data used.
   {
    mr <-
      manual_fract[, .(Mouse,
                       Colony_Prefix,
                       Genotype,
                       Gene_Name,
                       Sex,
                       Assay_Date)]
    if (isMatrix) {
      result_auto <-
        merge(mr, data.table(outside_auto, Mouse = rownames(outside_auto)), all.x = T)
      result_pre <-
        merge(mr,
              data.table(outside_pre_manual, Mouse = rownames(outside_pre_manual)),
              all.x = T)
      result_post <-
        merge(mr,
              data.table(outside_post_manual, Mouse = rownames(outside_post_manual)),
              all.x = T)
    } else{
      result_auto <- merge(mr, outside_auto, by = "Mouse", all.x = T)
      result_pre <-
        merge(mr, outside_pre_manual,  by = "Mouse", all.x = T)
      result_post <-
        merge(mr, outside_post_manual, by = "Mouse", all.x = T)
    }
    if (!all(colnames_for_RR %in% colnames(result_pre))) {
      newcols <- setdiff(colnames_for_RR, colnames(result_pre))
      add <-
        matrix(
          nrow = nrow(result_pre),
          ncol = length(newcols),
          data = NA
        )
      colnames(add) <- newcols
      result_pre <- cbind(result_pre, add)
    }
   # print(head(result_pre))
    
    ra <-
      as.matrix(result_auto[order(Mouse)][, colnames_for_RR, with = F])
    # print(dim(ra))
    rpre <-
      as.matrix(result_pre[order(Mouse)][, colnames_for_RR, with = F])
    # print(dim(rpre))
    rpost <-
      as.matrix(result_post[order(Mouse)][, colnames_for_RR, with = F])
    # print(dim(rpost))
    
    a = ifelse(is.na(rpre), ifelse(is.na(ra), rpost, ra), rpre)
    result <- data.table(mr[order(Mouse)], a)
    
    a = ifelse(is.na(rpre), ifelse(is.na(ra), NA, "auto"), "pre")
    b = ifelse(is.na(a), ifelse(is.na(rpost), NA, "post"), a)
    ref <- data.table(mr[order(Mouse)], b)
    return(list(res = result, ref = ref))
  }



fix_lack_of_entries <- function(df,
                                column_name,
                                suffixes){
  #Fixing lack of Colony Prefix/Gene name etc.
  #for a dt made form mergers, fixing the columns which had value only in one of the merged dt
  df <- data.table(df)
  df[ , (column_name):= NA_character_]
  #  print(df)
  for(suff in suffixes){
    print(suff)
    #    print(df[is.na(get(column_name))])
    df <- df[is.na(get(column_name)), (column_name):=ifelse(!is.na(get(paste0(column_name, suff))), get(paste0(column_name, suff)), NA_character_), by=Mouse]
  }
  return(df)
}


#simple RR, splits (by sex/timepoints); removes outlier days beforehand
prepare_RR_whole_dtWrapper <- function(WT_data,
                                       split_date=NULL,
                                       split_by_sex=TRUE,
                                       min_ind=100,
                                       params=NULL,
                                       assay_dates=NULL,
                                       percentile=0.025,
                                       outliers_dt=NULL){
  whole_dt_l <- lapply(params, function(single_par){
    #removing outliers if suggested
    if(!is.null(outliers_dt)){
      if(single_par%in%outliers_dt$par){
        #   print(outliers_dt[par==single_par])
        WT_data_corrected <-WT_data[!(Assay_Date%in%outliers_dt[par==single_par, Assay_Date])]
      }else{WT_data_corrected =WT_data}
    }else{WT_data_corrected =WT_data}
    #   print(single_par)
    one_par=prepare_RR_whole (WT_data=WT_data_corrected,
                              split_date=split_date,
                              split_by_sex=split_by_sex,
                              min_ind=min_ind,
                              single_par=single_par,
                              assay_dates=assay_dates,
                              percentile=percentile)[[1]]
    one_par$param <- single_par
    return(one_par)
  })
  
  whole_dt_t <- Reduce(rbind, whole_dt_l)
  return(whole_dt_t)
} 

prepare_RR_whole <- function (WT_data,
                              split_date=NULL,
                              split_by_sex=TRUE,
                              min_ind=100,
                              single_par=NULL,
                              assay_dates=NULL,
                              percentile=0.025){
  #this version takes ALL WT data as a dt. RR will be identical to all individuals of a given sex for a parameter (within split borders).
  WT_data <-data.table(WT_data)
  WT_data<- WT_data[order(Assay_Date),]
  
  if(is.null(assay_dates)){
    assay_dates <- unique(WT_data$Assay_Date)
  }
  assay_dates <- sort(as.Date(assay_dates))
  
  parameters <-list(WT_data_dim=dim(WT_data),
                    split_date=split_date,
                    sex_split=split_by_sex,
                    min_ind=min_ind,
                    single_par=single_par)
  
  if(!is.null(split_date)){
    print("split date detected")
    split_date <-as.Date(split_date)
    WT_data_after <- WT_data[Assay_Date >= split_date,]#print("before split")
    WT_data_before <- WT_data[Assay_Date < split_date,]
    assay_dates_after <-assay_dates[assay_dates>=split_date]
    assay_dates_before <-assay_dates[assay_dates<split_date]
    return(list(
      res=rbind(prepare_RR_whole(WT_data_before,split_date=NULL, split_by_sex=split_by_sex, min_ind=min_ind, single_par=single_par, assay_dates=assay_dates_before , percentile=percentile)[["res"]],
                data.frame(date=split_date,lower=NA, upper=NA, variable=c("Female", "Male")),
                prepare_RR_whole(WT_data_after,split_date=NULL, split_by_sex=split_by_sex, min_ind=min_ind, single_par=single_par, assay_dates=assay_dates_after, percentile=percentile)[["res"]]), parameters=parameters)
    )
  }
  
  if(split_by_sex){  #if RR is to be split by sex - divide data into male andfemale subset, run extract_data_for_local_RR separately and combine outputs
    data_f <- WT_data[Sex=="Female",]
    data_m <- WT_data[Sex=="Male",]
    local_RR_f <- prepare_RR_whole( data_f, split_date=NULL, split_by_sex=FALSE, min_ind=min_ind, single_par=single_par, percentile=percentile,assay_dates=assay_dates)
    #   print(local_RR_f)
    local_RR_m <- prepare_RR_whole( data_m, split_date=NULL, split_by_sex=FALSE, min_ind=min_ind, single_par=single_par, percentile=percentile,assay_dates=assay_dates)
    #   print(local_RR_m)
    res <- rbind(local_RR_f=data.table(local_RR_f$res, sex="Female"),local_RR_m=data.table(local_RR_m$res, sex="Male"))
    #   print(res)
    res <-res[order(as.Date(res$date)),]
    row.names(res)=NULL
    return(list(res=res, parameters=parameters))
  }
  
  #Here computes the RR
  available_datapoints <- WT_data[, sum(!is.na(get(single_par)))]
  if(available_datapoints >= min_ind){
    results <- get_percentiles_range_variable(WT_data[,get(single_par)], percentile)
  }else{
    print(paste(single_par, "not enough datapoints to get RR"))
    results <- c(NA, NA)
  }
  
  RRs_matrix <- matrix(rep(results, length(assay_dates) ), nrow=length(assay_dates) ,byrow = TRUE)
  colnames(RRs_matrix ) <- c("lower", "upper")
  RRs <- data.table(RRs_matrix )
  RRs$date <- assay_dates
  return(list(res=RRs, parameters=parameters))
  
}


identify_outliers <- function(dt,  coef_for_sd = 2.5, colns=NULL){
  if(is.null(colns)){
    colns <-grep("proc", colnames(dt), val=T)
  }
  outliers <- lapply(colns, function(x){
    a <-identify_outlying_days(take_only_WTs(dt), par=x, coef_for_sd = coef_for_sd )
    a[, wt_day:="YES"]
    
    print(x)
    b <-identify_outlying_days(dt = dt,  assay_dates = unique(dt[!(Genotype%in%c("WT","+/+", "+/Y")),Assay_Date]), par=x, coef_for_sd = coef_for_sd , extraCare = TRUE)
    b[, wt_day:="NO"]
    
    out <- rbindlist(list(a,b))
    out[, par:=x]
    return(out)
  })
  outliers_dt <-rbindlist(outliers)
  return(outliers_dt)
}

identify_outlying_days <- function(dt,
                                   par,
                                   coef_for_sd = 2,
                                   extraCare = FALSE,
                                   assay_dates = NULL,
                                   full_results = FALSE) {
  #extraCare: recommended when running on days without WT: will ty to get a subset of KO samples with a balanced set of genotypes/sex
  #also, will base general median on WT days
  #outputs only
  if (extraCare) {
    med_and_sd <-
      dt[Genotype %in% c("WT", "+/Y", "+/+")][, .(m = median(get(par), na.rm =
                                                               T), N = sum(!is.na(get(par)))), by = Assay_Date][N > 1, .(med = median(m, na.rm =
                                                                                                                                        T), sd = mad(m, na.rm = TRUE))]
    med <- unlist(med_and_sd[1, med])
    sd <- unlist(med_and_sd[1, sd])
    #  print(med_and_sd)
    balanced_subsets <-
      balance_draw_one_parameter(par, dt[Assay_Date %in% assay_dates], assay_dates =
                                   assay_dates)
    #  print(balanced_subsets)
    not_called <-
      assay_dates[unlist(lapply(balanced_subsets, is.null))]
    balanced_subsets_withoutNA <-
      balanced_subsets[sapply(balanced_subsets, function(x)
        ! is.null(x))]
    balanced_subsets <- Reduce(rbind, balanced_subsets_withoutNA)
    
    res <-
      balanced_subsets[, .(m = median(get(par), na.rm = T), N = sum(!is.na(get(par)))), by =
                         Assay_Date][, .(is_outside = abs(m - med) > coef_for_sd * sd, N), by = Assay_Date][is_outside ==
                                                                                                              TRUE & N > 2]
    
  } else{
    med_and_sd <-
      dt[, .(m = median(get(par), na.rm = T), N = sum(!is.na(get(par)))), by =
           Assay_Date][N > 1, .(med = median(m, na.rm = T), sd = mad(m))]
    
    #   print(med_and_sd )
    med <- unlist(med_and_sd[1, med])
    sd <- unlist(med_and_sd[1, sd])
    res <-
      dt[, .(m = median(get(par), na.rm = T), N = sum(!is.na(get(par)))), by =
           Assay_Date][, .(is_outside = abs(m - med) > coef_for_sd * sd, N), by = Assay_Date][is_outside ==
                                                                                                TRUE & N > 2]
    not_called <-
      dt[, .(m = median(get(par), na.rm = T), N = sum(!is.na(get(par)))), by =
           Assay_Date][, .(N), by = Assay_Date][N <= 2, Assay_Date]
    balanced_subsets <- dt[Assay_Date %in% assay_dates, get(par)]
  }
  if (full_results) {
    return(
      list(
        res = res,
        not_called = not_called,
        balanced_subsets = balanced_subsets,
        med_and_sd = med_and_sd
      )
    )
  } else{
    return(res)
  }
}

# spleen_data_auto_loess[!Assay_Date%in%spleen_data_auto_loess[Genotype%in%c("+/+", "+/Y"),Assay_Date,],Assay_Date]


balance_draw_one_parameter <- function(par, dt, assay_dates) {
  lapply(assay_dates, function(ud) {
    one_day_data <-
      dt[, c("Gene_Name", "Sex", par, "Assay_Date"), with = F][!is.na(get(par))][Assay_Date ==
                                                                                      ud]
    #print(one_day_data )
    one_day_data_balanced <- NULL
    if (nrow(one_day_data) > 2) {
      one_day_data_balanced <- get_a_balanced_subset(one_day_data)
    }
    else{
      one_day_data_balanced <- NULL
    }
    # return(list(one_day_data=one_day_data, one_day_data_balanced=one_day_data_balanced))
    return(one_day_data_balanced)
  })
}

get_a_balanced_subset <- function(one_day_data) {
  #takes only one animal per sex and strain
  one_day_data_collapsed <-
    one_day_data[, .SD[1], by = .(Gene_Name, Sex), .SDcols = colnames(one_day_data)]
  how_many_one_sex <-
    one_day_data_collapsed [, .SD[1], by = .(Gene_Name, Sex), .SDcols =
                              colnames(one_day_data)][, .N, by = Sex][, min(N)]
  if (length(unique(one_day_data_collapsed$Sex)) > 1) {
    which_sex <-
      one_day_data_collapsed [, .SD[1], by = .(Gene_Name, Sex), .SDcols =
                                colnames(one_day_data)][, .N, by = Sex][order(Sex)][, c("Female", "Male")[which.min(N)]]
    res = rbind(one_day_data_collapsed[Sex == which_sex], one_day_data_collapsed[!Sex ==
                                                                                      which_sex][sample(1:nrow(one_day_data_collapsed[!Sex == which_sex]), how_many_one_sex)])
    return(res)
  } else{
    return(NULL)
  }
}

find_perc_simple_dt <-function(WT_dt,  params,  dt,  cols_to_preserve=NULL,outliers_dt = NULL ){
  if(is.null(cols_to_preserve))
  { cols_to_preserve=setdiff(colnames(dt),params)}
  lisRes <- lapply(params, function(param){
    print(param)
    if(!is.null(outliers_dt)){
      if(param%in%outliers_dt$par){
        #   print(outliers_dt[par==single_par])
        WT_data_corrected <-WT_dt[!(Assay_Date%in%outliers_dt[par==param, Assay_Date])]
      }else{WT_data_corrected =WT_dt}
    }else{WT_data_corrected =WT_dt}
    
    find_perc_simple_par(WT_data_corrected,
                         param=param,
                         values=dt[,param, with=F],
                         sex_vals=dt[,Sex],
                         dt=dt,
                         other_cols_to_keep=cols_to_preserve)
  })
  mergef <-function(a, b){merge(a,b, by=cols_to_preserve, all=T)}
  Reduce(mergef, lisRes)
}


find_perc_simple_par <-function(WT_dt,  param,  values, sex_vals, dt, other_cols_to_keep){
  u=lapply(c("Male", "Female"), function(gender){
    a=find_perc_simple_gender_par(WT_dt, gender=gender, param,  value=unlist(values[sex_vals==gender]))
    b=data.table(dt[sex_vals==gender, other_cols_to_keep, with=F], a)
    colnames(b) <- c(other_cols_to_keep, param)
    return(b)
  })
  Reduce(rbind, u)
  
}

find_perc_simple_gender_par <-
  function(WT_dt, gender, param,  value){
    # works for a vector to
    # print(EndRR)
    db_of_parameter <- unlist(WT_dt[Sex==gender, param, with=F])
    # print(db_of_parameter 
    find_perc_value(value=value, db_of_parameter=db_of_parameter)
  }

find_perc_value <-function(value, db_of_parameter){
  Fn <- ecdf(db_of_parameter)
  return(Fn(value))
}


prepare_RR_byindnNum <-
  function(spleen_data_dt,
           WT_spleen_dt,
           minInd,
           colns = colns[proc_colns],
           tests1,
           samples_per_par_spleen,
           hard_split_date = NULL) {
    #apropriate when there is a temporal split
    #uses compute_RR_befaft
    #RR After splits:
    split_date = hard_split_date
    RR_spleen_81_ba_after <- lapply(colns , function(par) {
      print(par)
      if (is.null(hard_split_date)) {
        if (is.na(tests1$panel[which(tests1$subset == par)]) |
            tests1$panel[which(tests1$subset == par)] != "BCpanel") {
          split_date <- as.Date("2014-06-09")
        } else{
          split_date <- as.Date("2014-09-22")
        }
      } else{
        split_date = hard_split_date
      }
      
      print(split_date)
      days_without_WT = setdiff(spleen_data_dt[Assay_Date >= split_date, Assay_Date], WT_spleen_dt[Assay_Date >=
                                                                                                     split_date, Assay_Date])
      class(days_without_WT) <- "Date"
      compute_RR_befaft(
        WT_spleen_dt,
        split_by_sex = TRUE,
        minInd = minInd,
        par = par,
        percentile = 0.025,
        days_without_WT = days_without_WT,
        start = split_date,
        end = NULL
      )
    })
    
    names(RR_spleen_81_ba_after) <- colns
    
    if(split_date > min(spleen_data_dt[,Assay_Date])){ 
      RR_spleen_81_ba_before <- lapply(colns, function(par) {
        print(par)
        if (is.null(hard_split_date)) {
          if (is.na(tests1$panel[which(tests1$subset == par)]) |
              tests1$panel[which(tests1$subset == par)] != "BCpanel") {
            split_date <- as.Date("2014-06-09")
          } else{
            split_date <- as.Date("2014-09-22")
          }
        } else{
          split_date = hard_split_date
        }
        print(split_date)
        samples <- samples_per_par_spleen[par]
        
        days_without_WT = setdiff(spleen_data_dt[Assay_Date < split_date, Assay_Date],
                                  WT_spleen_dt[Assay_Date < split_date, Assay_Date])
        class(days_without_WT) <- "Date"
        if (samples >= minInd) {
          compute_RR_befaft (
            WT_spleen_dt,
            split_by_sex = TRUE,
            minInd = minInd,
            par = par,
            percentile = 0.025,
            end = split_date - 1,
            days_without_WT = days_without_WT
          )
        }
        else{
          compute_RR_simple (
            WT_spleen_dt[Assay_Date < split_date],
            split_by_sex = TRUE,
            par = par,
            percentile = 0.025,
            days_without_WT = days_without_WT
          )
        }
      })
      
      names(RR_spleen_81_ba_before) <- colns
    }
    
    RR_spleen_81_ba_dt_all <-
      lapply(names(RR_spleen_81_ba_after), function(par) {
        if(split_date >min(spleen_data_dt[,Assay_Date])){
          x = data.table(rbind(RR_spleen_81_ba_before[[par]], RR_spleen_81_ba_after[[par]]))
        }
        else{
          x = data.table(RR_spleen_81_ba_after[[par]])
        }
        x[, startRR := as.Date(startRR)]
        x[, endRR := as.Date(endRR)]
        x[, Assay_Date := as.Date(Assay_Date)]
        x$par <- par
        return(x)
      })
    names(RR_spleen_81_ba_dt_all) <- names(RR_spleen_81_ba_after)
    
    
    RR_spleen_81_ba_dt_all_dt <-
      unique(data.table(do.call("rbind", RR_spleen_81_ba_dt_all)))[order(Assay_Date)]
    setkeyv(RR_spleen_81_ba_dt_all_dt,
            cols = c("Assay_Date", "Sex", "par"))
    return(RR_spleen_81_ba_dt_all_dt)
  }




detect_hits_byparam_cp_1 <- function(param, male_genotype, data_to_extract_without_WT, RR, sex_column_RR="Sex", cutoff_in=5)
  {
  male_genotype <-as.character(male_genotype)
  if(grepl("(.*)\\/Y",male_genotype)){
    female_genotype <- sub("(.*)\\/Y$","\\1\\/\\1",male_genotype)
  }else{
    female_genotype <- male_genotype
  }
  genotype_to_pick<-unique(c(female_genotype, male_genotype))
  # print(genotype_to_pick) 
  setkeyv(data_to_extract_without_WT, cols=c("Genotype","Assay_Date", "Sex", param))
  data_param <-  data_to_extract_without_WT[Genotype%in%genotype_to_pick,c("Colony_Prefix", "Genotype", "Gene_Name", "Assay_Date", "Sex", "Mouse", param), with=F] 
  
  data_param [,par:=param] 
  data_param <-data_param [!is.na(param),]
  data_param[,below:= NA]
  data_param[,above:= NA]
  #  print(data_param)
  #  print(RR)
  if(nrow(data_param)<cutoff_in){  # <5????
    #  print(data_param)
    return(data.frame(param=param, data_param[1,.(Colony_Prefix,Gene_Name)] , sum_below=NA, sum_above=NA, male_genotype=male_genotype, female_genotype=female_genotype, samples=nrow(data_param), comm="not enough samples" ))}
  
  # print(data_param)  
  # other_data_the_same_days <-  data_to_extract_without_WT[Genotype%in%genotype_to_pick,c("Colony_Prefix", "Genotype", "Gene_Name", "Assay_Date", "Sex", "Mouse", param), with=F] 
  # print(RR)
  #  print(data_param)
  towork <- merge(data_param, RR, by=c("Assay_Date","Sex","par"))
  #  print(towork)
  towork[,below:=get(param) < lower]
  towork[,above:=get(param) > upper]
  
  sum_below <- sum(towork$below, na.rm=T)
  sum_above <- sum(towork$above, na.rm=T)
  samples <- nrow(towork)
  return(data.frame(param=param, data_param[1,.(Colony_Prefix,Gene_Name)], sum_below=sum_below,sum_above=sum_above, male_genotype=male_genotype, female_genotype=female_genotype, samples=samples, comm="OK" ))
}


#try the same but by matching with KOs fromother than this param background
match_by_factor_closest_n_eliminate_subset=
  #does the same as match_by_factor_closest_n, but excludes from matching a subset of individuals, sharing with the tested one a characteristi - for example, day of test, genotype etc.
  function(data_WT, data_one_ind_KO, parameter_to_match, closest_n=70, byGend=F, param_for_exclusion=NULL){
    #If tolerance_in_perc, it will take fraction tolerance_in_perc of the val_of_param_to_match
    val_of_param_to_match <- unlist(data_one_ind_KO[,parameter_to_match, with=F])
    
    #print(val_of_param_to_match)
    if(byGend){
      #   print(data_WT[Sex==unlist(data_one_ind_KO[,Sex])])
      return(match_by_factor_closest_n_eliminate_subset(data_WT[Sex==unlist(data_one_ind_KO[,Sex])],data_one_ind_KO,  parameter_to_match, closest_n=70, byGend=F,param_for_exclusion=param_for_exclusion)[1:closest_n]) 
    }else{
      if(!is.null(param_for_exclusion)){
        value_of_param_for_exclusion <- unlist(data_one_ind_KO[,param_for_exclusion, with=F])
        #   print(value_of_param_for_exclusion)
        data_for_match <- data_WT[get(param_for_exclusion)!=value_of_param_for_exclusion]
      }else{
        data_for_match <- data_WT
      }
      return(data_for_match[order(abs(get(parameter_to_match)-val_of_param_to_match))][1:closest_n])}
  }


match_and_call <-function(pars_for_testing,
                          parameter_to_match,
                          dss_data=dss_AUCt_histo[period=="C"],
                          min_number_of_animals = 3,
                          lower_cutoff = 0.025,
                          upper_cutoff = 0.975,
                          byGend=TRUE)
{
  print(paste0("Finding controls matching by sex and ",parameter_to_match,"..."))
  matched_WT <- lapply(1:nrow(dss_data[Genotype!="WT"]), function (a){
    list(tested=dss_data[Genotype!="WT"][a],
         wt=match_by_factor_closest_n_eliminate_subset(dss_data[Genotype=="WT"][Colony_Prefix=="CBSC"],
                                                       dss_data[Genotype!="WT"][a],closest_n =70, parameter_to_match = parameter_to_match, byGend = byGend))
  })
  
  
  print("Establishing location within matched WT for each KO animal...") 
  tests_sets <-
    lapply(pars_for_testing , function(parameter){ 
      print(parameter)
      sapply(1:nrow(dss_data[Genotype!="WT"]), function(x){
        locate_percentile <- ecdf(unlist(matched_WT[[x]]$wt[,parameter, with=F]))
        locate_percentile(matched_WT[[x]]$tested[, parameter, with=F])
      })
    })
  
  names(tests_sets ) <-pars_for_testing
  
  print("Heatmapping weight loss...")
  
  results_matched <- lapply(names(tests_sets), function(x){
    a=dss_data[Genotype!="WT"]
    a[,sigU:=tests_sets[[x]]>= upper_cutoff]
    a[,sigL:=tests_sets[[x]]<=lower_cutoff]
    a[,.(par_mean=round(mean(get(x), na.rm=T),2),
         age=round(mean(Age_In_Weeks, na.rm=T),1),
         .N,
         samples_above_cutoff=sum(sigU, na.rm=T),
         samples_below_cutoff=sum(sigL, na.rm=T), par=x),
      by=.(Colony_Prefix,Gene_Name,Genotype)][N>=min_number_of_animals]
    
  })%>%
    bind_rows()
  return(results_matched )
  
}

