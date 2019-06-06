require(data.table)
require(dplyr)

source("~/Documents/projects/3i/3i_paper_draft/NI_2018/after_reviewers/scripts_for_webpage/functions_heatmapping.R")

heatmap_tissue <-function(file_manual, file_automated){
  require(data.table)
  require(dplyr)
  print("Reading in files...")
  
  data_manual <- read.csv(file_manual, header=TRUE,  as.is=TRUE, check.names = FALSE)%>%
    data.table()%>%
    .[,Assay_Date:=as.Date(Assay_Date, format = "%Y-%m-%d")]
  
  
  data_automated <- read.csv(file_automated, header=TRUE,  as.is=TRUE, check.names = FALSE)%>%
    data.table()%>%
    .[,Assay_Date:=as.Date(Assay_Date, format = "%Y-%m-%d")]
  
  
  print("Dividing according to panel split and gating type...")
  
  #Divide manual gating results into old and new panel (new has parallel automated data)
  manual_data_presplit <- data_manual[Assay_Date < min(data_automated$Assay_Date)]
  manual_data_postsplit <- data_manual[Assay_Date >= min(data_automated$Assay_Date)]
  
  ########  Prepare Reference Ranges ######## 
  
  colnames_for_RR_data <- intersect(grep("%",colnames(data_automated), val=T), grep("%",colnames(data_manual), val=T))
  
  ##for automated gating (time-window RR)
  
  print("Preparing automated gating RR...")
  
  RR_data_auto <- prepare_RR(all_data_dt = data_automated, 
                           cols_for_RR  = colnames_for_RR_data,
                           minInd = 70)
  
  
  ##for manual gating (time-window RR)
  print("Preparing manual post-split RR...")
  
  RR_data_post_manual <- prepare_RR(all_data_dt = manual_data_postsplit,
                                  cols_for_RR =  colnames_for_RR_data,
                                  minInd = 70)
  
  print("Preparing manual pre-split RR...")
  
  ##for manual gating pre-panel change (all-samples RR)
  RR_data_pre_manual <- prepare_RR_simple_dt(all_data_dt = manual_data_presplit,
                                           cols_for_RR=colnames_for_RR_data  ) 
  
  ########  Check each mice whether it sits within its reference range ######## 
  print("Checking whether mice are in RR range...")
  
  outside_data_auto <- check_if_mouse_outside_dt (all_data_dt = data_automated ,
                                                RR_data_auto, pars=colnames_for_RR_data  )
  
  outside_data_post_manual <- check_if_mouse_outside_dt(all_data_dt = manual_data_postsplit ,
                                                      RR_data_post_manual,
                                                      pars=colnames_for_RR_data  ) 
  outside_data_pre_manual <- check_if_mouse_outside_dt(all_data_dt = manual_data_presplit ,
                                                     RR_data_pre_manual,
                                                     pars=colnames_for_RR_data  ) 
  
  BM_outside <- combine_outside(manual_fract = data_manual, 
                                outside_auto = outside_data_auto,
                                outside_post_manual = outside_data_post_manual,
                                outside_pre_manual = outside_data_pre_manual,
                                colnames_for_RR =colnames_for_RR_data,
                                isMatrix = TRUE)
  
  ########  Reformat and identify hits ######## 
  print("Identyfying hits...")
  
  
  add_zygosity(BM_outside$res, new_format = TRUE)
  melted_data <- melt(BM_outside$res[zygosity%in%c("Het","Hom")][Genotype!="WT"][,c("Mouse","Colony_Prefix","Genotype","zygosity", "Gene_Name","Sex",colnames_for_RR_data), with=F],
                    id.vars = c("Mouse","Genotype","zygosity", "Gene_Name","Sex","Colony_Prefix"))
  heatmapped_data  <- melted_data[zygosity%in%c("Het","Hom")][,tableh(value) ,by=.(Gene_Name, zygosity, Colony_Prefix, variable) ]
  
  hit_list_data <- heatmapped_data[(samples==3&(below==3|above==3))|
                                 (samples==4&(below>=3|above>=3))|
                                 (samples>=5&(below/samples>=0.6))|
                                 (samples>=5&(above/samples>=0.6))]
  
  return(list(hit_list_data=hit_list_data,
              heatmapped_data=heatmapped_data))
  
}

heatmap_ear <- function(ear_file) {
  ear_data <- read.csv(ear_file, header = TRUE, as.is = TRUE) %>%
    data.table()
  
  cols_to_analyse_ear <- colnames(ear_data)[8:33]
  
  print("Identyfying suspected samples - to eliminate from RR...")
  ear_outliers_raw <- lapply(cols_to_analyse_ear, function(par) {
    identify_outliers(ear_data[, c("Colony_Prefix",
                                   "Gene_Name",
                                   "Mouse" ,
                                   "Genotype",
                                   "Sex",
                                   "Assay_Date",
                                   par), with = F],
                      coef_for_sd = 4, colns = par)
  }) %>% bind_rows()
  
  ear_perc_raw <- find_perc_simple_dt(
    WT_dt = take_only_WTs(ear_data),
    dt = ear_data,
    params =  cols_to_analyse_ear,
    outliers_dt = ear_outliers_raw
  )
  
  ear_perc_raw <- add_zygosity(ear_perc_raw, new_format = TRUE)
  
  print("Identyfying mice outside RR...")
  
  ear_perc_raw[Gene_Name=="",Gene_Name:=gsub(":.*","",Genotype)]
  
  ear_heatmapping_raw <- lapply(cols_to_analyse_ear, function(spar) {
    print(spar)
    a = unique(ear_perc_raw[!Genotype %in% c("WT", "+/+", "+/Y")][, .(
      par = spar,
      all_above =
        sum(get(spar) >= 0.975, na.rm = T),
      all_below =
        sum(get(spar) <= 0.025, na.rm = T),
      N = sum(!is.na(get(spar)), na.rm =
                T)
    ),
    by = .(Gene_Name, zygosity)])
    #  print(a)
    c = NULL
    if (nrow(a) > 0) {
      if (spar %in% ear_outliers_raw$par) {
        suspected <- ear_outliers_raw[par == spar, Assay_Date]
        #      print(suspected)
        b <-
          unique(ear_perc_raw[!Genotype %in% c("WT", "+/+", "+/Y")][, .(
            above_susp = sum((get(spar) > 0.975) &
                               (Assay_Date %in% suspected), na.rm = T),
            below_susp =
              sum((get(spar) < 0.025) & (Assay_Date %in% suspected), na.rm = T),
            N_susp =
              sum(!is.na(
                get(spar) &
                  Assay_Date %in% suspected
              ), na.rm = T)
          ),  by = .(Gene_Name, zygosity)])
        c <- merge(a,
                   b,
                   by = c("Gene_Name", "zygosity"),
                   all.x = TRUE)
        #     print(c)
      }
    }
    
    if (is.null(c)) {
      c <- a
      c[, c("above_susp", "below_susp", "N_susp") := list(0, 0, 0)]
    }
    return(c)
  }) %>% bind_rows()
  
  return(ear_heatmapping_raw)
}

heatmap_blood <- function(blood_file) {
  print("Reading in the file...")
  flow16w_data_dt <-
    read.csv(blood_file, header = TRUE, as.is = TRUE) %>%
    data.table()
  
  colns <- colnames(flow16w_data_dt)
  data_colns <- (1:length(colns))%in%grep("proc|perc|number", colns,ignore.case = TRUE)
  WT_flow16w_dt <-take_only_WTs(flow16w_data_dt )
  WT_flow16w_dt[,Assay_Date:=as.Date(Assay_Date)]
  
  print("Preparing RR...")
  RR_flow16w_70_ba <-
    prepare_RR_byindnNum(
      flow16w_data_dt,
      WT_flow16w_dt,
      minInd = 70,
      colns = colns[data_colns],
      hard_split_date = min(flow16w_data_dt$Assay_Date)
    )
  
  print("Heatmapping both sexes together...")
  hit_calling_flow16w <-
    lapply(unique(RR_flow16w_70_ba$par), function(par) {
      print(par)
      x = do.call("rbind", lapply(unique(flow16w_data_dt[(Genotype != "+/+" &
                                                            Genotype != "+/Y" &
                                                            Genotype != "WT" & Sex == "Male"), Genotype]), function(x) {
                                                              detect_hits_byparam_cp_1(
                                                                param = par,
                                                                male_genotype = x,
                                                                data_to_extract_without_WT = flow16w_data_dt,
                                                                RR = RR_flow16w_70_ba,
                                                                sex_column_RR = "Sex",
                                                                cutoff_in = 4
                                                              )
                                                            }))
      return(x)
    }) %>%
    bind_rows() %>%
    data.table()
  
  
  print("Heatmapping males separately...")
  hit_calling_flow16w_male <-
    lapply(unique(RR_flow16w_70_ba$par), function(par) {
      print(par)
      x = do.call("rbind", lapply(unique(flow16w_data_dt[(Genotype != "+/+" &
                                                            Genotype != "+/Y" &
                                                            Genotype != "WT" & Sex == "Male"), Genotype]), function(x) {
                                                              detect_hits_byparam_cp_1(
                                                                param = par,
                                                                male_genotype = x,
                                                                data_to_extract_without_WT = flow16w_data_dt[Sex == "Male"],
                                                                RR = RR_flow16w_70_ba,
                                                                sex_column_RR = "Sex",
                                                                cutoff_in = 4
                                                              )
                                                            }))
      return(x)
    }) %>%
    bind_rows() %>%
    data.table()
  
  print("Heatmapping females separately...")
  hit_calling_flow16w_female <-
    lapply(unique(RR_flow16w_70_ba$par), function(par) {
      print(par)
      x = do.call("rbind", lapply(unique(flow16w_data_dt[(Genotype != "+/+" &
                                                            Genotype != "+/Y" & Sex == "Male"), Genotype]), function(x) {
                                                              detect_hits_byparam_cp_1(
                                                                param = par,
                                                                male_genotype = x,
                                                                data_to_extract_without_WT = flow16w_data_dt[Sex == "Female"],
                                                                RR = RR_flow16w_70_ba,
                                                                sex_column_RR = "Sex",
                                                                cutoff_in = 4
                                                              )
                                                            }))
      return(x)
    }) %>%
    bind_rows() %>%
    data.table()
  
  return(list(both_sexes=hit_calling_flow16w,
              female=hit_calling_flow16w_female,
              male=hit_calling_flow16w_male))
}