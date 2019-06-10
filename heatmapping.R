require(data.table)
require(dplyr)


source("./heatmapping_suite.R")

#Bone marrow
# Returns per-strain results ($heatmapped_data) as well as candidate hits ($hit_list_data).
# Final calls were made with carefull inspection of experiments which
#contained candidate genotypes by an experienced reseacher.

file_manual <- "data/3i_bulkDown_manualBM_190604.csv"
file_automated <- "data/3i_bulkDown_autoBM_190604.csv"
BM_results <- heatmap_tissue(file_manual, file_automated)  


#Spleen
# Returns per-strain results ($heatmapped_data) as well as candidate hits ($hit_list_data).
# Final calls were made with carefull inspection of experiments which
#contained candidate genotypes by an experienced reseacher.
file_manual <- "data/3i_bulkDown_manualspleen_190604.csv"
file_automated <- "data/3i_bulkDown_autospleen_190604.csv"
spleen_results <- heatmap_tissue(file_manual, file_automated)  



#Mesenteric Lymph Node
# Returns per-strain results ($heatmapped_data) as well as candidate hits ($hit_list_data).
# Final calls were made with carefull inspection of experiments which
#contained candidate genotypes by an experienced reseacher.

file_manual <- "data/3i_bulkDown_manualMLN_190604.csv"
file_automated <- "data/3i_bulkDown_autoMLN_190604.csv"
MLN_results <- heatmap_tissue(file_manual, file_automated)  


#Ear epidermis
# Final calls were made with carefull inspection of experiments which
#contained candidate genotypes by an experienced reseacher.
#In the output all_above, all_below,  N: number of samples above the RR, below the RR and tested
#   above_susp below_susp N_susp : number of samples above the RR on outlier days, below the RR on outlier days and tested for being  on outlier days
ear_file <- "data/3i_bulkDown_earEpidermis_190604.csv"
ear_results <- heatmap_ear(ear_file)


#Blood 
#Returns per-strain results for both sexes together, males-only and females-only testing
#Genotypes leading to any of (male-specific, female-specific, both sexes specific) calls were considered as candidates.
# Final calls were made with carefull inspection of experiments which
#contained candidate genotypes by an experienced reseacher.
blood_file <- "data/3i_bulkDown_blood_190604.csv"
blood_results <- heatmap_blood(blood_file)



#DSS There were 3 periods in this experiment, differing by DSS used and called
#A-C. Strength and effects of the chemical were different,as judged on WT data;
#hence the results from those periods cannot be compared directly. Only period C
#is heatmapped here, as periods A-B were analysed manually be researchers.
#
#Input data contains info about pre-day10 death, however this was tested
#independently. Weight values in input data were transformed into: difference
#from starting weight, as a  fraction of  starting weight (w_d7). Minimum of
#such transformed weights across alll measurements of a given mouse (weigth_min)
#Area under the curve for transformed weight versus experimental days, till day
#7 (auc_d7).
#
#Result contains histology part and weight loss part. Final calls were made
#based on those & premature death, with carefull inspection of experiments which
#contained candidate genotypes by an experienced reseacher.

dss_file <- "data/3i_bulkDown_DSS_190604.csv"
dss_results <- heatmap_dss(dss_file)
