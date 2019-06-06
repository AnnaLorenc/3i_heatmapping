require(data.table)
require(dplyr)

#change this line to your location of 'functions_heatmapping.R'

source("~/Documents/projects/3i/3i_paper_draft/NI_2018/after_reviewers/scripts_for_webpage/heatmapping_suite.R")

#Bone marrow
file_manual <- "data/3i_bulkDown_manualBM_190604.csv"
file_automated <- "data/3i_bulkDown_autoBM_190604.csv"
BM_results <- heatmap_tissue(file_manual, file_automated)  

#Spleen
file_manual <- "data/3i_bulkDown_manualspleen_190604.csv"
file_automated <- "data/3i_bulkDown_autospleen_190604.csv"
spleen_results <- heatmap_tissue(file_manual, file_automated)  

#Mesenteric Lymph Node
file_manual <- "data/3i_bulkDown_manualMLN_190604.csv"
file_automated <- "data/3i_bulkDown_autoMLN_190604.csv"
MLN_results <- heatmap_tissue(file_manual, file_automated)  


#Ear epidermis
ear_file <- "data/3i_bulkDown_earEpidermis_190604.csv"
ear_results <- heatmap_ear(ear_file)


#Blood 
blood_file <- "data/3i_bulkDown_blood_190604.csv"
blood_results <- heatmap_blood(blood_file)
