#'---
#' title: prediction_emb_tab
#' author: Xueqi Cao
#' wb:
#'  py:
#'   - |
#'    annotations = config["geneAnnotation"]
#'    inputDatasets = config["cohort"]["single_group"]
#'    outputDatasets = config["cohort"]["groups"]
#'  params:
#'    - projectPath: '`sm config["projectPath"]`'
#'    - annotations: '`sm annotations`'
#'    - gencode: '`sm config["gencode"]`'
#'    - sampAnno: '`sm config["sampleAnnotation"]`'
#'    - experimentDesign: '`sm config["experimentDesign"]`'
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - intogenDir: '`sm config["intogenDir"]`'
#'    - vep_path: '`sm config["vep_path"]`'   
#'    - mutsigCVdir: '`sm config["mutsigCVdir"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'    - manuscriptWording: '`sm config["manuscriptWording"]`'
#'  output:
#'    - prediction_emb_all_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/prediction_emb_all_tab.csv"`'
#'    - prediction_emb_diag_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/prediction_emb_diag_tab.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/prediction_emb_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/prediction_emb_tab.snakemake")
print("Snakemake saved") 



suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})
source("Scripts/manuscript/function.R")

intogen_dir <- snakemake@params$intogenDir
vep_dir <- snakemake@params$vep_path
mutsigcv_dir <- snakemake@params$mutsigCVdir
project_dir <- snakemake@params$projectPath
experiment_design <- fread(snakemake@params$experimentDesign)
experiment_design[is.na(experiment_design)] <- ''

single_group <- snakemake@params$inputDatasets
sample_groups <- snakemake@params$outputDatasets

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))

#### all sample #### 
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == single_group &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature != '', ]

exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)


# vale
if(exists("res_leu14")){remove(res_leu14)}

for (i in 1:nrow(exp_viz)) {
  res_temp <- fread(exp_viz[i, res_post_path])
  
  sample_group <- exp_viz[i, sample_group]
  sample_group_manuscript <- manuscript_wording[Study_group_during_analysis==sample_group, Study_group] %>% unique()
  res_temp[, StudyGroup := 'Pan-leukemia']
  res_temp[, ExternalFeature := exp_viz[i, coess_input_feature]]
  
  res_temp <- res_temp[, .( GeneID, GeneSymbol, StudyGroup, ExternalFeature, Label, Rank, Prediction, isCGC, RoleCGC, isIntOGen, RoleIntogen)]
  
  if(!exists("res_leu14")){
    res_leu14 <- res_temp
  }else{
    res_leu14 <- rbind(res_leu14, res_temp, fill=TRUE)
  }}

fwrite(res_leu14, snakemake@output$prediction_emb_all_tab) 


#### each group #### 
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group %in% sample_groups &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature != '', ]

exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)

# read in result for visualization
if(exists("res_groups")){remove(res_groups)}

for (i in 1:nrow(exp_viz)) {
  res_temp <- fread(exp_viz[i, res_post_path])
  
  sample_group <- exp_viz[i, sample_group]
  sample_group_manuscript <- manuscript_wording[Study_group_during_analysis==sample_group, Study_group] %>% unique()
  res_temp[, StudyGroup := sample_group_manuscript]
  res_temp[, ExternalFeature := exp_viz[i, coess_input_feature]]
  
  res_temp <- res_temp[, .( GeneID, GeneSymbol, StudyGroup, ExternalFeature, Label, Rank, Prediction, isCGC, RoleCGC, isIntOGen, RoleIntogen)]
  
  if(!exists("res_groups")){
    res_groups <- res_temp
  }else{
    res_groups <- rbind(res_groups, res_temp, fill=TRUE)
  }}

fwrite(res_groups, snakemake@output$prediction_emb_diag_tab) 