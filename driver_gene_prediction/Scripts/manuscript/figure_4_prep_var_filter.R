#'---
#' title: figure_4_prep_var_filter
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
#'  output:
#'    - mll_var_filter: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_var_filter.tsv"`'
#'    - exp_viz_var_filter: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/exp_viz_var_filter.tsv"`'
#'  type: script
#'  resources:
#'    - mem_mb: 8000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_4_prep_var_filter.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/figure_4_prep_var_filter.snakemake")
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



#### intOGen var filter #### 
project_dirs <- c(
  '/s/project/vale/driver_prediction_202401_0',
  '/s/project/vale/driver_prediction_202401_1',
  '/s/project/vale/driver_prediction_202401_2',
  '/s/project/vale/driver_prediction_202401_3',
  '/s/project/vale/driver_prediction_202401_4'
)

intogen_dirs <- c(
  '/s/project/mll/intogen/output_data/mll/leukemia_14groups/variants_wgs_filtered_0/output',
  '/s/project/mll/intogen/output_data/mll/leukemia_14groups/variants_wgs_filtered_1/output',
  '/s/project/mll/intogen/output_data/mll/leukemia_14groups/variants_wgs_filtered_2/output',
  '/s/project/mll/intogen/output_data/mll/leukemia_14groups/variants_wgs_filtered_3/output',
  '/s/project/mll/intogen/output_data/mll/leukemia_14groups/variants_wgs_filtered_4/output'
)

names(project_dirs) <- c(
  'driver_prediction_202401_0',
  'driver_prediction_202401_1',
  'driver_prediction_202401_2',
  'driver_prediction_202401_3',
  'driver_prediction_202401_4'
)


if(exists("exp_viz_var_filter")){remove(exp_viz_var_filter)}

for(i in seq_len(length(project_dirs))){
  experiment_design <- fread(paste0(project_dirs[i], '/experiment_design.tsv'))
  experiment_design[is.na(experiment_design)] <- ''
  
  exp_viz<- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                                model_method == 'rf' &
                                sample_group == 'leukemia_14group' &
                                intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                                outlier_input_feature == ''&
                                coess_input_feature == ''
  ]
  
  exp_viz <- add_result_paths(exp_viz, project_dirs[i], intogen_dirs[i], vep_dir)
  exp_viz[, intogen_input_dir := names(project_dirs)[i]]
  
  if(!exists("exp_viz_var_filter")){
    exp_viz_var_filter <- exp_viz
  }else{
    exp_viz_var_filter <- rbind(exp_viz_var_filter, exp_viz, fill=TRUE)
  }
}

fwrite(exp_viz_var_filter, snakemake@output$exp_viz_var_filter, sep='\t')

if(exists("res_viz_groups")){remove(res_viz_groups)}

for (i in 1:nrow(exp_viz_var_filter)) {
  print(i)
  res_temp <- fread(exp_viz_var_filter[i, res_post_path])
  res_temp[, Sample_group := exp_viz_var_filter[i, sample_group]]
  res_temp[, Prediction_post_scaled := Prediction_post/max(Prediction_post)]
  res_temp[, Method := exp_viz_var_filter[i, intogen_input_dir]]
  
  if(!exists("res_viz_groups")){
    res_viz_groups <- res_temp
  }else{
    res_viz_groups <- rbind(res_viz_groups, res_temp, fill=TRUE)
  }
}

res_vale_var_filter <- res_viz_groups[, .(GeneID, GeneSymbol, Prediction, Label, Sample_group, Method)]

fwrite(res_vale_var_filter, snakemake@output$mll_var_filter, sep='\t')
