#'---
#' title: figure_4_prep_benchmark
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
#'    - mll_benchmark: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_benchmark.tsv"`'
#'  type: script
#'  resources:
#'    - mem_mb: 16000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_4_prep_benchmark.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/figure_4_prep_benchmark.snakemake")
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



#### VALE: 7+4 #### 
exp_viz<- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                              model_method == 'rf' &
                              intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                              outlier_input_feature == 'or,ac,absplice,fr' &
                              coess_input_feature == '', ]

exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)


if(exists("res_viz_groups")){remove(res_viz_groups)}

for (i in 1:nrow(exp_viz)) {
  res_temp <- fread(exp_viz[i, res_post_path])
  res_temp[, Sample_group := exp_viz[i, sample_group]]
  res_temp[, Prediction_post_scaled := Prediction_post/max(Prediction_post)]
  
  if(!exists("res_viz_groups")){
    res_viz_groups <- res_temp
  }else{
    res_viz_groups <- rbind(res_viz_groups, res_temp, fill=TRUE)
  }
}

res_ref <- copy(res_viz_groups)
res_ref[, SAMPLE_GROUP := toupper(Sample_group)]

res_vale <- res_viz_groups[, .(GeneID, GeneSymbol, Prediction, Label, Sample_group)]
res_vale[, Method := 'vale']


#### VALE: 7 #### 
exp_viz<- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                              model_method == 'rf' &
                              intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                              outlier_input_feature == '' &
                              coess_input_feature == '', ]

exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)


if(exists("res_viz_groups")){remove(res_viz_groups)}

for (i in 1:nrow(exp_viz)) {
  res_temp <- fread(exp_viz[i, res_post_path])
  res_temp[, Sample_group := exp_viz[i, sample_group]]
  res_temp[, Prediction_post_scaled := Prediction_post/max(Prediction_post)]
  
  if(!exists("res_viz_groups")){
    res_viz_groups <- res_temp
  }else{
    res_viz_groups <- rbind(res_viz_groups, res_temp, fill=TRUE)
  }
}

res_vale_7tools <- res_viz_groups[, .(GeneID, GeneSymbol, Prediction, Label, Sample_group)]
res_vale_7tools[, Method := 'vale_7tools']


#### mutsigcv ####
mutsigcv_files <- list.files(path = mutsigcv_dir,
                             pattern = 'MLL_*',
                             full.names = TRUE, recursive = TRUE)
mutsigcv_files <- grep('sig_genes', mutsigcv_files, value = TRUE)


if(exists("res_viz_groups")){remove(res_viz_groups)}

for (i in mutsigcv_files) {
  res_temp <- fread(i)
  res_temp[, Sample_group := strsplit(i, 'MLL_|[.]sig_genes')[[1]][4]]
  
  if(!exists("res_viz_groups")){
    res_viz_groups <- res_temp
  }else{
    res_viz_groups <- rbind(res_viz_groups, res_temp, fill=TRUE)
  }
}

res_mutsigcv <- copy(res_viz_groups)
setnames(res_mutsigcv, 'gene', 'GeneSymbol')
res_mutsigcv[, Prediction := 1-q]
res_mutsigcv <- merge(res_mutsigcv, res_ref[, .(GeneID, GeneSymbol, Label, Sample_group)], 
                      by = c('GeneSymbol', 'Sample_group'), all.x = FALSE, all.y = TRUE)
res_mutsigcv <- res_mutsigcv[, .(GeneID, GeneSymbol, Prediction, Label, Sample_group)]
res_mutsigcv[, Method := 'mutsigcv']
# TODO: add leukemia 14group correctly


#### intogen wo post ####
intogen_comb_files <- list.files(paste0(intogen_dir, '/steps/combination/'), full.names = TRUE)

if(exists("res_viz_groups")){remove(res_viz_groups)}

for (i in intogen_comb_files) {
  res_temp <- fread(i)
  res_temp[, SAMPLE_GROUP := strsplit(i, 'MLL_WGS_MLL_|[.]')[[1]][2]]
  
  if(!exists("res_viz_groups")){
    res_viz_groups <- res_temp
  }else{
    res_viz_groups <- rbind(res_viz_groups, res_temp, fill=TRUE)
  }
}

res_intogen_raw <- copy(res_viz_groups)
setnames(res_intogen_raw, 'SYMBOL', 'GeneSymbol')
res_intogen_raw[, Prediction := 1-QVALUE_stouffer_w]
res_intogen_raw <- merge(res_intogen_raw, res_ref[, .(GeneID, GeneSymbol, Label, SAMPLE_GROUP, Sample_group)], 
                         by = c('GeneSymbol', 'SAMPLE_GROUP'), all.x = FALSE, all.y = TRUE)
res_intogen_raw <- res_intogen_raw[, .(GeneID, GeneSymbol, Prediction, Label, Sample_group)]
res_intogen_raw[, Method := 'intogen_wo_post']
res_intogen_raw <- res_intogen_raw[Sample_group!='leukemia_14group', ]


#### intogen ####
intogen_driver_files <- list.files(paste0(intogen_dir, '/steps/drivers/'), full.names = TRUE)

if(exists("res_viz_groups")){remove(res_viz_groups)}

for (i in intogen_driver_files) {
  res_temp <- fread(i)
  res_temp[, SAMPLE_GROUP := strsplit(i, 'MLL_WGS_MLL_|[.]')[[1]][2]]
  
  if(!exists("res_viz_groups")){
    res_viz_groups <- res_temp
  }else{
    res_viz_groups <- rbind(res_viz_groups, res_temp, fill=TRUE)
  }
}

res_intogen <- copy(res_viz_groups)
setnames(res_intogen, 'SYMBOL', 'GeneSymbol')
res_intogen[, Prediction := 1-QVALUE_COMBINATION]
res_intogen <- merge(res_intogen, res_ref[, .(GeneID, GeneSymbol, Label, SAMPLE_GROUP, Sample_group)], 
                     by = c('GeneSymbol', 'SAMPLE_GROUP'), all.x = FALSE, all.y = TRUE)
res_intogen <- res_intogen[, .(GeneID, GeneSymbol, Prediction, Label, Sample_group)]
res_intogen[, Method := 'intogen']
res_intogen <- res_intogen[Sample_group!='leukemia_14group', ]


#### 7 tools ####
res_7tools_paths <- list.files(paste0(project_dir, '/processed_data/intogen_feature/score/'), full.names = TRUE)
res_7tools_ls <- lapply(res_7tools_paths, function(i){
  res_7tools_temp <- fread(i)
  setnames(res_7tools_temp, 'gene_name', 'GeneSymbol')
  setnames(res_7tools_temp, 'group', 'SAMPLE_GROUP')
  res_7tools_temp[, Prediction := 1-exp(-score)]
  res_7tools_temp <- merge(res_7tools_temp, 
                           res_ref[, .(GeneID, GeneSymbol, Label, SAMPLE_GROUP, Sample_group)], 
                           by = c('GeneSymbol', 'SAMPLE_GROUP'), all.x = FALSE, all.y = TRUE)
  res_7tools_temp <- res_7tools_temp[, .(GeneID, GeneSymbol, Prediction, Label, Sample_group)]
  res_7tools_temp[, Method := tail(strsplit(i, '/|[.]')[[1]], n=2)[1]]
  res_7tools_temp <- res_7tools_temp[order(-Prediction),]
  return(res_7tools_temp)
})
res_7tools <- rbindlist(res_7tools_ls)
res_7tools <- res_7tools[Sample_group!='leukemia_14group', ]

res_benchmark <- rbind(res_vale, res_vale_7tools, res_mutsigcv, res_intogen, res_intogen_raw, res_7tools)
res_benchmark <- res_benchmark[!is.na(GeneID), ]
res_benchmark[is.na(Prediction), Prediction := 0]
res_benchmark[, Label := GeneSymbol %in% res_ref[Label==TRUE, GeneSymbol]]

fwrite(res_benchmark, snakemake@output$mll_benchmark, sep='\t')
