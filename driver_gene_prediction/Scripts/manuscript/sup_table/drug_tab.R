#'---
#' title: prediction_tab
#' author: Xueqi Cao
#' wb:
#'  py:
#'   - |
#'    inputDatasets = config["cohort"]["single_group"]
#'  params:
#'    - projectPath: '`sm config["projectPath"]`'
#'    - experimentDesign: '`sm config["experimentDesign"]`'
#'    - intogenDir: '`sm config["intogenDir"]`'
#'    - vep_path: '`sm config["vep_path"]`'   
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'    - inputDatasets: '`sm inputDatasets`'
#'    - drug_target: '`sm config["drug_target"]`'
#'  output:
#'    - drug_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/drug_tab.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/drug_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/drug_tab.snakemake")
print("Snakemake saved") 



suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(plyr)
})
source("Scripts/manuscript/function.R")

intogen_dir <- snakemake@params$intogenDir
vep_dir <- snakemake@params$vep_path
project_dir <- snakemake@params$projectPath
experiment_design <- fread(snakemake@params$experimentDesign)
experiment_design[is.na(experiment_design)] <- ''
single_group <- snakemake@params$inputDatasets
drug_target <- as.data.table(fread(snakemake@params$drug_target))

#### Formatting drug table ####
drug_target <- drug_target %>%
  rename(
    c(Status = "Manual_curation", `Approved Drugs` = "Approved_drugs", `MoA of the approved drugs` = "MoA_of_the_approved_drugs", Gene = "GeneSymbol"),
  )

#### all sample #### 
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == single_group &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == '', ]

exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)

# vale
res_leu14 <- fread(exp_viz[, res_post_path])
res_leu14[, StudyGroup := 'Pan-leukemia']
res_leu14 <- res_leu14[, .(GeneID, GeneSymbol, StudyGroup, Label, Rank, Prediction, isCGC, RoleCGC, isIntOGen, RoleIntogen)]

output <- join(drug_target, res_leu14[, .(GeneID, GeneSymbol, Prediction)])
output <- output[,c(7, 1, 8, seq(2, 6))]

fwrite(output, snakemake@output$drug_tab) 
