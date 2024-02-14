#'---
#' title: sample_summary_tab
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
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - manuscriptWording: '`sm config["manuscriptWording"]`'
#'  input:
#'    - ods_unfitted: '`sm expand(config["outriderDir"] +
#'                    "/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds",
#'                    annotation=annotations, dataset=inputDatasets)`'
#'  output:
#'    - sample_summary_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/sample_summary_tab.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/sample_summary_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/sample_summary_tab.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})




# get parameters
single_group <- snakemake@params$inputDatasets
sample_group <- snakemake@params$outputDatasets

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_14group <- samp_anno[grep(single_group, DROP_GROUP),]
samp_anno_14group[Age <= 18, .N]
samp_anno_14group[Age <= 10, .N]

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))

# create table
samp_anno_14group_gender <- samp_anno_14group[, .(.N, sum(Gender==1), sum(Gender==2)), by='Diag'] 
setnames(samp_anno_14group_gender, c("N", "V2", "V3"), c("Number_of_individual", "Number_of_female", "Number_of_male") )

samp_anno_14group[, Age_group := paste0("Number_of_individual_", 
                                        ceiling(as.numeric(Age)/10-1)*10, "-", 
                                        ceiling(as.numeric(Age)/10)*10, "_yo"),
                  by=rownames(samp_anno_14group)]
samp_anno_14group_age <- samp_anno_14group[, .N, by=c('Diag', 'Age_group')] 
samp_anno_14group_age <- dcast(samp_anno_14group_age, Diag ~ Age_group)
samp_anno_14group_age[is.na(samp_anno_14group_age)] <- 0

samp_anno_sum <- merge(samp_anno_14group_gender, samp_anno_14group_age, by='Diag')
samp_anno_sum <- merge(samp_anno_sum, 
                             manuscript_wording[, .(Cohort_abbreviation, Cohort_during_analysis)],
                             by.x='Diag', by.y='Cohort_during_analysis')
samp_anno_sum[, Diag := NULL]
setnames(samp_anno_sum, 'Cohort_abbreviation', 'DiseaseEntity')
setcolorder(samp_anno_sum, c("DiseaseEntity", head(colnames(samp_anno_sum), -1)))
samp_anno_sum <- samp_anno_sum[order(-Number_of_individual),]

fwrite(samp_anno_sum, snakemake@output$sample_summary_tab)