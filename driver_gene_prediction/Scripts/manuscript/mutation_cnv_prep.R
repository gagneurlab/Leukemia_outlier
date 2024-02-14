#'---
#' title: mutation_cnv_prep
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
#'    - predictedConsequence: '`sm config["predictedConsequence"]`'
#'    - AML_variants: '`sm config["AML_variants"]`'
#'    - AML_enformer: '`sm config["AML_enformer"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'  input:
#'    - mutation_gene_sample: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/mutation_gene_sample.csv"`'
#'    - mll_cnv: '`sm expand(config["projectPath"] + "/manuscript/cnv/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_cnv_full: '`sm expand(config["projectPath"] + "/manuscript/cnv_full/{dataset}.tsv",
#'                    dataset=outputDatasets)`'
#'  output:
#'    - mutation_cnv_gene_sample: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/mutation_cnv_gene_sample.csv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/mutation_cnv_prep.html"`'
#'  type: noindex
#'  threads: 1  
#'  resources:
#'    - mem_mb: 256000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/mutation_cnv_prep.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/mutation_cnv_prep.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})




### cnv ####
mutation_dt <- fread(snakemake@input$mutation_gene_sample)

for (x in snakemake@input$mll_cnv_full) {
  cnv_res_temp <- fread(x)
  cnv_res_temp[, LENGTH := END - START]
  
  if (exists("cnv_res")) {
    cnv_res <- rbind(cnv_res, cnv_res_temp, fill=TRUE)
  } else {
    cnv_res <- cnv_res_temp
  }
}

mutation_cnv <- merge(mutation_dt, cnv_res[, .(SYMBOL, ARRAY_ID, MEAN_LOG2_COPY_RATIO, MEAN_COPY_RATIO)], 
                  by.x=c('gene_name_orig', 'sampleID'),
                  by.y=c('SYMBOL', 'ARRAY_ID'),
                  all.x=TRUE, all.y=FALSE)
# mutation_cnv[, table(is.na(MEAN_COPY_RATIO))]
mutation_cnv[is.na(MEAN_COPY_RATIO), MEAN_COPY_RATIO := 1]
mutation_cnv[is.na(MEAN_LOG2_COPY_RATIO), MEAN_LOG2_COPY_RATIO := 0]
# mutation_cnv[, table(is.na(MEAN_COPY_RATIO))]
fwrite(mutation_cnv, snakemake@output$mutation_cnv_gene_sample)
