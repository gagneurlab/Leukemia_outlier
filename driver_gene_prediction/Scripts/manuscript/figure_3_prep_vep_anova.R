#'---
#' title: figure_3_prep_vep
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
#'    - mll_panel_genes: '`sm config["mll_panel_genes"]`'
#'    - AML_variants: '`sm config["AML_variants"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'    - CGC_leukemia_OCG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_OCG_list.tsv"`'
#'    - CGC_leukemia_TSG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_TSG_list.tsv"`'
#'    - anonymization_table_mll: '`sm config["anonymization_table_mll"]`'
#'    - abspliceDirNoHeader: '`sm config["abspliceDirNoHeader"]`'
#'  input:
#'    - vep_full: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_full.tsv"`'
#'  output:
#'    - vep_high: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_high.tsv"`'
#'  type: script
#'  resources:
#'    - mem_mb: 256000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_3.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_3.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})




# get parameters
set.seed(2023)

gencode <- fread(snakemake@params$gencode)

# read in vep 
vep_res <- fread(snakemake@input$vep_full)
vep_high <- vep_res[IMPACT=='HIGH']

fwrite(vep_high, snakemake@output$vep_high)

