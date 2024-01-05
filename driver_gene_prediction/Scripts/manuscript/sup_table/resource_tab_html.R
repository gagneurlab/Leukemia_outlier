#'---
#' title: resource_tab
#' author: xueqicao
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
#'    - resource_melt_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/resource_melt_tab.csv"`'         
#'  output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/resource_tab_html.html"`'
#'  type: noindex
#'  threads: 1  
#'  resources:
#'    - mem_mb: 32000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/resource_tab_html.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/resource_tab_html.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(DT)
})

resource_melt_tab <- fread(snakemake@input$resource_melt_tab)

datatable(resource_melt_tab)
