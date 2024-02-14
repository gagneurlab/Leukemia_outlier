#'---
#' title: fpkm_tab
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
#'    - fpkm_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/fpkm_tab.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/fpkm_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/fpkm_tab.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(rlist)
  library(OUTRIDER)
})




# get parameters
single_group <- snakemake@params$inputDatasets
sample_group <- snakemake@params$outputDatasets

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_14group <- samp_anno[grep(single_group, DROP_GROUP),]
diag_list <- samp_anno_14group[, table(Diag)] %>% sort(decreasing = TRUE) %>% names()

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]

ods_unfitted <- readRDS(snakemake@input$ods_unfitted)

fpkm_mtx <- rbind(fpkm(ods_unfitted))

fpkm_mean_ls <- lapply(diag_list, function(i){
  group_sample_id <- samp_anno_14group[Diag==i, ArrayID]
  fpkm_mtx_entity <- fpkm_mtx[, group_sample_id]
  fpkm_mean_group <- data.table(rowMeans(fpkm_mtx_entity))
  colnames(fpkm_mean_group) <- manuscript_wording[Cohort_during_analysis==i, Cohort_abbreviation]
  return(fpkm_mean_group)
})

fpkm_mean_dt <- list.cbind(fpkm_mean_ls)
fpkm_mean_dt[, gene_id_unique := rownames(ods_unfitted)]
fpkm_mean_dt <- merge(fpkm_mean_dt,
                      gencode[, .(gene_id_unique, gene_id, gene_name)], by='gene_id_unique', all.x=TRUE, all.y=FALSE)
setnames(fpkm_mean_dt, c("gene_id", "gene_name"), c('GeneID', 'GeneSymbol'))
fpkm_mean_dt[, gene_id_unique := NULL]

setcolorder(fpkm_mean_dt, c("GeneID", "GeneSymbol", head(colnames(fpkm_mean_dt), -2)))

fwrite(fpkm_mean_dt, snakemake@output$fpkm_tab)
