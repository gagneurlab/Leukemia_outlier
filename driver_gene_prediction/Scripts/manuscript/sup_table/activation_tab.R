#'---
#' title: generate activation_tab
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
#'  input:
#'    - ods: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods_filter_out.Rds",
#'             annotation=annotations, dataset=inputDatasets)`'
#'  output:
#'    - activation_tab: '`sm config["projectPath"] + 
#'                        "/manuscript/sup_table/activation_tab.csv"`'
#'  type: script
#'  threads: 24
#'  resources:
#'    - mem_mb: 28000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/activation_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/activation_tab.snakemake")




##### function #####
source("Scripts/preprocessing/outlier_feature/functions.R")

createPadjCutoffbyZScore <- function(ods, cutoffs, gene_dt, sample_group){
  padj <- assay(ods, 'padjust') 
  zscore <- zScore(ods) 
  
  padj_cutoffs_up <- lapply(cutoffs, nr_score_is_pvalCutoff_by_zScore, padj = padj, gene_dt = gene_dt, zScore=zscore, zScore_direction="up")
  padj_cutoffs_dn <- lapply(cutoffs, nr_score_is_pvalCutoff_by_zScore, padj = padj, gene_dt = gene_dt, zScore=zscore, zScore_direction="dn")
  
  padj_res <- rbindlist(padj_cutoffs_up)
  padj_res_wide <- dcast(padj_res, geneID ~ paste0("cutoff=", cutoff, "_up"),
                         value.var = "nrSignificant")
  padj_res_wide <- gene_dt[padj_res_wide, on="geneID"]
  padj_res_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  padj_res_wide_up <- copy(padj_res_wide)
  
  padj_res <- rbindlist(padj_cutoffs_dn)
  padj_res_wide <- dcast(padj_res, geneID ~ paste0("cutoff=", cutoff, "_dn"),
                         value.var = "nrSignificant")
  padj_res_wide <- gene_dt[padj_res_wide, on="geneID"]
  padj_res_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  padj_res_wide_dn <- copy(padj_res_wide)
  
  padj_res_wide <- merge(padj_res_wide_up, padj_res_wide_dn, by = c("geneSymbol", "geneID"))
  setnames(padj_res_wide, 'geneID', 'gene_id')
  padj_res_wide[, geneSymbol := NULL]
  
  colnames(padj_res_wide)[-1] <- paste(sample_group, "ac", colnames(padj_res_wide)[-1], sep="-")
  
  return(padj_res_wide)
}

createZscoreSig <- function(ods, cutoffs, gene_dt, sample_group){
  zScore_mtx <- zScore(ods)
  padj_mtx <- padj(ods)
  
  zScore_dt_list <- lapply(cutoffs, nr_score_is_zScore_by_padj, gene_dt=gene_dt, zScore_mtx=zScore_mtx, padj_mtx = padj_mtx)
  zScore_dt_long <- rbindlist(zScore_dt_list)
  zScore_dt <- dcast(zScore_dt_long, geneID+geneSymbol ~ variable, value.var = "value")
  
  setnames(zScore_dt, 'geneID', 'gene_id')
  zScore_dt[, geneSymbol := NULL]
  
  colnames(zScore_dt)[-1] <- paste(sample_group, "ac", colnames(zScore_dt)[-1], sep="-")
  
  return(zScore_dt)
}




##### process the data #####
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(OUTRIDER)
  library(BiocParallel)
})


#+ define paths
single_group <- snakemake@params$inputDatasets
gencode <- fread(snakemake@params$gencode)
sampAnno <- fread(snakemake@params$sampAnno)
ods_input_path <- snakemake@input$ods
activation_tab_path <- snakemake@output$activation_tab

sampAnno <- sampAnno[grep(single_group, DROP_GROUP)]
sampAnno[, Diag := gsub(" ", "_", Diag)]
sampAnno[, Diag := gsub("/", "_", Diag)]
sampAnno[, Diag := gsub("-", "_", Diag)]
diag_vec <- sampAnno[, unique(Diag)]

#+ write result
register(MulticoreParam(snakemake@threads))
message(date(), ": Running with ", bpworkers(), " workers ...")

cutoff_padj <- c(0.01, 0.05, 0.1)
cutoff_zscore <- c(2, 4, 6)

ods_input <- readRDS(ods_input_path)

padjCutoffs_list <- bplapply(diag_vec,
                             function(i){
                               
                               # read in 
                               diag_sampID <- sampAnno[Diag == i, ArrayID]
                               ods <- ods_input[, colnames(ods_input) %in% diag_sampID]
                               
                               gene_dt <- createGeneTable(ods, gencode)
                               rowData(ods)$geneID <- gene_dt[, geneID]
                               
                               # write
                               padjCutoffs_per_gene <- createPadjCutoffbyZScore(ods, cutoff_padj, gene_dt, i)
                               
                               return(padjCutoffs_per_gene)
                             })

padjCutoffs <- padjCutoffs_list[[1]]
for (num_diag in 2:length(diag_vec)) {
  padjCutoffs <- merge(padjCutoffs, padjCutoffs_list[[num_diag]], by='gene_id')
}


zScore_list <- bplapply(diag_vec,
                        function(i){
                          
                          # read in 
                          diag_sampID <- sampAnno[Diag == i, ArrayID]
                          ods <- ods_input[, colnames(ods_input) %in% diag_sampID]
                          
                          gene_dt <- createGeneTable(ods, gencode)
                          rowData(ods)$geneID <- gene_dt[, geneID]
                          
                          # write
                          zScore_per_gene <- createZscoreSig(ods, cutoff_zscore, gene_dt, i)
                          
                          return(zScore_per_gene)
                        })

zScore <- zScore_list[[1]]
for (num_diag in 2:length(diag_vec)) {
  zScore <- merge(zScore, zScore_list[[num_diag]], by='gene_id')
}

activation_tab <- merge(padjCutoffs, zScore, by='gene_id')
fwrite(activation_tab, activation_tab_path) 