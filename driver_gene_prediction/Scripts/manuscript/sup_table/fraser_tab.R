#'---
#' title: generate fraser_tab
#' author: xueqicao
#' wb:
#'  py:
#'   - |
#'    annotations = config["geneAnnotation"]
#'    datasets = config["cohort"]["groups"]
#'    singleDataset = config["cohort"]["single_group"]
#'  params:
#'    - projectPath: '`sm config["projectPath"]`'
#'    - annotations: '`sm annotations`'
#'    - gencode: '`sm config["gencode"]`'
#'    - sampAnno: '`sm config["sampleAnnotation"]`'
#'    - datasets: '`sm datasets`'
#'    - singleDataset: '`sm singleDataset`'
#'  input:
#'    - fds: '`sm expand(config["fraserDir"] +"/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/" + 
#'             "fds-object.RDS", annotation=annotations, dataset=datasets)`'
#'  output:
#'    - fraser_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/fraser_tab.csv"`'
#'  type: script
#'  threads: 14
#'  resources:
#'    - mem_mb: 48000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/fraser_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/fraser_tab.snakemake")




##### function #####
source("Scripts/preprocessing/outlier_feature/functions.R")

.libPaths("~/R/4.1/FRASER2")
# .libPaths()
# check with .libPaths() that indeed this directory shows up first
# devtools::install("../FraseR") # here, have the path to the folder in which the FRASER2 code from gitlab is

createPadjCutoffbyDeltaJac <- function(fds, cutoffs, gene_dt, gencode, sample_group){
  padj <- get_padj_gene(fds, gencode)
  deltaJac <- get_deltaJac_gene(fds, gencode)
  
  padj_cutoffs_up <- lapply(cutoffs, nr_score_is_pvalCutoff_by_zScore, padj = padj, gene_dt = gene_dt, zScore=deltaJac, zScore_direction="up")
  padj_cutoffs_dn <- lapply(cutoffs, nr_score_is_pvalCutoff_by_zScore, padj = padj, gene_dt = gene_dt, zScore=deltaJac, zScore_direction="dn")
  
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
  
  colnames(padj_res_wide)[-1] <- paste(sample_group, "fr", colnames(padj_res_wide)[-1], sep="-")
  
  return(padj_res_wide)
}

createDeltaJacSig <- function(fds, cutoffs, gene_dt, gencode, sample_group){
  padj_mtx <- get_padj_gene(fds, gencode)
  deltaJac_mtx <- get_deltaJac_gene(fds, gencode)
  
  zScore_dt_list <- lapply(cutoffs, nr_score_is_zScore_by_padj, gene_dt=gene_dt, zScore_mtx=deltaJac_mtx, padj_mtx = padj_mtx)
  zScore_dt_long <- rbindlist(zScore_dt_list)
  zScore_dt <- dcast(zScore_dt_long, geneID+geneSymbol ~ variable, value.var = "value")
  
  setnames(zScore_dt, 'geneID', 'gene_id')
  zScore_dt[, geneSymbol := NULL]
  
  colnames(zScore_dt)[-1] <- paste(sample_group, "fr", colnames(zScore_dt)[-1], sep="-")
  
  return(zScore_dt)
}




##### process the data #####
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(FRASER)
  library(BiocParallel)
})


#+ define paths
single_group <- snakemake@params$singleDataset
sample_groups <- snakemake@params$datasets
gencode <- fread(snakemake@params$gencode)
sampAnno <- fread(snakemake@params$sampAnno)
fds_paths <- snakemake@input$fds
names(fds_paths) <- sample_groups
fraser_tab_path <- snakemake@output$fraser_tab

sampAnno <- sampAnno[grep(single_group, DROP_GROUP)]
sampAnno[, Diag := gsub(" ", "_", Diag)]
sampAnno[, Diag := gsub("/", "_", Diag)]
sampAnno[, Diag := gsub("-", "_", Diag)]
diag_vec <- sampAnno[, unique(Diag)]

#+ write result
register(MulticoreParam(snakemake@threads))
message(date(), ": Running with ", bpworkers(), " workers ...")

cutoff_padj <- c(0.01, 0.05, 0.1)
cutoff_deltaJac <- c(0.1, 0.2, 0.3)
type <- 'jaccard'

fds_list <- bplapply(sample_groups,
                     function(i){
                       fds <- loadFraserDataSet(file=fds_paths[i])
                     })
names(fds_list) <- sample_groups

padjCutoffs_list <- vector(mode = "list", length = 1)

for (i in diag_vec) {
  message(i)
  sample_group <- sampAnno[Diag==i, unique(Cohort_group)]
  diag_sampID <- sampAnno[Diag == i, ArrayID]
  
  # read in 
  fds <- fds_list[[sample_group]]
  fds <- fds[, colnames(fds) %in% diag_sampID]
  currentType(fds) <- type
  gene_dt <- createGeneTableFds(fds, gencode)
  
  # write
  padjCutoffs_per_gene <- createPadjCutoffbyDeltaJac(fds, cutoff_padj, gene_dt, gencode, i)
  padjCutoffs_list[[i]] <- padjCutoffs_per_gene
}

padjCutoffs_list <- padjCutoffs_list[-1]
padjCutoffs <- padjCutoffs_list[[diag_vec[1]]]
for (j in diag_vec[-1]) {
  print(j)
  padjCutoffs <- merge(padjCutoffs, padjCutoffs_list[[j]], by='gene_id', all.x=TRUE, all.y=TRUE)
}

deltaJac_list <- vector(mode = "list", length = 1)

for (i in diag_vec) {
  message(i)
  sample_group <- sampAnno[Diag==i, unique(Cohort_group)]
  diag_sampID <- sampAnno[Diag == i, ArrayID]
  
  # read in 
  fds <- fds_list[[sample_group]]
  fds <- fds[, colnames(fds) %in% diag_sampID]
  currentType(fds) <- type
  gene_dt <- createGeneTableFds(fds, gencode)
  
  # write
  deltaJac_per_gene <- createDeltaJacSig(fds, cutoff_deltaJac, gene_dt, gencode, i)
  deltaJac_list[[i]] <- deltaJac_per_gene
}

deltaJac_list <- deltaJac_list[-1]
deltaJac <- deltaJac_list[[diag_vec[1]]]
for (j in diag_vec[-1]) {
  print(j)
  deltaJac <- merge(deltaJac, deltaJac_list[[j]], by='gene_id', all.x=TRUE, all.y=TRUE)
}

fraser_tab <- merge(padjCutoffs, deltaJac, by='gene_id', all.x=TRUE, all.y=TRUE)
fwrite(fraser_tab, fraser_tab_path) 