#'---
#' title: generate outrider features
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
#'    - nr_topK_per_gene: '`sm expand(config["projectPath"] + 
#'                          "/processed_data/outlier_feature/fr-{dataset}-topK.tsv", dataset=datasets)`'
#'    - padjCutoffs_per_gene: '`sm expand(config["projectPath"] + 
#'                                    "/processed_data/outlier_feature/fr-{dataset}-padjCutoffs.tsv", dataset=datasets)`'
#'    - deltaJac_per_gene: '`sm expand(config["projectPath"] + 
#'                                  "/processed_data/outlier_feature/fr-{dataset}-deltaJac.tsv", dataset=datasets)`'
#'    - nr_topK_per_gene_single: '`sm expand(config["projectPath"] + 
#'                          "/processed_data/outlier_feature/fr-{dataset}-topK.tsv", dataset=singleDataset)`'
#'    - padjCutoffs_per_gene_single: '`sm expand(config["projectPath"] + 
#'                                    "/processed_data/outlier_feature/fr-{dataset}-padjCutoffs.tsv", dataset=singleDataset)`'
#'    - deltaJac_per_gene_single: '`sm expand(config["projectPath"] + 
#'                                  "/processed_data/outlier_feature/fr-{dataset}-deltaJac.tsv", dataset=singleDataset)`'
#'  type: script
#'  threads: 14
#'  resources:
#'    - mem_mb: 64000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/generate_fraser_features_groupwise.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/generate_fraser_features_groupwise.snakemake")




##### function #####
source("Scripts/preprocessing/outlier_feature/functions.R")

.libPaths("~/R/4.1/FRASER2")
# .libPaths()
# check with .libPaths() that indeed this directory shows up first
# devtools::install("../FraseR") # here, have the path to the folder in which the FRASER2 code from gitlab is




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
gencode <- fread(snakemake@params$gencode)
sampAnno <- fread(snakemake@params$sampAnno)
datasets <- snakemake@params$datasets
singleDataset <- snakemake@params$singleDataset
fds_paths <- snakemake@input$fds
nr_topK_per_gene_paths <- snakemake@output$nr_topK_per_gene
padjCutoffs_per_gene_paths <- snakemake@output$padjCutoffs_per_gene
deltaJac_per_gene_paths <- snakemake@output$deltaJac_per_gene


#+ write result
register(MulticoreParam(snakemake@threads))
message(date(), ": Running with ", bpworkers(), " workers ...")

k <- c(1, 5, 10, 25, 50)
cutoff_padj <- c(0.01, 0.05, 0.1)
cutoff_deltaJac <- c(0.1, 0.2, 0.3)
type <- 'jaccard'

#+ 14 groups
res_lists <- bplapply(seq_len(length(nr_topK_per_gene_paths)),
                      function(i){
                        message(date(), ": ", datasets[i])
                        
                        # read in 
                        fds <- loadFraserDataSet(file=fds_paths[i])
                        currentType(fds) <- type
                        gene_dt <- createGeneTableFds(fds, gencode)
                        
                        # make dir
                        proResDir <- dirname(nr_topK_per_gene_paths[i])
                        if(!dir.exists(proResDir)){
                          dir.create(proResDir)
                        }
                        message("Output dir: ", proResDir)
                        
                        # write
                        nr_topK_per_gene <- createTopKbyDeltaJac(fds, k, gene_dt, gencode, nr_topK_per_gene_paths[i])
                        padjCutoffs_per_gene <- createPadjCutoffbyDeltaJac(fds, cutoff_padj, gene_dt, gencode, padjCutoffs_per_gene_paths[i])
                        deltaJac_per_gene <- createDeltaJacSig(fds, cutoff_deltaJac, gene_dt, gencode, deltaJac_per_gene_paths[i])
                        
                        nr_topK_per_gene_melt <- melt(nr_topK_per_gene, id.vars = 'gene_id')
                        nr_topK_per_gene_melt[, variable:=gsub(datasets[i], singleDataset, nr_topK_per_gene_melt[, variable])]
                        
                        padjCutoffs_per_gene_melt <- melt(padjCutoffs_per_gene, id.vars = 'gene_id')
                        padjCutoffs_per_gene_melt[, variable:=gsub(datasets[i], singleDataset, padjCutoffs_per_gene_melt[, variable])]
                        
                        deltaJac_per_gene_melt <- melt(deltaJac_per_gene, id.vars = 'gene_id')
                        deltaJac_per_gene_melt[, variable:=gsub(datasets[i], singleDataset, deltaJac_per_gene_melt[, variable])]
                        
                        res_list <- list(nr_topK_per_gene_melt, padjCutoffs_per_gene_melt, deltaJac_per_gene_melt)
                        names(res_list) <- c("nr_topK_per_gene_melt", "padjCutoffs_per_gene_melt", "deltaJac_per_gene_melt")
                        return(res_list)
                      })
names(res_lists) <- datasets
print( paste0("Total number of groups: ", length(res_lists)))




#+ single group
for (i in datasets) {
  
  if (exists('nr_topK_per_gene_melt')) {
    nr_topK_per_gene_melt <- rbind(nr_topK_per_gene_melt, res_lists[[i]][['nr_topK_per_gene_melt']])
  }else{
    nr_topK_per_gene_melt <- res_lists[[i]][['nr_topK_per_gene_melt']]
  }
  
  if (exists('padjCutoffs_per_gene_melt')) {
    padjCutoffs_per_gene_melt <- rbind(padjCutoffs_per_gene_melt, res_lists[[i]][['padjCutoffs_per_gene_melt']])
  }else{
    padjCutoffs_per_gene_melt <- res_lists[[i]][['padjCutoffs_per_gene_melt']]
  }
  
  if (exists('deltaJac_per_gene_melt')) {
    deltaJac_per_gene_melt <- rbind(deltaJac_per_gene_melt, res_lists[[i]][['deltaJac_per_gene_melt']])
  }else{
    deltaJac_per_gene_melt <- res_lists[[i]][['deltaJac_per_gene_melt']]
  }
}

nr_topK_per_gene <- dcast(nr_topK_per_gene_melt, gene_id ~ variable, fun.aggregate=sum)
padjCutoffs_per_gene <- dcast(padjCutoffs_per_gene_melt, gene_id ~ variable, fun.aggregate=sum)
deltaJac_per_gene <- dcast(deltaJac_per_gene_melt, gene_id ~ variable, fun.aggregate=sum)

fwrite(nr_topK_per_gene, snakemake@output$nr_topK_per_gene_single, sep = "\t", quote = F)
fwrite(padjCutoffs_per_gene, snakemake@output$padjCutoffs_per_gene_single, sep = "\t", quote = F)
fwrite(deltaJac_per_gene, snakemake@output$deltaJac_per_gene_single, sep = "\t", quote = F)