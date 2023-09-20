#'---
#' title: generate activation features
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
#'    - nr_topK_per_gene: '`sm expand(config["projectPath"] + 
#'                          "/processed_data/outlier_feature/ac-{dataset}-topK.tsv", dataset=outputDatasets)`'
#'    - padjCutoffs_per_gene: '`sm expand(config["projectPath"] + 
#'                                    "/processed_data/outlier_feature/ac-{dataset}-padjCutoffs.tsv", dataset=outputDatasets)`'
#'    - zScore_per_gene: '`sm expand(config["projectPath"] + 
#'                                  "/processed_data/outlier_feature/ac-{dataset}-zScore.tsv", dataset=outputDatasets)`'
#'  type: script
#'  threads: 14
#'  resources:
#'    - mem_mb: 16000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/generate_activation_features_groupwise.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202303/processed_data/snakemake/generate_activation_features_groupwise.snakemake")
 



##### function #####
source("Scripts/preprocessing/outlier_feature/functions.R")




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
gencode <- fread(snakemake@params$gencode)
sampAnno <- fread(snakemake@params$sampAnno)
ods_input_path <- snakemake@input$ods
outputDatasets <- snakemake@params$outputDatasets
nr_topK_per_gene_paths <- snakemake@output$nr_topK_per_gene
padjCutoffs_per_gene_paths <- snakemake@output$padjCutoffs_per_gene
zScore_per_gene_paths <- snakemake@output$zScore_per_gene


#+ write result
register(MulticoreParam(snakemake@threads))
message(date(), ": Running with ", bpworkers(), " workers ...")

k <- c(1, 5, 10, 25, 50)
cutoff_padj <- c(0.01, 0.05, 0.1)
cutoff_zscore <- c(2, 4, 6)

ods_input <- readRDS(ods_input_path)

res_lists <- bplapply(seq_len(length(nr_topK_per_gene_paths)),
                      function(i){
                        
                        # read in 
                        group <- outputDatasets[i]
                        sampAnno_sep <- separate_rows(sampAnno, 'DROP_GROUP', sep = ",") %>% as.data.table()
                        group_sampID <- sampAnno_sep[DROP_GROUP == group, unique(ArrayID)]
                        ods <- ods_input[, colnames(ods_input) %in% group_sampID]
                        
                        gene_dt <- createGeneTable(ods, gencode)
                        rowData(ods)$geneID <- gene_dt[, geneID]
                        
                        # make dir
                        proResDir <- dirname(nr_topK_per_gene_paths[i])
                        if(!dir.exists(proResDir)){
                          dir.create(proResDir)
                        }
                        message("Output dir: ", proResDir)
                        
                        # write
                        nr_topK_per_gene <- createTopKbyZScore(ods, k, gene_dt, nr_topK_per_gene_paths[i])
                        padjCutoffs_per_gene <- createPadjCutoffbyZScore(ods, cutoff_padj, gene_dt, padjCutoffs_per_gene_paths[i])
                        zScore_per_gene <- createZscoreSig(ods, cutoff_zscore, gene_dt, zScore_per_gene_paths[i])
                        
                        res_list <- list(nr_topK_per_gene, padjCutoffs_per_gene, zScore_per_gene)
                        names(res_list) <- c("nr_topK_per_gene", "padjCutoffs_per_gene", "zScore_per_gene")
                        return(res_list)
                      })
names(res_lists) <- outputDatasets 
print( paste0("Total number of groups: ", length(res_lists)))


