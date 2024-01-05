#'---
#' title: figure_3 preparation FRASER
#' author: xueqicao
#' wb:
#'  py:
#'   - |
#'    annotations = config["geneAnnotation"]
#'    inputDatasets = config["cohort"]["groups"]
#'    outputDatasets = config["cohort"]["groups"]
#'  params:
#'    - projectPath: '`sm config["projectPath"]`'
#'    - annotations: '`sm annotations`'
#'    - gencode: '`sm config["gencode"]`'
#'    - sampAnno: '`sm config["sampleAnnotation"]`'
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - CGC_cancer_gene_processed: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_processed.tsv"`'
#'  input:
#'    - fds: '`sm expand(config["fraserDir"] +
#'            "/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/fds-object.RDS",
#'             annotation=annotations, dataset=inputDatasets)`'
#'  output:
#'    - enrichmentTab: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_3/plot_data/enrichments_TopK_fr.csv"`'
#'  type: script
#'  threads: 84  
#'  resources:
#'    - mem_mb: 400000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_3_prep_fr.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_3_prep_fr.snakemake")
print("Snakemake saved") 




##### function #####
# define score function
createRankMeltFromFds <- function(fds, 
                                  gene_map,
                                  type = 'jaccard',
                                  sortBy = c('absDeltaPsi', 'padj'),
                                  signif_level=NULL){
  
  # extract value
  currentType(fds) <- type
  
  padj_dt <- padjVals(fds)
  absDeltaPsiVal_dt <- abs(deltaPsiValue(fds, type))
  
  # remove genes not found
  col_anno <- rowRanges(fds, type = type)
  
  padj_dt <- as.data.table(padj_dt)
  padj_dt[, hgnc_symbol := col_anno$hgnc_symbol]
  padj_dt <- padj_dt[!is.na(hgnc_symbol)] # remove na
  padj_dt <- padj_dt[hgnc_symbol %in% gene_map[, gene_name_orig]]
  
  absDeltaPsiVal_dt <- as.data.table(absDeltaPsiVal_dt)
  absDeltaPsiVal_dt[, hgnc_symbol := col_anno$hgnc_symbol]
  absDeltaPsiVal_dt <- absDeltaPsiVal_dt[!is.na(hgnc_symbol)] # remove na
  absDeltaPsiVal_dt <- absDeltaPsiVal_dt[hgnc_symbol %in% gene_map[, gene_name_orig]]
  
  # melt and filter by absDeltaPsi 0.1 and padj 0.05
  padj_melt_dt <- melt(padj_dt, id.vars = c("hgnc_symbol"))
  absDeltaPsiVal_melt_dt <- melt(absDeltaPsiVal_dt, id.vars = c("hgnc_symbol"))
  merge_dt <- cbind(padj_melt_dt, absDeltaPsiVal_melt_dt[, value])
  setnames(merge_dt, c("value", "V2"), c("padj", "absDeltaPsiVal"))
  merge_dt <- merge_dt[absDeltaPsiVal > 0.1 & padj < 0.05, ]
  
  
  rank <- switch(sortBy, 
                 "absDeltaPsi" = {
                   # dcast to gene-sample
                   absDeltaPsiVal_dcast_dt <- dcast(merge_dt[, .(hgnc_symbol, variable, absDeltaPsiVal)], hgnc_symbol ~ variable, fun.aggregate = max , fill = 0)
                   row_gene_symbol <- absDeltaPsiVal_dcast_dt[, hgnc_symbol]
                   absDeltaPsiVal_dcast_dt[, hgnc_symbol := NULL]
                   
                   # add missing column
                   missing_columns <- setdiff(colnames(fds), colnames(absDeltaPsiVal_dcast_dt))
                   absDeltaPsiVal_dcast_dt[, eval(missing_columns):= 0] # missing in missing sample with deltaPsi 0
                   
                   # convert to score
                   rc_score <- t(-absDeltaPsiVal_dcast_dt)
                   
                   # get the rank
                   rank <- apply(rc_score, 1,
                                 function(scores, ties.method){
                                   rank(scores, ties.method=ties.method)
                                 }, ties.method= "min")
                   
                   
                   # remove not meaningful
                   rank[absDeltaPsiVal_dcast_dt <= 0.3] <- NA
                   rank
                 },
                 "padj" = {
                   # dcast to gene-sample
                   padj_dcast_dt <- dcast(merge_dt[, .(hgnc_symbol, variable, padj)], hgnc_symbol ~ variable, fun.aggregate = min , fill = 1)
                   row_gene_symbol <- padj_dcast_dt[, hgnc_symbol]
                   padj_dcast_dt[, hgnc_symbol := NULL]
                   
                   # add missing column
                   missing_columns <- setdiff(colnames(fds), colnames(padj_dcast_dt))
                   padj_dcast_dt[, eval(missing_columns):= 1] # missing in missing sample with padj 1
                   
                   # convert to score
                   rc_score <- t(padj_dcast_dt)
                   
                   # get the rank 
                   rank <- apply(rc_score, 1, 
                                 function(scores, ties.method){
                                   rank(scores, ties.method=ties.method)
                                 }, ties.method= "min")
                   
                   
                   # remove not meaningful
                   rank[padj_dcast_dt >= 0.05] <- NA
                   rank
                 })
  
  
  # melt the rank
  rank <- as.data.table(rank)
  row_gene_id <- left_join(data.table(gene_name = row_gene_symbol),
                           gene_map[, .(gene_name, gene_id)], by = 'gene_name')$gene_id
  rank[, ENSGid := row_gene_id]
  
  rank_melt <- melt(rank, id.vars = "ENSGid")
  setnames(rank_melt, "value", "rank")
  rank_melt <- rank_melt[!is.na(rank), ]
  
  return(rank_melt)
}

# define enrichment function
computeEnrichmentFromRandomDriverList <- function(fds, type, rank, driver_genes, gene_map, n=99){
  set.seed(2020)
  
  # correct ENSGid
  fds_hgnc_symnol <- rowRanges(fds, type = type)$hgnc_symbol %>% unique()
  fds_hgnc_symnol <- fds_hgnc_symnol[!is.na(fds_hgnc_symnol)]
  fds_gene_id <- left_join(data.table(gene_name = fds_hgnc_symnol),
                           gene_map[, .(gene_name, gene_id)], by = 'gene_name')$gene_id
  ENSGid_short_fds <- sapply(fds_gene_id, function(x){
    strsplit(x, "[.]")[[1]][1]
  }) 
  
  driver_genes_short <- sapply(driver_genes, function(x){
    strsplit(x, "[.]")[[1]][1]
  })
  
  driver_genes_short_in_fds <- driver_genes_short[driver_genes_short %in% ENSGid_short_fds]
  ENSGid_short_rank <- sapply(rank[, ENSGid], function(x){
    strsplit(x, "[.]")[[1]][1]
  })
  rank[, ENSGid_short:= ENSGid_short_rank]
  
  # caculate drivers at K in cancer data
  k_case <- sapply(1:max(rank[, rank]), function(x){
    rank[rank==x, sum(ENSGid_short %in% driver_genes_short)]
  })
  
  # caculate drivers at K in random shuffle
  k_random <- sapply(1:n, function(i){
    driver_genes_random <- sample(ENSGid_short_fds, length(driver_genes_short_in_fds), replace = FALSE)
    k_random <- sapply(1:max(rank[, rank]), function(x){
      rank[rank==x, sum(ENSGid_short %in% driver_genes_random)]
    }) 
  })
  k_random_dt <- t(k_random) %>% as.data.table()
  # colnames(k_random_dt) <- paste0("k=", 1:max(rank[, rank]))
  
  # caculate enrichment
  rl <- sapply(1:max(rank[, rank]), function(i){
    sum(k_case[1:i]) 
  })
  rd <- sapply(1:max(rank[, rank]), function(i){
    mean(rowSums(k_random_dt[, 1:i]))
  })
  up <- sapply(1:max(rank[, rank]), function(i){
    quantile(rowSums(k_random_dt[, 1:i]), 0.95) 
  }) 
  dn <- sapply(1:max(rank[, rank]), function(i){
    quantile(rowSums(k_random_dt[, 1:i]), 0.05) 
  })   
  
  en_dt <- data.table(enrichment = rl/rd,
                      ci95 = up/rd,
                      ci05 = dn/rd,
                      real = rl,
                      random = rd,
                      random95 = up,
                      random05 = dn)
  en_dt[, k := 1:max(rank[, rank])]
  # en_dt <- melt(en_dt, id.vars = "k", variable.name = "term")
  
  return(en_dt)
}




##### process the data #####
.libPaths("~/R/4.1/FRASER2")
# check with .libPaths() that indeed this directory shows up first
# devtools::install("../FraseR") # here, have the path to the folder in which the FRASER2 code from gitlab is

suppressPackageStartupMessages({
  library(FRASER)
  library(data.table)
  library(dplyr)
  library(BiocParallel)
})

sessionInfo() # check that you got version 1.99


#+ get parameters
CGC_cancer_gene <- fread(snakemake@params$CGC_cancer_gene_processed)
samp_anno <- fread(snakemake@params$sampAnno)
gencode <- fread(snakemake@params$gencode)


#+ parallel setup
register(MulticoreParam(snakemake@threads))
message(date(), ": Running with ", bpworkers(), " workers ...")


#+ define paths
print(snakemake@input$fds)
sample_group <- snakemake@params$outputDatasets

fds_files <- snakemake@input$fds
names(fds_files) <- sample_group


#+ read in cancer gene lists
cancerDriverList <- list(all=CGC_cancer_gene,
                         leukemia=CGC_cancer_gene[grep("L", TissueType), ])

cancer_driver_lists <- lapply(cancerDriverList,
                              function(x){
                                tmp <- x 
                                tmp[,ensembleID:=ENSGid]
                                tmp[,ocg:=grepl("oncogene", RoleinCancer)]
                                tmp[,tsg:=grepl("TSG", RoleinCancer)]
                                tmp <- tmp[!is.na(ensembleID)]
                                res_list <- list(
                                  full=tmp[,ensembleID],
                                  tsg=tmp[tsg == TRUE,ensembleID],
                                  ocg=tmp[ocg == TRUE,ensembleID])
                                return(res_list)
                              })
cancer_driver_lists <- unlist(cancer_driver_lists, recursive=FALSE)
lapply(cancer_driver_lists, length)


#+ compute enrichment based on random shuffled driver gene list
type_to_plot <- "jaccard"
sortBy <- "padj"
n_rep <- 99 

combinations <- expand.grid(sample_group, names(cancer_driver_lists), type_to_plot)
message(date(), ": ", nrow(combinations), " combinations to plot")
enrichment_ls <- bplapply(seq_len(nrow(combinations)), 
                          function(i){ 
                            message(date(), ": ", combinations[i,'Var1'], "\t",
                                    combinations[i,'Var2'], "\t",
                                    combinations[i,'Var3'])
                            
                            # read in data
                            fds <- loadFraserDataSet(file = fds_files[combinations[i,'Var1']])
                            fds <- fds[rowData(fds)$hgnc_symbol %in% gencode[gene_type=='protein_coding', gene_name], ] 
                            
                            driver_genes <- cancer_driver_lists[[combinations[i,'Var2']]]
                            type <- as.character(combinations[i,'Var3'])
                            
                            rank <- createRankMeltFromFds(fds, 
                                                          gencode, 
                                                          type,
                                                          sortBy,
                                                          signif_level = signif_level)
                            
                            # calculate enrichment
                            en_dt <- computeEnrichmentFromRandomDriverList(fds, type, rank, driver_genes, gencode, n=n_rep) 
                            en_dt[, sample_subset := combinations[i,'Var1']]
                            en_dt[, cancer_driver_list := combinations[i,'Var2']]
                            en_dt[, type := type]
                            
                            return(en_dt)
                          })

#+ combine the list
enrichments <- rbindlist(enrichment_ls)
enrichments[, en95 := real/random95]
enrichments[, en05 := real/random05]


#+ save enrichments
fwrite(enrichments, snakemake@output$enrichmentTab)
