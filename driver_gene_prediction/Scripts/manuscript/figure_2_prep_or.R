#'---
#' title: figure_2 preparation OUTRIDER
#' author: Xueqi Cao, Ines Scheller
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
#'    - CGC_cancer_gene_processed: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_processed.tsv"`'
#'  input:
#'    - ods: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds",
#'             annotation=annotations, dataset=inputDatasets)`'
#'  output:
#'    - enrichmentTab: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_2/plot_data/enrichments_TopK_or.csv"`'
#'    - pvalTab: '`sm config["projectPath"] + 
#'                "/manuscript/figure_2/plot_data/pvals_TopK_or.csv"`'
#'  type: script
#'  threads: 90  
#'  resources:
#'    - mem_mb: 20000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_2_prep_or.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202207/processed_data/snakemake/figure_2_prep_or.snakemake")
print("Snakemake saved") 




##### function #####
# define score function
createRankMeltFromOds <- function(ods, selected_assay, 
                                  signif_level=NULL, 
                                  exp = c("all", "up", "down")){
  pvalues <- assay(ods, selected_assay)
  zSr <- zScore(ods)
  
  # make zScore not fulfilled as NA
  rc_score <- switch(exp, 
                     "all" = t(pvalues), 
                     "up" = {
                       pvalues[zSr < 0] <- NA
                       t(pvalues)
                     }, 
                     "down" = {
                       pvalues[zSr > 0] <- NA
                       t(pvalues)
                     })
  
  # get the rank 
  rank <- apply(rc_score, 1, 
                function(scores, ties.method){
                  rank(scores, ties.method=ties.method)
                }, ties.method= "min")
  
  # make the not fulfilled rank as NA
  rank <- switch(exp, 
                 "all" = rank, 
                 "up" = {
                   rank[zSr < 0] <- NA
                   rank
                 }, 
                 "down" = {
                   rank[zSr > 0] <- NA
                   rank
                 })
  
  if(!is.null(signif_level)){
    rank[pvalues > signif_level] <- NA
  }
  
  # melt the rank
  rank <- as.data.table(rank)
  rank[, ENSGid := rownames(rowData(ods))]
  
  rank_melt <- melt(rank, id.vars = "ENSGid")
  setnames(rank_melt, "value", "rank")
  rank_melt <- rank_melt[!is.na(rank), ]
  
  return(rank_melt)
}

# define enrichment function
computeEnrichmentFromRandomDriverList <- function(ods, rank, driver_genes, n=99){
  set.seed(2020)
  
  # correct ENSGid
  ENSGid_short_ods <- sapply(rownames(ods), function(x){
    strsplit(x, "[.]")[[1]][1]
  }) 
  
  driver_genes_short <- sapply(driver_genes, function(x){
    strsplit(x, "[.]")[[1]][1]
  })
  
  driver_genes_short_in_ods <- driver_genes_short[driver_genes_short %in% ENSGid_short_ods]
  
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
    driver_genes_random <- sample(ENSGid_short_ods, length(driver_genes_short_in_ods), replace = FALSE)
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
  
  # pval by permutation
  pval_dt <- sapply(1:max(rank[, rank]), function(i){
    occurance <- sum(sum(k_case[1:i]) < rowSums(k_random_dt[, 1:i]) ) + 1
    total <- n + 1
    pval <- occurance / total
  }) %>% as.data.table()
  colnames(pval_dt) <- "pval"
  pval_dt[, k := 1:max(rank[, rank])]
  
  return(list(en_dt, pval_dt))
}




##### process the data #####
suppressPackageStartupMessages({
  library(OUTRIDER)
  library(data.table)
  library(dplyr)
  library(BiocParallel)
})


#+ get parameters
CGC_cancer_gene <- fread(snakemake@params$CGC_cancer_gene_processed)
samp_anno <- fread(snakemake@params$sampAnno)


#+ parallel setup
register(MulticoreParam(snakemake@threads))
message(date(), ": Running with ", bpworkers(), " workers ...")


#+ define paths
print(snakemake@input$ods)
ods_all <- readRDS(snakemake@input$ods)
single_group <- snakemake@params$inputDatasets
sample_group <- snakemake@params$outputDatasets
  

#+ read in ods
odses <- lapply(sample_group, function(i){
  ods <- ods_all[, colnames(ods_all) %in% samp_anno[Cohort_group == i, ArrayID]]
})
names(odses) <- sample_group
odses <- append(odses, ods_all)
names(odses)[length(odses)] <- single_group
lapply(odses, dim)


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
sample_subsets <- names(odses)
signif_level <- 0.05
assay_to_plot <- "padjust"
n_rep <- 99

combinations <- expand.grid(sample_subsets, names(cancer_driver_lists))
enrichment_ls <- bplapply(seq_len(nrow(combinations)), 
                          function(i){ 
                            message(date(), ": ", combinations[i,'Var1'], "\t",
                                    combinations[i,'Var2'])
                            
                            # read in data
                            ods <- odses[[combinations[i,'Var1']]] 
                            driver_genes <- cancer_driver_lists[[combinations[i,'Var2']]]
                            exp_direction <- combinations[i,'Var2'] %>% as.character() %>% strsplit(., "[.]") %>% unlist() %>% tail(., 1)
                            
                            rank <- switch(exp_direction,
                                           "full" = createRankMeltFromOds(ods, assay_to_plot, signif_level, exp = "all"),
                                           "ocg" = createRankMeltFromOds(ods, assay_to_plot, signif_level, exp = "up"),
                                           "tsg" = createRankMeltFromOds(ods, assay_to_plot, signif_level, exp = "down")
                            )
                            
                            # calculate enrichment
                            res_ls <- computeEnrichmentFromRandomDriverList(ods, rank, driver_genes, n=n_rep) 
                            en_dt <- res_ls[[1]]
                            en_dt[, sample_subset := combinations[i,'Var1']]
                            en_dt[, cancer_driver_list := combinations[i,'Var2']]
                            
                            pval_dt <- res_ls[[2]]
                            pval_dt[, sample_subset := combinations[i,'Var1']]
                            pval_dt[, cancer_driver_list := combinations[i,'Var2']]
                            
                            return(list(en_dt, pval_dt))
                          }, BPPARAM = SerialParam())

#+ combine the list
enrichments <- rbindlist(lapply(enrichment_ls, "[[", 1))
enrichments[, en95 := real/random95]
enrichments[, en05 := real/random05]

pvals <- rbindlist(lapply(enrichment_ls, "[[", 2))


#+ save enrichments
fwrite(enrichments, snakemake@output$enrichmentTab)


#+ save pvals
fwrite(pvals, snakemake@output$pvalTab)

