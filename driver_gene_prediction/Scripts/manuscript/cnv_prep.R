#'---
#' title: cnv_prep
#' author: Xueqi Cao
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
#'    - mll_cnv: '`sm config["mll_cnv"]`'
#'  output:
#'    - mll_cnv: '`sm expand(config["projectPath"] + "/manuscript/cnv/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_cnv_full: '`sm expand(config["projectPath"] + "/manuscript/cnv_full/{dataset}.tsv",
#'                    dataset=outputDatasets)`'
#'  type: script
#'  resources:
#'    - mem_mb: 8000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/cnv_prep.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/cnv_prep.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(magrittr)
  # library(BiocParallel)
})

#+ parallel setup
# register(MulticoreParam(snakemake@threads))
# message(date(), ": Running with ", bpworkers(), " workers ...")

#+ read in 
sample_group <- snakemake@params$outputDatasets

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()

gencode <- fread(snakemake@params$gencode)
gencode[, seqnames := gsub("chr", "", seqnames)]



### cnv #### 
cnv <- fread(snakemake@params$mll_cnv)
cnv <- cnv[CALL!=0, ]
cnv[CONTIG=='NULL', CONTIG := 'M']
cnv[, MEAN_COPY_RATIO := 2^MEAN_LOG2_COPY_RATIO]

# get symbol for cnv
sapply(rev(sample_group), function(i){
  message(paste0(i, " start"))
  start_time <- Sys.time()
  
  cnv_cohort <- cnv[ARRAY_ID %in% samp_anno_exp[DROP_GROUP==i, ArrayID], ]
  
  overlap_symbol_ls <- sapply(1:nrow(cnv_cohort), function(x){
    gencode_temp <- copy(gencode)
    gencode_temp <- gencode_temp[seqnames==cnv_cohort[x, CONTIG], ]
    non_overlap <- gencode_temp[, end] < cnv_cohort[x, START] | gencode_temp[, start] > cnv_cohort[x, END] 
    overlap_symbol <- gencode_temp[!non_overlap, paste(gene_name, collapse = ",")]
  })

  cnv_cohort[, SYMBOL := unlist(overlap_symbol_ls)]
  cnv_cohort <- separate_rows(cnv_cohort, 'SYMBOL', sep=',') %>% as.data.table()
  
  fwrite(cnv_cohort, snakemake@output$mll_cnv[sample_group==i], sep = '\t')
  
  end_time <- Sys.time()
  message(end_time - start_time)
  message(paste0(i, " done"))
})


### cnv_full #### 
cnv <- fread(snakemake@params$mll_cnv)
cnv[CONTIG=='NULL', CONTIG := 'M']
cnv[, MEAN_COPY_RATIO := 2^MEAN_LOG2_COPY_RATIO]

# get symbol for cnv
sapply(rev(sample_group), function(i){
  message(paste0(i, " start"))
  start_time <- Sys.time()
  
  cnv_cohort <- cnv[ARRAY_ID %in% samp_anno_exp[DROP_GROUP==i, ArrayID], ]
  
  overlap_symbol_ls <- sapply(1:nrow(cnv_cohort), function(x){
    gencode_temp <- copy(gencode)
    gencode_temp <- gencode_temp[seqnames==cnv_cohort[x, CONTIG], ]
    non_overlap <- gencode_temp[, end] < cnv_cohort[x, START] | gencode_temp[, start] > cnv_cohort[x, END] 
    overlap_symbol <- gencode_temp[!non_overlap, paste(gene_name, collapse = ",")]
  })

  cnv_cohort[, SYMBOL := unlist(overlap_symbol_ls)]
  cnv_cohort <- separate_rows(cnv_cohort, 'SYMBOL', sep=',') %>% as.data.table()
  
  fwrite(cnv_cohort, snakemake@output$mll_cnv_full[sample_group==i], sep = '\t')
  
  end_time <- Sys.time()
  message(end_time - start_time)
  message(paste0(i, " done"))
})

