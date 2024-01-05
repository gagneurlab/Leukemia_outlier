#'---
#' title: cnv_tet2
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
#'    - cnv_tet2: '`sm config["projectPath"] + "/manuscript/cnv/cnv_tet2.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/cnv_tet2.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/cnv_tet2.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(magrittr)
})

#+ read in 
sample_group <- snakemake@params$outputDatasets

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()

gencode <- fread(snakemake@params$gencode)
gencode[, seqnames := gsub("chr", "", seqnames)]

cnv <- fread(snakemake@params$mll_cnv)
cnv[CONTIG=='NULL', CONTIG := 'M']
cnv[, MEAN_COPY_RATIO := 2^MEAN_LOG2_COPY_RATIO]
cnv <- cnv[CONTIG==4, ]
non_overlap <- gencode[gene_name=='TET2', end] < cnv[, START] | gencode[gene_name=='TET2', start] > cnv[, END]
table(non_overlap)
cnv <- cnv[!non_overlap,]

# get symbol for cnv
overlap_symbol_ls <- sapply(1:nrow(cnv), function(x){
  gencode_temp <- copy(gencode)
  gencode_temp <- gencode_temp[gene_name=='TET2', ]
  non_overlap <- gencode_temp[, end] < cnv[x, START] | gencode_temp[, start] > cnv[x, END] 
  overlap_symbol <- gencode_temp[!non_overlap, paste(gene_name, collapse = ",")]
})

cnv[, SYMBOL := unlist(overlap_symbol_ls)]
cnv <- separate_rows(cnv, 'SYMBOL', sep=',') %>% as.data.table()

fwrite(cnv, snakemake@output$cnv_tet2)
