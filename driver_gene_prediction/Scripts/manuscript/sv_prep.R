#'---
#' title: sv_prep
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
#'    - mll_manta_sv: '`sm config["mll_manta_sv"]`'
#'  output:
#'    - mll_manta_sv: '`sm expand(config["projectPath"] + "/manuscript/manta_sv/{dataset}.tsv",
#'                      dataset=outputDatasets)`'
#'  type: script
#'  threads: 48
#'  resources:
#'    - mem_mb: 48000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/sv_prep.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/sv_prep.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(magrittr)
  library(BiocParallel)
})

#+ parallel setup
register(MulticoreParam(snakemake@threads))
message(date(), ": Running with ", bpworkers(), " workers ...")




# get parameters
set.seed(2023)

sample_group <- snakemake@params$outputDatasets

gencode <- fread(snakemake@params$gencode)

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno[, arriba_file := paste0(snakemake@params$mll_arriba_dir, '/', ArrayID, '.arriba.tsv')]
samp_anno[, star_sv_file := paste0(snakemake@params$mll_star_sv_dir, '/', ArrayID, '.star-fusion.tsv')]
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()

sv_manta <- fread(snakemake@params$mll_manta_sv)

# merge by cohort
sapply(sample_group, function(i){

  sv_manta <- sv_manta[from_gene!='' & to_gene!='', ]
  sv_cohort_manta <- sv_manta[array_id %in% samp_anno_exp[DROP_GROUP==i, ArrayID], ]
  sv_cohort_manta <- separate(data=sv_cohort_manta, remove=FALSE,
                                  col='from_gene', sep="[(]|[)]", into=c('gene1', "gene1_strand")) %>% as.data.table()
  sv_cohort_manta <- separate(data=sv_cohort_manta, remove=FALSE,
                                  col='to_gene', sep="[(]|[)]", into=c('gene2', "gene2_strand")) %>% as.data.table()
  
  fwrite(sv_cohort_manta, snakemake@output$mll_manta_sv[sample_group==i], sep = '\t')
})

