#'---
#' title: fusion_prep
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
#'    - mll_manta: '`sm config["mll_manta"]`'
#'    - mll_arriba_dir: '`sm config["mll_arriba_dir"]`'
#'    - mll_star_fusion_dir: '`sm config["mll_star_fusion_dir"]`'
#'  output:
#'    - mll_arriba: '`sm expand(config["projectPath"] + "/manuscript/arriba/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_star_fusion: '`sm expand(config["projectPath"] + "/manuscript/star_fusion/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_manta: '`sm expand(config["projectPath"] + "/manuscript/manta/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'  type: script
#'  threads: 48
#'  resources:
#'    - mem_mb: 48000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/fusion_prep.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/fusion_prep.snakemake")
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
samp_anno[, star_fusion_file := paste0(snakemake@params$mll_star_fusion_dir, '/', ArrayID, '.star-fusion.tsv')]
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()

fusion_manta <- fread(snakemake@params$mll_manta)

# merge by cohort
sapply(sample_group, function(i){
  # read in arriba
  file_paths <- samp_anno_exp[DROP_GROUP==i, unique(arriba_file)]
  fusion_cohort_ls <- bplapply(file_paths, function(j){ 
    fusion_file <- fread(j) 
    fusion_file[, array_id := samp_anno[arriba_file==j, ArrayID]]
  }, BPPARAM = SerialParam())
  fusion_cohort_arriba <- rbindlist(fusion_cohort_ls)
  setnames(fusion_cohort_arriba, '#gene1', 'gene1')
  fusion_cohort_arriba[, fusion_id_1 := paste(array_id, gene1, breakpoint1, gene2, breakpoint2, sep='--')]
  fusion_cohort_arriba[, fusion_id_2 := paste(array_id, gene2, breakpoint2, gene1, breakpoint1, sep='--')]
  
  # read in star
  file_paths <- samp_anno_exp[DROP_GROUP==i, unique(star_fusion_file)]
  fusion_cohort_ls <- bplapply(file_paths, function(j){ 
    fusion_file <- fread(j) 
    fusion_file[, array_id := samp_anno[star_fusion_file==j, ArrayID]]
  }, BPPARAM = SerialParam())
  fusion_cohort_star <- rbindlist(fusion_cohort_ls)
  setnames(fusion_cohort_star, '#FusionName', 'FusionName')
  fusion_cohort_star <- separate(data=fusion_cohort_star, 
                                 col='FusionName', sep='--', into=c('gene1', "gene2")) %>% as.data.table()
  fusion_cohort_star <- separate(data=fusion_cohort_star, 
                                 col='LeftBreakpoint', sep=':', 
                                 into=c('LeftBreakpoint_chr', "LeftBreakpoint_pos", "LeftBreakpoint_strand"), 
                                 remove = FALSE) %>% as.data.table()
  fusion_cohort_star <- separate(data=fusion_cohort_star, 
                                 col='RightBreakpoint', sep=':', 
                                 into=c('RightBreakpoint_chr', "RightBreakpoint_pos", "RightBreakpoint_strand"), 
                                 remove = FALSE) %>% as.data.table()
  fusion_cohort_star[, LeftBreakpoint_chr := gsub("chr", "", LeftBreakpoint_chr)]
  fusion_cohort_star[, breakpoint1 := paste(LeftBreakpoint_chr, LeftBreakpoint_pos, sep=':')]
  fusion_cohort_star[, RightBreakpoint_chr := gsub("chr", "", RightBreakpoint_chr)]
  fusion_cohort_star[, breakpoint2 := paste(RightBreakpoint_chr, RightBreakpoint_pos, sep=':')]
  
  fusion_cohort_star[, fusion_id_1 := paste(array_id, gene1, breakpoint1, gene2, breakpoint2, sep='--')]
  fusion_cohort_star[, fusion_id_2 := paste(array_id, gene2, breakpoint2, gene1, breakpoint1, sep='--')]
  
  # read in manta
  fusion_cohort_manta <- fusion_manta[array_id %in% samp_anno_exp[DROP_GROUP==i, ArrayID], ]
  setnames(fusion_cohort_manta, c('Gene1', 'Gene2'), c('gene1', 'gene2'))
  fusion_cohort_manta[, Chr1 := gsub("chr", "", Chr1)]
  fusion_cohort_manta[, breakpoint1 := paste(Chr1, Pos1, sep=':')]
  fusion_cohort_manta[, Chr2 := gsub("chr", "", Chr2)]
  fusion_cohort_manta[, breakpoint2 := paste(Chr2, Pos2, sep=':')]
  
  fusion_cohort_manta[, fusion_id_1 := paste(array_id, gene1, breakpoint1, gene2, breakpoint2, sep='--')]
  fusion_cohort_manta[, fusion_id_2 := paste(array_id, gene2, breakpoint2, gene1, breakpoint1, sep='--')]
  
  
  # get fusion that confirmed by at least 2 tools 
  fusion_ids <- c(unique(c(fusion_cohort_arriba[, fusion_id_1], fusion_cohort_arriba[, fusion_id_2])),
                  unique(c(fusion_cohort_star[, fusion_id_1], fusion_cohort_star[, fusion_id_2])),
                  unique(c(fusion_cohort_manta[, fusion_id_1], fusion_cohort_manta[, fusion_id_2])))
  dup_fusion_ids <- fusion_ids[duplicated(fusion_ids)]
  
  fusion_cohort_arriba_dup <- fusion_cohort_arriba[fusion_id_1 %in% dup_fusion_ids | fusion_id_2 %in% dup_fusion_ids, ]
  fusion_cohort_star_dup <- fusion_cohort_star[fusion_id_1 %in% dup_fusion_ids | fusion_id_2 %in% dup_fusion_ids, ]
  fusion_cohort_manta_dup <- fusion_cohort_manta[fusion_id_1 %in% dup_fusion_ids | fusion_id_2 %in% dup_fusion_ids, ]
  
  fwrite(fusion_cohort_arriba_dup, snakemake@output$mll_arriba[sample_group==i], sep = '\t')
  fwrite(fusion_cohort_star_dup, snakemake@output$mll_star_fusion[sample_group==i], sep = '\t')
  fwrite(fusion_cohort_manta_dup, snakemake@output$mll_manta[sample_group==i], sep = '\t')
}) 

