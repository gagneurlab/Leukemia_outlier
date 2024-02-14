#'---
#' title: mutation_prep
#' author: Xueqi Cao
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
#'    - predictedConsequence: '`sm config["predictedConsequence"]`'
#'    - AML_variants: '`sm config["AML_variants"]`'
#'    - AML_enformer: '`sm config["AML_enformer"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'  input:
#'    - mll_cnv: '`sm expand(config["projectPath"] + "/manuscript/cnv/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_cnv_full: '`sm expand(config["projectPath"] + "/manuscript/cnv_full/{dataset}.tsv",
#'                    dataset=outputDatasets)`'
#'    - ods: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds",
#'             annotation=annotations, dataset=inputDatasets)`'
#'    - ods_filter_out: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods_filter_out.Rds",
#'                        annotation=annotations, dataset=inputDatasets)`'
#'    - fraserRes: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv",
#'                    annotation=annotations, dataset=outputDatasets)`'
#'    - fraserResJunc: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv",
#'                    annotation=annotations, dataset=outputDatasets)`'
#'    - fds: '`sm expand(config["fraserDir"] +
#'            "/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/fds-object.RDS",
#'             annotation=annotations, dataset=outputDatasets)`'
#'    - abspliceRes: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice.tsv"`'
#'    - abspliceResVar: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice_var.tsv"`'
#'    - mutation_fr: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/mutation_fr.csv"`'
#'  output:
#'    - fr_mut_en: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/fr_mut_en.csv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/mutation_prep.html"`'
#'  type: noindex
#'  threads: 1  
#'  resources:
#'    - mem_mb: 256000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/mutation_prep.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/mutation_prep.snakemake")

.libPaths("~/R/4.1/FRASER2")
suppressPackageStartupMessages({
  library(OUTRIDER)
  library(FRASER)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]
gencode <- gencode[gene_type == 'protein_coding', ]

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()
samp_anno_14groups <- samp_anno_exp[DROP_GROUP %in% snakemake@params$inputDatasets, ] %>% unique()
single_group_id <- samp_anno_14groups[, ArrayID]
  
mutation_fr <- fread(snakemake@input$mutation_fr)

mutation_type_list <- c("promoter_variant", "frameshift_variant", "stop_gained",
                        "splice_related_variant", "structural_variant", "absplice_variant",
                        "copy_number_gain", "copy_number_loss")           

fr_all_mut_en_list <- lapply(mutation_type_list, function(x){
  x_val <- mutation_fr[, ..x]
  x_tab <- table(data.table(x = x_val,
                            aberrant = mutation_fr[, aberrant]))
  x_fisher <- fisher.test(x_tab)
  
  fr_all_mut_en <- data.table(
    mutation_type = x,
    odds_ratio = x_fisher$estimate,
    odds_ratio_ci_low = x_fisher$conf.int[1],
    odds_ratio_ci_up = x_fisher$conf.int[2],
    p_val = x_fisher$p.value,
    total = x_tab[2,2]
  )
  return(fr_all_mut_en)
})

fr_all_mut_en <- rbindlist(fr_all_mut_en_list)
fr_all_mut_en[, cutoff :='All']

fr_top3_mut_en_list <- lapply(mutation_type_list, function(x){
  x_val <- mutation_fr[, ..x]
  x_tab <- table(data.table(x = x_val,
                            aberrant = mutation_fr[, top3]))
  x_fisher <- fisher.test(x_tab)
  
  fr_top3_mut_en <- data.table(
    mutation_type = x,
    odds_ratio = x_fisher$estimate,
    odds_ratio_ci_low = x_fisher$conf.int[1],
    odds_ratio_ci_up = x_fisher$conf.int[2],
    p_val = x_fisher$p.value,
    total = x_tab[2,2]
  )
  return(fr_top3_mut_en)
})

fr_top3_mut_en <- rbindlist(fr_top3_mut_en_list)
fr_top3_mut_en[, cutoff :='At most three']

fr_mut_en <- rbind(fr_all_mut_en, fr_top3_mut_en)
fwrite(fr_mut_en, snakemake@output$fr_mut_en)