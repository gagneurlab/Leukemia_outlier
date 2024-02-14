#'---
#' title: figure_2_or_anova_gene
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
#'    - vep_high: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_high.tsv"`'
#'    - vep_splice: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_splice.tsv"`'
#'    - mll_cnv: '`sm expand(config["projectPath"] + "/manuscript/cnv/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_cnv_full: '`sm expand(config["projectPath"] + "/manuscript/cnv_full/{dataset}.tsv",
#'                    dataset=outputDatasets)`'
#'    - mll_arriba: '`sm expand(config["projectPath"] + "/manuscript/arriba/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_star_fusion: '`sm expand(config["projectPath"] + "/manuscript/star_fusion/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_manta: '`sm expand(config["projectPath"] + "/manuscript/manta/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_manta_sv: '`sm expand(config["projectPath"] + "/manuscript/manta_sv/{dataset}.tsv",
#'                      dataset=outputDatasets)`'
#'    - vepRes: '`sm expand(config["vep_path"] +
#'              "/MLL_{dataset}.vep.tsv.gz", dataset=outputDatasets)`'
#'    - ods: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds",
#'             annotation=annotations, dataset=inputDatasets)`'
#'    - ods_filter_out: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods_filter_out.Rds",
#'                        annotation=annotations, dataset=inputDatasets)`'
#'    - mutation_gene_sample: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/mutation_gene_sample.csv"`'
#'    - mutation_cnv_gene_sample: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/mutation_cnv_gene_sample.csv"`'
#'  output:
#'    - or_anova_gene: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_anova_gene.csv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_2_or_anova_gene.html"`'
#'  type: noindex
#'  threads: 1  
#'  resources:
#'    - mem_mb: 128000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_2_or_anova_gene.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2_or_anova_gene.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(OUTRIDER)
})

source("Scripts/manuscript/manuscript_theme.R")

vep_high <- fread(snakemake@input$vep_high)
vep_high_sum <- vep_high[, .N, by=c('SYMBOL', 'array_id')]
setnames(vep_high_sum, 'N', 'vep_high_count')

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]

mutation_cnv <- fread(snakemake@input$mutation_cnv_gene_sample)
or_lm_dt <- mutation_cnv[tool=='outrider', ]
or_lm_dt[, aberrant := padjust < 0.05]
or_lm_dt <- merge(or_lm_dt, gencode[, .(gene_id_unique, gene_id)] %>% unique(),
                  by='gene_id', all.x=TRUE, all.y=FALSE)

ods <- readRDS(snakemake@input$ods)

or_aeCount_zscore_dt <- zScore(ods) %>% as.data.table()
or_aeCount_zscore_dt[, gene_id_unique := rownames(ods)] 
or_aeCount_zscore_dt <- melt(or_aeCount_zscore_dt, id.vars = 'gene_id_unique', 
                             variable.name = "sampleID", value.name = "aeCount_zscore")

or_aeCount_dt <- counts(ods, normalized=TRUE) %>% as.data.table()
or_aeCount_dt[, gene_id_unique := rownames(ods)] 
or_aeCount_dt <- melt(or_aeCount_dt, id.vars = 'gene_id_unique', 
                      variable.name = "sampleID", value.name = "aeCount")

or_sfCount_dt <- counts(ods, normalized = FALSE) %>% as.data.table()
size_factor <- sizeFactors(ods)
or_sfCount_dt <- t(t(or_sfCount_dt)/size_factor) %>% as.data.table()
or_sfCount_dt[, gene_id_unique := rownames(ods)] 
or_sfCount_dt <- melt(or_sfCount_dt, id.vars = 'gene_id_unique', 
                      variable.name = "sampleID", value.name = "sfCount")

or_lm_dt <- merge(or_lm_dt, or_aeCount_zscore_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# or_lm_dt[, table(is.na(zScore))]

or_lm_dt <- merge(or_lm_dt, or_aeCount_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# or_lm_dt[, table(is.na(aeCount))]

or_lm_dt <- merge(or_lm_dt, or_sfCount_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# or_lm_dt[, table(is.na(sfCount))]

or_lm_dt[, vep_high := samp_symbol %in% vep_high[, samp_symbol]]


or_lm_dn_dt <- or_lm_dt[aeCount_zscore < 0, ]
or_lm_up_dt <- or_lm_dt[aeCount_zscore > 0, ]



### aeCount_zscore dn ####
anova_zscore_dn_list <- lapply(or_lm_dn_dt[, unique(gene_id)], function(gene_lm){
  or_lm_gene_dt <- or_lm_dn_dt[gene_id == gene_lm]
  or_lm_gene <- lm(aeCount_zscore ~ log(MEAN_COPY_RATIO) + promoter_variant + vep_high + structural_variant,
                   data=or_lm_gene_dt)
  anova_gene <- anova(or_lm_gene)
  anova_gene_dt <- data.table(var = rownames(anova_gene),
                              sum_sq = anova_gene$`Sum Sq`,
                              mean_sq = anova_gene$`Mean Sq`,
                              gene_id = gene_lm)
  return(anova_gene_dt)
})

anova_zscore <- rbindlist(anova_zscore_dn_list) %>% as.data.table()
anova_zscore <- anova_zscore[var!='Residuals', ]
anova_zscore[, sum_sq_ratio := sum_sq/sum(sum_sq), by='gene_id']
anova_zscore_dn_plot <- anova_zscore[, mean(sum_sq_ratio), by='var']
setnames(anova_zscore_dn_plot, 'V1', 'mean_variance_explained_per_gene')
anova_zscore_dn_plot[, category := 'Underexpression outlier']
anova_zscore_dn_plot[, response := 'Autoencoder corrected zScore']


### sfCount dn ####
anova_sf_dn_list <- lapply(or_lm_dn_dt[, unique(gene_id)], function(gene_lm){
  or_lm_gene_dt <- or_lm_dn_dt[gene_id == gene_lm]
  or_lm_gene <- lm(log(sfCount+1) ~ log(MEAN_COPY_RATIO) + promoter_variant + vep_high + structural_variant,
                   data=or_lm_gene_dt)
  anova_gene <- anova(or_lm_gene)
  anova_gene_dt <- data.table(var = rownames(anova_gene),
                              sum_sq = anova_gene$`Sum Sq`,
                              mean_sq = anova_gene$`Mean Sq`,
                              gene_id = gene_lm)
  return(anova_gene_dt)
})

anova_sf <- rbindlist(anova_sf_dn_list) %>% as.data.table()
anova_sf <- anova_sf[var!='Residuals', ]
anova_sf[, sum_sq_ratio := sum_sq/sum(sum_sq), by='gene_id']
anova_sf_dn_plot <- anova_sf[, mean(sum_sq_ratio), by='var']
setnames(anova_sf_dn_plot, 'V1', 'mean_variance_explained_per_gene')
anova_sf_dn_plot[, category := 'Underexpression outlier']
anova_sf_dn_plot[, response := 'Size factor corrected expression level']


### aeCount dn ####
anova_ae_dn_list <- lapply(or_lm_dn_dt[, unique(gene_id)], function(gene_lm){
  or_lm_gene_dt <- or_lm_dn_dt[gene_id == gene_lm]
  or_lm_gene <- lm(log(aeCount+1) ~ log(MEAN_COPY_RATIO) + promoter_variant + vep_high + structural_variant,
                   data=or_lm_gene_dt)
  anova_gene <- anova(or_lm_gene)
  anova_gene_dt <- data.table(var = rownames(anova_gene),
                              sum_sq = anova_gene$`Sum Sq`,
                              mean_sq = anova_gene$`Mean Sq`,
                              gene_id = gene_lm)
  return(anova_gene_dt)
})

anova_ae <- rbindlist(anova_ae_dn_list) %>% as.data.table()
anova_ae <- anova_ae[var!='Residuals', ]
anova_ae[, sum_sq_ratio := sum_sq/sum(sum_sq), by='gene_id']
anova_ae_dn_plot <- anova_ae[, mean(sum_sq_ratio), by='var']
setnames(anova_ae_dn_plot, 'V1', 'mean_variance_explained_per_gene')
anova_ae_dn_plot[, category := 'Underexpression outlier']
anova_ae_dn_plot[, response := 'Autoencoder corrected expression level']



### aeCount_zscore up ####
anova_zscore_up_list <- lapply(or_lm_up_dt[, unique(gene_id)], function(gene_lm){
  or_lm_gene_dt <- or_lm_up_dt[gene_id == gene_lm]
  or_lm_gene <- lm(aeCount_zscore ~ log(MEAN_COPY_RATIO) + promoter_variant + vep_high + structural_variant,
                   data=or_lm_gene_dt)
  anova_gene <- anova(or_lm_gene)
  anova_gene_dt <- data.table(var = rownames(anova_gene),
                              sum_sq = anova_gene$`Sum Sq`,
                              mean_sq = anova_gene$`Mean Sq`,
                              gene_id = gene_lm)
  return(anova_gene_dt)
})

anova_zscore <- rbindlist(anova_zscore_up_list) %>% as.data.table()
anova_zscore <- anova_zscore[var!='Residuals', ]
anova_zscore[, sum_sq_ratio := sum_sq/sum(sum_sq), by='gene_id']
anova_zscore_up_plot <- anova_zscore[, mean(sum_sq_ratio), by='var']
setnames(anova_zscore_up_plot, 'V1', 'mean_variance_explained_per_gene')
anova_zscore_up_plot[, category := 'Overexpression outlier']
anova_zscore_up_plot[, response := 'Autoencoder corrected zScore']



### sfCount up ####
anova_sf_up_list <- lapply(or_lm_up_dt[, unique(gene_id)], function(gene_lm){
  or_lm_gene_dt <- or_lm_up_dt[gene_id == gene_lm]
  or_lm_gene <- lm(log(sfCount+1) ~ log(MEAN_COPY_RATIO) + promoter_variant + vep_high + structural_variant,
                   data=or_lm_gene_dt)
  anova_gene <- anova(or_lm_gene)
  anova_gene_dt <- data.table(var = rownames(anova_gene),
                              sum_sq = anova_gene$`Sum Sq`,
                              mean_sq = anova_gene$`Mean Sq`,
                              gene_id = gene_lm)
  return(anova_gene_dt)
})

anova_sf <- rbindlist(anova_sf_up_list) %>% as.data.table()
anova_sf <- anova_sf[var!='Residuals', ]
anova_sf[, sum_sq_ratio := sum_sq/sum(sum_sq), by='gene_id']
anova_sf_up_plot <- anova_sf[, mean(sum_sq_ratio), by='var']
setnames(anova_sf_up_plot, 'V1', 'mean_variance_explained_per_gene')
anova_sf_up_plot[, category := 'Overexpression outlier']
anova_sf_up_plot[, response := 'Size factor corrected expression level']


### aeCount up ####
anova_ae_up_list <- lapply(or_lm_up_dt[, unique(gene_id)], function(gene_lm){
  or_lm_gene_dt <- or_lm_up_dt[gene_id == gene_lm]
  or_lm_gene <- lm(log(aeCount+1) ~ log(MEAN_COPY_RATIO) + promoter_variant + vep_high + structural_variant,
                   data=or_lm_gene_dt)
  anova_gene <- anova(or_lm_gene)
  anova_gene_dt <- data.table(var = rownames(anova_gene),
                              sum_sq = anova_gene$`Sum Sq`,
                              mean_sq = anova_gene$`Mean Sq`,
                              gene_id = gene_lm)
  return(anova_gene_dt)
})

anova_ae <- rbindlist(anova_ae_up_list) %>% as.data.table()
anova_ae <- anova_ae[var!='Residuals', ]
anova_ae[, sum_sq_ratio := sum_sq/sum(sum_sq), by='gene_id']
anova_ae_up_plot <- anova_ae[, mean(sum_sq_ratio), by='var']
setnames(anova_ae_up_plot, 'V1', 'mean_variance_explained_per_gene')
anova_ae_up_plot[, category := 'Overexpression outlier']
anova_ae_up_plot[, response := 'Autoencoder corrected expression level']

anova_or_plot <- rbind(anova_sf_dn_plot,
                       anova_sf_up_plot,
                       anova_ae_dn_plot,
                       anova_ae_up_plot,
                       anova_zscore_dn_plot,
                       anova_zscore_up_plot)

fwrite(anova_or_plot, snakemake@output$or_anova_gene)
