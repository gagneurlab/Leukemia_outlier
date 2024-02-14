#'---
#' title: figure_2_or_anova
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
#'    - or_anova: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_anova.csv"`'
#'    - or_anova_res: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_anova_res.csv"`'
#'    - or_lm_cnv_ae: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_lm_cnv_ae.Rds"`'
#'    - or_lm_cnv_exist_ae: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_lm_cnv_exist_ae.Rds"`'
#'    - or_lm_cnv_sf: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_lm_cnv_sf.Rds"`'
#'    - or_lm_cnv_exist_sf: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_lm_cnv_exist_sf.Rds"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_2_or_anova.html"`'
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
                             "/processed_data/snakemake/figure_2_or_anova.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2_or_anova.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("Scripts/manuscript/manuscript_theme.R")

vep_high <- fread(snakemake@input$vep_high)
vep_high_sum <- vep_high[, .N, by=c('SYMBOL', 'array_id')]
setnames(vep_high_sum, 'N', 'vep_high_count')

vep_splice <- fread(snakemake@input$vep_splice)
vep_splice_sum <- vep_splice[, .N, by=c('SYMBOL', 'array_id')]
setnames(vep_splice_sum, 'N', 'vep_splice_count')

vep_stop <- vep_high[!Consequence %in% c('splice_acceptor_variant', 'splice_donor_variant'), ]
vep_stop_sum <- vep_stop[, .N, by=c('SYMBOL', 'array_id')]
setnames(vep_stop_sum, 'N', 'vep_stop_count')

library(OUTRIDER)

ods <- readRDS(snakemake@input$ods)

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]

for (x in snakemake@input$mll_cnv_full) {
  cnv_res_temp <- fread(x)
  cnv_res_temp[, LENGTH := END - START]
  
  if (exists("cnv_res")) {
    cnv_res <- rbind(cnv_res, cnv_res_temp, fill=TRUE)
  } else {
    cnv_res <- cnv_res_temp
  }
}

or_zscore_dt <- zScore(ods) %>% as.data.table()
or_zscore_dt[, gene_id_unique := rownames(ods)] 
or_zscore_dt <- melt(or_zscore_dt, id.vars = 'gene_id_unique', 
                     variable.name = "sampleID", value.name = "zScore")

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

or_lm_dt <- aberrant(ods) %>% as.data.table()
or_lm_dt[, gene_id_unique := rownames(ods)] 
or_lm_dt <- merge(or_lm_dt, 
                  gencode[, .(gene_id_unique, gene_id, gene_name_orig)], 
                  by='gene_id_unique', all.x=TRUE, all.y=FALSE)
or_lm_dt <- melt(or_lm_dt, id.vars = c('gene_id_unique', 'gene_id', 'gene_name_orig'),
                 variable.name = "sampleID", value.name = "aberrant")

or_lm_dt <- merge(or_lm_dt, or_zscore_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
or_lm_dt[, table(is.na(zScore))]

or_lm_dt <- merge(or_lm_dt, or_aeCount_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
or_lm_dt[, table(is.na(aeCount))]

or_lm_dt <- merge(or_lm_dt, or_sfCount_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
or_lm_dt[, table(is.na(sfCount))]

or_lm_dt <- merge(or_lm_dt, cnv_res[, .(SYMBOL, ARRAY_ID, MEAN_LOG2_COPY_RATIO, MEAN_COPY_RATIO)], 
                  by.x=c('gene_name_orig', 'sampleID'),
                  by.y=c('SYMBOL', 'ARRAY_ID'),
                  all.x=TRUE, all.y=FALSE)
or_lm_dt[, table(is.na(MEAN_COPY_RATIO))]
or_lm_dt[is.na(MEAN_COPY_RATIO), MEAN_COPY_RATIO := 1]
or_lm_dt[is.na(MEAN_LOG2_COPY_RATIO), MEAN_LOG2_COPY_RATIO := 0]
or_lm_dt[, table(is.na(MEAN_COPY_RATIO))]

or_lm_dt <- merge(or_lm_dt, vep_high_sum, 
                  by.x=c('gene_name_orig', 'sampleID'),
                  by.y=c('SYMBOL', 'array_id'),
                  all.x=TRUE, all.y=FALSE)
or_lm_dt[, table(is.na(vep_high_count))]
or_lm_dt[is.na(vep_high_count), vep_high_count := 0]
or_lm_dt[, vep_high_exist := as.factor(as.numeric(vep_high_count >0))]
or_lm_dt[, table(vep_high_exist)]

or_lm_dt <- merge(or_lm_dt, vep_stop_sum, 
                  by.x=c('gene_name_orig', 'sampleID'),
                  by.y=c('SYMBOL', 'array_id'),
                  all.x=TRUE, all.y=FALSE)
or_lm_dt[, table(is.na(vep_stop_count))]
or_lm_dt[is.na(vep_stop_count), vep_stop_count := 0]
or_lm_dt[, vep_stop_exist := as.factor(as.numeric(vep_stop_count >0))]
or_lm_dt[, table(vep_stop_exist)]

or_lm_dt <- merge(or_lm_dt, vep_splice_sum, 
                  by.x=c('gene_name_orig', 'sampleID'),
                  by.y=c('SYMBOL', 'array_id'),
                  all.x=TRUE, all.y=FALSE)
or_lm_dt[, table(is.na(vep_splice_count))]
or_lm_dt[is.na(vep_splice_count), vep_splice_count := 0]
or_lm_dt[, vep_splice_exist := as.factor(as.numeric(vep_splice_count >0))]
or_lm_dt[, table(vep_splice_exist)]

fwrite(or_lm_dt, snakemake@output$or_anova)
# or_lm_dt <- fread(snakemake@output$or_anova)

or_lm_dn_dt <- or_lm_dt[zScore < 0, ]
or_lm_up_dt <- or_lm_dt[zScore > 0, ]


### aberrant - downregulated ####
# or_lm_dn_cnv_count <- lm(aberrant ~ MEAN_COPY_RATIO + vep_stop_count + vep_splice_count,
#                     data=or_lm_dn_dt)
# or_lm_dn_cnv <- lm(aberrant ~ MEAN_COPY_RATIO,
#                        data=or_lm_dn_dt)
# 
# summary(or_lm_dn_cnv)
# summary(or_lm_dn_cnv_count)
# anova(or_lm_dn_cnv, or_lm_dn_cnv_count)



### aberrant - not overexp outlier ####
# or_lm_no_up_dt <- or_lm_dt[!(zScore > 0 & aberrant == 0), ]
# 
# or_lm_no_up_cnv_exist <- lm(aberrant ~ MEAN_COPY_RATIO + vep_high_exist,
#                      data=or_lm_no_up_dt)
# or_lm_no_up_cnv <- lm(aberrant ~ MEAN_COPY_RATIO,
#                      data=or_lm_no_up_dt)
# 
# summary(or_lm_no_up_cnv)
# summary(or_lm_no_up_cnv_exist)
# anova(or_lm_no_up_cnv, or_lm_no_up_cnv_exist)



### aeCount - underexp outlier ####
or_lm_dn_ab_dt <- or_lm_dt[zScore < 0 & aberrant == 1, ]

or_lm_dn_ab_cnv_exist_ae <- lm(log(aeCount+1) ~ log(MEAN_COPY_RATIO) + vep_high_exist,
                               data=or_lm_dn_ab_dt)
or_lm_dn_ab_cnv_ae <- lm(log(aeCount+1) ~ log(MEAN_COPY_RATIO),
                         data=or_lm_dn_ab_dt)


summary(or_lm_dn_ab_cnv_ae)
summary(or_lm_dn_ab_cnv_exist_ae)
anova(or_lm_dn_ab_cnv_ae, or_lm_dn_ab_cnv_exist_ae)



### aeCount - overexp outlier ####
or_lm_up_ab_dt <- or_lm_dt[zScore > 0 & aberrant == 1, ]

or_lm_up_ab_cnv_exist_ae <- lm(log(aeCount+1) ~ log(MEAN_COPY_RATIO) + vep_high_exist,
                               data=or_lm_up_ab_dt)
or_lm_up_ab_cnv_ae <- lm(log(aeCount+1) ~ log(MEAN_COPY_RATIO),
                         data=or_lm_up_ab_dt)


summary(or_lm_up_ab_cnv_ae)
summary(or_lm_up_ab_cnv_exist_ae)
anova(or_lm_up_ab_cnv_ae, or_lm_up_ab_cnv_exist_ae)


### aeCount - all ####
or_lm_cnv_exist_ae <- lm(log(aeCount+1) ~ log(MEAN_COPY_RATIO) + vep_high_exist,
                         data=or_lm_dt)
saveRDS(or_lm_cnv_exist_ae, snakemake@output$or_lm_cnv_exist_ae)
# or_lm_cnv_exist_ae <- readRDS(snakemake@output$or_lm_cnv_exist_ae)
or_lm_cnv_ae <- lm(log(aeCount+1) ~ log(MEAN_COPY_RATIO),
                   data=or_lm_dt)
saveRDS(or_lm_cnv_ae, snakemake@output$or_lm_cnv_ae)
# or_lm_cnv_ae <- readRDS(snakemake@output$or_lm_cnv_ae)

summary(or_lm_cnv_ae)
summary(or_lm_cnv_exist_ae)
anova(or_lm_cnv_ae, or_lm_cnv_exist_ae)



### sfCount - underexp outlier ####
or_lm_dn_ab_dt <- or_lm_dt[zScore < 0 & aberrant == 1, ]

or_lm_dn_ab_cnv_exist_sf <- lm(log(sfCount+1) ~ log(MEAN_COPY_RATIO) + vep_high_exist,
                               data=or_lm_dn_ab_dt)
or_lm_dn_ab_cnv_sf <- lm(log(sfCount+1) ~ log(MEAN_COPY_RATIO),
                         data=or_lm_dn_ab_dt)


summary(or_lm_dn_ab_cnv_sf)
summary(or_lm_dn_ab_cnv_exist_sf)
anova(or_lm_dn_ab_cnv_exist_sf)
anova(or_lm_dn_ab_cnv_sf, or_lm_dn_ab_cnv_exist_sf)

or_lm_dn_ab_cnv_exist_sf_anova <- anova(or_lm_dn_ab_cnv_exist_sf)
or_lm_dn_ab_cnv_exist_sf_anova_plot <- data.table(var = c('Copy ratio', 'VEP HIGH IMPACT variant', 'Residuals'),
                                                  sum_sq = or_lm_dn_ab_cnv_exist_sf_anova$`Sum Sq`)
or_lm_dn_ab_cnv_exist_sf_anova_plot[, sum_sq_ratio := sum_sq/sum(sum_sq)]
or_lm_dn_ab_cnv_exist_sf_anova_plot[, Category := 'Underexpression outliers']

p_a <- ggplot(or_lm_dn_ab_cnv_exist_sf_anova_plot[var!="Residuals", ], 
              aes(x=var, y=sum_sq_ratio)) +
  geom_bar(stat="identity", width=0.7) +
  xlab('Mutation type') +
  ylab('Ratio of variance explained') +
  ggtitle('Underexpression outliers') +
  theme_vale 

p_a


### sfCount - overexp outlier ####
or_lm_up_ab_dt <- or_lm_dt[zScore > 0 & aberrant == 1, ]

or_lm_up_ab_cnv_exist_sf <- lm(log(sfCount+1) ~ log(MEAN_COPY_RATIO) + vep_high_exist,
                               data=or_lm_up_ab_dt)
or_lm_up_ab_cnv_sf <- lm(log(sfCount+1) ~ log(MEAN_COPY_RATIO),
                         data=or_lm_up_ab_dt)


summary(or_lm_up_ab_cnv_sf)
summary(or_lm_up_ab_cnv_exist_sf)
anova(or_lm_up_ab_cnv_sf, or_lm_up_ab_cnv_exist_sf)

or_lm_up_ab_cnv_exist_sf_anova <- anova(or_lm_up_ab_cnv_exist_sf)
or_lm_up_ab_cnv_exist_sf_anova_plot <- data.table(var = c('Copy ratio', 'VEP HIGH IMPACT variant', 'Residuals'),
                                                  sum_sq = or_lm_up_ab_cnv_exist_sf_anova$`Sum Sq`)
or_lm_up_ab_cnv_exist_sf_anova_plot[, sum_sq_ratio := sum_sq/sum(sum_sq)]
or_lm_up_ab_cnv_exist_sf_anova_plot[, Category := 'Overexpression outliers']

or_anova_res <- rbind(or_lm_up_ab_cnv_exist_sf_anova_plot, or_lm_dn_ab_cnv_exist_sf_anova_plot)

fwrite(or_anova_res, snakemake@output$or_anova_res)

p_b <- ggplot(or_lm_up_ab_cnv_exist_sf_anova_plot[var!="Residuals", ], 
              aes(x=var, y=sum_sq_ratio)) +
  geom_bar(stat="identity", width=0.7) +
  xlab('Mutation type') +
  ylab('Ratio of variance explained') +
  ggtitle('Overexpression outliers') +
  theme_vale 

p_b



### sfCount - all ####
or_lm_cnv_exist_sf <- lm(log(sfCount+1) ~ log(MEAN_COPY_RATIO) + vep_high_exist,
                         data=or_lm_dt)
saveRDS(or_lm_cnv_exist_sf, snakemake@output$or_lm_cnv_exist_sf)
# or_lm_cnv_exist_sf <- readRDS(snakemake@output$or_lm_cnv_exist_sf)
or_lm_cnv_sf <- lm(log(sfCount+1) ~ log(MEAN_COPY_RATIO),
                   data=or_lm_dt)
saveRDS(or_lm_cnv_sf, snakemake@output$or_lm_cnv_sf)
# or_lm_cnv_sf <- readRDS(snakemake@output$or_lm_cnv_sf)

summary(or_lm_cnv_sf)
summary(or_lm_cnv_exist_sf)
anova(or_lm_cnv_sf, or_lm_cnv_exist_sf)
