#'---
#' title: figure_2_anova
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
#'    - or_anova_res: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_anova_res.csv"`'
#'    - ac_anova_res: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/ac_anova_res.csv"`'
#'  output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_2_anova.html"`'
#'  type: noindex
#'  threads: 1  
#'  resources:
#'    - mem_mb: 8000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_2_anova.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2_anova.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("Scripts/manuscript/manuscript_theme.R")

or_anova_res <- fread(snakemake@input$or_anova_res)
ac_anova_res <- fread(snakemake@input$ac_anova_res)
# ac_anova_res <- fread('/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/ac_anova_res.csv')
# or_anova_res <- fread('/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/or_anova_res.csv')

anova_res <- rbind(ac_anova_res, or_anova_res)
anova_res[, Category := factor(Category, 
                               levels = c('Underexpression outliers', 'Overexpression outliers', 'Activation outliers'))]
anova_res[, var := factor(var, 
                          levels =  c('Copy ratio', 'VEP HIGH IMPACT variant', 'Residuals'),
                          labels =  c('Copy ratio', 'VEP HIGH\nIMPACT variant', 'Residuals'))]

p <- ggplot(anova_res[var!="Residuals", ], 
            aes(x=var, y=sum_sq_ratio)) +
  facet_wrap('Category') + 
  geom_bar(stat="identity", width=0.7) +
  xlab('Mutation type') +
  ylab('Ratio of variance explained') +
  theme_vale 
# theme(axis.text.x = element_text(angle = 45, hjust=1))

p
