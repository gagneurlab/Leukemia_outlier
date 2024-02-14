#'---
#' title: figure_2_anova_gene
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
#'    - or_anova_gene: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/or_anova_gene.csv"`'
#'    - ac_anova_gene: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/ac_anova_gene.csv"`'
#'  output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_2_anova_gene.html"`'
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
                             "/processed_data/snakemake/figure_2_anova_gene.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2_anova_gene.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("Scripts/manuscript/manuscript_theme.R")

or_anova_gene <- fread(snakemake@input$or_anova_gene)
ac_anova_gene <- fread(snakemake@input$ac_anova_gene)

anova_gene <- rbind(or_anova_gene, ac_anova_gene)

### plot ####
anova_gene[, var := factor(var, 
                           levels = c(
                             "log(MEAN_COPY_RATIO)",   
                             "vep_high",
                             "promoter_variant",
                             "structural_variant",
                             "Residuals"
                           ),
                           labels = c(
                             "Copy ratio",   
                             "VEP HIGH\nIMPACT",
                             "Promoter\n(TSS±2Kbp)",
                             "Structural",
                             "Residuals"
                           ))]

# anova_gene[, response := factor(response, 
#                                 levels=c('Size factor corrected expression level', 
#                                          'Autoencoder corrected expression level'),
#                                 labels=c('Size factor corrected\nexpression level', 
#                                          'Autoencoder corrected\nexpression level'
#                                          
#                                 ))]

anova_gene[, category := factor(category, levels = c(
  'Underexpression outlier', 'Overexpression outlier', 'Activation outlier'
))]

p_a <- ggplot(anova_gene[response=="Autoencoder corrected zScore", ], 
              aes(x=var, y=mean_variance_explained_per_gene)) +
  facet_wrap('category') +
  geom_bar(stat="identity", width=0.7) +
  xlab('Mutation type') +
  ylab('Mean variance\nexplained per gene (%)') +
  scale_y_continuous(labels = scales::percent) +
  theme_vale +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

p_a

p_a <- ggplot(anova_gene[response=="Autoencoder corrected expression level", ], 
              aes(x=var, y=mean_variance_explained_per_gene)) +
  facet_wrap('category') +
  geom_bar(stat="identity", width=0.7) +
  xlab('Mutation type') +
  ylab('Mean variance\nexplained per gene (%)') +
  scale_y_continuous(labels = scales::percent) +
  theme_vale +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

p_a
