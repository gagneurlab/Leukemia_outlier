#'---
#' title: figure_4_sup
#' author: Xueqi Cao, Ata Jadid Ahari
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
#'    - experimentDesign: '`sm config["experimentDesign"]`'
#'    - manuscriptWording: '`sm config["manuscriptWording"]`'
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - intogenDir: '`sm config["intogenDir"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'    - MLL_CGC_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/MLL_CGC_leukemia_gene_list.tsv"`'
#'  input:
#'    - mll_prc: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_prc.tsv"`'
#'    - mll_ap: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_ap.tsv"`'
#'    - mll_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_ap_full.tsv"`'
#'    - mll_var_filter_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_var_filter_prc.tsv"`'
#'    - mll_var_filter_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_var_filter_ap.tsv"`'
#'    - mll_var_filter_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_var_filter_ap_full.tsv"`'
#'    - mll_emb_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_emb_prc.tsv"`'
#'    - mll_emb_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_emb_ap.tsv"`'
#'    - mll_emb_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_emb_ap_full.tsv"`'
#'    - mll_benchmark_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_benchmark_prc.tsv"`'
#'    - mll_benchmark_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_benchmark_ap.tsv"`'
#'    - p_c: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_c.Rds"`' 
#'    - p_d: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_d.Rds"`' 
#'    - p_c_emb: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_c_emb.Rds"`' 
#'    - p_d_emb: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_d_emb.Rds"`'
#'    - p_c_emb_exp: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_c_emb_exp.Rds"`' 
#'    - p_d_emb_exp: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_d_emb_exp.Rds"`' 
#'    - p_c_coess: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_c_coess.Rds"`' 
#'    - p_d_coess: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_d_coess.Rds"`' 
#'  output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_4_sup.html"`'
#'  type: noindex
#'  resources:
#'    - mem_mb: 8000 
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_4_sup.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/figure_4_sup.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

p_c <- readRDS(snakemake@input$p_c)
p_d <- readRDS(snakemake@input$p_d)
p_c_emb <- readRDS(snakemake@input$p_c_emb)
p_d_emb <- readRDS(snakemake@input$p_d_emb)
p_c_emb_exp <- readRDS(snakemake@input$p_c_emb_exp)
p_d_emb_exp <- readRDS(snakemake@input$p_d_emb_exp)
p_c_coess <- readRDS(snakemake@input$p_c_coess)
p_d_coess <- readRDS(snakemake@input$p_d_coess)




### Supplement #####
#' ## Sup
#' ### s6 raw
#+ plot s6 raw, fig.width=10, fig.height=6
p_s6_raw <- ggarrange(p_c, p_c_emb, nrow = 1, labels = c('A', 'B'))
p_s6_raw
#' ### s6 annotated
#+ plot s6, fig.width=10, fig.height=6
p_s6 <- p_s6_raw
p_s6

#' ### s7 raw
#+ plot s7 raw, fig.width=10, fig.height=14
p_s7_raw <- ggarrange(p_d, p_d_emb, nrow = 1, labels = c('A', 'B'))
p_s7_raw
#' ### s7 annotated
#+ plot s7, fig.width=10, fig.height=14
p_s7 <- p_s7_raw
p_s7

#' ### s8 raw
#+ plot s8 raw, fig.width=10, fig.height=6
p_s8_raw <- ggarrange(p_c, p_c_emb_exp, nrow = 1, labels = c('A', 'B'))
p_s8_raw
#' ### s8 annotated
#+ plot s8, fig.width=10, fig.height=6
p_s8 <- p_s8_raw
p_s8

#' ### s9 raw
#+ plot s9 raw, fig.width=10, fig.height=14
p_s9_raw <- ggarrange(p_d, p_d_emb_exp, nrow = 1, labels = c('A', 'B'))
p_s9_raw
#' ### s9 annotated
#+ plot s9, fig.width=10, fig.height=14
p_s9 <- p_s9_raw
p_s9
