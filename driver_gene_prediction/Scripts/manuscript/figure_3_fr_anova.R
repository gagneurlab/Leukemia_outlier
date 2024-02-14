#'---
#' title: figure_3_fr_anova
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
#'  output:
#'    - fr_anova: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/anova/fr_anova.csv"`'
#'    - fr_lm_ab: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/anova/fr_lm_ab.RData"`'
#'    - fr_lm_deltaPsi: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/anova/fr_lm_deltaPsi.RData"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_3_fr_anova.html"`'
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
                             "/processed_data/snakemake/figure_3_fr_anova.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_3_fr_anova.snakemake")

.libPaths("~/R/4.1/FRASER2")
suppressPackageStartupMessages({
  library(FRASER)
  library(data.table)
  library(dplyr)
  library(tidyr)
})


vep_high <- fread(snakemake@input$vep_high)
vep_high_sum <- vep_high[, .N, by=c('SYMBOL', 'array_id')]
setnames(vep_high_sum, 'N', 'vep_high_count')
vep_high_sum[, samp_symbol := paste0(array_id, "-", SYMBOL)]

vep_splice <- fread(snakemake@input$vep_splice)
vep_splice_sum <- vep_splice[, .N, by=c('SYMBOL', 'array_id')]
setnames(vep_splice_sum, 'N', 'vep_splice_count')
vep_splice_sum[, samp_symbol := paste0(array_id, "-", SYMBOL)]

vep_stop <- vep_high[!Consequence %in% c('splice_acceptfr_variant', 'splice_donfr_variant'), ]
vep_stop_sum <- vep_stop[, .N, by=c('SYMBOL', 'array_id')]
setnames(vep_stop_sum, 'N', 'vep_stop_count')
vep_stop_sum[, samp_symbol := paste0(array_id, "-", SYMBOL)]

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]
gencode <- gencode[gene_type == 'protein_coding', ]

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()
samp_anno_14groups <- samp_anno_exp[DROP_GROUP %in% snakemake@params$inputDatasets, ] %>% unique()

for (x in snakemake@input$fraserRes) {
  fr_res_temp <- fread(x)
  fr_res_temp <- merge(fr_res_temp, gencode[, .(gene_name, gene_type)], 
                       by.x='hgncSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  fr_res_temp <- fr_res_temp[gene_type=='protein_coding', ]
  
  if (exists("fr_res")) {
    fr_res <- rbind(fr_res, fr_res_temp, fill=TRUE)
  } else {
    fr_res <- fr_res_temp
  }
}
fr_res[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]

for (x in snakemake@input$mll_cnv_full) {
  # for (x in snakemake@input$mll_cnv) {
  cnv_res_temp <- fread(x)
  cnv_res_temp[, LENGTH := END - START]
  
  if (exists("cnv_res")) {
    cnv_res <- rbind(cnv_res, cnv_res_temp, fill=TRUE)
  } else {
    cnv_res <- cnv_res_temp
  }
}

absplice_res <- fread(snakemake@input$abspliceRes)
absplice_res[, samp_symbol := paste0(sampleID, "-", gene_name)]
absplice_res <- separate(absplice_res, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_res <- separate(absplice_res, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()
absplice_res <- absplice_res[AbSplice_DNA>=0.2, ]

fr_lm_dt <- matrix(data=FALSE, nrow = nrow(gencode), ncol = nrow(samp_anno_14groups), 
                   dimnames = list(gencode[,gene_id_unique], samp_anno_14groups[, ArrayID])) %>% as.data.table()
fr_lm_dt[, gene_id_unique := gencode[,gene_id_unique]] 
fr_lm_dt <- merge(fr_lm_dt, 
                  gencode[, .(gene_id_unique, gene_id, gene_name_orig)], 
                  by='gene_id_unique', all.x=TRUE, all.y=FALSE)
fr_lm_dt <- melt(fr_lm_dt, id.vars = c('gene_id_unique', 'gene_id', 'gene_name_orig'),
                 variable.name = "sampleID", value.name = "aberrant")
fr_lm_dt[, samp_symbol := paste0(sampleID, "-", gene_name_orig)]

fr_lm_dt <- merge(fr_lm_dt, 
                  fr_res[, .(sampleID, hgncSymbol, deltaPsi)], 
                  by.x=c('sampleID', 'gene_name_orig'),
                  by.y=c('sampleID', 'hgncSymbol'),
                  all.x=TRUE, all.y=FALSE)
fr_lm_dt[is.na(deltaPsi), deltaPsi := 0]

fr_lm_dt[samp_symbol %in% fr_res[, samp_symbol], aberrant := TRUE]
fr_lm_dt[, table(aberrant)]

fr_lm_dt[, absplice := samp_symbol %in% absplice_res[, samp_symbol]]
fr_lm_dt[, table(absplice)]

fr_lm_dt[, vep_splice := samp_symbol %in% vep_splice[, samp_symbol]]
fr_lm_dt[, table(vep_splice)]

fr_lm_dt[, vep_high := samp_symbol %in% vep_high[, samp_symbol]]
fr_lm_dt[, table(vep_high)]

fr_lm_dt[, vep_stop := samp_symbol %in% vep_stop[, samp_symbol]]
fr_lm_dt[, table(vep_stop)]

fr_lm_dt <- merge(fr_lm_dt, cnv_res[, .(SYMBOL, ARRAY_ID, MEAN_LOG2_COPY_RATIO, MEAN_COPY_RATIO)], 
                  by.x=c('gene_name_orig', 'sampleID'),
                  by.y=c('SYMBOL', 'ARRAY_ID'),
                  all.x=TRUE, all.y=FALSE)

fr_lm_dt[, table(is.na(MEAN_COPY_RATIO))]
fr_lm_dt[is.na(MEAN_COPY_RATIO), MEAN_COPY_RATIO := 1]
fr_lm_dt[is.na(MEAN_LOG2_COPY_RATIO), MEAN_LOG2_COPY_RATIO := 0]
fr_lm_dt[, table(is.na(MEAN_COPY_RATIO))]

fwrite(fr_lm_dt, snakemake@output$fr_anova)
# fr_lm_dt <- fread(snakemake@output$fr_anova)




### aberrant - all ####
fr_lm_ab_cnv_vep_absplice <- lm(aberrant ~ MEAN_COPY_RATIO + absplice + vep_splice + vep_stop,
                             data=fr_lm_dt)

summary(fr_lm_ab_cnv_vep_absplice)
anova(fr_lm_ab_cnv_vep_absplice)

fr_lm_ab_vep_absplice <- lm(aberrant ~ absplice + vep_splice,
                         data=fr_lm_dt)
fr_lm_ab_vep <- lm(aberrant ~ vep_splice,
                data=fr_lm_dt)
fr_lm_ab_absplice <- lm(aberrant ~ absplice,
                     data=fr_lm_dt)
save(fr_lm_ab_cnv_vep_absplice, fr_lm_ab_vep_absplice, fr_lm_ab_vep, fr_lm_ab_absplice, 
     file=snakemake@output$fr_lm_ab)
# load(snakemake@output$fr_lm_ab)

anova(fr_lm_ab_vep_absplice)
anova(fr_lm_ab_vep, fr_lm_ab_vep_absplice)
anova(fr_lm_ab_absplice, fr_lm_ab_vep_absplice)




### deltaPsi - outlier ####
fr_lm_deltaPsi_outlier_dt <- fr_lm_dt[aberrant == TRUE, ]

fr_lm_deltaPsi_outlier_vep_absplice_psi <- lm(deltaPsi ~ absplice + vep_splice,
                               data=fr_lm_deltaPsi_outlier_dt)
fr_lm_deltaPsi_outlier_absplice_psi <- lm(deltaPsi ~ absplice,
                               data=fr_lm_deltaPsi_outlier_dt)
fr_lm_deltaPsi_outlier_vep_psi <- lm(deltaPsi ~ vep_splice,
                               data=fr_lm_deltaPsi_outlier_dt)

summary(fr_lm_deltaPsi_outlier_vep_psi)
summary(fr_lm_deltaPsi_outlier_absplice_psi)
anova(fr_lm_deltaPsi_outlier_vep_psi, fr_lm_deltaPsi_outlier_vep_absplice_psi)
anova(fr_lm_deltaPsi_outlier_absplice_psi, fr_lm_deltaPsi_outlier_vep_absplice_psi)


### deltaPsi - all ####
fr_lm_deltaPsi_vep_absplice_psi <- lm(deltaPsi ~ absplice + vep_splice,
                             data=fr_lm_dt)
fr_lm_deltaPsi_absplice_psi <- lm(deltaPsi ~ absplice,
                         data=fr_lm_dt)
fr_lm_deltaPsi_vep_psi <- lm(deltaPsi ~ vep_splice,
                    data=fr_lm_dt)
save(fr_lm_deltaPsi_vep_absplice_psi, fr_lm_deltaPsi_absplice_psi, fr_lm_deltaPsi_vep_psi, 
     file=snakemake@output$fr_lm_deltaPsi) 
# load(snakemake@output$fr_lm_deltaPsi)

summary(fr_lm_deltaPsi_vep_psi)
summary(fr_lm_deltaPsi_absplice_psi)
anova(fr_lm_deltaPsi_vep_psi, fr_lm_deltaPsi_vep_absplice_psi)
anova(fr_lm_deltaPsi_absplice_psi, fr_lm_deltaPsi_vep_absplice_psi)
