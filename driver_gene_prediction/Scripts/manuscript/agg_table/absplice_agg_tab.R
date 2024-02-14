#'---
#' title: generate absplice_agg_tab
#' author: xueqicao
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
#'    - manuscriptWording: '`sm config["manuscriptWording"]`'
#'  input:
#'    - abspliceResVar: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice_var.tsv"`'
#'    - n_var_gene_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/n_var_gene_tab.csv"`'
#'  output:
#'    - absplice_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/absplice_agg_tab.csv"`'
#'    - absplice_ratio_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/absplice_ratio_tab.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/absplice_agg_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/absplice_agg_tab.snakemake")




##### process the data #####
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(tidyr)
})


#+ define paths
single_group <- snakemake@params$inputDatasets

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_14group <- samp_anno[grep(single_group, DROP_GROUP),]
samp_anno_14group[, Diag_sum := .N, by='Diag']

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))
diag_order <- manuscript_wording[order(-Number_of_samples_per_cohort), unique(Cohort_abbreviation)]

ab_res <- fread(snakemake@input$abspliceResVar)
ab_res <- ab_res[AbSplice_DNA >= 0.2, ]
ab_res <- merge(ab_res, samp_anno_14group[, .(ArrayID, Diag)], by.x='sampleID', by.y='ArrayID')

ab_gene <- ab_res[, .N, by=c('gene_id', 'Diag')] 
ab_gene[, n_var_gene := sum(N), by='gene_id']

ab_gene <- merge(ab_gene, 
                 manuscript_wording[, .(Cohort_abbreviation, Cohort_during_analysis)],
                 by.x='Diag', by.y='Cohort_during_analysis')
ab_gene <- merge(ab_gene,
                 gencode[, .(gene_id, gene_name)], 
                 by='gene_id', all.x=TRUE, all.y=FALSE)

ab_sum <- dcast(ab_gene, gene_id + gene_name ~ Cohort_abbreviation, value.var='N')
ab_sum[is.na(ab_sum)] <- 0 

ab_sum <- merge(ab_sum, ab_gene[, .(gene_id, n_var_gene)] %>% unique(), by='gene_id')
ab_sum <- ab_sum[order(-n_var_gene),]

setnames(ab_sum, c("gene_id", "gene_name", "n_var_gene"), c('GeneID', 'GeneSymbol', "Total"))
setcolorder(ab_sum, c("GeneID", "GeneSymbol", diag_order, "Total"))

#+ write result
fwrite(ab_sum, snakemake@output$absplice_agg_tab)

#+ create ratio
n_var_gene_entity_dt <- fread(snakemake@input$n_var_gene_tab)
n_var_gene_entity <- copy(n_var_gene_entity_dt)
n_var_gene_entity <- melt(n_var_gene_entity[, Total:=NULL], 
                          id.vars=c("GeneID", "GeneSymbol"), 
                          variable.name = "entity", value.name = "n_var_total",)

ab_var_gene_entity <- ab_gene[, .(gene_id, gene_name, Cohort_abbreviation, N)]
colnames(ab_var_gene_entity) <- c("GeneID", "GeneSymbol", "entity", "n_var_ab")

var_gene_entity <- merge(n_var_gene_entity, ab_var_gene_entity, by=c("GeneID", "GeneSymbol", "entity"),all=TRUE)
var_gene_entity[, ab_ratio_entity := n_var_ab/n_var_total]

ratio_gene_entity <- dcast(var_gene_entity, GeneID + GeneSymbol ~ entity, value.var='ab_ratio_entity')
ratio_gene_entity[is.na(ratio_gene_entity)] <- 0 


ab_var_gene_total <- ab_gene[, .(gene_id, gene_name, n_var_gene)] %>% unique()
colnames(ab_var_gene_total) <- c("GeneID", "GeneSymbol", "n_var_ab")

var_gene_total <- merge(n_var_gene_entity_dt[, .(GeneID, GeneSymbol, Total)], 
                        ab_var_gene_total, by=c("GeneID", "GeneSymbol"), all=TRUE)
var_gene_total[, ab_ratio_total := n_var_ab/Total]
var_gene_total[is.na(var_gene_total)] <- 0

ratio_gene_entity <- merge(ratio_gene_entity, var_gene_total[, .(GeneID, GeneSymbol, ab_ratio_total)], 
                           by=c('GeneID', 'GeneSymbol'))
ratio_gene_entity <- ratio_gene_entity[order(-ab_ratio_total),]

setnames(ratio_gene_entity, c("ab_ratio_total"), c("All_samples"))
setcolorder(ratio_gene_entity, c("GeneID", "GeneSymbol", diag_order, "All_samples"))


fwrite(ratio_gene_entity, snakemake@output$absplice_ratio_tab)