#'---
#' title: generate activation_agg_tab
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
#'    - activtionRes: '`sm expand(config["outriderDir"] +
#'                    "/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/res_filter_out.tsv",
#'                    annotation=annotations, dataset=inputDatasets)`'
#'  output:
#'    - activation_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/activation_agg_tab.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/activation_agg_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/activation_agg_tab.snakemake")




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

ac_res <- fread(snakemake@input$activtionRes)
ac_res <- merge(ac_res, samp_anno_14group[, .(ArrayID, Diag)], by.x='sampleID', by.y='ArrayID')

ac_gene <- ac_res[, .N, by=c('geneID', 'Diag')] 
ac_gene[, n_var_gene := sum(N), by='geneID']

ac_gene <- merge(ac_gene, 
                 manuscript_wording[, .(Cohort_abbreviation, Cohort_during_analysis)],
                 by.x='Diag', by.y='Cohort_during_analysis')
ac_gene <- merge(ac_gene,
                 gencode[, .(gene_id_unique, gene_id, gene_name)], 
                 by.x='geneID', by.y='gene_id_unique', all.x=TRUE, all.y=FALSE)
setnames(ac_gene, c("gene_id", "gene_name"), c('GeneID', 'GeneSymbol'))

ac_sum <- dcast(ac_gene, GeneID + GeneSymbol + geneID ~ Cohort_abbreviation, value.var='N')
ac_sum[is.na(ac_sum)] <- 0 

ac_sum <- merge(ac_sum, ac_gene[, .(geneID, n_var_gene)] %>% unique(), by='geneID')
ac_sum <- ac_sum[order(-n_var_gene),]
ac_sum[, geneID := NULL]

setnames(ac_sum, "n_var_gene", "Total")
setcolorder(ac_sum, c("GeneID", "GeneSymbol", diag_order, "Total"))


#+ write result
fwrite(ac_sum, snakemake@output$activation_agg_tab)
