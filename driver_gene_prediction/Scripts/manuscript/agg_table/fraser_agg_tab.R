#'---
#' title: generate fraser_agg_tab
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
#'    - fraserRes: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv",
#'                    annotation=annotations, dataset=outputDatasets)`'
#'  output:
#'    - fraser_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/fraser_agg_tab.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/fraser_agg_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/fraser_agg_tab.snakemake")




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

fr_gene <- fr_res[, .N, by=c('hgncSymbol', 'Diag')] 
fr_gene[, n_var_gene := sum(N), by='hgncSymbol']

fr_gene <- merge(fr_gene, 
                 manuscript_wording[, .(Cohort_abbreviation, Cohort_during_analysis)],
                 by.x='Diag', by.y='Cohort_during_analysis')
fr_gene <- merge(fr_gene,
                 gencode[, .(gene_id_unique, gene_id, gene_name)], 
                 by.x='hgncSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)

fr_sum <- dcast(fr_gene, gene_id + hgncSymbol ~ Cohort_abbreviation, value.var='N')
fr_sum[is.na(fr_sum)] <- 0 

fr_sum <- merge(fr_sum, fr_gene[, .(hgncSymbol, n_var_gene)] %>% unique(), by='hgncSymbol')
fr_sum <- fr_sum[order(-n_var_gene),]

setnames(fr_sum, c("gene_id", "hgncSymbol", "n_var_gene"), c('GeneID', 'GeneSymbol', "Total"))
setcolorder(fr_sum, c("GeneID", "GeneSymbol", diag_order, "Total"))


#+ write result
fwrite(fr_sum, snakemake@output$fraser_agg_tab)
