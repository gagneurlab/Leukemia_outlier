#'---
#' title: generate fusion_agg_tab
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
#'    - mll_arriba: '`sm expand(config["projectPath"] + "/manuscript/arriba/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_star_fusion: '`sm expand(config["projectPath"] + "/manuscript/star_fusion/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_manta: '`sm expand(config["projectPath"] + "/manuscript/manta/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'  output:
#'    - fusion_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/fusion_agg_tab.csv"`'
#'  type: script
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/fusion_agg_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/fusion_agg_tab.snakemake")




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

#+ read in 
for (x in snakemake@input$mll_manta) {
  fusion_manta_temp <- fread(x)[, .(array_id, gene1, gene2)]
  
  if (exists("fusion_manta")) {
    fusion_manta <- rbind(fusion_manta, fusion_manta_temp)
  } else {
    fusion_manta <- fusion_manta_temp
  }
}

for (x in snakemake@input$mll_arriba) {
  fusion_arriba_temp <- fread(x)[, .(array_id, gene1, gene2)]
  
  if (exists("fusion_arriba")) {
    fusion_arriba <- rbind(fusion_arriba, fusion_arriba_temp)
  } else {
    fusion_arriba <- fusion_arriba_temp
  }
}

for (x in snakemake@input$mll_star_fusion) {
  fusion_star_fusion_temp <- fread(x)[, .(array_id, gene1, gene2)]
  
  if (exists("fusion_star_fusion")) {
    fusion_star_fusion <- rbind(fusion_star_fusion, fusion_star_fusion_temp)
  } else {
    fusion_star_fusion <- fusion_star_fusion_temp
  }
}

fusion_res <- rbind(fusion_manta, fusion_arriba, fusion_star_fusion) %>% unique()

dup_fusion <- intersect(fusion_res[, paste(array_id, gene1, gene2)],
                        fusion_res[, paste(array_id, gene2, gene1)])
dup_fusion_half <- fusion_res[paste(array_id, gene1, gene2) %in% dup_fusion & gene1<gene2, paste(array_id, gene1, gene2)]
fusion_res <- merge(fusion_res, samp_anno_14group[, .(ArrayID, Diag)], by.x='array_id', by.y='ArrayID')
fusion_res <- fusion_res[!(paste(array_id, gene1, gene2) %in% dup_fusion_half), ]

fusion_res[, Gene_pair := paste0(gene1,"--", gene2) ]
fusion_gene_pair <- fusion_res[, .N, by=c('Gene_pair', 'Diag')] 
fusion_gene_pair[, n_gene_pair := sum(N), by='Gene_pair']

fusion_gene_pair <- merge(fusion_gene_pair, 
                          manuscript_wording[, .(Cohort_abbreviation, Cohort_during_analysis)],
                          by.x='Diag', by.y='Cohort_during_analysis')

fusion_sum <- dcast(fusion_gene_pair, Gene_pair ~ Cohort_abbreviation, value.var='N')
fusion_sum[is.na(fusion_sum)] <- 0 

fusion_sum <- fusion_sum %>% separate('Gene_pair', 
                                      c("GeneSymbol_1", "GeneSymbol_2"), 
                                      sep='--', remove = FALSE) %>% as.data.table()

fusion_sum <- merge(fusion_sum,
                    gencode[, .(gene_id, gene_name)], 
                    by.x='GeneSymbol_1', by.y='gene_name', all.x=TRUE, all.y=FALSE)
setnames(fusion_sum, 'gene_id', 'GeneID_1')
fusion_sum <- merge(fusion_sum,
                    gencode[, .(gene_id, gene_name)], 
                    by.x='GeneSymbol_1', by.y='gene_name', all.x=TRUE, all.y=FALSE)
setnames(fusion_sum, 'gene_id', 'GeneID_2')


fusion_sum <- merge(fusion_sum, fusion_gene_pair[, .(Gene_pair, n_gene_pair)] %>% unique(), by='Gene_pair')
fusion_sum <- fusion_sum[order(-n_gene_pair),]

setnames(fusion_sum, "n_gene_pair", "Total")
setcolorder(fusion_sum, c("Gene_pair", 
                          "GeneID_1", "GeneSymbol_1",
                          "GeneID_2", "GeneSymbol_2",
                          diag_order, "Total"))


#+ write result
fwrite(fusion_sum, snakemake@output$fusion_agg_tab)
