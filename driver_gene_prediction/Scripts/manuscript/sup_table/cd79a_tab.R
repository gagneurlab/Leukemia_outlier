#'---
#' title: cd79a_tab
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
#'    - anonymization_table_mll: '`sm config["anonymization_table_mll"]`' 
#'  input:
#'    - cd79a_tab_raw: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/cd79a_tab_raw.csv"`'
#'  output:
#'    - cd79a_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/cd79a_tab.csv"`'
#'  type: script
#'  threads: 1
#'  resources:
#'    - mem_mb: 8000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/cd79a_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/cd79a_tab.snakemake")




##### process the data #####
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(tidyr)
})

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))
manuscript_wording[, Diag := Cohort_during_analysis]
manuscript_wording[, Diag := gsub(" ", "_", Diag)]
manuscript_wording[, Diag := gsub("/", "_", Diag)]
manuscript_wording[, Diag := gsub("-", "_", Diag)]
manuscript_wording[, Cohort_abbr_filename := Cohort_abbreviation]
manuscript_wording[, Cohort_abbr_filename := gsub(" ", "_", Cohort_abbr_filename)]
manuscript_wording[, Cohort_abbr_filename := gsub("/", "_", Cohort_abbr_filename)]
manuscript_wording[, Cohort_abbr_filename := gsub("-", "_", Cohort_abbr_filename)]

manuscript_wording_sub <- manuscript_wording[, .(Cohort_abbreviation, Cohort_during_analysis)]

gencode <- fread(snakemake@params$gencode)

anonymization_table_mll <- fread(snakemake@params$anonymization_table_mll)

#+ define paths
cd79a_tab_raw <- fread(snakemake@input$cd79a_tab_raw)

cd79a_tab <- merge(cd79a_tab_raw, manuscript_wording_sub, by.x='Diag', by.y='Cohort_during_analysis')
cd79a_tab <- merge(cd79a_tab, manuscript_wording[, .(Study_group_during_analysis, Study_group)] %>% unique(), 
                   by.x='DROP_GROUP.y', by.y='Study_group_during_analysis')
cd79a_tab <- merge(cd79a_tab, gencode[, .(gene_name, gene_id)], by.x='hgncSymbol', by.y='gene_name')
cd79a_tab <- merge(cd79a_tab, anonymization_table_mll, by='ArrayID')

cd79a_tab <- setnames(cd79a_tab, 
                      c('AnonamizedID', 'Cohort_abbreviation', 'Study_group', 'hgncSymbol', 'gene_id'),
                      c('sampleID', 'EntityAbbreviation', 'StudyGroup', 'GeneSymbol', 'GeneID')
)

cd79a_tab[, GeneID := sapply(cd79a_tab[, GeneID], function(x){strsplit(x, '[.]')[[1]][1]})]

cd79a_tab <- cd79a_tab[, .(GeneSymbol, GeneID, sampleID, EntityAbbreviation, StudyGroup,
                           seqnames, start, end, deltaPsi, counts, totalCounts, 
                           variant, AbSplice_DNA)]
cd79a_tab <- cd79a_tab[order(deltaPsi), ]

fwrite(cd79a_tab, snakemake@output$cd79a_tab) 

