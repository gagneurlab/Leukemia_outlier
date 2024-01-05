#'---
#' title: tet2_tab
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
#'    - tet2_tab_raw: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/tet2_tab_raw.csv"`'
#'  output:
#'    - tet2_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/tet2_tab.csv"`'
#'  type: script
#'  threads: 1
#'  resources:
#'    - mem_mb: 8000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/tet2_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/tet2_tab.snakemake")




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

anonymization_table_mll <- fread(snakemake@params$anonymization_table_mll)

#+ define paths
tet2_tab_raw <- fread(snakemake@input$tet2_tab_raw)

tet2_tab <- merge(tet2_tab_raw, manuscript_wording_sub, by.x='Diag', by.y='Cohort_during_analysis')
tet2_tab <- merge(tet2_tab, anonymization_table_mll, by='ArrayID')

tet2_tab <- setnames(tet2_tab, 
                     c('AnonamizedID', 'Cohort_abbreviation', 'hgncSymbol', 'geneID'),
                     c('sampleID', 'EntityAbbreviation', 'GeneSymbol', 'GeneID')
)

tet2_tab[, GeneID := sapply(tet2_tab[, GeneID], function(x){strsplit(x, '[.]')[[1]][1]})]

tet2_tab <- tet2_tab[, .(GeneSymbol, GeneID, sampleID, EntityAbbreviation, 
                         pValue, padjust, zScore, rawcounts, expected_counts)]
tet2_tab <- tet2_tab[order(zScore), ]

fwrite(tet2_tab, snakemake@output$tet2_tab) 

