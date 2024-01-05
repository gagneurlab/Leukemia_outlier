#'---
#' title: association_tab
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
#'    - diagFisher: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/diag_fisher.tsv"`'
#'  output:
#'    - association_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/association_tab.csv"`'
#'  type: script
#'  threads: 1
#'  resources:
#'    - mem_mb: 8000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/association_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/association_tab.snakemake")




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


#+ define paths
diag_fisher <- fread(snakemake@input$diagFisher)
diag_fisher_sig <- diag_fisher[p_greater_adjust < 0.05, ]
diag_fisher_sig <- full_join(diag_fisher_sig, manuscript_wording_sub, by =c("Diag" = "Cohort_during_analysis"))
diag_fisher_sig <- subset(diag_fisher_sig, select = -Diag )
diag_fisher_sig <- diag_fisher_sig[, c(1:3, 10, 4:9)] 

setnames(diag_fisher_sig, 
         c('gene_id', 'Cohort_abbreviation', 'oddsr', 'p_greater', 'p_greater_adjust', 'geneSymbol', 'input_res'),
         c('GeneID', 'DiseaseEntity', 'OddsRatio', 'Pvalue', 'FDR', 'GeneSymbol', 'Method')
         )

diag_fisher_sig[, table(Method)]
diag_fisher_sig[Method=='absplice', Method:='AbSplice']
diag_fisher_sig[Method=='ac', Method:='Activation']
diag_fisher_sig[Method=='fr', Method:='FRASER']
diag_fisher_sig[Method=='or_dn', Method:='Underexpression']
diag_fisher_sig[Method=='or_up', Method:='Overexpression']
diag_fisher_sig[, table(Method)]

diag_fisher_sig <- diag_fisher_sig[, .(Method, GeneID, GeneSymbol, DiseaseEntity, OddsRatio, Pvalue, FDR, 
                    isCGC, RoleinCancer, `TumourTypes(Somatic)`)]

fwrite(diag_fisher_sig, snakemake@output$association_tab) 

