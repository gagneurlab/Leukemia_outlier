#'---
#' title: resource_tab
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
#'    - outrider_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/outrider_tab.csv"`'
#'    - activation_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/activation_tab.csv"`'
#'    - fraser_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/fraser_tab.csv"`'
#'    - absplice_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/absplice_tab.csv"`'
#'    - intogen_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/intogen_tab.csv"`'
#'    - prediction_all_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/prediction_all_tab.csv"`'
#'  output:
#'    - resource_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/resource_tab.csv"`'
#'    - col_anno: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/resource_diag/col_anno.csv"`'
#'    - resource_melt_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/resource_melt_tab.csv"`'              
#'  type: script
#'  threads: 1
#'  resources:
#'    - mem_mb: 8000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/resource_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/resource_tab.snakemake")




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
gencode <- gencode[gene_type=='protein_coding', ]
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := sapply(gencode[, gene_id_unique], function(x){strsplit(x, "[.]")[[1]][1]})]

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))
manuscript_wording[, Diag := Cohort_during_analysis]
manuscript_wording[, Diag := gsub(" ", "_", Diag)]
manuscript_wording[, Diag := gsub("/", "_", Diag)]
manuscript_wording[, Diag := gsub("-", "_", Diag)]
manuscript_wording[, Study_group_underscore := Study_group]
manuscript_wording[, Study_group_underscore := gsub(" ", "_", Study_group_underscore)]
manuscript_wording[, Study_group_underscore := gsub("/", "_", Study_group_underscore)]
manuscript_wording[, Study_group_underscore := gsub("-", "_", Study_group_underscore)]
manuscript_wording[, Cohort_abbr_filename := Cohort_abbreviation]
manuscript_wording[, Cohort_abbr_filename := gsub(" ", "_", Cohort_abbr_filename)]
manuscript_wording[, Cohort_abbr_filename := gsub("/", "_", Cohort_abbr_filename)]
manuscript_wording[, Cohort_abbr_filename := gsub("-", "_", Cohort_abbr_filename)]

sampAnno <- fread(snakemake@params$sampAnno)
sampAnno <- sampAnno[grep(single_group, DROP_GROUP)]
sampAnno[, Cohort_during_analysis := Diag]
sampAnno[, Diag := gsub(" ", "_", Diag)]
sampAnno[, Diag := gsub("/", "_", Diag)]
sampAnno[, Diag := gsub("-", "_", Diag)]
diag_vec <- sampAnno[, unique(Diag)]

outrider_tab <- fread(snakemake@input$outrider_tab)
activation_tab <- fread(snakemake@input$activation_tab)
fraser_tab <- fread(snakemake@input$fraser_tab)
absplice_tab <- fread(snakemake@input$absplice_tab)
intogen_tab <- fread(snakemake@input$intogen_tab)
prediction_all_tab <- fread(snakemake@input$prediction_all_tab)

##### old version #####
resource_tab <- merge(outrider_tab, activation_tab, by='gene_id', all.x=TRUE, all.y=TRUE)
resource_tab <- merge(resource_tab, fraser_tab, by='gene_id', all.x=TRUE, all.y=TRUE)
resource_tab <- merge(resource_tab, absplice_tab, by='gene_id', all.x=TRUE, all.y=TRUE)
resource_tab[is.na(resource_tab)] <- 0
resource_tab <- merge(gencode[, .(gene_id, gene_name)], resource_tab, by='gene_id', all.x=TRUE, all.y=FALSE)
resource_tab <- merge(resource_tab, intogen_tab, by='gene_name', all.x=TRUE, all.y=FALSE)

setnames(resource_tab, 
         c("gene_id", "gene_name"),
         c("GeneID", "GeneSymbol")
)
colnames(resource_tab) <- gsub("fr-zScore", "fr-deltaPsi", colnames(resource_tab))

fwrite(resource_tab, snakemake@output$resource_tab) 

resource_table_diag_dir <- dirname(snakemake@output$col_anno)

for (i in diag_vec) {
  diag_col <- grep(paste0(i, "-"), colnames(resource_tab), value = TRUE)
  diag_col <- c("GeneID", "GeneSymbol", diag_col)
  resource_tab_diag <- resource_tab[, ..diag_col]
  cohort_abbr_eng <- manuscript_wording[Diag==i, Cohort_abbr_filename]
  colnames(resource_tab_diag) <- gsub(i, cohort_abbr_eng, colnames(resource_tab_diag))
  fwrite(resource_tab_diag, paste0(resource_table_diag_dir, "/resource_tab_", cohort_abbr_eng, ".csv")) 
}

col_anno <- data.table(Colname = grep(paste0("AML-"), colnames(resource_tab), value = TRUE)[c(-1, -2)])
col_anno <- col_anno %>% separate(Colname, c('Entity', 'Method', 'MetricAbbreviation'), 
                                  sep = "-", remove = FALSE) %>% as.data.table()

col_anno[Method=='or', Method := "OUTRIDER"]
col_anno[Method=='ac', Method := "NB-act"]
col_anno[Method=='fr', Method := "FRASER2"]
col_anno[Method=='ab', Method := "AbSplice"]
col_anno[Method=='hotmaps', Method := "HotMAPS"]
col_anno[Method=='clustl', Method := "OncodriveCLUSTL"]
col_anno[Method=='smregions', Method := "smRegions"]
col_anno[Method=='fml', Method := "OncodriveFML"]
col_anno[Method=='mutpanning', Method := "MutPanning"]
col_anno[Method=='dndscv', Method := "dNdScv"]
col_anno[Method=='cbase', Method := "CBaSE"]

col_anno[MetricAbbreviation=='cutoff=0.01_up', Metric := "Number of samples predicted when filtering for padjust<0.01 and zScore>0"]
col_anno[MetricAbbreviation=='cutoff=0.05_up', Metric := "Number of samples predicted when filtering for padjust<0.05 and zScore>0"]
col_anno[MetricAbbreviation=='cutoff=0.1_up', Metric := "Number of samples predicted when filtering for padjust<0.1 and zScore>0"]
col_anno[MetricAbbreviation=='cutoff=0.01_dn', Metric := "Number of samples predicted when filtering for padjust<0.01 and zScore<0"]
col_anno[MetricAbbreviation=='cutoff=0.05_dn', Metric := "Number of samples predicted when filtering for padjust<0.05 and zScore<0"]
col_anno[MetricAbbreviation=='cutoff=0.1_dn', Metric := "Number of samples predicted when filtering for padjust<0.1 and zScore<0"]
col_anno[MetricAbbreviation=='zScore_2_up', Metric := "Number of samples predicted when filtering for zScore>2"]
col_anno[MetricAbbreviation=='zScore_4_up', Metric := "Number of samples predicted when filtering for zScore>4"]
col_anno[MetricAbbreviation=='zScore_6_up', Metric := "Number of samples predicted when filtering for zScore>6"]
col_anno[MetricAbbreviation=='zScore_2_dn', Metric := "Number of samples predicted when filtering for zScore<2"]
col_anno[MetricAbbreviation=='zScore_4_dn', Metric := "Number of samples predicted when filtering for zScore<4"]
col_anno[MetricAbbreviation=='zScore_6_dn', Metric := "Number of samples predicted when filtering for zScore<6"]
col_anno[MetricAbbreviation=='deltaPsi_0.1_up', Metric := "Number of samples predicted when filtering for deltaPsi>0.1"]
col_anno[MetricAbbreviation=='deltaPsi_0.2_up', Metric := "Number of samples predicted when filtering for deltaPsi>0.2"]
col_anno[MetricAbbreviation=='deltaPsi_0.3_up', Metric := "Number of samples predicted when filtering for deltaPsi>0.3"]
col_anno[MetricAbbreviation=='deltaPsi_0.1_dn', Metric := "Number of samples predicted when filtering for deltaPsi<0.1"]
col_anno[MetricAbbreviation=='deltaPsi_0.2_dn', Metric := "Number of samples predicted when filtering for deltaPsi<0.2"]
col_anno[MetricAbbreviation=='deltaPsi_0.3_dn', Metric := "Number of samples predicted when filtering for deltaPsi<0.3"]

col_anno[MetricAbbreviation=='max_0.2', Metric := "Max prediction value when filtering for AbSplice_DNA>0.2"]
col_anno[MetricAbbreviation=='mean_0.2', Metric := "Mean prediction value when filtering for AbSplice_DNA>0.2"]
col_anno[MetricAbbreviation=='freq_0.2', Metric := "Number of samples predicted when filtering for AbSplice_DNA>0.2"]
col_anno[MetricAbbreviation=='max_0.05', Metric := "Max prediction value when filtering for AbSplice_DNA>0.05"]
col_anno[MetricAbbreviation=='mean_0.05', Metric := "Mean prediction value when filtering for AbSplice_DNA>0.05"]
col_anno[MetricAbbreviation=='freq_0.05', Metric := "Number of samples predicted when filtering for AbSplice_DNA>0.05"]
col_anno[MetricAbbreviation=='max_0.01', Metric := "Max prediction value when filtering for AbSplice_DNA>0.01"]
col_anno[MetricAbbreviation=='mean_0.01', Metric := "Mean prediction value when filtering for AbSplice_DNA>0.01"]
col_anno[MetricAbbreviation=='freq_0.01', Metric := "Number of samples predicted when filtering for AbSplice_DNA>0.01"]

col_anno[MetricAbbreviation=='q_value', Metric := "q_value from HotMAPS if available"]
col_anno[MetricAbbreviation=='SCORE', Metric := "SCORE from OncodriveCLUSTL if available"]
col_anno[MetricAbbreviation=='U', Metric := "U from smRegions if available"]
col_anno[MetricAbbreviation=='Q_VALUE', Metric := "Q_VALUE from OncodriveFML if available"]
col_anno[MetricAbbreviation=='FDR', Metric := "FDR from MutPanning if available"]
col_anno[MetricAbbreviation=='qallsubs_cv', Metric := "qallsubs_cv from dNdScv if available"]
col_anno[MetricAbbreviation=='q_pos', Metric := "q_pos from CBaSE if available"]

fwrite(col_anno, snakemake@output$col_anno)


##### new aggregated version #####
resource_melt_tab <- melt(resource_tab, id.vars=c('GeneSymbol', 'GeneID'))

resource_melt_tab_1 <- resource_melt_tab[grep('-or-|-ac-|-fr-|-ab-', variable), ]
resource_melt_tab_2 <- resource_melt_tab[-grep('-or-|-ac-|-fr-|-ab-', variable), ]
resource_melt_tab <- rbind(resource_melt_tab_1[value !=0, ],
                           resource_melt_tab_2[!is.na(value), ])

resource_melt_tab <- resource_melt_tab %>% separate(variable, c('Entity', 'Method', 'MetricAbbreviation'), 
                                                    sep = "-", remove = FALSE) %>% as.data.table()

resource_melt_tab[Method=='or', Method := "OUTRIDER"]
resource_melt_tab[Method=='ac', Method := "NB-act"]
resource_melt_tab[Method=='fr', Method := "FRASER2"]
resource_melt_tab[Method=='ab', Method := "AbSplice"]
resource_melt_tab[Method=='hotmaps', Method := "HotMAPS"]
resource_melt_tab[Method=='clustl', Method := "OncodriveCLUSTL"]
resource_melt_tab[Method=='smregions', Method := "smRegions"]
resource_melt_tab[Method=='fml', Method := "OncodriveFML"]
resource_melt_tab[Method=='mutpanning', Method := "MutPanning"]
resource_melt_tab[Method=='dndscv', Method := "dNdScv"]
resource_melt_tab[Method=='cbase', Method := "CBaSE"]

resource_melt_tab[MetricAbbreviation=='cutoff=0.01_up', Metric := "Number of samples predicted when filtering for padjust<0.01 and zScore>0"]
resource_melt_tab[MetricAbbreviation=='cutoff=0.05_up', Metric := "Number of samples predicted when filtering for padjust<0.05 and zScore>0"]
resource_melt_tab[MetricAbbreviation=='cutoff=0.1_up', Metric := "Number of samples predicted when filtering for padjust<0.1 and zScore>0"]
resource_melt_tab[MetricAbbreviation=='cutoff=0.01_dn', Metric := "Number of samples predicted when filtering for padjust<0.01 and zScore<0"]
resource_melt_tab[MetricAbbreviation=='cutoff=0.05_dn', Metric := "Number of samples predicted when filtering for padjust<0.05 and zScore<0"]
resource_melt_tab[MetricAbbreviation=='cutoff=0.1_dn', Metric := "Number of samples predicted when filtering for padjust<0.1 and zScore<0"]
resource_melt_tab[MetricAbbreviation=='zScore_2_up', Metric := "Number of samples predicted when filtering for zScore>2"]
resource_melt_tab[MetricAbbreviation=='zScore_4_up', Metric := "Number of samples predicted when filtering for zScore>4"]
resource_melt_tab[MetricAbbreviation=='zScore_6_up', Metric := "Number of samples predicted when filtering for zScore>6"]
resource_melt_tab[MetricAbbreviation=='zScore_2_dn', Metric := "Number of samples predicted when filtering for zScore<2"]
resource_melt_tab[MetricAbbreviation=='zScore_4_dn', Metric := "Number of samples predicted when filtering for zScore<4"]
resource_melt_tab[MetricAbbreviation=='zScore_6_dn', Metric := "Number of samples predicted when filtering for zScore<6"]
resource_melt_tab[MetricAbbreviation=='deltaPsi_0.1_up', Metric := "Number of samples predicted when filtering for deltaPsi>0.1"]
resource_melt_tab[MetricAbbreviation=='deltaPsi_0.2_up', Metric := "Number of samples predicted when filtering for deltaPsi>0.2"]
resource_melt_tab[MetricAbbreviation=='deltaPsi_0.3_up', Metric := "Number of samples predicted when filtering for deltaPsi>0.3"]
resource_melt_tab[MetricAbbreviation=='deltaPsi_0.1_dn', Metric := "Number of samples predicted when filtering for deltaPsi<0.1"]
resource_melt_tab[MetricAbbreviation=='deltaPsi_0.2_dn', Metric := "Number of samples predicted when filtering for deltaPsi<0.2"]
resource_melt_tab[MetricAbbreviation=='deltaPsi_0.3_dn', Metric := "Number of samples predicted when filtering for deltaPsi<0.3"]

resource_melt_tab[MetricAbbreviation=='max_0.2', Metric := "Max prediction value when filtering for AbSplice_DNA>0.2"]
resource_melt_tab[MetricAbbreviation=='mean_0.2', Metric := "Mean prediction value when filtering for AbSplice_DNA>0.2"]
resource_melt_tab[MetricAbbreviation=='freq_0.2', Metric := "Number of samples predicted when filtering for AbSplice_DNA>0.2"]
resource_melt_tab[MetricAbbreviation=='max_0.05', Metric := "Max prediction value when filtering for AbSplice_DNA>0.05"]
resource_melt_tab[MetricAbbreviation=='mean_0.05', Metric := "Mean prediction value when filtering for AbSplice_DNA>0.05"]
resource_melt_tab[MetricAbbreviation=='freq_0.05', Metric := "Number of samples predicted when filtering for AbSplice_DNA>0.05"]
resource_melt_tab[MetricAbbreviation=='max_0.01', Metric := "Max prediction value when filtering for AbSplice_DNA>0.01"]
resource_melt_tab[MetricAbbreviation=='mean_0.01', Metric := "Mean prediction value when filtering for AbSplice_DNA>0.01"]
resource_melt_tab[MetricAbbreviation=='freq_0.01', Metric := "Number of samples predicted when filtering for AbSplice_DNA>0.01"]

resource_melt_tab[MetricAbbreviation=='q_value', Metric := "q_value from HotMAPS if available"]
resource_melt_tab[MetricAbbreviation=='SCORE', Metric := "SCORE from OncodriveCLUSTL if available"]
resource_melt_tab[MetricAbbreviation=='U', Metric := "U from smRegions if available"]
resource_melt_tab[MetricAbbreviation=='Q_VALUE', Metric := "Q_VALUE from OncodriveFML if available"]
resource_melt_tab[MetricAbbreviation=='FDR', Metric := "FDR from MutPanning if available"]
resource_melt_tab[MetricAbbreviation=='qallsubs_cv', Metric := "qallsubs_cv from dNdScv if available"]
resource_melt_tab[MetricAbbreviation=='q_pos', Metric := "q_pos from CBaSE if available"]


resource_melt_tab_1 <- merge(resource_melt_tab[Method %in% c("OUTRIDER", "NB-act", "FRASER2", "AbSplice"), ], 
                             manuscript_wording[, .(Cohort_abbreviation, Diag)],
                             by.x='Entity', by.y='Diag')

resource_melt_tab_2 <- merge(resource_melt_tab[! Method %in% c("OUTRIDER", "NB-act", "FRASER2", "AbSplice"), ], 
                             manuscript_wording[, .(Study_group, Study_group_underscore)] %>% unique(),
                             by.x='Entity', by.y='Study_group_underscore')

resource_melt_tab <- rbind(resource_melt_tab_1, resource_melt_tab_2, fill=TRUE)

resource_melt_tab <- resource_melt_tab[, .(GeneSymbol, GeneID, Cohort_abbreviation, Study_group, Method, Metric, value)]
setnames(resource_melt_tab, 
         c('Cohort_abbreviation', 'Study_group', 'value'), 
         c('Entity', 'Study group', 'Value'))

resource_melt_tab <- merge(prediction_all_tab[, .(GeneID, GeneSymbol)], 
                           resource_melt_tab, 
                           by=c('GeneID', 'GeneSymbol'), sort=FALSE)

fwrite(resource_melt_tab, snakemake@output$resource_melt_tab)
