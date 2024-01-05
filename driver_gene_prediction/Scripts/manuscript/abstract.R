#'---
#' title: abstract
#' author: Xueqi Cao
#' wb:
#'  py:
#'   - |
#'    annotations = config["geneAnnotation"]
#'    singleGroupDatasets = config["cohort"]["single_group"]
#'    groupsDatasets = config["cohort"]["groups"]
#'  params:
#'    - projectPath: '`sm config["projectPath"]`'
#'    - annotations: '`sm annotations`'
#'    - gencode: '`sm config["gencode"]`'
#'    - sampAnno: '`sm config["sampleAnnotation"]`'
#'    - singleGroupDatasets: '`sm singleGroupDatasets`'
#'    - groupsDatasets: '`sm groupsDatasets`'
#'    - CGC_cancer_gene_processed: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_processed.tsv"`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'  input:
#'    - outriderRes: '`sm expand(config["outriderDir"] +
#'                    "/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv",
#'                    annotation=annotations, dataset=singleGroupDatasets)`'
#'    - activtionRes: '`sm expand(config["outriderDir"] +
#'                    "/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/res_filter_out.tsv",
#'                    annotation=annotations, dataset=singleGroupDatasets)`'
#'    - fraserRes: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv",
#'                    annotation=annotations, dataset=groupsDatasets)`'
#'    - fraserResJunc: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv",
#'                    annotation=annotations, dataset=groupsDatasets)`'
#'    - abspliceRes: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice.tsv"`'
#'    - abspliceResVar: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice_var.tsv"`'
#'  output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/abstract.html"`'
#'  type: noindex
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/abstract.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/abstract.snakemake")
print("Snakemake saved") 

.libPaths("~/R/4.1/FRASER2.2")
# check with .libPaths() that indeed this directory shows up first
# devtools::install("../FraseR") # here, have the path to the folder in which the FRASER2 code from gitlab is

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(magrittr)
  library(UpSetR)
  library(FRASER)
})

theme_vale <- theme_bw()
options(bitmapType='cairo')




# get parameters
set.seed(2023)

output_dir <- snakemake@params$htmlOutputPath

single_group <- snakemake@params$singleGroupDatasets
groups <- snakemake@params$groupsDatasets

sample_group <- snakemake@params$groupsDatasets
sample_group <- append(sample_group, single_group)

samp_anno <- fread(snakemake@params$sampAnno)
# fwrite(samp_anno[grep(single_group, DROP_GROUP)], "final_cohort_for_stephan.csv")
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()
samp_anno_aml <- samp_anno_exp[DROP_GROUP=='AML',]

cohort_dt <- samp_anno[grep(single_group, DROP_GROUP), .N, by = "Cohort_group"]
cohort_dt <- cohort_dt[order(-N),]
cohort_dt <- rbind(data.table(Cohort_group=single_group, 
                              N=samp_anno[grep(single_group, DROP_GROUP), .N]),
                   cohort_dt)

cgc_cancer_gene <- fread(snakemake@params$CGC_cancer_gene_processed)

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_short := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]
gencode[, isCGC := gene_id_short %in% cgc_cancer_gene[, ENSGid]]



##### a. expression outlier ##### 
or_res <- fread(snakemake@input$outriderRes)
ac_res <- fread(snakemake@input$activtionRes)

exp_res <- rbind(or_res[, .(geneID, sampleID)], ac_res[, .(geneID, sampleID)])
exp_sum <- exp_res[, .N, by='sampleID']
setnames(exp_sum, 'N', 'aberrantBySample')

exp_sum <- merge(exp_sum, samp_anno[grep(single_group, DROP_GROUP), .(ArrayID)],
                 by.x='sampleID', by.y='ArrayID', all.x=TRUE, all.y=TRUE)
exp_sum[is.na(aberrantBySample), aberrantBySample := 0]

#' OUTRIDER2 expression outlier gene 
exp_sum[, sum(aberrantBySample)]

#' OUTRIDER2 expression outlier gene median
exp_sum[, median(aberrantBySample)]




##### b. splicing outlier ##### 
fr_sum <- data.table(sampleID=character(),
                     N=numeric())
for (x in groups) {
  fr_res <- fread(snakemake@input$fraserResJunc[grep(paste0("/", x, "/"), snakemake@input$fraserResJunc)])
  fr_res <- merge(fr_res, gencode[, .(gene_name, isCGC, gene_type)], 
                  by.x='hgncSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  fr_res <- fr_res[gene_type=='protein_coding', ]
  # print(fr_res[, table(gene_type)])
  
  fr_sum <- rbind(fr_sum, fr_res[, .N, by='sampleID'])
  
  cohort_dt[Cohort_group==x, n_fr_cgc := fr_res[, sum(isCGC, na.rm = TRUE)]]
  cohort_dt[Cohort_group==x, n_fr := nrow(fr_res)]
}
setnames(fr_sum, 'N', 'aberrantBySample')

fr_sum <- merge(fr_sum, samp_anno[grep(single_group, DROP_GROUP), .(ArrayID)],
                by.x='sampleID', by.y='ArrayID', all.x=TRUE, all.y=TRUE)
fr_sum[is.na(aberrantBySample), aberrantBySample := 0]

#' FRASER2 splicing outlier junction
fr_sum[, sum(aberrantBySample)]

#' FRASER2 splicing outlier junction median
fr_sum[, median(aberrantBySample)]


fr_sum <- data.table(sampleID=character(),
                     N=numeric())
for (x in groups) {
  fr_res <- fread(snakemake@input$fraserRes[grep(paste0("/", x, "/"), snakemake@input$fraserRes)])
  fr_res <- merge(fr_res, gencode[, .(gene_name, isCGC, gene_type)], 
                  by.x='hgncSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  fr_res <- fr_res[gene_type=='protein_coding', ]
  # print(fr_res[, table(gene_type)])
  
  fr_sum <- rbind(fr_sum, fr_res[, .N, by='sampleID'])
  
  cohort_dt[Cohort_group==x, n_fr_cgc := fr_res[, sum(isCGC, na.rm = TRUE)]]
  cohort_dt[Cohort_group==x, n_fr := nrow(fr_res)]
}
setnames(fr_sum, 'N', 'aberrantBySample')

fr_sum <- merge(fr_sum, samp_anno[grep(single_group, DROP_GROUP), .(ArrayID)],
                by.x='sampleID', by.y='ArrayID', all.x=TRUE, all.y=TRUE)
fr_sum[is.na(aberrantBySample), aberrantBySample := 0]

#' FRASER2 splicing outlier gene
fr_sum[, sum(aberrantBySample)]

#' FRASER2 splicing outlier gene median
fr_sum[, median(aberrantBySample)]




##### c. splicing variant ##### 
absplice_res <- fread(snakemake@input$abspliceResVar)
absplice_res <- absplice_res[AbSplice_DNA >= 0.2, ]
absplice_res <- absplice_res[sampleID %in% samp_anno[grep(single_group, DROP_GROUP), ArrayID], ]

absplice_sum <- absplice_res[, .N, by='sampleID']
setnames(absplice_sum, 'N', 'aberrantBySample')

absplice_sum <- merge(absplice_sum, samp_anno[grep(single_group, DROP_GROUP), .(ArrayID)],
                      by.x='sampleID', by.y='ArrayID', all.x=TRUE, all.y=TRUE)
absplice_sum[is.na(aberrantBySample), aberrantBySample := 0]

#' abSplice splicing variant 
absplice_sum[, sum(aberrantBySample)]

#' abSplice splicing variant median
absplice_sum[, median(aberrantBySample)]

absplice_res <- fread(snakemake@input$abspliceRes)
absplice_res <- absplice_res[AbSplice_DNA >= 0.2, ]
absplice_res <- absplice_res[sampleID %in% samp_anno[grep(single_group, DROP_GROUP), ArrayID], ]

absplice_sum <- absplice_res[, .N, by='sampleID']
setnames(absplice_sum, 'N', 'aberrantBySample')

absplice_sum <- merge(absplice_sum, samp_anno[grep(single_group, DROP_GROUP), .(ArrayID)],
                      by.x='sampleID', by.y='ArrayID', all.x=TRUE, all.y=TRUE)
absplice_sum[is.na(aberrantBySample), aberrantBySample := 0]

#' abSplice splicing gene 
absplice_sum[, sum(aberrantBySample)]

#' abSplice splicing gene median
absplice_sum[, median(aberrantBySample)]

