#'---
#' title: figure_2_vep
#' author: Xueqi Cao, Ata Jadid Ahari
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
#'    - mll_panel_genes: '`sm config["mll_panel_genes"]`'
#'    - AML_variants: '`sm config["AML_variants"]`'
#'    - AML_enformer: '`sm config["AML_enformer"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'    - cnv_MLL_29041: '`sm config["cnv_MLL_29041"]`'
#'    - cnv_MLL_12658: '`sm config["cnv_MLL_12658"]`'
#'    - cnv_MLL_14744: '`sm config["cnv_MLL_14744"]`'
#'    - cnv_MLL_29369: '`sm config["cnv_MLL_29369"]`'
#'    - CGC_leukemia_OCG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_OCG_list.tsv"`'
#'    - CGC_leukemia_TSG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_TSG_list.tsv"`'
#'    - anonymization_table_mll: '`sm config["anonymization_table_mll"]`' 
#'  input:
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
#'    - cnv_tet2: '`sm config["projectPath"] + "/manuscript/cnv/cnv_tet2.csv"`'
#'    - vepRes: '`sm expand(config["vep_path"] +
#'              "/MLL_{dataset}.vep.tsv.gz", dataset=outputDatasets)`'
#'    - enrichmentTabOr: '`sm config["projectPath"] + 
#'                        "/manuscript/figure_2/plot_data/enrichments_TopK_or.csv"`'
#'    - enrichmentTabAc: '`sm config["projectPath"] + 
#'                        "/manuscript/figure_2/plot_data/enrichments_TopK_ac.csv"`'
#'    - outriderRes: '`sm expand(config["outriderDir"] +
#'                    "/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv",
#'                    annotation=annotations, dataset=inputDatasets)`'
#'    - activtionRes: '`sm expand(config["outriderDir"] +
#'                    "/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/res_filter_out.tsv",
#'                    annotation=annotations, dataset=inputDatasets)`'
#'    - ods: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds",
#'             annotation=annotations, dataset=inputDatasets)`'
#'    - ods_filter_out: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods_filter_out.Rds",
#'                        annotation=annotations, dataset=inputDatasets)`'
#'  output:
#'    - tet2_tab_raw: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/tet2_tab_raw.csv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_2_vep.html"`'
#'  type: noindex
#'  threads: 1  
#'  resources:
#'    - mem_mb: 256000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_2_vep.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/figure_2_vep.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggplot2)
  library(ggbio)
  library(ggpubr)
  library(ggrepel)
  library(UpSetR)
  library(grid)
  library(gridExtra)
  library(tidyverse)
  library(RColorBrewer)
  library(rstatix)
  library(OUTRIDER)
})
source("Scripts/manuscript/function.R")
source("Scripts/manuscript/manuscript_theme.R")

options(bitmapType='cairo')




# get parameters
output_dir <- snakemake@params$htmlOutputPath

single_group <- snakemake@params$inputDatasets
sample_group <- snakemake@params$outputDatasets
sample_group <- append(sample_group, single_group)

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()
samp_anno_aml <- samp_anno_exp[DROP_GROUP=='AML',]
single_group_id <- samp_anno_exp[DROP_GROUP==single_group, ArrayID]

cohort_dt <- samp_anno[grep(single_group, DROP_GROUP), .N, by = "Cohort_group"]
cohort_dt <- cohort_dt[order(-N),]
cohort_dt <- rbind(data.table(Cohort_group=single_group, 
                              N=samp_anno[grep(single_group, DROP_GROUP), .N]),
                   cohort_dt)

gencode <- fread(snakemake@params$gencode)

en_or_dt <- fread(snakemake@input$enrichmentTabOr)
en_ac_dt <- fread(snakemake@input$enrichmentTabAc)

predictedConsequence <- fread(snakemake@params$predictedConsequence)
predictedConsequence[, Display_term := tolower(Display_term)]
predictedConsequence[, Display_term := gsub(' ', '_', Display_term)]

or_res <- fread(snakemake@input$outriderRes)
or_res[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]
or_res[, geneID_short := strsplit(geneID, "[.]")[[1]][1], by=rownames(or_res)]

ac_res <- fread(snakemake@input$activtionRes)
ac_res[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]
ac_res[, geneID_short := strsplit(geneID, "[.]")[[1]][1], by=rownames(ac_res)]

leu_ocg <- fread(snakemake@params$CGC_leukemia_OCG_list)
leu_tsg <- fread(snakemake@params$CGC_leukemia_TSG_list)
mll_panel_genes <- fread(snakemake@params$mll_panel_genes, header = FALSE)$V1

AML_enformer <- fread(snakemake@params$AML_enformer)

AML_variants <- read_xlsx(snakemake@params$AML_variants) %>% as.data.table()
AML_variants <- AML_variants[array_id %in% samp_anno_aml[, ArrayID]]
AML_variants <- separate_rows(AML_variants, 'consequence', sep = ",") %>% as.data.table()
AML_variants <- merge(AML_variants, predictedConsequence[, .(Display_term, IMPACT)], 
                      by.x='consequence', by.y='Display_term',
                      all.x=TRUE, all.y=FALSE)
AML_variants[, samp_symbol := paste0(array_id, "-", symbol)]

for (x in snakemake@input$vepRes) {
  vep_res_temp <- fread(x)
  vep_res_temp <- merge(vep_res_temp, gencode[, .(gene_name, gene_type)], 
                        by.x='SYMBOL', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  vep_res_temp <- vep_res_temp[gene_type=='protein_coding', ]
  
  if (exists("vep_res")) {
    vep_res <- rbind(vep_res, vep_res_temp, fill=TRUE)
  } else {
    vep_res <- vep_res_temp
  }
}
vep_res <- vep_res %>% separate('#Uploaded_variation', 
                                c("col1", "array_id", "alt", "ref", "col2"), 
                                sep='__') %>% as.data.table()
vep_res[, samp_symbol := paste0(array_id, "-", SYMBOL)]
vep_res_aml <- vep_res[array_id %in% samp_anno_aml[, ArrayID], ]

for (x in snakemake@input$mll_cnv) {
  cnv_res_temp <- fread(x)
  cnv_res_temp[, LENGTH := END - START]
  
  if (exists("cnv_res")) {
    cnv_res <- rbind(cnv_res, cnv_res_temp, fill=TRUE)
  } else {
    cnv_res <- cnv_res_temp
  }
}
cnv_res[, samp_symbol := paste0(ARRAY_ID, "-", SYMBOL)]
cnv_res <- cnv_res[abs(MEAN_LOG2_COPY_RATIO)>0.3, ]
cnv_res <- cnv_res[LENGTH >= 1000000, ]
cnv_res_aml <- cnv_res[ARRAY_ID %in% samp_anno_aml[, ArrayID], ]

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

ods <- readRDS(snakemake@input$ods)
ods_filter_out <- readRDS(snakemake@input$ods_filter_out)

ods_symbol <- ods
ods_symbol_rownames <- merge(data.table(gene_id = rownames(ods)),
                             gencode[, .(gene_id, gene_name)], by='gene_id', all.x=TRUE, all.y=FALSE)
rownames(ods_symbol) <- ods_symbol_rownames[, gene_name]


cnv_MLL_29041 <- read_delim(snakemake@params$cnv_MLL_29041,
                            delim = "\t",
                            comment = "@",
                            col_types = "ciid") %>% as.data.table()

cnv_MLL_12658 <- read_delim(snakemake@params$cnv_MLL_12658,
                            delim = "\t",
                            comment = "@",
                            col_types = "ciid") %>% as.data.table()

cnv_MLL_14744 <- read_delim(snakemake@params$cnv_MLL_14744,
                            delim = "\t",
                            comment = "@",
                            col_types = "ciid") %>% as.data.table()

cnv_MLL_29369 <- read_delim(snakemake@params$cnv_MLL_29369,
                            delim = "\t",
                            comment = "@",
                            col_types = "ciid") %>% as.data.table()

anonymization_table_mll <- fread(snakemake@params$anonymization_table_mll)




### link or to curated variants in AML panel genes #####
or_res_dn_aml_panel <- or_res[zScore<0 & 
                                sampleID %in% samp_anno_aml[, ArrayID] & 
                                hgncSymbol %in% mll_panel_genes, ]
curated_variant_aml_panel <- AML_variants[(consequence=='stop_hained' | 
                                             consequence=='frameshift_variant' |
                                             consequence=='splice_acceptor_variant' |
                                             consequence=='splice_donor_variant' ) & 
                                            symbol %in% mll_panel_genes, ]

or_curated_variant_aml_panel <- merge(or_res_dn_aml_panel, AML_variants, 
                                      by='samp_symbol', all.x=TRUE, all.y=FALSE)
or_curated_variant_aml_panel[, table(!is.na(array_id))]

or_curated_variant_stop_aml_panel <- merge(or_res_dn_aml_panel, curated_variant_aml_panel, 
                                           by='samp_symbol', all.x=TRUE, all.y=FALSE)
or_curated_variant_stop_aml_panel[, table(!is.na(array_id))]

curated_variant_stop_or_aml_panel <- merge(curated_variant_aml_panel, or_res_dn_aml_panel, 
                                           by='samp_symbol', all.x=TRUE, all.y=FALSE)
curated_variant_stop_or_aml_panel[, table(!is.na(sampleID))]




### link or to vep in AML panel genes #####
vep_res_aml_stop <- vep_res_aml[(Consequence=='stop_hained' | 
                                   Consequence=='frameshift_variant' |
                                   Consequence=='splice_acceptor_variant' |
                                   Consequence=='splice_donor_variant' ) & 
                                  SYMBOL %in% mll_panel_genes, ]

or_vep_aml_panel <- merge(or_res_dn_aml_panel, vep_res_aml, 
                          by='samp_symbol', all.x=TRUE, all.y=FALSE)
or_vep_aml_panel[, table(!is.na(array_id))]

or_vep_stop_aml_panel <- merge(or_res_dn_aml_panel, vep_res_aml_stop, 
                               by='samp_symbol', all.x=TRUE, all.y=FALSE)
or_vep_stop_aml_panel[, table(!is.na(array_id))]

vep_stop_or_aml_panel <- merge(vep_res_aml_stop, or_res_dn_aml_panel, 
                               by='samp_symbol', all.x=TRUE, all.y=FALSE)
vep_stop_or_aml_panel[, table(!is.na(sampleID))]




### link TET2 outlier to cnv, vep, sv, fusion #####
symbol_to_analyze <- 'TET2'
sample_to_analyze <- or_res[hgncSymbol=='TET2', sampleID]
samp_anno[ArrayID %in% sample_to_analyze, table(Cohort_group)]
samp_anno[ArrayID %in% sample_to_analyze, table(Diag)]

samp_anno_tet2_or <- samp_anno[ArrayID %in% sample_to_analyze, 
                               .(ArrayID, Diag)]
samp_anno_tet2_or <- merge(samp_anno_tet2_or, 
                           or_res[hgncSymbol=='TET2', .(hgncSymbol, geneID, sampleID, pValue, padjust, zScore, 
                                                        rawcounts, expected_counts)], 
                           by.x='ArrayID', by.y='sampleID')

fwrite(samp_anno_tet2_or, snakemake@output$tet2_tab_raw)

# sample_to_analyze <- c('MLL_29041', 'MLL_14744')

# cnv
cnv_res[SYMBOL==symbol_to_analyze, table(CALL)]
cnv_samp <- cnv_res[ARRAY_ID %in% sample_to_analyze, ]
cnv_samp[SYMBOL==symbol_to_analyze, ]
#' copy number loss found in MLL_29041

# vep
vep_res_tet2 <- vep_res[SYMBOL==symbol_to_analyze, ]
vep_res_tet2[, table(Consequence)]
vep_res_tet2[array_id %in% sample_to_analyze, ]
#' nothing found for vep

# sv
sv_manta_bcp_all <- fread(grep('BCP_ALL', snakemake@input$mll_manta_sv, value = TRUE))
# sv_manta_bcp_all[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
sv_manta_samp <- sv_manta_bcp_all[array_id %in% sample_to_analyze, ]
sv_manta_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ) & array_id %in% sample_to_analyze, ]

sv_manta_aml <- fread(grep('AML', snakemake@input$mll_manta_sv, value = TRUE))
# sv_manta_aml[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
sv_manta_samp <- sv_manta_aml[array_id %in% sample_to_analyze, ]
sv_manta_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ) & array_id %in% sample_to_analyze, ]
#' nothing found for sv

# fusion
fusion_manta_bcp_all <- fread(grep('BCP_ALL', snakemake@input$mll_manta, value = TRUE))
# fusion_manta_bcp_all[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_manta_samp <- fusion_manta_bcp_all[array_id %in% sample_to_analyze, ]
fusion_manta_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_manta_aml <- fread(grep('AML', snakemake@input$mll_manta, value = TRUE))
# fusion_manta_aml[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_manta_samp <- fusion_manta_aml[array_id %in% sample_to_analyze, ]
fusion_manta_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_manta_MDS <- fread(grep('MDS[.]', snakemake@input$mll_manta, value = TRUE))
# fusion_manta_MDS[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_manta_samp <- fusion_manta_MDS[array_id %in% sample_to_analyze, ]
fusion_manta_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_manta_TL_group <- fread(grep('TL_group', snakemake@input$mll_manta, value = TRUE))
# fusion_manta_TL_group[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_manta_samp <- fusion_manta_TL_group[array_id %in% sample_to_analyze, ]
fusion_manta_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
#' nothing found for manta

fusion_arriba_bcp_all <- fread(grep('BCP_ALL', snakemake@input$mll_arriba, value = TRUE))
# fusion_arriba_bcp_all[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_arriba_samp <- fusion_arriba_bcp_all[array_id %in% sample_to_analyze, ]
fusion_arriba_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_arriba_aml <- fread(grep('AML', snakemake@input$mll_arriba, value = TRUE))
# fusion_arriba_aml[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_arriba_samp <- fusion_arriba_aml[array_id %in% sample_to_analyze, ]
fusion_arriba_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_arriba_MDS <- fread(grep('MDS[.]', snakemake@input$mll_arriba, value = TRUE))
# fusion_arriba_MDS[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_arriba_samp <- fusion_arriba_MDS[array_id %in% sample_to_analyze, ]
fusion_arriba_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_arriba_TL_group <- fread(grep('TL_group', snakemake@input$mll_arriba, value = TRUE))
# fusion_arriba_TL_group[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_arriba_samp <- fusion_arriba_TL_group[array_id %in% sample_to_analyze, ]
fusion_arriba_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
#' nothing found for arriba

fusion_star_fusion_bcp_all <- fread(grep('BCP_ALL', snakemake@input$mll_star_fusion, value = TRUE))
# fusion_star_fusion_bcp_all[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_star_fusion_samp <- fusion_star_fusion_bcp_all[array_id %in% sample_to_analyze, ]
fusion_star_fusion_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_star_fusion_aml <- fread(grep('AML', snakemake@input$mll_star_fusion, value = TRUE))
# fusion_star_fusion_aml[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_star_fusion_samp <- fusion_star_fusion_aml[array_id %in% sample_to_analyze, ]
fusion_star_fusion_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_star_fusion_MDS <- fread(grep('MDS[.]', snakemake@input$mll_star_fusion, value = TRUE))
# fusion_star_fusion_MDS[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_star_fusion_samp <- fusion_star_fusion_MDS[array_id %in% sample_to_analyze, ]
fusion_star_fusion_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]

fusion_star_fusion_TL_group <- fread(grep('TL_group', snakemake@input$mll_star_fusion, value = TRUE))
# fusion_star_fusion_TL_group[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
fusion_star_fusion_samp <- fusion_star_fusion_TL_group[array_id %in% sample_to_analyze, ]
fusion_star_fusion_samp[( gene1==symbol_to_analyze | gene2==symbol_to_analyze ), ]
#' nothing found for star_fusion




### TET2 expression vs. VEP #####
geneID_TET2 <- gencode[gene_name == 'TET2', gene_id]

counts_tet2 <- data.table(sampleID = colnames(ods),
                          raw_count = OUTRIDER::counts(ods, normalized = FALSE)[geneID_TET2,],
                          expected_count = OUTRIDER::normalizationFactors(ods)[geneID_TET2,],
                          autoencoder_norm_count = OUTRIDER::counts(ods, normalized = TRUE)[geneID_TET2,],
                          size_factor = OUTRIDER::sizeFactors(ods))
counts_tet2 <- merge(counts_tet2, 
                     vep_res[SYMBOL==symbol_to_analyze, .(array_id, Consequence, Location, ref, alt)] %>% unique(), 
                     by.x='sampleID', by.y='array_id', all.x=TRUE, all.y=FALSE)

cnv_tet2 <- fread(snakemake@input$cnv_tet2)
counts_tet2 <- merge(counts_tet2, 
                     cnv_tet2[SYMBOL=='TET2', .(ARRAY_ID, MEAN_LOG2_COPY_RATIO, CALL)] %>% unique(), 
                     by.x='sampleID', by.y='ARRAY_ID', all.x=TRUE, all.y=FALSE)
counts_tet2[is.na(MEAN_LOG2_COPY_RATIO), MEAN_LOG2_COPY_RATIO := 0]
counts_tet2[, size_factor_norm_count := (raw_count + 1)/size_factor]
counts_tet2[, outlier := sampleID %in% or_res[hgncSymbol=='TET2', sampleID]]
counts_tet2[, table(Consequence)]

stat.test <- ggpubr::compare_means(formula = size_factor_norm_count ~ Consequence, data = counts_tet2, method = "wilcox.test")

ggplot(counts_tet2, aes(x=Consequence, y=size_factor_norm_count+1)) +
  geom_violin(aes(colour=Consequence, fill=Consequence), alpha=0.3) + 
  geom_boxplot(fill=NA, width=0.25) + 
  stat_pvalue_manual(stat.test, label="p.signif", y.position = 5) + 
  scale_y_log10()

stat.test <- ggpubr::compare_means(formula = autoencoder_norm_count ~ Consequence, data = counts_tet2, method = "wilcox.test")

ggplot(counts_tet2, aes(x=Consequence, y=autoencoder_norm_count+1)) +
  geom_violin(aes(colour=Consequence, fill=Consequence), alpha=0.3) + 
  geom_boxplot(fill=NA, width=0.25) + 
  stat_pvalue_manual(stat.test, label="p.signif", y.position = 4.5) + 
  scale_y_log10()

stat.test <- ggpubr::compare_means(formula = expected_count ~ Consequence, data = counts_tet2, method = "wilcox.test")

ggplot(counts_tet2, aes(x=Consequence, y=expected_count+1)) +
  geom_violin(aes(colour=Consequence, fill=Consequence), alpha=0.3) + 
  geom_boxplot(fill=NA, width=0.25) + 
  stat_pvalue_manual(stat.test, label="p.signif", y.position = 4.5) + 
  scale_y_log10()

ggplot(counts_tet2, aes(x=autoencoder_norm_count+1, y=size_factor_norm_count+1)) +
  geom_point(aes(color=outlier), alpha=0.3) +
  geom_smooth(method='lm', formula= y~x) + 
  scale_x_log10() + 
  scale_y_log10()

ggplot(counts_tet2, aes(x=size_factor_norm_count+1, y=raw_count+1)) +
  geom_point(aes(color=outlier), alpha=0.3) +
  geom_smooth(method='lm', formula= y~x) + 
  scale_x_log10() + 
  scale_y_log10()

ggplot(counts_tet2, aes(x=expected_count+1, y=raw_count+1)) +
  geom_point(aes(color=outlier), alpha=0.3) +
  geom_smooth(method='lm', formula= y~x) + 
  scale_x_log10() + 
  scale_y_log10()

ggplot(counts_tet2, aes(x=expected_count+1, y=autoencoder_norm_count+1)) +
  geom_point(aes(color=outlier), alpha=0.3) +
  geom_smooth(method='lm', formula= y~x) + 
  scale_x_log10() + 
  scale_y_log10()

ggplot(counts_tet2, aes(x=MEAN_LOG2_COPY_RATIO, y=size_factor_norm_count+1)) +
  geom_point(aes(color=outlier), alpha=0.3) +
  geom_vline(xintercept=c(-0.3, 0.3), linetype="dashed") + 
  geom_smooth(method='lm', formula= y~x) + 
  scale_y_log10()

ggplot() +
  geom_point(data=counts_tet2[is.na(Consequence), ], 
             aes(x=autoencoder_norm_count+1, y=size_factor_norm_count+1, color=Consequence), alpha=0.3) +
  geom_point(data=counts_tet2[!is.na(Consequence), ], 
             aes(x=autoencoder_norm_count+1, y=size_factor_norm_count+1, color=Consequence), alpha=0.7) +
  geom_smooth(data=counts_tet2, aes(x=autoencoder_norm_count+1, y=size_factor_norm_count+1), 
              method='lm', formula= y~x) + 
  scale_x_log10() + 
  scale_y_log10()
