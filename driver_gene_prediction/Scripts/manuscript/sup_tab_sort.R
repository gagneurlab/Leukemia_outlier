#'---
#' title: sup_table_sort
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
#'    - prediction_all_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/prediction_all_tab.csv"`'
#'    - prediction_diag_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/prediction_diag_tab.csv"`'
#'    - prediction_emb_all_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/prediction_emb_all_tab.csv"`'
#'    - prediction_emb_diag_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/prediction_emb_diag_tab.csv"`'
#'    - association_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/association_tab.csv"`'
#'    - curation_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/curation_tab.csv"`'
#'    - drug_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/sup_table/drug_tab.csv"`'
#'    - n_var_samp_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/n_var_samp_tab.csv"`'
#'    - n_var_gene_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/n_var_gene_tab.csv"`'
#'    - n_var_vep_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/n_var_vep_tab.csv"`'
#'    - fpkm_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/fpkm_tab.csv"`'
#'    - sample_summary_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/sample_summary_tab.csv"`'
#'    - or_dn_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/or_dn_agg_tab.csv"`'
#'    - or_up_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/or_up_agg_tab.csv"`'
#'    - activation_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/activation_agg_tab.csv"`'
#'    - fraser_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/fraser_agg_tab.csv"`'
#'    - absplice_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/absplice_agg_tab.csv"`'
#'    - absplice_ratio_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/absplice_ratio_tab.csv"`'
#'    - fusion_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table/fusion_agg_tab.csv"`'
#'    - resource_melt_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/resource_melt_tab.csv"`' 
#'  output:
#'    - sample_summary_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S1_sample_summary.csv"`'
#'    - n_var_samp_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S2_number_variant_sample.csv"`'
#'    - n_var_gene_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S3_number_variant_entity_gene.csv"`'
#'    - n_var_vep_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S4_number_variant_entity_gene_vep.csv"`'
#'    - or_dn_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S5_underexpression_entity_gene.csv"`'
#'    - or_up_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S6_overexpression_entity_gene.csv"`'
#'    - activation_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S7_activation_entity_gene.csv"`'
#'    - fraser_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S8_fraser_entity_gene.csv"`'
#'    - absplice_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S9_absplice_entity_gene.csv"`'
#'    - absplice_ratio_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S10_absplice_ratio_entity_gene.csv"`'
#'    - prediction_emb_all_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/supplementary_table/S11_prediction_complete_dataset_embedding.csv"`'
#'    - prediction_emb_diag_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/supplementary_table/S12_prediction_study_groups_embedding.csv"`'
#'    - prediction_all_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/supplementary_table/S13_prediction_complete_dataset.csv"`'
#'    - drug_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/supplementary_table/S14_drug.csv"`'
#'    - prediction_diag_tab: '`sm config["projectPath"] + 
#'                           "/manuscript/supplementary_table/S15_prediction_study_groups.csv"`'   
#'    - association_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S16_association.csv"`'
#'    - curation_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_table/S17_curation.csv"`'
#'    - resource_melt_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/Catalog_of_transcriptomic_and_genomic_aberrations_of_24_hematologic_malignancy_entities.csv"`' 
#'    - fpkm_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_file/F2_fpkm_entity_gene.csv"`'
#'    - fusion_agg_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/supplementary_file/F3_fusion_entity_gene.csv"`'
#'  type: script
#'  threads: 1
#'  resources:
#'    - mem_mb: 8000 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/sup_table_sort.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/sup_table_sort.snakemake")




##### process the data #####
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(dplyr)
  library(tidyr)
})


tab_file <- fread(snakemake@input$prediction_all_tab)
fwrite(tab_file, snakemake@output$prediction_all_tab)

tab_file <- fread(snakemake@input$prediction_diag_tab)
fwrite(tab_file, snakemake@output$prediction_diag_tab)

tab_file <- fread(snakemake@input$prediction_emb_all_tab)
fwrite(tab_file, snakemake@output$prediction_emb_all_tab)

tab_file <- fread(snakemake@input$prediction_emb_diag_tab)
fwrite(tab_file, snakemake@output$prediction_emb_diag_tab)

tab_file <- fread(snakemake@input$association_tab)
fwrite(tab_file, snakemake@output$association_tab)

tab_file <- fread(snakemake@input$curation_tab)
fwrite(tab_file, snakemake@output$curation_tab)

tab_file <- fread(snakemake@input$drug_tab)
fwrite(tab_file, snakemake@output$drug_tab)

tab_file <- fread(snakemake@input$n_var_samp_tab)
fwrite(tab_file, snakemake@output$n_var_samp_tab)

tab_file <- fread(snakemake@input$n_var_gene_tab)
fwrite(tab_file, snakemake@output$n_var_gene_tab)

tab_file <- fread(snakemake@input$n_var_vep_tab)
fwrite(tab_file, snakemake@output$n_var_vep_tab)

tab_file <- fread(snakemake@input$fpkm_tab)
fwrite(tab_file, snakemake@output$fpkm_tab)

tab_file <- fread(snakemake@input$sample_summary_tab)
fwrite(tab_file, snakemake@output$sample_summary_tab)

tab_file <- fread(snakemake@input$or_dn_agg_tab)
fwrite(tab_file, snakemake@output$or_dn_agg_tab)

tab_file <- fread(snakemake@input$or_up_agg_tab)
fwrite(tab_file, snakemake@output$or_up_agg_tab)

tab_file <- fread(snakemake@input$activation_agg_tab)
fwrite(tab_file, snakemake@output$activation_agg_tab)

tab_file <- fread(snakemake@input$fraser_agg_tab)
fwrite(tab_file, snakemake@output$fraser_agg_tab)

tab_file <- fread(snakemake@input$absplice_ratio_tab)
fwrite(tab_file, snakemake@output$absplice_ratio_tab)

tab_file <- fread(snakemake@input$absplice_agg_tab)
fwrite(tab_file, snakemake@output$absplice_agg_tab)

tab_file <- fread(snakemake@input$fusion_agg_tab)
fwrite(tab_file, snakemake@output$fusion_agg_tab)

tab_file <- fread(snakemake@input$resource_melt_tab)
fwrite(tab_file, snakemake@output$resource_melt_tab)

