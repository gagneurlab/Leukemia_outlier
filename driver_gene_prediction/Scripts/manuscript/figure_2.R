#'---
#' title: figure_2
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
#'    - mll_arriba: '`sm expand(config["projectPath"] + "/manuscript/arriba/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_star_fusion: '`sm expand(config["projectPath"] + "/manuscript/star_fusion/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_manta: '`sm expand(config["projectPath"] + "/manuscript/manta/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_manta_sv: '`sm expand(config["projectPath"] + "/manuscript/manta_sv/{dataset}.tsv",
#'                      dataset=outputDatasets)`'
#'    - cnv_tet2: '`sm config["projectPath"] + "/manuscript/cnv/cnv_tet2.csv"`'
#'    - vepRes: '`sm expand(config["intogenDir"] +
#'              "/steps/vep/MLL_WGS_MLL_{dataset}.tsv.gz", dataset=[x.upper() for x in outputDatasets])`'
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
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_2.html"`'
#'  type: noindex
#'  threads: 1  
#'  resources:
#'    - mem_mb: 32000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_2.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/figure_2.snakemake")
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




### Figure #### 
#### a. OUTRIDER, leu13, leukemia tsg #### 
#' OUTRIDER number of input genes & samples
dim(ods)
#' OUTRIDER number of genes detected with outliers
or_res[, unique(geneID)] %>% length()
#' OUTRIDER number of underexpression outliers
or_res[zScore<0, .N]
#' OUTRIDER number of underexpression outliers per sample - median
or_res_dn <- or_res[zScore<0, ]
or_dn_samp_freq <- data.table(sampleID = single_group_id)
or_dn_samp_freq <- merge(or_dn_samp_freq, or_res_dn[, .N, by='sampleID'], all=TRUE)
or_dn_samp_freq[is.na(N), N := 0]
or_dn_samp_freq[, sample_rank := rank(N, ties.method = 'random')]
median(or_dn_samp_freq[, N])

exp_gene <- rownames(ods)
exp_gene <- sapply(exp_gene, function(x){strsplit(x, "[.]")[[1]][1]})

# all 
or_dn <- or_res[zScore < 0, ]
or_dn_freq <- or_dn[, .N, by='geneID_short']
or_dn_gene <- or_dn_freq[, geneID_short]
or_dn_gene_0 <- exp_gene[!exp_gene %in% or_dn_freq[, geneID_short]]
or_dn_gene_1 <- or_dn_freq[N==1, geneID_short]
or_dn_gene_2 <- or_dn_freq[N>1 & N<5, geneID_short]
or_dn_gene_5 <- or_dn_freq[N>=5, geneID_short]

#' enrichment overall
fisher_test(or_dn_gene, exp_gene, leu_tsg[, ENSGid])$estimate
fisher_test(or_dn_gene, exp_gene, leu_tsg[, ENSGid])$p.value


or_dn_all_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(or_dn_gene_0, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(or_dn_gene_1, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(or_dn_gene_2, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(or_dn_gene_5, exp_gene, leu_tsg[, ENSGid])$estimate
  ),
  p_val = c(
    fisher_test(or_dn_gene_0, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(or_dn_gene_1, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(or_dn_gene_2, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(or_dn_gene_5, exp_gene, leu_tsg[, ENSGid])$p.value
  ),
  total = c(
    length(or_dn_gene_0),
    length(or_dn_gene_1),
    length(or_dn_gene_2),
    length(or_dn_gene_5)
  ))

or_dn_all_en <- add_significance(or_dn_all_en, "p_val")
or_dn_all_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
or_dn_all_en[, cutoff := 'All']

# top 3
or_dn <- merge(or_dn, or_dn[, .(geneID_short, rank(padjust)), by='sampleID'], 
               by = c('geneID_short', 'sampleID'))
setnames(or_dn, 'V2', 'padjustRankInSample')
or_dn_top3 <- or_dn[padjustRankInSample <= 3, ]
or_dn_top3_freq <- or_dn_top3[, .N, by='geneID_short']
or_dn_top3_gene_0 <- exp_gene[!exp_gene %in% or_dn_top3_freq[, geneID_short]]
or_dn_top3_gene_1 <- or_dn_top3_freq[N==1, geneID_short]
or_dn_top3_gene_2 <- or_dn_top3_freq[N>1 & N<5, geneID_short]
or_dn_top3_gene_5 <- or_dn_top3_freq[N>=5, geneID_short]

or_dn_top3_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(or_dn_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(or_dn_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(or_dn_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(or_dn_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$estimate
  ),
  p_val = c(
    fisher_test(or_dn_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(or_dn_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(or_dn_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(or_dn_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$p.value
  ),
  total = c(
    length(or_dn_top3_gene_0),
    length(or_dn_top3_gene_1),
    length(or_dn_top3_gene_2),
    length(or_dn_top3_gene_5)
  ))

or_dn_top3_en <- add_significance(or_dn_top3_en, "p_val")
or_dn_top3_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
or_dn_top3_en[, cutoff := 'At most three']

or_dn_en <- rbind(or_dn_all_en, or_dn_top3_en)
or_dn_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=rownames(or_dn_en)]
or_dn_en[total < 1000, total_label := as.character(total)]
or_dn_en[category=='0', total_pos := 1]
or_dn_en[category=='1', total_pos := 2]
or_dn_en[category=='2-4', total_pos := 3]
or_dn_en[category=="\u22655", total_pos := 4]
or_dn_en[cutoff=='All', total_pos := total_pos-0.15]
or_dn_en[cutoff=='At most three', total_pos := total_pos+0.15]

signif_height <- 5
total_height <- 8.5

p_a_raw <- ggplot(or_dn_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_text(data = or_dn_en[cutoff=='All'], aes(y=signif_height, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = or_dn_en[cutoff=='At most three'], aes(y=signif_height, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=or_dn_en[, total_pos], y=total_height, label=or_dn_en[, total_label], size = 3) +
  annotate("text", x=0.25, y=total_height, label='n =', size = 3) +
  scale_y_log10(
    breaks = c(0.3, 0.5, 1, 3, 5), 
    minor_breaks = c(c(3:9)/10, c(1:5))
  ) +
  coord_cartesian(ylim=c(0.3, 6), xlim=c(1, 4), clip="off")

# p_a_raw

p_a <- p_a_raw +
  xlab('Samples per gene') +
  ylab('Enrichment for hematologic\ntumor supressor genes (Odds ratio)') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey')) +
  guides(fill=guide_legend('Underexpression outliers (OUTRIDER)', 
                           title.position = "left", reverse = TRUE)) + 
  theme_vale + 
  theme(
    legend.position = c(0.5, 1.25),
    legend.direction = "vertical",
    plot.margin = margin(67, 14, 14, 20, "points"), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

# p_a




#### b. OUTRIDER, leu13, leukemia ocg #### 
#' OUTRIDER number of overexpression outliers
or_res[zScore>0, .N]
#' OUTRIDER number of overexpression outliers per sample - median
or_res_up <- or_res[zScore>0, ]
or_up_samp_freq <- data.table(sampleID = single_group_id)
or_up_samp_freq <- merge(or_up_samp_freq, or_res_up[, .N, by='sampleID'], all=TRUE)
or_up_samp_freq[is.na(N), N := 0]
or_up_samp_freq[, sample_rank := rank(N, ties.method = 'random')]
median(or_up_samp_freq[, N])

exp_gene <- rownames(ods)
exp_gene <- sapply(exp_gene, function(x){strsplit(x, "[.]")[[1]][1]})

# all 
or_up <- or_res[zScore > 0, ]
or_up_freq <- or_up[, .N, by='geneID_short']
or_up_gene <- or_up_freq[, geneID_short]
or_up_gene_0 <- exp_gene[!exp_gene %in% or_up_freq[, geneID_short]]
or_up_gene_1 <- or_up_freq[N==1, geneID_short]
or_up_gene_2 <- or_up_freq[N>1 & N<5, geneID_short]
or_up_gene_5 <- or_up_freq[N>=5, geneID_short]

#' enrichment overall
fisher_test(or_up_gene, exp_gene, leu_ocg[, ENSGid])$estimate
fisher_test(or_up_gene, exp_gene, leu_ocg[, ENSGid])$p.value %>% format(scientific = TRUE)

or_up_all_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(or_up_gene_0, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(or_up_gene_1, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(or_up_gene_2, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(or_up_gene_5, exp_gene, leu_ocg[, ENSGid])$estimate
  ),
  p_val = c(
    fisher_test(or_up_gene_0, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(or_up_gene_1, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(or_up_gene_2, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(or_up_gene_5, exp_gene, leu_ocg[, ENSGid])$p.value
  ),
  total = c(
    length(or_up_gene_0),
    length(or_up_gene_1),
    length(or_up_gene_2),
    length(or_up_gene_5)
  ))

or_up_all_en <- add_significance(or_up_all_en, "p_val")
or_up_all_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
or_up_all_en[, cutoff := 'All']

# top 3
or_up <- merge(or_up, or_up[, .(geneID_short, rank(padjust)), by='sampleID'], 
               by = c('geneID_short', 'sampleID'))
setnames(or_up, 'V2', 'padjustRankInSample')
or_up_top3 <- or_up[padjustRankInSample <= 3, ]
or_up_top3_freq <- or_up_top3[, .N, by='geneID_short']
or_up_top3_gene_0 <- exp_gene[!exp_gene %in% or_up_top3_freq[, geneID_short]]
or_up_top3_gene_1 <- or_up_top3_freq[N==1, geneID_short]
or_up_top3_gene_2 <- or_up_top3_freq[N>1 & N<5, geneID_short]
or_up_top3_gene_5 <- or_up_top3_freq[N>=5, geneID_short]

or_up_top3_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(or_up_top3_gene_0, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(or_up_top3_gene_1, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(or_up_top3_gene_2, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(or_up_top3_gene_5, exp_gene, leu_ocg[, ENSGid])$estimate
  ),
  p_val = c(
    fisher_test(or_up_top3_gene_0, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(or_up_top3_gene_1, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(or_up_top3_gene_2, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(or_up_top3_gene_5, exp_gene, leu_ocg[, ENSGid])$p.value
  ),
  total = c(
    length(or_up_top3_gene_0),
    length(or_up_top3_gene_1),
    length(or_up_top3_gene_2),
    length(or_up_top3_gene_5)
  ))

or_up_top3_en <- add_significance(or_up_top3_en, "p_val")
or_up_top3_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
or_up_top3_en[, cutoff := 'At most three']

or_up_en <- rbind(or_up_all_en, or_up_top3_en)
or_up_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=rownames(or_up_en)]
or_up_en[total < 1000, total_label := as.character(total)]
or_up_en[category=='0', total_pos := 1]
or_up_en[category=='1', total_pos := 2]
or_up_en[category=='2-4', total_pos := 3]
or_up_en[category=="\u22655", total_pos := 4]
or_up_en[cutoff=='All', total_pos := total_pos-0.15]
or_up_en[cutoff=='At most three', total_pos := total_pos+0.15]

signif_height <- 5
total_height <- 8.5

p_b_raw <- ggplot(or_up_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_text(data = or_up_en[cutoff=='All'], aes(y=signif_height, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = or_up_en[cutoff=='At most three'], aes(y=signif_height, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=or_up_en[, total_pos], y=total_height, label=or_up_en[, total_label], size = 3) +
  annotate("text", x=0.25, y=total_height, label='n =', size = 3) +
  scale_y_log10(
    breaks = c(0.3, 0.5, 1, 3, 5), 
    minor_breaks = c(c(3:9)/10, c(1:5))
  ) +
  coord_cartesian(ylim=c(0.3, 6), xlim=c(1, 4), clip="off")

# p_b_raw

p_b <- p_b_raw +
  xlab('Samples per gene') +
  ylab('Enrichment for hematologic\noncogenes (Odds ratio)') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey')) +
  guides(fill=guide_legend('Overexpression outliers (OUTRIDER)', 
                           title.position = "left", reverse = TRUE)) + 
  theme_vale + 
  theme(
    legend.position = c(0.5, 1.25),
    legend.direction = "vertical",
    plot.margin = margin(67, 14, 14, 20, "points"), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

# p_b




#### c. Activation, leu13, leukemia ocg ####
#' Activation number of input genes & samples
dim(ods_filter_out)
#' Activation number of genes detected with outliers
ac_res[, unique(geneID)] %>% length()
#' Activation number of activation outliers
ac_res[, .N]
#' Activation number of activation outliers per sample - median
ac_samp_freq <- data.table(sampleID = single_group_id)
ac_samp_freq <- merge(ac_samp_freq, ac_res[, .N, by='sampleID'], all=TRUE)
ac_samp_freq[is.na(N), N := 0]
ac_samp_freq[, sample_rank := rank(N, ties.method = 'random')]
median(ac_samp_freq[, N])

exp_gene <- rownames(ods_filter_out)
exp_gene <- sapply(exp_gene, function(x){strsplit(x, "[.]")[[1]][1]})

# all 
ac_freq <- ac_res[, .N, by='geneID_short']
ac_gene <- ac_freq[, geneID_short]
ac_gene_0 <- exp_gene[!exp_gene %in% ac_freq[, geneID_short]]
ac_gene_1 <- ac_freq[N==1, geneID_short]
ac_gene_2 <- ac_freq[N>1 & N<5, geneID_short]
ac_gene_5 <- ac_freq[N>=5, geneID_short]

#' enrichment overall
fisher_test(ac_gene, exp_gene, leu_ocg[, ENSGid])$estimate
fisher_test(ac_gene, exp_gene, leu_ocg[, ENSGid])$p.value %>% format(scientific = TRUE)

ac_all_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(ac_gene_0, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(ac_gene_1, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(ac_gene_2, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(ac_gene_5, exp_gene, leu_ocg[, ENSGid])$estimate
  ),
  p_val = c(
    fisher_test(ac_gene_0, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(ac_gene_1, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(ac_gene_2, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(ac_gene_5, exp_gene, leu_ocg[, ENSGid])$p.value
  ),
  total = c(
    length(ac_gene_0),
    length(ac_gene_1),
    length(ac_gene_2),
    length(ac_gene_5)
  ))

ac_all_en <- add_significance(ac_all_en, "p_val")
ac_all_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
ac_all_en[, cutoff := 'All']

# top 3
ac_res <- merge(ac_res, ac_res[, .(geneID_short, rank(padjust)), by='sampleID'], 
                by = c('geneID_short', 'sampleID'))
setnames(ac_res, 'V2', 'padjustRankInSample')
ac_top3 <- ac_res[padjustRankInSample <= 3, ]
ac_top3_freq <- ac_top3[, .N, by='geneID_short']
ac_top3_gene_0 <- exp_gene[!exp_gene %in% ac_top3_freq[, geneID_short]]
ac_top3_gene_1 <- ac_top3_freq[N==1, geneID_short]
ac_top3_gene_2 <- ac_top3_freq[N>1 & N<5, geneID_short]
ac_top3_gene_5 <- ac_top3_freq[N>=5, geneID_short]

ac_top3_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(ac_top3_gene_0, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(ac_top3_gene_1, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(ac_top3_gene_2, exp_gene, leu_ocg[, ENSGid])$estimate,
    fisher_test(ac_top3_gene_5, exp_gene, leu_ocg[, ENSGid])$estimate
  ),
  p_val = c(
    fisher_test(ac_top3_gene_0, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(ac_top3_gene_1, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(ac_top3_gene_2, exp_gene, leu_ocg[, ENSGid])$p.value,
    fisher_test(ac_top3_gene_5, exp_gene, leu_ocg[, ENSGid])$p.value
  ),
  total = c(
    length(ac_top3_gene_0),
    length(ac_top3_gene_1),
    length(ac_top3_gene_2),
    length(ac_top3_gene_5)
  ))

ac_top3_en <- add_significance(ac_top3_en, "p_val")
ac_top3_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
ac_top3_en[, cutoff := 'At most three']

ac_en <- rbind(ac_all_en, ac_top3_en) 
ac_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=rownames(ac_en)]
ac_en[total < 1000, total_label := as.character(total)]
ac_en[category=='0', total_pos := 1]
ac_en[category=='1', total_pos := 2]
ac_en[category=='2-4', total_pos := 3]
ac_en[category=="\u22655", total_pos := 4]
ac_en[cutoff=='All', total_pos := total_pos-0.15]
ac_en[cutoff=='At most three', total_pos := total_pos+0.15]

signif_height <- 5
total_height <- 8.5

p_c_raw <- ggplot(ac_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_text(data = ac_en[cutoff=='All'], aes(y=signif_height, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = ac_en[cutoff=='At most three'], aes(y=signif_height, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=ac_en[, total_pos], y=total_height, label=ac_en[, total_label], size = 3) +
  annotate("text", x=0.25, y=total_height, label='n =', size = 3) +
  scale_y_log10(
    breaks = c(0.3, 0.5, 1, 3, 5), 
    minor_breaks = c(c(3:9)/10, c(1:5))
  ) +
  coord_cartesian(ylim=c(0.3, 6), xlim=c(1, 4), clip="off")

# p_c_raw

p_c <- p_c_raw +
  xlab('Samples per gene') +
  ylab('Enrichment for hematologic\noncogenes (Odds ratio)') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey')) +
  guides(fill=guide_legend('Activation outliers (NB-act)', 
                           title.position = "left", reverse = TRUE)) + 
  theme_vale + 
  theme(
    legend.position = c(0.5, 1.25),
    legend.direction = "vertical",
    plot.margin = margin(67, 14, 14, 20, "points"), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

# p_c




#### d. samp strat MLL panel genes #####
AML_variants_unique <- unique(AML_variants[, .(array_id, symbol)])
AML_variants_unique_panel <- AML_variants_unique[symbol %in% mll_panel_genes]

cnv_aml_panel <- cnv_res_aml[SYMBOL %in% mll_panel_genes, ]
cnv_aml_panel[, LENGTH := END - START]

fusion_aml <- rbind(
  fread(grep('AML', snakemake@input$mll_manta, value = TRUE))[, .(array_id, gene1, gene2)],
  fread(grep('AML', snakemake@input$mll_arriba, value = TRUE))[, .(array_id, gene1, gene2)],
  fread(grep('AML', snakemake@input$mll_star_fusion, value = TRUE))[, .(array_id, gene1, gene2)]
)
fusion_aml <- melt(fusion_aml, id.vars = 'array_id', value.name = "symbol")
fusion_aml_panel <- fusion_aml[symbol %in% mll_panel_genes, ]

fusion_cnv_union <- rbind(unique(fusion_aml_panel[, .(array_id, symbol)]),
                          unique(cnv_aml_panel[, .(ARRAY_ID, SYMBOL)]),
                          use.names=FALSE)
fusion_cnv_union <- unique(fusion_cnv_union)

fusion_cnv_variants_union <- rbind(unique(fusion_aml_panel[, .(array_id, symbol)]),
                                   unique(cnv_aml_panel[, .(ARRAY_ID, SYMBOL)]),
                                   AML_variants_unique_panel, 
                                   use.names=FALSE)
fusion_cnv_variants_union <- unique(fusion_cnv_variants_union)

fusion_cnv_variants_or_ac_union <- rbind(unique(fusion_aml_panel[, .(array_id, symbol)]),
                                         unique(cnv_aml_panel[, .(ARRAY_ID, SYMBOL)]),
                                         AML_variants_unique_panel, 
                                         or_res[sampleID %in% samp_anno_aml$ArrayID, .(sampleID, hgncSymbol)], 
                                         ac_res[sampleID %in% samp_anno_aml$ArrayID, .(sampleID, hgncSymbol)], 
                                         use.names=FALSE)
fusion_cnv_variants_or_ac_union <- unique(fusion_cnv_variants_or_ac_union[symbol %in% mll_panel_genes])

fusion_cnv_union <- as.data.table(table(sort(table(fusion_cnv_union[, array_id]))))
colnames(fusion_cnv_union) <- c("n", "fusion_cnv_union")

fusion_cnv_variants_union <- as.data.table(table(sort(table(fusion_cnv_variants_union[, array_id]))))
colnames(fusion_cnv_variants_union) <- c("n", "fusion_cnv_variants_union")

fusion_cnv_variants_or_ac_union <- as.data.table(table(sort(table(fusion_cnv_variants_or_ac_union[, array_id]))))
colnames(fusion_cnv_variants_or_ac_union) <- c("n", "fusion_cnv_variants_or_ac_union")

df_merge <- merge(fusion_cnv_union, fusion_cnv_variants_union, by="n", all.x = TRUE, all.y=TRUE)
df_merge <- merge(df_merge, fusion_cnv_variants_or_ac_union, by="n", all.x = TRUE, all.y=TRUE)

df_merge$n <- as.integer(df_merge$n)
df_merge <- df_merge[order(df_merge$n),]
num_of_aml_samples <- length(unique(samp_anno_aml$ArrayID))

df_merge <- rbindlist(list(data.table("n" = 0, 
                                      "fusion_cnv_union" = num_of_aml_samples - sum(df_merge$fusion_cnv_union,na.rm=TRUE ),
                                      "fusion_cnv_variants_union" = num_of_aml_samples - sum(df_merge$fusion_cnv_variants_union,na.rm=TRUE),
                                      "fusion_cnv_variants_or_ac_union" = num_of_aml_samples - sum(df_merge$fusion_cnv_variants_or_ac_union,na.rm=TRUE )),
                           df_merge))
dt_merge <- df_merge %>% as.data.table()
#' number of samples
colSums(dt_merge[n>=5, ], na.rm = TRUE)

subset_n <- 10
df_merge[n==subset_n, "fusion_cnv_variants_or_ac_union"] <- sum(df_merge[n >= subset_n, "fusion_cnv_variants_or_ac_union"], na.rm = TRUE)
df_merge[n==subset_n, "fusion_cnv_union"] <- sum(df_merge[n >= subset_n, "fusion_cnv_union"], na.rm = TRUE)
df_merge[n==subset_n, "fusion_cnv_variants_union"] <- sum(df_merge[n >= subset_n, "fusion_cnv_variants_union"], na.rm = TRUE)
df_merge <- df_merge[n <= subset_n , ]
df_merge$n <- as.character(df_merge$n)
df_merge[n==subset_n, "n"] <- paste0("\u2265", subset_n)
df_merge$n <- factor(df_merge$n, levels = df_merge$n)

df_merge_melted <- melt(df_merge, id.vars = c("n"),
                        measure.vars = c("fusion_cnv_union", "fusion_cnv_variants_union", "fusion_cnv_variants_or_ac_union" ))

p_d_raw <- ggplot(df_merge_melted, aes(factor(n), value, fill=variable))+
  geom_bar(stat="identity", position="dodge") 

p_d <- p_d_raw +
  xlab("Predicted aberrant hematologic panel genes per sample") + 
  ylab("Samples in AML") +
  scale_fill_manual('Stratified based on',
                    labels=c('fusion_cnv_union'='Fusions or CNV ',
                             'fusion_cnv_variants_union'='Fusions or CNV \nor Curated variants',
                             'fusion_cnv_variants_or_ac_union'='Fusions or CNV \nor Curated variants \nor Expression outliers'
                    ),
                    values=c(RColorBrewer::brewer.pal(12, "Paired")[1],
                             RColorBrewer::brewer.pal(12, "Paired")[3],
                             RColorBrewer::brewer.pal(12, "Paired")[5]),) +
  theme_vale + 
  theme(
    legend.position = 'top',
    axis.title.x = element_text(hjust=0.9)
  ) +
  guides(fill=guide_legend(title.position = "top"))  

# p_d




#### s4. validated by cnv percentage #####
or_dn <- or_res[zScore < 0, .(sampleID, hgncSymbol, samp_symbol)] 
or_up <- or_res[zScore > 0, .(sampleID, hgncSymbol, samp_symbol)] 
ac <- ac_res[, .(sampleID, hgncSymbol, samp_symbol)]

cnv_up <- cnv_res[CALL=="+",]
cnv_dn <- cnv_res[CALL=="-",]

or_dn <- or_res[zScore < 0, .(sampleID, hgncSymbol, samp_symbol)] 
or_up <- or_res[zScore > 0, .(sampleID, hgncSymbol, samp_symbol)] 
ac <- ac_res[, .(sampleID, hgncSymbol, samp_symbol)]

cnv_up <- cnv_res[CALL=="+",]
cnv_dn <- cnv_res[CALL=="-",]

#' number of overlap
or_dn[, samp_symbol %in% cnv_dn[, samp_symbol]] %>% table()
nrow(or_dn)
or_up[, samp_symbol %in% cnv_up[, samp_symbol]] %>% table()
nrow(or_up)
ac[, samp_symbol %in% cnv_up[, samp_symbol]] %>% table()
nrow(ac)

prop_dt <- data.table(
  or_dn[, samp_symbol %in% cnv_dn[, samp_symbol]] %>% table() %>% prop.table() %>% as.numeric(),
  or_up[, samp_symbol %in% cnv_up[, samp_symbol]] %>% table() %>% prop.table() %>% as.numeric(),
  ac[, samp_symbol %in% cnv_up[, samp_symbol]] %>% table() %>% prop.table() %>% as.numeric()
) %>% t()

colnames(prop_dt) <- c(FALSE, TRUE)
rownames(prop_dt) <- c('Underexpression\noutlier', 'Overexpression\noutlier', 'Activation\noutlier')
prop_dt <- melt(prop_dt) %>% as.data.table()
colnames(prop_dt) <- c('Category', 'Validated', 'Percentage')

num_dt <- data.table(
  or_dn[, samp_symbol %in% cnv_dn[, samp_symbol]] %>% table() %>% as.numeric(),
  or_up[, samp_symbol %in% cnv_up[, samp_symbol]] %>% table() %>% as.numeric(),
  ac[, samp_symbol %in% cnv_up[, samp_symbol]] %>% table() %>% as.numeric()
) %>% t() %>% as.data.table()

num_dt[, Sum := sum(V1, V2), by=rownames(num_dt)]
colnames(num_dt) <- c('Num_not_val', 'Num_val', 'Sum')
num_dt[, Category := factor(c('Underexpression\noutlier', 'Overexpression\noutlier', 'Activation\noutlier'))]

prop_dt <- merge(prop_dt, num_dt, by='Category')
prop_dt[, text_label := paste0( label_percent(accuracy=0.1)(Percentage), "\n",
                                "(", Num_val, "/", Sum, ")")]
prop_dt[Validated==FALSE, text_label := '']

p_s4_raw <- ggplot(prop_dt, aes(x=Category, y=Percentage, fill=Validated)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=text_label), vjust=-0.5)

# p_s4_raw

p_s4 <- p_s4_raw +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits = c(0,1)) +
  scale_fill_manual(values=c('grey', RColorBrewer::brewer.pal(8, "Paired")[4]), name='Validated by CNV') +
  # ylab(expression(paste(bold("Ratio of "), bolditalic("LRP1B"), bold("-activated samples"))))+
  xlab('') +
  theme_vale 

# p_s4




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




#### e. TET2 raw counts vs. predicted count #####
#' TET2 number of outliers
or_res[hgncSymbol=='TET2', .N]
or_res[hgncSymbol=='TET2', sort(zScore)]
or_res_tet2 <- or_res[hgncSymbol=='TET2', .(zScore, foldChange, sampleID)]

plotExpressionRank(ods_symbol, symbol_to_analyze, basePlot = TRUE)

plotAberrantPerSample(ods_symbol, title='')

p_e_raw <- plotExpectedVsObservedCounts(ods_symbol, symbol_to_analyze, basePlot = TRUE)
p_e_raw$layers[[3]] <- NULL

anonym_id_dict <- c(
  'MLL_29041' = anonymization_table_mll[anonymization_table_mll$ArrayID == "MLL_29041", AnonamizedID],
  'MLL_14744' = anonymization_table_mll[anonymization_table_mll$ArrayID == "MLL_14744", AnonamizedID]
)

p_e <- p_e_raw +
  geom_text_repel(aes(label=ifelse(sampleID %in% c("MLL_14744", "MLL_29041"), anonym_id_dict[sampleID], '')),
                  nudge_x = 0.5,
                  color="firebrick", fontface = "bold") +
  ggtitle('TET2') +
  xlab('Expected counts + 1') +
  ylab('Raw counts + 1') +
  scale_color_manual(guide = guide_legend(reverse=TRUE), 
                     values = c('darkgrey', 'firebrick'),
                     labels = c("Non-outlier", "Outlier")) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x)(c(1, 1e6)),
    labels = trans_format("log10", math_format(10^.x)),
    minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9))
  ) +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x)(c(1, 1e6)),
    labels = trans_format("log10", math_format(10^.x)),
    minor_breaks = rep(1:9, 4)*(10^rep(0:3, each = 9))
  ) +
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = c(.05, .95),
    legend.justification = c("left", "top"),
    legend.margin = margin(0, 8, 4, 4),
    plot.title = element_text(face = "bold.italic"),
    legend.title = element_blank()
  )

# p_e




#### f. TET2 OUTRIDER zScore posi ctrl #####
sample_to_curate <- 'MLL_29041'
or_res[sampleID==sample_to_curate, .N]
or_res[sampleID==sample_to_curate, table(zScore>0)]

p_s8_raw <- plotVolcano(ods_symbol, sample_to_curate, basePlot = TRUE)

p_f_raw <- plotChromOUTRIDERzScore(ods_symbol, sample_to_curate, 'zScore', chrom_to_plot = '4')[[1]]
p_f_raw_data <- plotChromOUTRIDERzScore(ods_symbol, sample_to_curate, 'zScore', chrom_to_plot = '4')[[2]]
vline_pos <- p_f_raw_data[aberrant==TRUE, START_mbp]

p_f <- p_f_raw +
  scale_color_manual(guide = guide_legend(reverse=TRUE), 
                     values = c('grey40', 'firebrick'),
                     labels = c("Non-outlier", "Outlier")) +
  geom_text_repel(aes(label=ifelse(aberrant, symbol, '')),
                  nudge_x = 0.1, nudge_y = 0.1, 
                  color="firebrick", fontface = "bold.italic") +
  xlab('Genomic position on chromosome 4 (MB)') +
  ylab('Gene expression z-score') +
  ggtitle( anonym_id_dict[sample_to_curate]) +
  theme_vale  +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = c(.05, .95),
    legend.justification = c("left", "top"),
    legend.margin = margin(0, 8, 4, 4),
    legend.title = element_blank()
  )

# p_f




#### g. TET2 CNV Cumulative posi ctrl #####
cnv_MLL_29041[, COPY_RATIO := 2^LOG2_COPY_RATIO]

p_g_raw <- plotChromCNV(cnv_MLL_29041, chrom_to_plot = '4')

p_g <- p_g_raw +
  geom_vline(xintercept=vline_pos, alpha=0.5, color='firebrick') + 
  scale_y_log10(
    limits = c(0.1, 10),
    breaks = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10),
    minor_breaks = c(c(1:9)/10, c(1:10))
  ) +
  geom_hline(yintercept=2^0.3, linetype="dashed", color = "black") +
  geom_hline(yintercept=2^(-0.3), linetype="dashed", color = "black") +
  xlab('Genomic position on chromosome 4 (MB)') +
  ylab('Copy ratio') +
  ggtitle(anonym_id_dict[sample_to_curate]) +
  theme_vale

# p_g




#### h. TET2 OUTRIDER zScore posi ctrl #####
# sample_to_curate <- 'MLL_12658'
# sample_to_curate <- 'MLL_29369'
sample_to_curate <- 'MLL_14744'
or_res[sampleID==sample_to_curate, .N]
or_res[sampleID==sample_to_curate, table(zScore>0)]

p_s9_raw <- plotVolcano(ods_symbol, sample_to_curate, basePlot = TRUE)

p_h_raw <- plotChromOUTRIDERzScore(ods_symbol, sample_to_curate, 'zScore', chrom_to_plot = '4')[[1]]
p_h_raw_data <- plotChromOUTRIDERzScore(ods_symbol, sample_to_curate, 'zScore', chrom_to_plot = '4')[[2]]
vline_pos <- p_h_raw_data[aberrant==TRUE, START_mbp]

p_h <- p_h_raw +
  scale_color_manual(guide = guide_legend(reverse=TRUE), 
                     values = c('grey40', 'firebrick'),
                     labels = c("Non-outlier", "Outlier")) +
  geom_text_repel(aes(label=ifelse(aberrant, symbol, '')),
                  nudge_x = 0.1, nudge_y = 0.1, 
                  color="firebrick", fontface = "bold.italic") +
  xlab('Genomic position on chromosome 4 (MB)') +
  ylab('Gene expression z-score') +
  ggtitle(anonym_id_dict[sample_to_curate]) +
  theme_vale  +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = c(.05, .95),
    legend.justification = c("left", "top"),
    legend.margin = margin(0, 8, 4, 4),
    legend.title = element_blank()
  )

# p_h




#### i. TET2 CNV Cumulative posi ctrl #####
cnv_MLL_14744[, COPY_RATIO := 2^LOG2_COPY_RATIO]

p_i_raw <- plotChromCNV(cnv_MLL_14744, chrom_to_plot = '4')

p_i <- p_i_raw +
  geom_vline(xintercept=vline_pos, alpha=0.5, color='firebrick') + 
  scale_y_log10(
    limits = c(0.1, 10),
    breaks = c(0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10),
    minor_breaks = c(c(1:9)/10, c(1:10))
  ) +
  geom_hline(yintercept=2^0.3, linetype="dashed", color = "black") +
  geom_hline(yintercept=2^(-0.3), linetype="dashed", color = "black") +
  xlab('Genomic position on chromosome 4 (MB)') +
  ylab('Copy ratio') +
  ggtitle(anonym_id_dict[sample_to_curate]) +
  theme_vale

# p_i




### arrange plot #### 
#' ### raw
#+ plot p2 raw, fig.width=10, fig.height=14
p2_raw_ab <- grid.arrange(p_a_raw, p_b_raw, nrow = 1)
p2_raw_ce <- grid.arrange(p_c_raw, p_e_raw, nrow = 1)
p2_raw_fghi <- grid.arrange(p_f_raw, p_h_raw, p_g_raw, p_i_raw, nrow = 2)
p2_raw <- grid.arrange(p2_raw_ab, p2_raw_ce, p2_raw_fghi, nrow = 3, heights = c(1,1,2))

p2_raw


#' ### annotated 
#+ plot p2, fig.width=10, fig.height=14
p2_ab <- grid.arrange(p_a, p_b, nrow = 1)
p2_ce <- grid.arrange(p_c, p_e, nrow = 1)
p2_fghi <- grid.arrange(p_f, p_h, p_g, p_i, nrow = 2)
p2 <- grid.arrange(p2_ab, p2_ce, p2_fghi, nrow = 3, heights = c(1,1,2))

p2


# png(paste0(output_dir, "/figure_2_poster.png"), width = 11*0.8, height = 7*0.8, units = "in", res = 600)
# grid.arrange(p_f, p_h, p_g, p_i, nrow = 2)
# dev.off() 


pdf(paste0(output_dir, "/figure_2.pdf"), width = 10, height=14)
p2_ab <- grid.arrange(p_a, p_b, nrow = 1)
p2_ce <- grid.arrange(p_c, p_e, nrow = 1)
p2_fghi <- grid.arrange(p_f, p_h, p_g, p_i, nrow = 2)
p2 <- grid.arrange(p2_ab, p2_ce, p2_fghi, nrow = 3, heights = c(1,1,2))
p2

grid.text("A", x = 0.01, y=0.99, gp=gpar(fontsize=12, fontface='bold'))
grid.text("B", x = 0.51, y=0.99, gp=gpar(fontsize=12, fontface='bold'))
grid.text("C", x = 0.01, y=0.74, gp=gpar(fontsize=12, fontface='bold'))
grid.text("D", x = 0.51, y=0.74, gp=gpar(fontsize=12, fontface='bold'))
grid.text("E", x = 0.01, y=0.50, gp=gpar(fontsize=12, fontface='bold'))
grid.text("F", x = 0.51, y=0.50, gp=gpar(fontsize=12, fontface='bold'))
grid.text("G", x = 0.01, y=0.24, gp=gpar(fontsize=12, fontface='bold'))
grid.text("H", x = 0.51, y=0.24, gp=gpar(fontsize=12, fontface='bold'))
dev.off() 


png(paste0(output_dir, "/figure_2.png"), width = 10, height = 14, units = "in", res = 600)
p2_ab <- grid.arrange(p_a, p_b, nrow = 1)
p2_ce <- grid.arrange(p_c, p_e, nrow = 1)
p2_fghi <- grid.arrange(p_f, p_h, p_g, p_i, nrow = 2)
p2 <- grid.arrange(p2_ab, p2_ce, p2_fghi, nrow = 3, heights = c(1,1,2))
p2

grid.text("A", x = 0.01, y=0.99, gp=gpar(fontsize=12, fontface='bold'))
grid.text("B", x = 0.51, y=0.99, gp=gpar(fontsize=12, fontface='bold'))
grid.text("C", x = 0.01, y=0.74, gp=gpar(fontsize=12, fontface='bold'))
grid.text("D", x = 0.51, y=0.74, gp=gpar(fontsize=12, fontface='bold'))
grid.text("E", x = 0.01, y=0.50, gp=gpar(fontsize=12, fontface='bold'))
grid.text("F", x = 0.51, y=0.50, gp=gpar(fontsize=12, fontface='bold'))
grid.text("G", x = 0.01, y=0.24, gp=gpar(fontsize=12, fontface='bold'))
grid.text("H", x = 0.51, y=0.24, gp=gpar(fontsize=12, fontface='bold'))
dev.off() 




#' ## Sup
### Supplement #####
#### s1. all cohorts, or dn, leukemia tsg #####
cancer_driver_list_selected <- 'leukemia.tsg'

en_or_tsg_dt <- en_or_dt[cancer_driver_list==cancer_driver_list_selected, ] 
en_or_tsg_dt[, sample_subset := factor(sample_subset, levels=sample_group)]
en_or_tsg_dt[, k_max := max(k), by='sample_subset'] 
en_or_tsg_dt <- en_or_tsg_dt[k==3|k==k_max , ]
en_or_tsg_dt[, x_lab := as.character(k)]
en_or_tsg_dt[k == k_max, x_lab := 'max']
en_or_tsg_dt[, x_lab := factor(x_lab, levels=c('max', '3'))]

#' ### s1 raw
p_s1_raw <- ggplot(en_or_tsg_dt, aes(x=x_lab, y=enrichment)) +
  geom_bar(stat="identity", width = 0.3, fill='lightgrey') + 
  geom_errorbar(aes(ymin=en05, ymax=en95), width=0.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  
  facet_wrap('sample_subset', nrow=2)  +
  scale_y_log10()

#+ plot s1 raw, fig.width=10, fig.height=3
p_s1_raw

#' ### s1 annotated
p_s1 <- p_s1_raw +
  xlab('') +
  ylab('Fold of enrichment for\nhematologic TSG') +
  theme_vale +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#+ plot s1, fig.width=10, fig.height=3
p_s1




#### s2. all cohorts, or up, leukemia ocg #### 
cancer_driver_list_selected <- 'leukemia.ocg'

en_or_ocg_dt <- en_or_dt[cancer_driver_list==cancer_driver_list_selected, ] 
en_or_ocg_dt[, k_max := max(k), by='sample_subset'] 
en_or_ocg_dt <- en_or_ocg_dt[k==3|k==k_max , ]
en_or_ocg_dt[, x_lab := as.character(k)]
en_or_ocg_dt[k == k_max, x_lab := 'max']
en_or_ocg_dt[, x_lab := factor(x_lab, levels=c('max', '3'))]
en_or_ocg_dt[, sample_subset := factor(sample_subset, levels=sample_group)]

#' ### s2 raw
p_s2_raw <- ggplot(en_or_ocg_dt, aes(x=x_lab, y=enrichment)) +
  geom_bar(stat="identity", width = 0.3, fill='lightgrey') + 
  geom_errorbar(aes(ymin=en05, ymax=en95), width=0.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  
  facet_wrap('sample_subset', nrow=2) +
  scale_y_log10()

#+ plot s2 raw, fig.width=10, fig.height=3
p_s2_raw  

#' ### s2 annotated 
p_s2 <- p_s2_raw +  
  xlab('') +
  ylab('Fold of enrichment for\nhematologic oncogene') +
  theme_vale +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#+ plot s2, fig.width=10, fig.height=3
p_s2




#### s3. all cohorts, ac, leukemia ocg #### 
cancer_driver_list_selected <- 'leukemia.ocg'

en_ac_ocg_dt <- en_ac_dt[cancer_driver_list==cancer_driver_list_selected, ] 
en_ac_ocg_dt[, k_max := max(k), by='sample_subset'] 
en_ac_ocg_dt <- en_ac_ocg_dt[k==3|k==k_max , ]
en_ac_ocg_dt[, x_lab := as.character(k)]
en_ac_ocg_dt[k == k_max, x_lab := 'max']
en_ac_ocg_dt[, x_lab := factor(x_lab, levels=c('max', '3'))]
en_ac_ocg_dt[, sample_subset := factor(sample_subset, levels=sample_group)]

#' ### s3 raw
p_s3_raw <- ggplot(en_ac_ocg_dt, aes(x=x_lab, y=enrichment)) +
  geom_bar(stat="identity", width = 0.3, fill='lightgrey') + 
  geom_errorbar(aes(ymin=en05, ymax=en95), width=0.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  
  facet_wrap('sample_subset', nrow=2) +
  scale_y_log10()

#+ plot s3 raw, fig.width=10, fig.height=3
p_s3_raw  

#' ### s3 annotated 
p_s3 <- p_s3_raw +  
  xlab('') +
  ylab('Fold of enrichment for\nhematologic oncogene') +
  theme_vale +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#+ plot s3, fig.width=10, fig.height=3
p_s3




#### s4. validated by cnv percentage #####
#' ### s4 raw
#+ plot s4 raw, fig.width=6, fig.height=4
p_s4_raw

#' ### s4 annotated
#+ plot s4, fig.width=6, fig.height=4
p_s4




#### s5. sample rank underexpression #####
#' ### s5 raw
#+ plot s5 raw, fig.width=6, fig.height=4
hlines <- quantile(or_dn_samp_freq[, N], c(0.5, 0.9)) + 1

p_s5_raw <- ggplot(or_dn_samp_freq, aes(x=sample_rank, y=N+1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=hlines) + 
  scale_y_log10(limits = c(1, or_dn_samp_freq[, 10^ceiling(log10(max(N)))]+1),
                breaks = c(1, 2, 4, 11, 31, 101, 301, 1001), 
                labels = c(0, 1, 3, 10, 30, 100, 300, 1000)) +
  annotate("text", label=c("Median", "90^th ~ percentile"), 
           x=1, y=hlines*1.2, hjust=0, parse=TRUE)

p_s5_raw

#' ### s5 annotated
#+ plot s5, fig.width=6, fig.height=4
p_s5 <- p_s5_raw +
  xlab('Sample rank') +
  ylab('Underexpressed genes (OUTRIDER)') +
  theme_vale 

p_s5




#### s6. sample rank overexpression #####
#' ### s6 raw
#+ plot s6 raw, fig.width=6, fig.height=4
hlines <- quantile(or_up_samp_freq[, N], c(0.5, 0.9)) + 1

p_s6_raw <- ggplot(or_up_samp_freq, aes(x=sample_rank, y=N+1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=hlines) + 
  scale_y_log10(limits = c(1, or_up_samp_freq[, 10^ceiling(log10(max(N)))]+1),
                breaks = c(1, 2, 4, 11, 31, 101, 301, 1001), 
                labels = c(0, 1, 3, 10, 30, 100, 300, 1000)) +
  annotate("text", label=c("Median", "90^th ~ percentile"), 
           x=1, y=hlines*1.2, hjust=0, parse=TRUE)

p_s6_raw

#' ### s6 annotated
#+ plot s6, fig.width=6, fig.height=4
p_s6 <- p_s6_raw +
  xlab('Sample rank') +
  ylab('Overexpressed genes (OUTRIDER)') +
  theme_vale 

p_s6




#### s7. sample rank activation #####
#' ### s7 raw
#+ plot s7 raw, fig.width=6, fig.height=4
hlines <- quantile(ac_samp_freq[, N], 0.9) + 1

p_s7_raw <- ggplot(ac_samp_freq, aes(x=sample_rank, y=N+1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=hlines) + 
  scale_y_log10(limits = c(1, ac_samp_freq[, 10^ceiling(log10(max(N)))]+1),
                breaks = c(1, 2, 4, 11, 31, 101, 301, 1001), 
                labels = c(0, 1, 3, 10, 30, 100, 300, 1000)) +
  annotate("text", label="90^th ~ percentile", 
           x=1, y=hlines*1.2, hjust=0, parse=TRUE)

p_s7_raw

#' ### s7 annotated
#+ plot s7, fig.width=6, fig.height=4
p_s7 <- p_s7_raw +
  xlab('Sample rank') +
  ylab('Activated genes (NB-act)') +
  theme_vale 

p_s7




#### s8. volcano #####
#' ### s8 raw
#+ plot s8 raw, fig.width=6, fig.height=4
p_s8_raw

#' ### s8 annotated
#+ plot s8, fig.width=6, fig.height=4
sample_to_curate <- "MLL_29041"
p_s8 <- p_s8_raw + 
  scale_color_manual(labels=c('Outlier', 'Non-outlier'), values=c('firebrick', 'gray')) +
  xlab('Z-score') +
  theme_vale +
  guides(color=guide_legend('Expression outlier\n(OUTRIDER)')) + 
  ggtitle(paste0("Volcano plot: ", anonym_id_dict[sample_to_curate]))

p_s8




#### s9. volcano #####
#' ### s9 raw
#+ plot s9 raw, fig.width=6, fig.height=4
p_s9_raw

#' ### s9 annotated
#+ plot s9, fig.width=6, fig.height=4
sample_to_curate <- "MLL_14744"
p_s9 <- p_s9_raw + 
  scale_color_manual(labels=c('Outlier', 'Non-outlier'), values=c('firebrick', 'gray')) +
  xlab('Z-score') +
  theme_vale +
  guides(color=guide_legend('Expression outlier\n(OUTRIDER)')) + 
  ggtitle(paste0("Volcano plot: ", anonym_id_dict[sample_to_curate]))

p_s9



### thesis ####
png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "or_sample_rank.png"), 
    width = 6, height = 8, units = "in", res = 600)
ggarrange(p_s5, p_s6, ncol = 1, labels = c('A', 'B'))
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "ac_sample_rank.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_s7
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "or_enrich.png"), 
    width = 6, height = 8, units = "in", res = 600)
ggarrange(p_a, p_b, ncol = 1, labels = c('A', 'B'))
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "ac_enrich.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_c
dev.off() 

p_s4_1 <- ggplot(prop_dt[Category != 'Activation\noutlier'], aes(x=Category, y=Percentage, fill=Validated)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=text_label), vjust=-0.5) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits = c(0,1)) +
  scale_fill_manual(values=c('grey', RColorBrewer::brewer.pal(8, "Paired")[4]), name='Validated by CNV') +
  # ylab(expression(paste(bold("Ratio of "), bolditalic("LRP1B"), bold("-activated samples"))))+
  xlab('') +
  theme_vale 

# p_s4_1

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "or_cnv.png"), 
    width = 5, height = 4, units = "in", res = 600)
p_s4_1
dev.off() 

p_s4_2 <- ggplot(prop_dt[Category == 'Activation\noutlier'], aes(x=Category, y=Percentage, fill=Validated)) + 
  geom_bar(stat="identity") +
  geom_text(aes(label=text_label), vjust=-0.5) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits = c(0,1)) +
  scale_fill_manual(values=c('grey', RColorBrewer::brewer.pal(8, "Paired")[4]), name='Validated by CNV') +
  # ylab(expression(paste(bold("Ratio of "), bolditalic("LRP1B"), bold("-activated samples"))))+
  xlab('') +
  theme_vale 

# p_s4_2

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "ac_cnv.png"), 
    width = 4, height = 4, units = "in", res = 600)
p_s4_2
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "tet2_count.png"), 
    width = 4, height = 4, units = "in", res = 600)
p_e
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "tet2_vol_1.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_s8
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "tet2_vol_2.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_s9
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "tet2_man_1.png"), 
    width = 6, height = 8, units = "in", res = 600)
ggarrange(p_f, p_g, ncol = 1, labels = c('A', 'B'))
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "tet2_man_2.png"), 
    width = 6, height = 8, units = "in", res = 600)
ggarrange(p_h, p_i, ncol = 1, labels = c('A', 'B'))
dev.off() 
