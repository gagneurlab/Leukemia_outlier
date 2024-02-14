#'---
#' title: figure_3
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
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
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
#'    - enrichmentTabFr: '`sm config["projectPath"] + 
#'                        "/manuscript/figure_3/plot_data/enrichments_TopK_fr.csv"`'
#'    - enrichmentTabAb: '`sm config["projectPath"] + 
#'                        "/manuscript/figure_3/plot_data/enrichments_absplice.csv"`'
#'    - vepRes: '`sm expand(config["vep_path"] +
#'              "/MLL_{dataset}.vep.tsv.gz", dataset=outputDatasets)`'
#'    - fraserRes: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv",
#'                    annotation=annotations, dataset=outputDatasets)`'
#'    - fraserResJunc: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv",
#'                    annotation=annotations, dataset=outputDatasets)`'
#'    - fds: '`sm expand(config["fraserDir"] +
#'            "/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/fds-object.RDS",
#'             annotation=annotations, dataset=outputDatasets)`'
#'    - abspliceRes: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice.tsv"`'
#'    - abspliceResVar: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice_var.tsv"`'
#'    - venn_list: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_3/plot_data/venn_list.Rds"`'
#'    - vep_splice: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_splice.tsv"`'
#'    - vep_res_abSplice: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_res_abSplice"`'
#'    - vep_full: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_full.tsv"`'
#'    - prop_mutation_fr: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/prop_mutation_fr.csv"`'
#'    - fr_mut_en: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/fr_mut_en.csv"`'
#'  output:
#'    - cd79a_tab_raw: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/cd79a_tab_raw.csv"`'
#'    - figure_3a: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/figure_3a.csv"`'
#'    - figure_3b: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/figure_3b.csv"`'
#'    - figure_3d: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/figure_3d.csv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_3.html"`'
#'  type: noindex
#'  resources:
#'    - mem_mb: 96000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_3.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_3.snakemake")
print("Snakemake saved") 

.libPaths("~/R/4.1/FRASER2")
# check with .libPaths() that indeed this directory shows up first
# devtools::install("../FraseR") # here, have the path to the folder in which the FRASER2 code from gitlab is

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(magrittr)
  library(UpSetR)
  library(FRASER)
  library(ggvenn)
  library(scales)
  library(rstatix)
  library(rtracklayer)
  library(ggrepel)
})
source("Scripts/manuscript/function.R")
source("Scripts/manuscript/manuscript_theme.R")

options(bitmapType='cairo')




# get parameters
set.seed(2023)

output_dir <- snakemake@params$htmlOutputPath

single_group <- snakemake@params$inputDatasets

sample_group <- snakemake@params$outputDatasets

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()
samp_anno_aml <- samp_anno_exp[DROP_GROUP=='AML',]
samp_anno_14groups <- samp_anno_exp[DROP_GROUP %in% sample_group, ] %>% unique()
single_group_id <- samp_anno_14groups[, ArrayID]


cohort_dt <- samp_anno[grep(single_group, DROP_GROUP), .N, by = "Cohort_group"]
cohort_dt <- cohort_dt[order(-N),]
cohort_dt <- rbind(data.table(Cohort_group=single_group, 
                              N=samp_anno[grep(single_group, DROP_GROUP), .N]),
                   cohort_dt)

gencode <- fread(snakemake@params$gencode)

en_fr_dt <- fread(snakemake@input$enrichmentTabFr)
en_ab_dt <- fread(snakemake@input$enrichmentTabAb)

predictedConsequence <- fread(snakemake@params$predictedConsequence)
predictedConsequence[, Display_term := tolower(Display_term)]
predictedConsequence[, Display_term := gsub(' ', '_', Display_term)]

leu_ocg <- fread(snakemake@params$CGC_leukemia_OCG_list)
leu_tsg <- fread(snakemake@params$CGC_leukemia_TSG_list)
mll_panel_genes <- fread(snakemake@params$mll_panel_genes, header = FALSE)$V1

AML_variants <- read_xlsx(snakemake@params$AML_variants) %>% as.data.table()
AML_variants <- AML_variants[array_id %in% samp_anno_aml[, ArrayID]]
AML_variants <- separate_rows(AML_variants, 'consequence', sep = ",") %>% as.data.table()
AML_variants <- merge(AML_variants, predictedConsequence[, .(Display_term, IMPACT)], 
                      by.x='consequence', by.y='Display_term',
                      all.x=TRUE, all.y=FALSE)
AML_variants[, samp_symbol := paste0(array_id, "-", symbol)]

absplice_res <- fread(snakemake@input$abspliceRes)
absplice_res[, samp_symbol := paste0(sampleID, "-", gene_name)]
absplice_res <- separate(absplice_res, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_res <- separate(absplice_res, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()
absplice_res_sig <- absplice_res[AbSplice_DNA >= 0.2, ]

absplice_res_var <- fread(snakemake@input$abspliceResVar)
absplice_res_var[, samp_symbol := paste0(sampleID, "-", gene_name)]
absplice_res_var <- separate(absplice_res_var, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_res_var <- separate(absplice_res_var, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()

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
fr_res[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]

for (x in snakemake@input$fraserResJunc) {
  fr_res_temp <- fread(x)
  fr_res_temp <- merge(fr_res_temp, gencode[, .(gene_name, gene_type)], 
                       by.x='hgncSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  fr_res_temp <- fr_res_temp[gene_type=='protein_coding', ]
  
  if (exists("fr_res_junc")) {
    fr_res_junc <- rbind(fr_res_junc, fr_res_temp)
  } else {
    fr_res_junc <- fr_res_temp
  }
}
fr_res_junc[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]
fr_res_junc[, junction_id := paste0("Intron ", seqnames, ": ", 
                                    format(start, nsmall=1, big.mark=","), "-", 
                                    format(end, nsmall=1, big.mark=",")), 
            by = rownames(fr_res_junc)]

# vep_res <- fread(snakemake@input$vep_full)
vep_splice <- fread(snakemake@input$vep_splice)
vep_res_abSplice <- fread(snakemake@input$vep_res_abSplice)
vep_lym <- fread(grep('Lym_group', snakemake@input$vepRes, value = TRUE))
vep_lym <- vep_lym %>% unique()
vep_lym <- vep_lym %>% separate('#Uploaded_variation', 
                                c("col1", "array_id", "alt", "ref", "pos"), 
                                sep='__', remove = FALSE) %>% as.data.table()
vep_lym[, samp_symbol := paste0(array_id, "-", SYMBOL)]

anonymization_table_mll <- fread(snakemake@params$anonymization_table_mll)




### Figure #### 
#### a. FRASER, leu13, leukemia tsg #### 
#' FRASER total splice outlier junctions
fr_res_junc[, .N]
#' FRASER number samples junction
fr_res_junc[, length(unique(sampleID))]
#' FRASER number genes 
fr_res_junc[, length(unique(hgncSymbol))]
#' FRASER total splice outlier genes
fr_res[, .N]
#' FRASER number samples genes
fr_res[, length(unique(sampleID))]
#' FRASER number of splicing outlier genes per sample - median
fr_samp_freq <- data.table(sampleID = single_group_id)
fr_samp_freq <- merge(fr_samp_freq, fr_res[, .N, by='sampleID'], all=TRUE)
fr_samp_freq[is.na(N), N := 0]
fr_samp_freq[, sample_rank := rank(N, ties.method = 'random')]
median(fr_samp_freq[, N])

exp_gene_fr <- c()
for (fds_path in snakemake@input$fds) {
  print(fds_path)
  fds_temp <- loadFraserDataSet(file = fds_path)
  exp_gene_fr_temp <- rowData(fds_temp)$hgnc_symbol
  exp_gene_fr <- c(exp_gene_fr, exp_gene_fr_temp)
}
exp_gene_fr <- exp_gene_fr[!is.na(exp_gene_fr)] %>% unique()
exp_gene_fr <- sapply(exp_gene_fr, function(x){strsplit(x, ";")[[1]]}) %>% unlist()
names(exp_gene_fr) <- NULL
exp_gene_fr <- unique(exp_gene_fr)

# all 
fr_freq <- fr_res[, .N, by='hgncSymbol']
fr_gene <- fr_freq[, hgncSymbol]
fr_gene_0 <- exp_gene_fr[!exp_gene_fr %in% fr_freq[, hgncSymbol]]
fr_gene_1 <- fr_freq[N==1, hgncSymbol]
fr_gene_2 <- fr_freq[N>1 & N<5, hgncSymbol]
fr_gene_5 <- fr_freq[N>=5, hgncSymbol]

#' enrichment overall
fisher_test(fr_gene, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate
fisher_test(fr_gene, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value

fr_all_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(fr_gene_0, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate,
    fisher_test(fr_gene_1, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate,
    fisher_test(fr_gene_2, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate,
    fisher_test(fr_gene_5, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate
  ),
  odds_ratio_ci_low = c(
    fisher_test(fr_gene_0, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[1],
    fisher_test(fr_gene_1, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[1],
    fisher_test(fr_gene_2, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[1],
    fisher_test(fr_gene_5, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(fr_gene_0, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[2],
    fisher_test(fr_gene_1, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[2],
    fisher_test(fr_gene_2, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[2],
    fisher_test(fr_gene_5, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[2]
  ),
  p_val = c(
    fisher_test(fr_gene_0, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value,
    fisher_test(fr_gene_1, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value,
    fisher_test(fr_gene_2, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value,
    fisher_test(fr_gene_5, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value
  ),
  total = c(
    length(fr_gene_0),
    length(fr_gene_1),
    length(fr_gene_2),
    length(fr_gene_5)
  ))

fr_all_en <- add_significance(fr_all_en, "p_val")
fr_all_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
fr_all_en[, cutoff := 'All']

# top 3
fr_res <- merge(fr_res, fr_res[, .(hgncSymbol, rank(padjust)), by='sampleID'], 
                by = c('hgncSymbol', 'sampleID'))
setnames(fr_res, 'V2', 'padjustRankInSample')
fr_top3 <- fr_res[padjustRankInSample <= 3, ]
fr_top3_freq <- fr_top3[, .N, by='hgncSymbol']
fr_top3_gene_0 <- exp_gene_fr[!exp_gene_fr %in% fr_top3_freq[, hgncSymbol]]
fr_top3_gene_1 <- fr_top3_freq[N==1, hgncSymbol]
fr_top3_gene_2 <- fr_top3_freq[N>1 & N<5, hgncSymbol]
fr_top3_gene_5 <- fr_top3_freq[N>=5, hgncSymbol]

fr_top3_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(fr_top3_gene_0, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate,
    fisher_test(fr_top3_gene_1, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate,
    fisher_test(fr_top3_gene_2, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate,
    fisher_test(fr_top3_gene_5, exp_gene_fr, leu_tsg[, GeneSymbol])$estimate
  ),
  odds_ratio_ci_low = c(
    fisher_test(fr_top3_gene_0, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[1],
    fisher_test(fr_top3_gene_1, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[1],
    fisher_test(fr_top3_gene_2, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[1],
    fisher_test(fr_top3_gene_5, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(fr_top3_gene_0, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[2],
    fisher_test(fr_top3_gene_1, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[2],
    fisher_test(fr_top3_gene_2, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[2],
    fisher_test(fr_top3_gene_5, exp_gene_fr, leu_tsg[, GeneSymbol])$conf.int[2]
  ),
  p_val = c(
    fisher_test(fr_top3_gene_0, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value,
    fisher_test(fr_top3_gene_1, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value,
    fisher_test(fr_top3_gene_2, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value,
    fisher_test(fr_top3_gene_5, exp_gene_fr, leu_tsg[, GeneSymbol])$p.value
  ),
  total = c(
    length(fr_top3_gene_0),
    length(fr_top3_gene_1),
    length(fr_top3_gene_2),
    length(fr_top3_gene_5)
  ))

fr_top3_en <- add_significance(fr_top3_en, "p_val")
fr_top3_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
fr_top3_en[, cutoff := 'At most three']

fr_en <- rbind(fr_all_en, fr_top3_en)
fr_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=list(rownames(fr_en))]
fr_en[total < 1000, total_label := as.character(total)]
fr_en[category=='0', total_pos := 1]
fr_en[category=='1', total_pos := 2]
fr_en[category=='2-4', total_pos := 3]
fr_en[category=="\u22655", total_pos := 4]
fr_en[cutoff=='All', total_pos := total_pos-0.2]
fr_en[cutoff=='At most three', total_pos := total_pos+0.2]

fr_en[p_val.signif=='ns', odds_ratio_ci_low := NA]
fr_en[p_val.signif=='ns', odds_ratio_ci_up := NA]
fr_en[, max(odds_ratio_ci_up, na.rm=TRUE)]

signif_height_a <- 25
total_height_a <- 50

p_a_raw <- ggplot(fr_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=odds_ratio_ci_low, ymax=odds_ratio_ci_up), 
                width=0.1, position=position_dodge(width=0.7)) +
  geom_text(data = fr_en[cutoff=='All'], aes(y=signif_height_a, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = fr_en[cutoff=='At most three'], aes(y=signif_height_a, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=fr_en[, total_pos], y=total_height_a, label=fr_en[, total_label], size = 2) +
  annotate("text", x=0.25, y=total_height_a, label='n =', size = 2) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  scale_y_log10(
    breaks = c(0.1, 0.3, 1, 3, 10),
    minor_breaks = c(c(1:9)/10, c(1:10))
  ) +
  coord_cartesian(ylim=c(0.1, 30), xlim=c(1, 4), clip="off")

# p_a_raw

p_a <- p_a_raw +
  xlab('Samples per gene') +
  ylab('Enrichment for hematologic\ntumor supressor genes (Odds ratio)') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey')) +
  guides(fill=guide_legend('Called splicing outliers\n(FRASER)', 
                           title.position = "left", reverse = TRUE)) + 
  theme_vale + 
  theme(
    legend.position = c(0.5, 1.25),
    legend.direction = "vertical",
    plot.margin = margin(67, 14, 14, 20, "points"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

# p_a

figure_3a_dt <- data.table(exp_gene_fr = exp_gene_fr)
figure_3a_dt[, fr_gene_0 := exp_gene_fr %in% fr_gene_0]
figure_3a_dt[, fr_gene_1 := exp_gene_fr %in% fr_gene_1]
figure_3a_dt[, fr_gene_2 := exp_gene_fr %in% fr_gene_2]
figure_3a_dt[, fr_gene_5 := exp_gene_fr %in% fr_gene_5]
figure_3a_dt[, fr_top3_gene_0 := exp_gene_fr %in% fr_top3_gene_0]
figure_3a_dt[, fr_top3_gene_1 := exp_gene_fr %in% fr_top3_gene_1]
figure_3a_dt[, fr_top3_gene_2 := exp_gene_fr %in% fr_top3_gene_2]
figure_3a_dt[, fr_top3_gene_5 := exp_gene_fr %in% fr_top3_gene_5]

fwrite(figure_3a_dt, snakemake@output$figure_3a)




#### b. abSplice, leu13, leukemia tsg ####
#' abSplice number of splice-affecting variant predicted
absplice_res_var[AbSplice_DNA>=0.2, .N]
#' abSplice number of samples predicted with splice-affecting variant
absplice_res_var[AbSplice_DNA>=0.2, unique(sampleID)] %>% length()
#' abSplice number of splice-affecting variant predicted
absplice_res[AbSplice_DNA>=0.2, .N]
#' abSplice number of samples predicted with splice-affecting variant
absplice_res[AbSplice_DNA>=0.2, unique(sampleID)] %>% length()
#' abSplice number of samples predicted with splice-affecting gene - median
absplice_samp_freq <- data.table(sampleID = single_group_id)
absplice_samp_freq <- merge(absplice_samp_freq, absplice_res[AbSplice_DNA>=0.2, .N, by='sampleID'], all=TRUE)
absplice_samp_freq[is.na(N), N := 0]
absplice_samp_freq[, sample_rank := rank(N, ties.method = 'random')]
median(absplice_samp_freq[, N])

exp_gene <- gencode[gene_type=='protein_coding', gene_id]
exp_gene <- sapply(exp_gene, function(x){strsplit(x, "[.]")[[1]][1]})

# all 
absplice_sig <- absplice_res[AbSplice_DNA >= 0.2 , ]
absplice_freq <- absplice_sig[, .N, by='gene_id']
absplice_gene <- absplice_freq[, gene_id]
absplice_gene_0 <- exp_gene[!exp_gene %in% absplice_freq[, gene_id]]
absplice_gene_1 <- absplice_freq[N==1, gene_id]
absplice_gene_2 <- absplice_freq[N>1 & N<5, gene_id]
absplice_gene_5 <- absplice_freq[N>=5, gene_id]

#' enrichment overall
fisher_test(absplice_gene, exp_gene, leu_tsg[, ENSGid])$estimate
fisher_test(absplice_gene, exp_gene, leu_tsg[, ENSGid])$p.value %>% format(scientific = TRUE)

absplice_all_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(absplice_gene_0, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_gene_1, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_gene_2, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_gene_5, exp_gene, leu_tsg[, ENSGid])$estimate
  ),
  odds_ratio_ci_low = c(
    fisher_test(absplice_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(absplice_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[2]
  ),
  p_val = c(
    fisher_test(absplice_gene_0, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_gene_1, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_gene_2, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_gene_5, exp_gene, leu_tsg[, ENSGid])$p.value
  ),
  total = c(
    length(absplice_gene_0),
    length(absplice_gene_1),
    length(absplice_gene_2),
    length(absplice_gene_5)
  ))

absplice_all_en <- add_significance(absplice_all_en, "p_val")
absplice_all_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
absplice_all_en[, cutoff := 'All']

# top 3
absplice_sig <- merge(absplice_sig, absplice_sig[, .(gene_id, rank(-AbSplice_DNA)), by='sampleID'], 
                      by = c('gene_id', 'sampleID'))
setnames(absplice_sig, 'V2', 'padjustRankInSample')
absplice_top3 <- absplice_sig[padjustRankInSample <= 3, ]
absplice_top3_freq <- absplice_top3[, .N, by='gene_id']
absplice_top3_gene_0 <- exp_gene[!exp_gene %in% absplice_top3_freq[, gene_id]]
absplice_top3_gene_1 <- absplice_top3_freq[N==1, gene_id]
absplice_top3_gene_2 <- absplice_top3_freq[N>1 & N<5, gene_id]
absplice_top3_gene_5 <- absplice_top3_freq[N>=5, gene_id]

absplice_top3_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(absplice_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$estimate
  ),
  odds_ratio_ci_low = c(
    fisher_test(absplice_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(absplice_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[2]
  ),
  p_val = c(
    fisher_test(absplice_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$p.value
  ),
  total = c(
    length(absplice_top3_gene_0),
    length(absplice_top3_gene_1),
    length(absplice_top3_gene_2),
    length(absplice_top3_gene_5)
  ))

absplice_top3_en <- add_significance(absplice_top3_en, "p_val")
absplice_top3_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
absplice_top3_en[, cutoff := 'At most three']

absplice_en <- rbind(absplice_all_en, absplice_top3_en)
absplice_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=list(rownames(absplice_en))]
absplice_en[total < 1000, total_label := as.character(total)]
absplice_en[category=='0', total_pos := 1]
absplice_en[category=='1', total_pos := 2]
absplice_en[category=='2-4', total_pos := 3]
absplice_en[category=="\u22655", total_pos := 4]
absplice_en[cutoff=='All', total_pos := total_pos-0.2]
absplice_en[cutoff=='At most three', total_pos := total_pos+0.2]

absplice_en[p_val.signif=='ns', odds_ratio_ci_low := NA]
absplice_en[p_val.signif=='ns', odds_ratio_ci_up := NA]
absplice_en[, max(odds_ratio_ci_up, na.rm=TRUE)]

signif_height_b <- 30
total_height_b <- 60

p_b_raw <- ggplot(absplice_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=odds_ratio_ci_low, ymax=odds_ratio_ci_up), 
                width=0.1, position=position_dodge(width=0.7)) +
  geom_text(data = absplice_en[cutoff=='All'], aes(y=signif_height_b, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = absplice_en[cutoff=='At most three'], aes(y=signif_height_b, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=absplice_en[, total_pos], y=total_height_b, label=absplice_en[, total_label], size = 2) +
  annotate("text", x=0.25, y=total_height_b, label='n =', size = 2) +
  scale_y_log10(
    breaks = c(0.1, 0.3, 1, 3, 10),
    minor_breaks = c(c(1:9)/10, c(1:10))
  ) +
  coord_cartesian(ylim=c(0.1, 35), xlim=c(1, 4), clip="off")

# p_b_raw

p_b <- p_b_raw +
  xlab('Samples per gene') +
  ylab('Enrichment for hematologic\ntumor supressor genes (Odds ratio)') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey')) +
  guides(fill=guide_legend('Predicted splice-affecting\nvariants (AbSplice)', 
                           title.position = "left", reverse = TRUE)) + 
  theme_vale + 
  theme(
    legend.position = c(0.5, 1.25),
    legend.direction = "vertical",
    plot.margin = margin(67, 14, 14, 20, "points"), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

# p_b

figure_3b_dt <- data.table(exp_gene = exp_gene)
figure_3b_dt[, absplice_gene_0 := exp_gene %in% absplice_gene_0]
figure_3b_dt[, absplice_gene_1 := exp_gene %in% absplice_gene_1]
figure_3b_dt[, absplice_gene_2 := exp_gene %in% absplice_gene_2]
figure_3b_dt[, absplice_gene_5 := exp_gene %in% absplice_gene_5]
figure_3b_dt[, absplice_top3_gene_0 := exp_gene %in% absplice_top3_gene_0]
figure_3b_dt[, absplice_top3_gene_1 := exp_gene %in% absplice_top3_gene_1]
figure_3b_dt[, absplice_top3_gene_2 := exp_gene %in% absplice_top3_gene_2]
figure_3b_dt[, absplice_top3_gene_5 := exp_gene %in% absplice_top3_gene_5]

fwrite(figure_3b_dt, snakemake@output$figure_3b)




#### s8 s9. validated by mutation #####
prop_mutation_fr <- fread(snakemake@input$prop_mutation_fr)

prop_mutation_fr[, Category := factor(Category, levels = rev(c('Splicing outlier', 'None outlier')))]
prop_mutation_fr[, mutation_type := factor(mutation_type, 
                                           levels = rev(c(
                                             "copy_number_loss",   
                                             "copy_number_gain",
                                             "structural_variant",
                                             
                                             "promoter_variant",
                                             "frameshift_variant",
                                             "stop_gained",
                                             
                                             "absplice_vep",
                                             "vep_only",
                                             "absplice_only",
                                             
                                             "multiple_type", 
                                             "none"
                                           )),
                                           labels = rev(c(
                                             "Copy number loss",   
                                             "Copy number gain",
                                             "Structural",
                                             
                                             "Promoter (TSS±2Kbp)",
                                             "VEP frameshift",
                                             "VEP stop-gained",
                                             
                                             "AbSplice and VEP splice-related",
                                             "VEP splice-related exclusive",
                                             "AbSplice exclusive",
                                             
                                             "Multiple",
                                             "None"
                                           )))]

p_s8_raw <- ggplot(prop_mutation_fr, aes(x=Category, y=Percentage, fill=mutation_type)) + 
  geom_bar(stat="identity",position = "stack")
# p_s8_raw

p_s8 <- p_s8_raw +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits = c(0,1), breaks = c(0:10)/10) +
  scale_fill_manual(values=rev(c(
    RColorBrewer::brewer.pal(8, "YlOrRd")[c(8,6,4)],
    RColorBrewer::brewer.pal(8, "YlGnBu")[c(2,6,8,3,4,5)],
    RColorBrewer::brewer.pal(8, "BuPu")[c(7)],
    'lightgrey'
  )), name='Mutation type',  guide = guide_legend(reverse = TRUE)) +
  xlab('') +
  theme_vale +
  coord_flip()

# p_s8

p_s9_raw <- ggplot(prop_mutation_fr[mutation_type != 'none'], aes(x=Category, y=Percentage, fill=mutation_type)) + 
  geom_bar(stat="identity",position = "stack")
# p_s9_raw

p_s9 <- p_s9_raw +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits = c(0,0.2), breaks = c(0:10)/10) +
  scale_fill_manual(values=rev(c(
    RColorBrewer::brewer.pal(8, "YlOrRd")[c(8,6,4)],
    RColorBrewer::brewer.pal(8, "YlGnBu")[c(2,6,8,3,4,5)],
    RColorBrewer::brewer.pal(8, "BuPu")[c(7)]
  )), name='Mutation type',  guide = guide_legend(reverse = TRUE)) +
  xlab('') +
  theme_vale +
  coord_flip()

# p_s9




#### s10. validated by mutaiton ####
fr_mut_en <- fread(snakemake@input$fr_mut_en)
fr_mut_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=list(rownames(fr_mut_en))]
fr_mut_en[total < 1000, total_label := as.character(total)]
fr_mut_en <- add_significance(fr_mut_en, "p_val")

fr_mut_en[, total_pos := 0]
fr_mut_en[mutation_type=='promoter_variant', total_pos := 1]
fr_mut_en[mutation_type=='structural_variant', total_pos := 2]
fr_mut_en[mutation_type=='splice_related_variant', total_pos := 3]
fr_mut_en[mutation_type=="absplice_variant", total_pos := 4]
fr_mut_en[mutation_type=="copy_number_gain", total_pos := 5]
fr_mut_en[mutation_type=="copy_number_loss", total_pos := 6]
fr_mut_en[mutation_type=="frameshift_variant", total_pos := 7]
fr_mut_en[mutation_type=="stop_gained", total_pos := 8]
fr_mut_en[cutoff=='All', total_pos := total_pos-0.2]
fr_mut_en[cutoff=='At most three', total_pos := total_pos+0.2]

fr_mut_en[, mutation_type := factor(mutation_type, 
                                    levels = c(
                                      "promoter_variant",
                                      "structural_variant",
                                      "splice_related_variant",
                                      "absplice_variant",
                                      "copy_number_gain",
                                      "copy_number_loss",   
                                      "frameshift_variant",
                                      "stop_gained"
                                    ),
                                    labels = c(
                                      "Promoter (TSS±2Kbp)",
                                      "Structural",
                                      "VEP splice-related",
                                      "AbSplice",
                                      "Copy number gain",
                                      "Copy number loss",   
                                      "VEP frameshift",
                                      "VEP stop-gained"
                                    ))]

signif_height_b <- 1700
total_height_b <- 4000

p_s10_raw <- ggplot(fr_mut_en, aes(x=mutation_type, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=odds_ratio_ci_low, ymax=odds_ratio_ci_up), 
                width=0.1, position=position_dodge(width=0.7)) +
  geom_text(data = fr_mut_en[cutoff=='All'], aes(y=signif_height_b, label = p_val.signif),
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = fr_mut_en[cutoff=='At most three'], aes(y=signif_height_b, label = p_val.signif),
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=fr_mut_en[, total_pos], y=total_height_b, label=fr_mut_en[, total_label], size = 3) +
  annotate("text", x=0.25, y=total_height_b, label='n =', size = 3) +
  scale_y_log10(
    breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000)
  ) +
  coord_cartesian(ylim=c(1, 2000), xlim=c(1, 8), clip="off")

# p_s10_raw

p_s10 <- p_s10_raw +
  xlab('Mutation type') +
  ylab('Enrichment\nfor specific mutations (Odds ratio)') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey')) +
  guides(fill=guide_legend('Splicing outliers (FRASER)', 
                           title.position = "left", reverse = TRUE)) + 
  theme_vale + 
  theme(
    legend.position = c(0.5, 1.25),
    legend.direction = "vertical",
    plot.margin = margin(67, 14, 14, 20, "points"), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

# p_s10




#### i. FRASER, abSplice, leu13, leukemia tsg ####
exp_gene <- gencode[gene_type=='protein_coding', gene_id]
exp_gene <- sapply(exp_gene, function(x){strsplit(x, "[.]")[[1]][1]})

# all 
absplice_sig <- absplice_res[AbSplice_DNA >= 0.2 , ]
absplice_fr <- absplice_sig[samp_symbol %in% fr_res_junc[, samp_symbol], ]

absplice_fr_freq <- absplice_fr[, .N, by='gene_id']
absplice_fr_gene <- absplice_fr_freq[, gene_id]
absplice_fr_gene_0 <- exp_gene[!exp_gene %in% absplice_fr_freq[, gene_id]]
absplice_fr_gene_1 <- absplice_fr_freq[N==1, gene_id]
absplice_fr_gene_2 <- absplice_fr_freq[N>1 & N<5, gene_id]
absplice_fr_gene_5 <- absplice_fr_freq[N>=5, gene_id]

#' enrichment overall
fisher_test(absplice_fr_gene, exp_gene, leu_tsg[, ENSGid])$estimate
fisher_test(absplice_fr_gene, exp_gene, leu_tsg[, ENSGid])$p.value %>% format(scientific = TRUE)

absplice_fr_all_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(absplice_fr_gene_0, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_fr_gene_1, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_fr_gene_2, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_fr_gene_5, exp_gene, leu_tsg[, ENSGid])$estimate
  ),
  odds_ratio_ci_low = c(
    fisher_test(absplice_fr_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_fr_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_fr_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_fr_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(absplice_fr_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_fr_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_fr_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_fr_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[2]
  ),
  p_val = c(
    fisher_test(absplice_fr_gene_0, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_fr_gene_1, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_fr_gene_2, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_fr_gene_5, exp_gene, leu_tsg[, ENSGid])$p.value
  ),
  total = c(
    length(absplice_fr_gene_0),
    length(absplice_fr_gene_1),
    length(absplice_fr_gene_2),
    length(absplice_fr_gene_5)
  ))

absplice_fr_all_en <- add_significance(absplice_fr_all_en, "p_val")
absplice_fr_all_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
absplice_fr_all_en[, cutoff := 'All']

# top 3
absplice_fr <- merge(absplice_fr, absplice_fr[, .(gene_id, rank(-AbSplice_DNA)), by='sampleID'], 
                     by = c('gene_id', 'sampleID'))
setnames(absplice_fr, 'V2', 'padjustRankInSample')
absplice_fr_top3 <- absplice_fr[padjustRankInSample <= 3, ]
absplice_fr_top3_freq <- absplice_fr_top3[, .N, by='gene_id']
absplice_fr_top3_gene_0 <- exp_gene[!exp_gene %in% absplice_fr_top3_freq[, gene_id]]
absplice_fr_top3_gene_1 <- absplice_fr_top3_freq[N==1, gene_id]
absplice_fr_top3_gene_2 <- absplice_fr_top3_freq[N>1 & N<5, gene_id]
absplice_fr_top3_gene_5 <- absplice_fr_top3_freq[N>=5, gene_id]

absplice_fr_top3_en <- data.table(
  category = c(
    "0", "1", "2-4", "\u22655"
  ),
  odds_ratio = c(
    fisher_test(absplice_fr_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_fr_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_fr_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$estimate,
    fisher_test(absplice_fr_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$estimate
  ),
  odds_ratio_ci_low = c(
    fisher_test(absplice_fr_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_fr_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_fr_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(absplice_fr_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(absplice_fr_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_fr_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_fr_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(absplice_fr_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[2]
  ),
  p_val = c(
    fisher_test(absplice_fr_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_fr_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_fr_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$p.value,
    fisher_test(absplice_fr_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$p.value
  ),
  total = c(
    length(absplice_fr_top3_gene_0),
    length(absplice_fr_top3_gene_1),
    length(absplice_fr_top3_gene_2),
    length(absplice_fr_top3_gene_5)
  ))

absplice_fr_top3_en <- add_significance(absplice_fr_top3_en, "p_val")
absplice_fr_top3_en[, category := factor(category, levels = c("0", "1", "2-4", "\u22655"))]
absplice_fr_top3_en[, cutoff := 'At most three']

absplice_fr_en <- rbind(absplice_fr_all_en, absplice_fr_top3_en)
absplice_fr_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=list(rownames(absplice_fr_en))]
absplice_fr_en[total < 1000, total_label := as.character(total)]
absplice_fr_en[category=='0', total_pos := 1]
absplice_fr_en[category=='1', total_pos := 2]
absplice_fr_en[category=='2-4', total_pos := 3]
absplice_fr_en[category=="\u22655", total_pos := 4]
absplice_fr_en[cutoff=='All', total_pos := total_pos-0.2]
absplice_fr_en[cutoff=='At most three', total_pos := total_pos+0.2]

absplice_fr_en[p_val.signif=='ns', odds_ratio_ci_low := NA]
absplice_fr_en[p_val.signif=='ns', odds_ratio_ci_up := NA]
absplice_fr_en[, max(odds_ratio_ci_up, na.rm=TRUE)]

signif_height_i <- 400
total_height_i <- 1000

p_i_raw <- ggplot(absplice_fr_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=odds_ratio_ci_low, ymax=odds_ratio_ci_up), 
                width=0.1, position=position_dodge(width=0.7)) +
  geom_text(data = absplice_fr_en[cutoff=='All'], aes(y=signif_height_i, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = absplice_fr_en[cutoff=='At most three'], aes(y=signif_height_i, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=absplice_fr_en[, total_pos], y=total_height_i, label=absplice_fr_en[, total_label], size = 2) +
  annotate("text", x=0.25, y=total_height_i, label='n =', size = 2) +
  # scale_y_log10(breaks=c(0.1, 1, 10, 100)) +
  scale_y_log10(
    breaks = c(0.1, 1, 10, 100),
    minor_breaks = c(rep(1:9, 3)*(10^rep(0:2, each = 9))/10)
  ) +
  coord_cartesian(ylim=c(0.1, 470), xlim=c(1, 4), clip="off")

# p_i_raw

p_i <- p_i_raw +
  xlab('Samples per gene') +
  ylab('Enrichment for hematologic\ntumor supressor genes (Odds ratio)') +
  scale_fill_manual(values=c('lightgrey', 'darkgrey')) +
  guides(fill=guide_legend('Overlap of\nsplicing outliers (FRASER) and\nsplice-affecting variants (AbSplice)', 
                           title.position = "left", reverse = TRUE)) + 
  theme_vale + 
  theme(
    plot.title = element_text(hjust = -0.45, vjust=2.12),
    legend.position = c(0.5, 1.25),
    legend.direction = "vertical",
    plot.margin = margin(67, 14, 14, 20, "points"), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

# p_i

figure_3d_dt <- data.table(exp_gene = exp_gene)
figure_3d_dt[, absplice_fr_gene_0 := exp_gene %in% absplice_fr_gene_0]
figure_3d_dt[, absplice_fr_gene_1 := exp_gene %in% absplice_fr_gene_1]
figure_3d_dt[, absplice_fr_gene_2 := exp_gene %in% absplice_fr_gene_2]
figure_3d_dt[, absplice_fr_gene_5 := exp_gene %in% absplice_fr_gene_5]
figure_3d_dt[, absplice_fr_top3_gene_0 := exp_gene %in% absplice_fr_top3_gene_0]
figure_3d_dt[, absplice_fr_top3_gene_1 := exp_gene %in% absplice_fr_top3_gene_1]
figure_3d_dt[, absplice_fr_top3_gene_2 := exp_gene %in% absplice_fr_top3_gene_2]
figure_3d_dt[, absplice_fr_top3_gene_5 := exp_gene %in% absplice_fr_top3_gene_5]

fwrite(figure_3d_dt, snakemake@output$figure_3d)




### link fr to curated/vep/absplice in AML and panel genes ####
fr_aml_panel <- fr_res[sampleID %in% samp_anno_aml[, ArrayID] &
                         hgncSymbol %in% mll_panel_genes, ]
curated_splice_aml_panel <- AML_variants[(consequence=='splice_acceptor_variant' |
                                            consequence=='splice_donor_variant' ), ]
vep_splice_aml_panel <- vep_splice[array_id %in% samp_anno_aml[, ArrayID] &
                                     SYMBOL %in% mll_panel_genes, ]
absplice_aml_panel <- absplice_res[tissue=='AML' & gene_name %in% mll_panel_genes,]
absplice_high_aml_panel <- absplice_res[tissue=='AML' & gene_name %in% mll_panel_genes & AbSplice_DNA >= 0.2,]


fr_curated_aml_panel <- merge(fr_aml_panel, curated_splice_aml_panel, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_curated_aml_panel[, table(!is.na(array_id))]

fr_vep_aml_panel <- merge(fr_aml_panel, vep_splice_aml_panel, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_vep_aml_panel[, table(!is.na(array_id))]

fr_absplice_aml_panel <- merge(fr_aml_panel, absplice_aml_panel, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_absplice_aml_panel[, table(!is.na(sampleID.y))]

fr_absplice_high_aml_panel <- merge(fr_aml_panel, absplice_high_aml_panel, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_absplice_high_aml_panel[, table(!is.na(sampleID.y))]


curated_splice_fr_aml_panel <- merge(curated_splice_aml_panel, fr_aml_panel, by='samp_symbol', all.x=TRUE, all.y=FALSE)
curated_splice_fr_aml_panel[, table(!is.na(sampleID))]

vep_splice_fr_aml_panel <- merge(vep_splice_aml_panel, fr_aml_panel, by='samp_symbol', all.x=TRUE, all.y=FALSE)
vep_splice_fr_aml_panel[, table(!is.na(sampleID))]

absplice_fr_aml_panel <- merge(absplice_aml_panel, fr_aml_panel, by='samp_symbol', all.x=TRUE, all.y=FALSE)
absplice_fr_aml_panel[, table(!is.na(sampleID.y))]

absplice_high_fr_aml_panel <- merge(absplice_high_aml_panel, fr_aml_panel, by='samp_symbol', all.x=TRUE, all.y=FALSE)
absplice_high_fr_aml_panel[, table(!is.na(sampleID.y))]




### link fr to vep/absplice in AML and all genes ####
fr_aml <- fr_res[sampleID %in% samp_anno_aml[, ArrayID], ]
vep_splice_aml <- vep_splice[, ]
absplice_aml <- absplice_res[tissue=='AML',]
absplice_high_aml <- absplice_res[tissue=='AML' & AbSplice_DNA >= 0.2,]


fr_vep_aml <- merge(fr_aml, vep_splice_aml, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_vep_aml[, table(!is.na(array_id))]

fr_absplice_aml <- merge(fr_aml, absplice_aml, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_absplice_aml[, table(!is.na(sampleID.y))]

fr_absplice_high_aml <- merge(fr_aml, absplice_high_aml, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_absplice_high_aml[, table(!is.na(sampleID.y))]


vep_splice_fr_aml <- merge(vep_splice_aml, fr_aml, by='samp_symbol', all.x=TRUE, all.y=FALSE)
vep_splice_fr_aml[, table(!is.na(sampleID))]

absplice_fr_aml <- merge(absplice_aml, fr_aml, by='samp_symbol', all.x=TRUE, all.y=FALSE)
absplice_fr_aml[, table(!is.na(sampleID.y))]

absplice_high_fr_aml <- merge(absplice_high_aml, fr_aml, by='samp_symbol', all.x=TRUE, all.y=FALSE)
absplice_high_fr_aml[, table(!is.na(sampleID.y))]




### link fr to absplice/vep in all samples and all genes ####
fr <- copy(fr_res)
absplice <- copy(absplice_res)
absplice_high <- absplice[AbSplice_DNA >= 0.2,]


fr_vep <- merge(fr, vep_splice, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_vep[, table(!is.na(array_id))]

fr_absplice <- merge(fr, absplice, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_absplice[, table(!is.na(sampleID.y))]

fr_absplice_high <- merge(fr, absplice_high, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_absplice_high[, table(!is.na(sampleID.y))]


vep_splice_fr <- merge(vep_splice, fr, by='samp_symbol', all.x=TRUE, all.y=FALSE)
vep_splice_fr[, table(!is.na(sampleID))]

absplice_fr <- merge(absplice, fr, by='samp_symbol', all.x=TRUE, all.y=FALSE)
absplice_fr[, table(!is.na(sampleID.y))]

absplice_high_fr <- merge(absplice_high, fr, by='samp_symbol', all.x=TRUE, all.y=FALSE)
absplice_high_fr[, table(!is.na(sampleID.y))]




#### c. overlap venn ####
vep_splice[, table(Consequence)]

x <- readRDS(snakemake@input$venn_list)

p_c_raw <- ggvenn(
  x, 
  fill_color = c(RColorBrewer::brewer.pal(12, "Paired")[2],
                 RColorBrewer::brewer.pal(12, "Paired")[6],
                 RColorBrewer::brewer.pal(12, "Paired")[11]
  ),
  stroke_size = 0.5, set_name_size = 0, show_percentage = FALSE, text_size = 3.5
)

p_c <- p_c_raw +
  annotate("text", x=-1, y=2, label='Splice-affecting variants\n(AbSplice)', size = 3.5, fontface = "bold") +
  annotate("text", x=1.2, y=2, label='Splicing outliers\n(FRASER)', size = 3.5, fontface = "bold") +
  annotate("text", x=0, y=-2, label='Splice-related variants\n(VEP)', size = 3.5, fontface = "bold") +
  theme(plot.margin = margin(14, 20, 20, 20, "points"))
# p_c




### case study searching ####
fr_junc_absplice_var <- merge(fr_res_junc, absplice_res_var, by='samp_symbol', all.x=TRUE, all.y=FALSE)
fr_junc_absplice_var[, region_igv := paste0(seqnames, ":", start, "-", end)]
fr_junc_absplice_var[, region_igv100 := paste0(seqnames, ":", start-100, "-", end+100)]
fr_junc_absplice_var <- fr_junc_absplice_var[hgncSymbol %in% mll_panel_genes, ]
fr_junc_absplice_var <- merge(fr_junc_absplice_var, samp_anno_14groups[, .(ArrayID, DROP_GROUP)],
                              by.x='sampleID.x', by.y='ArrayID', all.x=TRUE, all.y=FALSE)
curate_sub <- fr_junc_absplice_var[order(-abs(deltaPsi)), .(samp_symbol, DROP_GROUP.y, seqnames, start, end, deltaPsi, counts, region_igv100,
                                                            tissue, variant, chr, pos, ref, alt, AbSplice_DNA)]
curate_sub[1:10, ]

# samp_symbol_p <- 'MLL_19546-KMT2D' # outlier messy, coverage low, intron retention
# samp_symbol_p <- 'MLL_10583-GATA2' # coverage low, whole transcript missing, expression outlier more likely
# samp_symbol_p <- 'MLL_11493-GATA1' # outlier messy, coverage low, intron retention
samp_symbol_p <- 'MLL_19282-RB1' # exon skipping
res_p <- fr_res_junc[samp_symbol==samp_symbol_p,][1]
# fds_p <- loadFraserDataSet(file = snakemake@input$fds[
#   grep(paste0('/', fr_junc_absplice_var[samp_symbol==samp_symbol_p, DROP_GROUP.y][1], '/'), snakemake@input$fraserRes)])
# plotExpression(fds_p, result = res_p)
# plotBamCoverageFromResultTable(fds_p, result = res_p)




### case study searching not in VEP ####
absplice_var_sig <- absplice_res_var[AbSplice_DNA >= 0.2, ] 
absplice_var_sig <- absplice_var_sig[order(-AbSplice_DNA), ] 
absplice_var_sig <- absplice_var_sig[isCGC == TRUE, ]
absplice_var_sig_fraser_junction <- absplice_var_sig[samp_symbol %in% fr_res_junc[,samp_symbol],]

chain <- import.chain("/s/project/vale/Resource/hg19ToHg38.over.chain")

absplice_hg38 <- data.table()
for (i in  1:nrow(absplice_var_sig_fraser_junction)){
  
  variant <- absplice_var_sig_fraser_junction[i, ]
  chromosome <- paste0("chr", variant[,chr])
  position <- as.numeric(variant[,pos])
  ref <- variant[, ref]
  alt <- variant[, alt]
  
  # Create a GRanges object with the input variant
  input_gr <- GRanges(
    seqnames = Rle(chromosome),
    ranges = IRanges(start = position, end = position),
    strand = "*",
    ref = ref,
    alt = alt
  )
  
  output_gr <- as.data.table(liftOver(input_gr, chain))
  vep_location <- paste0(sub("^chr", "", output_gr$seqnames), ":",  output_gr$start)
  lifted_over_variant <- data.table("gene_name" = variant$gene_name, 
                                    "samp_symbol" = variant$samp_symbol, 
                                    "Location_hg38" = vep_location, 
                                    "variant" = variant$variant)
  absplice_hg38 <- rbindlist(list(absplice_hg38, lifted_over_variant))  
}

absplice_in_vep <- merge(absplice_hg38, vep_splice, 
                         by.x = "Location_hg38", by.y='Location',
                         all.x = TRUE, all.y = FALSE)

absplice_not_in_vep <- absplice_in_vep[is.na(SYMBOL),  ]
absplice_not_in_vep <- absplice_not_in_vep[, .(gene_name, samp_symbol.x, variant, Location_hg38)]
setnames(absplice_not_in_vep, 'samp_symbol.x', 'samp_symbol')

absplice_not_in_vep <- merge(absplice_not_in_vep, 
                             absplice_res_var[, .(samp_symbol, variant, AbSplice_DNA)], 
                             by=c('samp_symbol', 'variant'))
absplice_not_in_vep <- separate(absplice_not_in_vep, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_not_in_vep <- separate(absplice_not_in_vep, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()
absplice_not_in_vep[, pos := as.numeric(pos)]

absplice_not_in_vep <- merge(absplice_not_in_vep, fr_res_junc, by='samp_symbol')         

absplice_not_in_vep[, min_abs_dist := min(abs(pos-start), abs(pos-end)), by=list(rownames(absplice_not_in_vep))]
absplice_not_in_vep[, dist_pos_start := pos-start, by=list(rownames(absplice_not_in_vep))]
absplice_not_in_vep[, dist_pos_end := pos-end, by=list(rownames(absplice_not_in_vep))]
absplice_not_in_vep[, .(variant, junction_id, min_abs_dist)]      
hist(absplice_not_in_vep[, min_abs_dist], 100)

absplice_not_in_vep[, region_igv := paste0(seqnames, ":", start, "-", end)]
absplice_not_in_vep[, region_igv100 := paste0(seqnames, ":", start-100, "-", end+100)]
absplice_not_in_vep[, junction_length := end - start]
absplice_not_in_vep <- absplice_not_in_vep[order(-AbSplice_DNA), ]

# no good recurrent variant example found
absplice_not_in_vep_dup <- absplice_not_in_vep[duplicated(variant), 
                                               .(samp_symbol, 
                                                 variant, Location_hg38, AbSplice_DNA, 
                                                 deltaPsi, min_abs_dist, dist_pos_start, dist_pos_end, region_igv, region_igv100)]
absplice_not_in_vep_dup <- absplice_not_in_vep_dup[!duplicated(samp_symbol), ]

# less than 1000 bp junction
absplice_not_in_vep_short <- absplice_not_in_vep[junction_length < 1000, ]
# absplice_not_in_vep_short <- absplice_not_in_vep_short[, c(1: 22, 73:79) ]
absplice_not_in_vep[variant=='19:42384927:ATCCCAGGGCCTGAACCTGGACGACTGCTCCATGTATGAGGACATC>A', 
                    .(samp_symbol, variant, AbSplice_DNA, 
                      padjust, deltaPsi, junction_id, dist_pos_start, dist_pos_end, region_igv100)]
# samp_symbol                                                      variant AbSplice_DNA    padjust
# 1: MLL_22178-CD79A 19:42384927:ATCCCAGGGCCTGAACCTGGACGACTGCTCCATGTATGAGGACATC>A     0.355067 2.2529e-03
# 2: MLL_60467-CD79A 19:42384927:ATCCCAGGGCCTGAACCTGGACGACTGCTCCATGTATGAGGACATC>A     0.355067 1.9822e-57
# deltaPsi                         junction_id dist_pos_start dist_pos_end           region_igv100
# 1:    -0.35 Intron chr19: 42,384,806-42,384,933            121           -6 chr19:42384706-42385033
# 2:    -0.96 Intron chr19: 42,384,806-42,384,933            121           -6 chr19:42384706-42385033




### case study searching not in VEP first run (before intogen) ####

# checking to see if CD79A variant was in vep before intogen
sample_vcf <- fread("/s/project/mll/preprocess_202302/Somatic_Analysis_abSplice/NovaSeq_Male_MLL_60467-M060_G1_P1.somatic.abSplice.txt")
variant <- sample_vcf[pos == 42384927]
variant[, INFO]
# MLL_60467-CD79A has splicing annotation in vep

mll_22179_vcf <- fread("/s/project/mll/preprocess_202302/Somatic_Analysis_abSplice/-M069_comb-MLL_22178_G1_P1.somatic.abSplice.txt")
variant <- mll_22179_vcf[pos == 42384927]
variant[, INFO]
# we look for other candidates
absplice_var_sig <- absplice_res_var[AbSplice_DNA >= 0.2, ] 
absplice_var_sig <- absplice_var_sig[order(-AbSplice_DNA), ] 
absplice_var_sig <- absplice_var_sig[isCGC == TRUE, ]
absplice_var_sig_fraser_junction <- absplice_var_sig[samp_symbol %in% fr_res_junc[,samp_symbol],]
absplice_var_sig_fraser_junction$pos = as.integer(absplice_var_sig_fraser_junction$pos)


absplice_in_vep <- merge(absplice_var_sig_fraser_junction, vep_res_abSplice, 
                         by='pos',
                         all.x = TRUE, all.y = FALSE)

absplice_not_in_vep <- absplice_in_vep[is.na(sample_id),  ]
absplice_not_in_vep <- absplice_not_in_vep[, .(gene_name, samp_symbol, variant, pos)]

absplice_not_in_vep <- merge(absplice_not_in_vep, 
                             absplice_res_var[, .(samp_symbol, variant, AbSplice_DNA)], 
                             by=c('samp_symbol', 'variant'))
absplice_not_in_vep <- separate(absplice_not_in_vep, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_not_in_vep <- separate(absplice_not_in_vep, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()
absplice_not_in_vep[, pos := as.numeric(pos)]

absplice_not_in_vep <- merge(absplice_not_in_vep, fr_res_junc, by='samp_symbol')         

absplice_not_in_vep[, min_abs_dist := min(abs(pos-start), abs(pos-end)), by=list(rownames(absplice_not_in_vep))]
absplice_not_in_vep[, dist_pos_start := pos-start, by=list(rownames(absplice_not_in_vep))]
absplice_not_in_vep[, dist_pos_end := pos-end, by=list(rownames(absplice_not_in_vep))]
absplice_not_in_vep[, .(variant, junction_id, min_abs_dist)]      
hist(absplice_not_in_vep[, min_abs_dist], 100)

absplice_not_in_vep[, region_igv := paste0(seqnames, ":", start, "-", end)]
absplice_not_in_vep[, region_igv100 := paste0(seqnames, ":", start-100, "-", end+100)]
absplice_not_in_vep[, junction_length := end - start]
absplice_not_in_vep <- absplice_not_in_vep[order(-AbSplice_DNA), ]

absplice_subset <- absplice_not_in_vep[, .(samp_symbol, variant,region_igv, min_abs_dist, pos, start, end, dist_pos_end, dist_pos_start, deltaPsi) ]

absplice_subset_dist <- absplice_subset[dist_pos_end < 0 & dist_pos_start > 0 &
                                          min_abs_dist !=2 &
                                          abs(dist_pos_start) < 1000 &
                                          abs(dist_pos_end) < 1000, ]

absplice_subset_dist <- absplice_subset[dist_pos_end < 0 & dist_pos_start > 0 &
                                          min_abs_dist !=2, ]
## found 3 candidates, no duplicate variant and no duplicate junction
# MLL_24594-FUBP1  1:78430665:A>C   chr1:78430654-78430752
# MLL_22145-UBR5 8:103270939:G>A chr8:103269946-103270944 --> find a good control 
# MLL_22145-UBR5 8:103270939:G>A chr8:103269946-103271212

# checking vep effect of found variants
# TODO: fix path
MLL_24594_FUBP1_vcf = fread("/s/project/mll/preprocess_202302/Somatic_Analysis_abSplice/-M058_comb-MLL_24594_G1_P1.somatic.abSplice.txt")
variant <- MLL_24594_FUBP1_vcf[pos == 78430665]
variant[,CSQT ]

MLL_22145_UBR5_vcf = fread("/s/project/mll/preprocess_202302/Somatic_Analysis_abSplice/comb-MaleCon_comb-MLL_22145_G1_P1.somatic.abSplice.txt")
variant <- MLL_22145_UBR5_vcf[pos == 103270939]
variant[,CSQT ]




#### MLL_60467-CD79A vep/CNV/fusion ####
samp_symbol_p_d <- 'MLL_60467-CD79A'

vep_lym[samp_symbol == samp_symbol_p_d, ]
vep_lym[Location == "19:41880860", ]

absplice_res_var[variant=='19:42384927:ATCCCAGGGCCTGAACCTGGACGACTGCTCCATGTATGAGGACATC>A', ]
fr_res_junc[junction_id == 'Intron chr19: 42,384,806-42,384,933', ]
#' this splicing outlier junction is still unique to this variant

# absplice_res_var[gene_name=='CD79A', ]
# absplice_res_var[samp_symbol=='MLL_13378-CD79A', ]
# vep_lym[SYMBOL=='CD79A', ]
# vep_lym[samp_symbol=='MLL_13378-CD79A', ]
#' no variant identified that would cause this change in MLL_13378

# cnv
cnv_lym <- fread(grep('Lym_group', snakemake@input$mll_cnv, value = TRUE))
cnv_lym <- cnv_lym[abs(MEAN_LOG2_COPY_RATIO)>0.3, ]
cnv_lym[, LENGTH := END - START]
cnv_lym <- cnv_lym[LENGTH >= 1000000, ]
gencode[gene_name=='CD79A', .(start, end, seqnames)]
cnv_lym[ARRAY_ID=='MLL_60467' & CONTIG=='19', .(START, END, CALL, MEAN_COPY_RATIO, SYMBOL)] %>% unique()
#' no cnv found explainable

# sv
sv_manta_lym <- fread(grep('Lym_group', snakemake@input$mll_manta_sv, value = TRUE))
sv_manta_lym[( gene1=='CD79A' | gene2=='CD79A' ), ]
sv_manta_samp <- sv_manta_lym[array_id %in% 'MLL_60467', ]
sv_manta_samp[( gene1=='CD79A' | gene2=='CD79A' ) & array_id %in% 'MLL_60467', ]
#' no sv found explainable

# fusion
fusion_manta_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_manta, value = TRUE))
fusion_manta_Lym_group[( gene1=='CD79A' | gene2=='CD79A' ), ]
fusion_manta_samp <- fusion_manta_Lym_group[array_id %in% 'MLL_60467', ]
fusion_manta_samp[( gene1=='CD79A' | gene2=='CD79A' ), ]

fusion_arriba_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_arriba, value = TRUE))
fusion_arriba_Lym_group[( gene1=='CD79A' | gene2=='CD79A' ), ]
fusion_arriba_samp <- fusion_arriba_Lym_group[array_id %in% 'MLL_60467', ]
fusion_arriba_samp[( gene1=='CD79A' | gene2=='CD79A' ), ]

fusion_star_fusion_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_star_fusion, value = TRUE))
fusion_star_fusion_Lym_group[( gene1=='CD79A' | gene2=='CD79A' ), ]
fusion_star_fusion_samp <- fusion_star_fusion_Lym_group[array_id %in% 'MLL_60467', ]
fusion_star_fusion_samp[( gene1=='CD79A' | gene2=='CD79A' ), ]
#' nothing found 




#### MLL_22178-CD79A vep/CNV/fusion ####
samp_symbol_p_d <- 'MLL_22178-CD79A'

vep_lym[samp_symbol == samp_symbol_p_d, ]
vep_lym[Location == "19:41880860", ]

absplice_res_var[variant=='19:42384927:ATCCCAGGGCCTGAACCTGGACGACTGCTCCATGTATGAGGACATC>A', ]
fr_res_junc[junction_id == 'Intron chr19: 42,384,806-42,384,933', ]
#' this splicing outlier junction is still unique to this variant

# cnv
cnv_lym <- fread(grep('Lym_group', snakemake@input$mll_cnv, value = TRUE))
cnv_lym <- cnv_lym[abs(MEAN_LOG2_COPY_RATIO)>0.3, ]
cnv_lym[, LENGTH := END - START]
cnv_lym <- cnv_lym[LENGTH >= 1000000, ]
gencode[gene_name=='CD79A', .(start, end, seqnames)]
cnv_lym[ARRAY_ID=='MLL_22178' & CONTIG=='19', .(START, END, CALL, MEAN_COPY_RATIO, SYMBOL)] %>% unique()
#' cnv gain found, but doesn't seem explainable

# sv
sv_manta_lym <- fread(grep('Lym_group', snakemake@input$mll_manta_sv, value = TRUE))
sv_manta_lym[( gene1=='CD79A' | gene2=='CD79A' ), ]
sv_manta_samp <- sv_manta_lym[array_id %in% 'MLL_22178', ]
sv_manta_samp[( gene1=='CD79A' | gene2=='CD79A' ) & array_id %in% 'MLL_22178', ]
#' no sv found explainable

# fusion
fusion_manta_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_manta, value = TRUE))
fusion_manta_Lym_group[( gene1=='CD79A' | gene2=='CD79A' ), ]
fusion_manta_samp <- fusion_manta_Lym_group[array_id %in% 'MLL_22178', ]
fusion_manta_samp[( gene1=='CD79A' | gene2=='CD79A' ), ]

fusion_arriba_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_arriba, value = TRUE))
fusion_arriba_Lym_group[( gene1=='CD79A' | gene2=='CD79A' ), ]
fusion_arriba_samp <- fusion_arriba_Lym_group[array_id %in% 'MLL_22178', ]
fusion_arriba_samp[( gene1=='CD79A' | gene2=='CD79A' ), ]

fusion_star_fusion_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_star_fusion, value = TRUE))
fusion_star_fusion_Lym_group[( gene1=='CD79A' | gene2=='CD79A' ), ]
fusion_star_fusion_samp <- fusion_star_fusion_Lym_group[array_id %in% 'MLL_22178', ]
fusion_star_fusion_samp[( gene1=='CD79A' | gene2=='CD79A' ), ]
#' nothing found 




#### d. n vs. k, absplice vs. FRASER: CD79A ####
samp_symbol_p_d <- 'MLL_60467-CD79A'
samp_symbol_p_d_2 <- c('MLL_60467-CD79A', 'MLL_22178-CD79A')
samp_p_d_2 <- c('MLL_60467', 'MLL_22178')
fr_junc_absplice_var[samp_symbol==samp_symbol_p_d, 
                     .(samp_symbol, DROP_GROUP.y, seqnames, start, end, deltaPsi, counts, region_igv100,
                       tissue, variant, chr, pos, ref, alt, AbSplice_DNA)]
fr_junc_absplice_var[samp_symbol%in%samp_symbol_p_d_2, 
                     .(samp_symbol, DROP_GROUP.y, seqnames, start, end, deltaPsi, counts, region_igv100,
                       tissue, variant, chr, pos, ref, alt, AbSplice_DNA)]
samp_anno_cd79a <- samp_anno[ArrayID %in% samp_p_d_2, 
                             .(ArrayID, Diag)]
samp_anno_cd79a_fr <- merge(samp_anno_cd79a, 
                            fr_junc_absplice_var[samp_symbol%in%samp_symbol_p_d_2, 
                                                 .(sampleID.x, hgncSymbol, DROP_GROUP.y, 
                                                   seqnames, start, end, deltaPsi, counts, totalCounts, 
                                                   variant, AbSplice_DNA)], 
                            by.x='ArrayID', by.y='sampleID.x')
fwrite(samp_anno_cd79a_fr, snakemake@output$cd79a_tab_raw)


res_p_d <- fr_res_junc[samp_symbol==samp_symbol_p_d,]
fds_p_d_path <- snakemake@input$fds[
  grep(paste0('/', unique(fr_junc_absplice_var[samp_symbol==samp_symbol_p_d, DROP_GROUP.y]), '/'), 
       snakemake@input$fraserRes)]
fds_p_d <- loadFraserDataSet(file = fds_p_d_path)
res_p_d[, junction_id]
p_d_title <- expression(paste(
  bold("Intron chr19: 42,384,806-42,384,933"), 
  bold(", "), 
  bolditalic("CD79A")))

p_d_raw <- plotExpression(fds_p_d, result = res_p_d)
p_d_raw$layers[[3]] <- NULL


anonym_id_dict <- c(
  'MLL_22178' = anonymization_table_mll[anonymization_table_mll$ArrayID == "MLL_22178", AnonamizedID],
  'MLL_60467' = anonymization_table_mll[anonymization_table_mll$ArrayID == "MLL_60467", AnonamizedID], 
  'MLL_19282' = anonymization_table_mll[anonymization_table_mll$ArrayID == "MLL_19282", AnonamizedID]
)
options(ggrepel.max.overlaps = Inf)

p_d <- p_d_raw + 
  
  geom_text_repel(aes(label=ifelse(sampleID %in% c("MLL_22178", "MLL_60467"), anonym_id_dict[as.character(sampleID)], '')),
                  
                  color="firebrick", fontface = "bold") +
  
  scale_color_manual(values = c('firebrick', 'darkgrey'),
                     labels = c("Outlier", "Non-outlier")) +
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
  ggtitle(p_d_title) +
  xlab('Total junction coverage + 2') +
  ylab('Junction counts + 1') +
  labs(color='Called as outlier')  +
  
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = c(.05, .95),
    legend.justification = c("left", "top"),
    legend.margin = margin(0, 8, 4, 4),
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),  
    legend.key=element_rect(fill=NA)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

# p_d




#### f1. bam cov, absplice vs. FRASER: CD79A ####
#+ plot p_f_1, fig.width=8, fig.height=8
samp_symbol_p_f_1 <- 'MLL_60467-CD79A'
res_p_f_1 <- fr_res_junc[samp_symbol==samp_symbol_p_d,]

# p_f_1 <- plotBamCoverageFromResultTable(fds_p_d, result = res_p_f_1)




#### f2. bam cov, absplice vs. FRASER: CD79A ####
#+ plot p_f_2, fig.width=8, fig.height=8
samp_symbol_p_f_2 <- 'MLL_22178-CD79A'
res_p_f_2 <- fr_res_junc[samp_symbol==samp_symbol_p_d,]

# p_f_2 <- plotBamCoverageFromResultTable(fds_p_d, result = res_p_d)




#### MLL_19282-RB1 vep/CNV/fusion ####
samp_symbol_p_e <- 'MLL_19282-RB1'
vep_lym[samp_symbol == samp_symbol_p_e, ]

# cnv
cnv_lym <- fread(grep('Lym_group', snakemake@input$mll_cnv, value = TRUE))
cnv_lym <- cnv_lym[abs(MEAN_LOG2_COPY_RATIO)>0.3, ]
cnv_lym[, LENGTH := END - START]
cnv_lym <- cnv_lym[LENGTH >= 1000000, ]
cnv_lym[ARRAY_ID=='MLL_19282' & CONTIG=='13', .(START, END, CALL, MEAN_COPY_RATIO, SYMBOL)] %>% unique()
gencode[gene_name=='RB1', .(start, end)]
#' cnv loss found, but doesn't seem explainable

# sv
sv_manta_lym <- fread(grep('Lym_group', snakemake@input$mll_manta_sv, value = TRUE))
sv_manta_lym[( gene1=='RB1' | gene2=='RB1' ), ]
sv_manta_samp <- sv_manta_lym[array_id %in% 'MLL_19282', ]
sv_manta_samp[( gene1=='RB1' | gene2=='RB1' ) & array_id %in% 'MLL_19282', ]
#' breakends found, but doesn't seem explainable

# fusion
fusion_manta_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_manta, value = TRUE))
fusion_manta_Lym_group[( gene1=='RB1' | gene2=='RB1' ), ]
fusion_manta_samp <- fusion_manta_Lym_group[array_id %in% 'MLL_19282', ]
fusion_manta_samp[( gene1=='RB1' | gene2=='RB1' ), ]

fusion_arriba_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_arriba, value = TRUE))
fusion_arriba_Lym_group[( gene1=='RB1' | gene2=='RB1' ), ]
fusion_arriba_samp <- fusion_arriba_Lym_group[array_id %in% 'MLL_19282', ]
fusion_arriba_samp[( gene1=='RB1' | gene2=='RB1' ), ]

fusion_star_fusion_Lym_group <- fread(grep('Lym_group', snakemake@input$mll_star_fusion, value = TRUE))
fusion_star_fusion_Lym_group[( gene1=='RB1' | gene2=='RB1' ), ]
fusion_star_fusion_samp <- fusion_star_fusion_Lym_group[array_id %in% 'MLL_19282', ]
fusion_star_fusion_samp[( gene1=='RB1' | gene2=='RB1' ), ]
#' nothing found 




#### e. n vs. k, FRASER ####
samp_symbol_p_e <- 'MLL_19282-RB1'
fr_junc_absplice_var[samp_symbol==samp_symbol_p_e, .(samp_symbol, DROP_GROUP.y, seqnames, start, end, deltaPsi, counts, region_igv100,
                                                     tissue, variant, chr, pos, ref, alt, AbSplice_DNA)]
fr_res_junc[junction_id == 'Intron chr13: 49,030,486-49,033,823', ]

res_p_e <- fr_res_junc[samp_symbol==samp_symbol_p_e,][2]
fds_p_e_path <- snakemake@input$fds[
  grep(paste0('/', unique(fr_junc_absplice_var[samp_symbol==samp_symbol_p_e, DROP_GROUP.y]), '/'), 
       snakemake@input$fraserRes)]
fds_p_e <- loadFraserDataSet(file = fds_p_e_path)
res_p_e[, junction_id]
p_e_title <- expression(paste(
  bold("Intron chr13: 49,033,970-49,037,866"), 
  bold(", "), 
  bolditalic("RB1")))

#+ plot e raw, fig.width=6, fig.height=6
p_e_raw <- plotExpression(fds_p_e, result = res_p_e)
p_e_raw$layers[[3]] <- NULL
p_e_raw

#+ plot e, fig.width=6, fig.height=6


p_e <- p_e_raw + 
  
  geom_text_repel(aes(label=ifelse(sampleID %in% c("MLL_19282"), anonym_id_dict[as.character(sampleID)], '')),
                  
                  color="firebrick", fontface = "bold") +
  scale_color_manual(values = c('firebrick', 'darkgrey'),
                     labels = c("Outlier", "Non-outlier")) +
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
  ggtitle(p_e_title) +
  xlab('Total junction coverage + 2') +
  ylab('Junction count + 1') +
  labs(color='Called as outlier')  +
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", size=0.3),
    legend.position = c(.05, .95),
    legend.justification = c("left", "top"),
    legend.margin = margin(0, 8, 4, 4),
    legend.title = element_blank()
  ) + 
  guides(colour = guide_legend(override.aes = list(size = 2)))

p_e




#### g. bam cov, FRASER ####
#+ plot p_g 1, fig.width=8, fig.height=8
# p_g <- plotBamCoverageFromResultTable(fds_p_e, result = res_p_e)




#### h. AML, MLL panel genes strat by num total ####
fr_res_aml_panel <- fr_res[sampleID %in% samp_anno_aml[, ArrayID] &
                             hgncSymbol %in% mll_panel_genes, ]
absplice_res_aml_panel <- absplice_res[tissue=='AML' & gene_name %in% mll_panel_genes,]
absplice_res_aml_panel <- absplice_res_aml_panel[AbSplice_DNA >= 0.2]

AML_variants_unique <- unique(AML_variants[, .(array_id, symbol)])
AML_variants_unique_panel <- AML_variants_unique[symbol %in% mll_panel_genes]

cnv_aml <- fread(grep('AML', snakemake@input$mll_cnv, value = TRUE))
cnv_aml_panel <- cnv_aml[SYMBOL %in% mll_panel_genes, ]
cnv_aml_panel <- cnv_aml_panel[abs(MEAN_LOG2_COPY_RATIO)>0.3, ]
cnv_aml_panel[, LENGTH := END - START]
cnv_aml_panel <- cnv_aml_panel[LENGTH >= 1000000, ]

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

fusion_cnv_variants_fr_absplice_union <- rbind(unique(fusion_aml_panel[, .(array_id, symbol)]),
                                               unique(cnv_aml_panel[, .(ARRAY_ID, SYMBOL)]),
                                               AML_variants_unique_panel, 
                                               fr_res_aml_panel[, .(sampleID, hgncSymbol)], 
                                               absplice_res_aml_panel[, .(sampleID, gene_name)], 
                                               use.names=FALSE)
fusion_cnv_variants_fr_absplice_union <- unique(fusion_cnv_variants_fr_absplice_union[symbol %in% mll_panel_genes])

fusion_cnv_union <- as.data.table(table(sort(table(fusion_cnv_union[, array_id]))))
colnames(fusion_cnv_union) <- c("n", "fusion_cnv_union")

fusion_cnv_variants_union <- as.data.table(table(sort(table(fusion_cnv_variants_union[, array_id]))))
colnames(fusion_cnv_variants_union) <- c("n", "fusion_cnv_variants_union")

fusion_cnv_variants_fr_absplice_union <- as.data.table(table(sort(table(fusion_cnv_variants_fr_absplice_union[, array_id]))))
colnames(fusion_cnv_variants_fr_absplice_union) <- c("n", "fusion_cnv_variants_fr_absplice_union")

df_merge <- merge(fusion_cnv_union, fusion_cnv_variants_union, by="n", all.x = TRUE, all.y=TRUE)
df_merge <- merge(df_merge, fusion_cnv_variants_fr_absplice_union, by="n", all.x = TRUE, all.y=TRUE)

df_merge$n <- as.integer(df_merge$n)
df_merge <- df_merge[order(df_merge$n),]
num_of_aml_samples <- length(unique(samp_anno_aml$ArrayID))

df_merge <- rbindlist(list(data.table("n" = 0, 
                                      "fusion_cnv_union" = num_of_aml_samples - sum(df_merge$fusion_cnv_union,na.rm=TRUE ),
                                      "fusion_cnv_variants_union" = num_of_aml_samples - sum(df_merge$fusion_cnv_variants_union,na.rm=TRUE),
                                      "fusion_cnv_variants_fr_absplice_union" = num_of_aml_samples - sum(df_merge$fusion_cnv_variants_fr_absplice_union,na.rm=TRUE )),
                           df_merge))
dt_merge <- df_merge %>% as.data.table()
#' number of samples
colSums(dt_merge[n>=5, ], na.rm = TRUE)

subset_n <- 10
df_merge[n==subset_n, "fusion_cnv_variants_fr_absplice_union"] <- sum(df_merge[n >= subset_n, "fusion_cnv_variants_fr_absplice_union"], na.rm = TRUE)
df_merge[n==subset_n, "fusion_cnv_union"] <- sum(df_merge[n >= subset_n, "fusion_cnv_union"], na.rm = TRUE)
df_merge[n==subset_n, "fusion_cnv_variants_union"] <- sum(df_merge[n >= subset_n, "fusion_cnv_variants_union"], na.rm = TRUE)
df_merge <- df_merge[n <= subset_n , ]
df_merge$n <- as.character(df_merge$n)
df_merge[n==subset_n, "n"] <- paste0("\u2265", subset_n)
df_merge$n <- factor(df_merge$n, levels = df_merge$n)

df_merge_melted <- melt(df_merge, id.vars = c("n"),
                        measure.vars = c("fusion_cnv_union", "fusion_cnv_variants_union", "fusion_cnv_variants_fr_absplice_union" ))

p_h_raw <- ggplot(df_merge_melted, aes(factor(n), value, fill=variable))+
  geom_bar(stat="identity", position="dodge")

p_h <- p_h_raw +
  xlab("Number of predicted aberrant hematologic panel genes per sample") + 
  ylab("Number of samples in AML") +
  scale_fill_manual('Stratified based on',
                    labels=c('fusion_cnv_union'='Fusions or CNV ',
                             'fusion_cnv_variants_union'='Fusions or CNV \nor Curated variant',
                             'fusion_cnv_variants_fr_absplice_union'='Fusions or CNV \nor Curated variants \nor Splicing aberrations'
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

# p_h




#### s5. upserR ####
listInput <- list("Splice variants (VEP)" = vep_splice[, samp_symbol], 
                  "Splicing outlier junctions (FRASER)" = fr_res_junc[, samp_symbol], 
                  "Splice-affecting variants (AbSplice)" = absplice_res_var[AbSplice_DNA>=0.2, samp_symbol])
p_s5_raw <- upset(fromList(listInput), order.by = "freq")

p_s5 <- p_s5_raw
# p_s5




### arrange plot #### 
#' ### raw
#+ plot p3 raw, fig.width=10, fig.height=12
p3_top_raw <- grid.arrange(p_a_raw, p_b_raw, 
                           p_c_raw, p_i_raw, 
                           nrow = 2)
p3_bottom_raw <- grid.arrange(p_d_raw, ggplot() + theme(panel.background = element_blank()), 
                              nrow = 1)

p3_raw <- grid.arrange(p3_top_raw, p3_bottom_raw, 
                       heights = c(8, 4), nrow = 2)
p3_raw

#' ### annotated
#+ plot p3, fig.width=10, fig.height=12
p3_top <- grid.arrange(p_a, p_b, 
                       p_c, p_i, 
                       nrow = 2)
p3_bottom <- grid.arrange(p_d, ggplot() + theme(panel.background = element_blank()),
                          nrow = 1)

p3 <- grid.arrange(p3_top, p3_bottom, 
                   heights = c(8, 4), nrow = 2)
p3

pdf(paste0(output_dir, "/figure_3.pdf"), width = 10, height = 12)

p3_top <- grid.arrange(p_a, p_b, 
                       p_c, p_i, 
                       nrow = 2)
p3_bottom <- grid.arrange(p_d, ggplot() + theme(panel.background = element_blank()),
                          nrow = 1)

p3 <- grid.arrange(p3_top, p3_bottom, 
                   heights = c(8, 4), nrow = 2)
p3

grid.text("A", x = 0.01, y=0.99, gp=gpar(fontsize=12, fontface='bold'))
grid.text("B", x = 0.51, y=0.99, gp=gpar(fontsize=12, fontface='bold'))
grid.text("C", x = 0.01, y=0.67, gp=gpar(fontsize=12, fontface='bold'))
grid.text("D", x = 0.51, y=0.67, gp=gpar(fontsize=12, fontface='bold'))
grid.text("E", x = 0.01, y=0.34, gp=gpar(fontsize=12, fontface='bold'))
grid.text("F", x = 0.51, y=0.34, gp=gpar(fontsize=12, fontface='bold'))
dev.off() 



png(paste0(output_dir, "/figure_3.png"), width = 10, height = 12, units = "in", res = 600)

p3_top <- grid.arrange(p_a, p_b, 
                       p_c, p_i, 
                       nrow = 2)
p3_bottom <- grid.arrange(p_d, ggplot() + theme(panel.background = element_blank()),
                          nrow = 1)

p3 <- grid.arrange(p3_top, p3_bottom, 
                   heights = c(8, 4), nrow = 2)
p3

grid.text("A", x = 0.01, y=0.99, gp=gpar(fontsize=12, fontface='bold'))
grid.text("B", x = 0.51, y=0.99, gp=gpar(fontsize=12, fontface='bold'))
grid.text("C", x = 0.01, y=0.67, gp=gpar(fontsize=12, fontface='bold'))
grid.text("D", x = 0.51, y=0.67, gp=gpar(fontsize=12, fontface='bold'))
grid.text("E", x = 0.01, y=0.34, gp=gpar(fontsize=12, fontface='bold'))
grid.text("F", x = 0.51, y=0.34, gp=gpar(fontsize=12, fontface='bold'))
dev.off() 




#' ## Sup
### Supplement #### cancer
#### s1. abSplice, leu13, enrichment ####
#' ### s1 raw
en_ab_dt[, num_samp_abs_per_gene := factor(num_samp_abs_per_gene, levels = c('>=8', as.character(7:0)))]
p_s1_raw <- ggplot(en_ab_dt[gene_list %in% c('leu', 'CGCleu'), ], aes(x=num_samp_abs_per_gene, y=enrichment)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=en05, ymax=en95), width=.2,
                position=position_dodge(.9)) +
  facet_wrap('gene_list', ncol=1) +
  geom_hline(yintercept = 1) +
  scale_y_log10()

#+ plot s1 raw, fig.width=6, fig.height=4
p_s1_raw

#' ### s1 annotated
p_s1 <- p_s1_raw +
  xlab('Number of samples predicted of the gene') +
  ylab('Enrichment') +
  theme_vale 

#+ plot s1, fig.width=6, fig.height=4
p_s1




#### s2. FRASER, all cohorts, ocg/tsg #### 
#' ### s2 raw
en_fr_leu13_dt <- en_fr_dt
en_fr_leu13_dt[, x_lab := as.character(k)]
en_fr_leu13_dt[, sample_subset_max_k := max(k), by='sample_subset']
en_fr_leu13_dt <- en_fr_leu13_dt[k==3|k==sample_subset_max_k , ]
en_fr_leu13_dt[k==sample_subset_max_k, x_lab := 'max']
en_fr_leu13_dt[, x_lab := factor(x_lab, levels=c('3', 'max'))]

p_s2_raw <- ggplot(en_fr_leu13_dt, aes(x=x_lab, y=enrichment)) +
  geom_bar(stat="identity", width = 0.3, fill="lightgrey") +
  geom_errorbar(aes(ymin=en05, ymax=en95), width=0.1) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  
  facet_grid(cancer_driver_list~sample_subset)  +
  
  scale_y_log10()

#+ plot s2 raw, fig.width=10, fig.height=7
p_s2_raw

#' ### s2 annotated
p_s2 <- p_s2_raw +
  scale_x_discrete(labels = c('top 3 splicing \n outliers in each sample',
                              'all significant \n splicing outliers')) +
  
  # ggtitle(paste0(sample_subset_selected, ' (n=', cohort_dt[Cohort_group==sample_subset_selected, N] ,')')) +
  xlab('') +
  ylab('Enrichmnet') +
  theme_vale +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#+ plot s2, fig.width=10, fig.height=7
p_s2




#### s5. upserR ####
#' ### s5 raw
#+ plot s5 raw, fig.width=6, fig.height=4
p_s5_raw

#' ### s5 annotated
#+ plot s5, fig.width=6, fig.height=4
p_s5




#### s6. sample rank fraser gene #####
#' ### s6 raw
#+ plot s6 raw, fig.width=6, fig.height=4
hlines <- quantile(fr_samp_freq[, N], c(0.5, 0.9)) + 1

p_s6_raw <- ggplot(fr_samp_freq, aes(x=sample_rank, y=N+1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=hlines) + 
  scale_y_log10(limits = c(1, fr_samp_freq[, 10^ceiling(log10(max(N)))]+1),
                breaks = c(1, 2, 4, 11, 31, 101, 301, 1001), 
                labels = c(0, 1, 3, 10, 30, 100, 300, 1000)) +
  annotate("text", label=c("Median", "90^th ~ percentile"), 
           x=1, y=hlines*1.2, hjust=0, parse=TRUE)

p_s6_raw

#' ### s6 annotated
#+ plot s6, fig.width=6, fig.height=4
p_s6 <- p_s6_raw +
  xlab('Sample rank') +
  ylab('Called aberrantly spliced genes (FRASER)') +
  theme_vale 

p_s6




#### s7. sample rank absplice gene #####
#' ### s7 raw
#+ plot s7 raw, fig.width=6, fig.height=4
hlines <- quantile(absplice_samp_freq[, N], c(0.5, 0.9)) + 1

p_s7_raw <- ggplot(absplice_samp_freq, aes(x=sample_rank, y=N+1)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=hlines) + 
  scale_y_log10(limits = c(1, 301),
                breaks = c(1, 2, 4, 11, 31, 101, 301), 
                labels = c(0, 1, 3, 10, 30, 100, 300)) +
  annotate("text", label=c("Median", "90^th ~ percentile"), 
           x=1, y=hlines*1.2, hjust=0, parse=TRUE)

p_s7_raw

#' ### s7 annotated
#+ plot s7, fig.width=6, fig.height=4
p_s7 <- p_s7_raw +
  xlab('Sample rank') +
  ylab('Predicted splice-affected genes (AbSplice)') +
  theme_vale 

p_s7



#### s8. upserR ####
#' ### s8 raw
#+ plot s8 raw, fig.width=12, fig.height=4
p_s8_raw

#' ### s8 annotated
#+ plot s8, fig.width=12, fig.height=4
p_s8


#### s9. upserR ####
#' ### s9 raw
#+ plot s9 raw, fig.width=12, fig.height=4
p_s9_raw

#' ### s9 annotated
#+ plot s9, fig.width=12, fig.height=4
p_s9



#### s10. upserR ####
#' ### s10 raw
#+ plot s10 raw, fig.width=12, fig.height=8
p_s10_raw

#' ### s10 annotated
#+ plot s10, fig.width=12, fig.height=8
p_s10



### thesis ####
png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "fr_sample_rank.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_s6
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "fr_enrich.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_a
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "rb1_count.png"), 
    width = 5, height = 5, units = "in", res = 600)
p_e
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "ab_sample_rank.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_s7
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "ab_enrich.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_b
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "splice_venn.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_c
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "fr_ab_enrich.png"), 
    width = 6, height = 4, units = "in", res = 600)
p_i
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "cd79a_count.png"), 
    width = 5, height = 5, units = "in", res = 600)
p_d
dev.off() 