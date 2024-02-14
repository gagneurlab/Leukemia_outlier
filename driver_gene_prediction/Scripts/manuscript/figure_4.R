#'---
#' title: figure_4
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
#'    - experimentDesign: '`sm config["experimentDesign"]`'
#'    - manuscriptWording: '`sm config["manuscriptWording"]`'
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - intogenDir: '`sm config["intogenDir"]`'
#'    - vep_path: '`sm config["vep_path"]`'   
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'    - CGC_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_gene_list.tsv"`'
#'    - MLL_CGC_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/MLL_CGC_leukemia_gene_list.tsv"`'
#'    - drug_target: '`sm config["drug_target"]`'
#'    - figure_4_heatmap: '`sm config["figure_4_heatmap"]`'
#'  input:
#'    - mll_prc: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_prc.tsv"`'
#'    - mll_ap: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_ap.tsv"`'
#'    - mll_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_ap_full.tsv"`'
#'    - mll_var_filter_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_var_filter_prc.tsv"`'
#'    - mll_var_filter_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_var_filter_ap.tsv"`'
#'    - mll_var_filter_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_var_filter_ap_full.tsv"`'
#'    - mll_emb_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_emb_prc.tsv"`'
#'    - mll_emb_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_emb_ap.tsv"`'
#'    - mll_emb_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_emb_ap_full.tsv"`'
#'    - mll_method_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_method_prc.tsv"`'
#'    - mll_method_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_method_ap.tsv"`'
#'    - mll_method_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_4/plot_data/mll_method_ap_full.tsv"`'
#'    - mll_benchmark_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_benchmark_prc.tsv"`'
#'    - mll_benchmark_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_benchmark_ap.tsv"`'
#'    - mll_benchmark_rp: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_benchmark_rp.tsv"`'
#'  output:
#'    - p_c: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_c.Rds"`' 
#'    - p_d: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_d.Rds"`' 
#'    - CGC_leukemia_gene: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/CGC_leukemia_gene.csv"`'
#'    - MLL_CGC_leukemia_gene: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/MLL_CGC_leukemia_gene.csv"`'
#'    - gencode: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/gencode.csv"`'
#'    - mll_prc: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/mll_prc.tsv"`'
#'    - mll_ap: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/mll_ap.tsv"`'
#'    - mll_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/mll_ap_full.tsv"`'
#'    - rcl_list: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/rcl_list.Rds"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_4.html"`'
#'  type: noindex
#'  resources:
#'    - mem_mb: 8000 
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_4.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_4.snakemake")
print("Snakemake saved") 



suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(pheatmap)
  library(ComplexHeatmap)
  library(ggnewscale)
})
source("Scripts/manuscript/function.R")
source("Scripts/manuscript/manuscript_theme.R")

options(bitmapType='cairo')

benchmark_color <- c(RColorBrewer::brewer.pal(8, "Set2")[6],
                     RColorBrewer::brewer.pal(12, "Paired")[c(1:7, 10:12, 8:9)])
names(benchmark_color) <- c("mutsigcv", 
                            "cbase", 
                            "clustl", 
                            "dndscv",  
                            "fml", 
                            "hotmaps", 
                            "mutpanning", 
                            "smregions",
                            "intogen_wo_post", 
                            "intogen", 
                            "oncovar", 
                            "vale_7tools",
                            "vale")
benchmark_label_dt <- data.table(Method = names(benchmark_color), 
                                 benchmark_label = c("MutSigCV", 
                                                     "CBaSE", 
                                                     "OncodriveCLUSTL", 
                                                     "dNdScv", 
                                                     "OncodriveFML", 
                                                     "HotMAPS", 
                                                     "MutPanning", 
                                                     "smRegions",
                                                     "intOGen: 7 tools", 
                                                     "intOGen", 
                                                     "oncoVar",
                                                     "Integrative model: 7 tools",
                                                     "Integrative model"))
benchmark_label_dt[, feature_source := 'Genomic']
benchmark_label_dt[Method %in% c("vale") , feature_source := 'Transcriptomic']

feat_label_dt <- data.table(input_feature = c("clustl",
                                              "hotmaps",
                                              "smregions",
                                              "fml",
                                              "cbase",
                                              "mutpanning",
                                              "dndscv",
                                              
                                              "clustl,hotmaps",
                                              "clustl,hotmaps,smregions",
                                              "clustl,hotmaps,smregions,fml",
                                              "clustl,hotmaps,smregions,fml,cbase",
                                              "clustl,hotmaps,smregions,fml,cbase,mutpanning",
                                              "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv",
                                              
                                              "absplice",
                                              "fr",
                                              "or,ac",
                                              
                                              "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice",
                                              "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice,fr",
                                              "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice,fr"
),
feat_label = c("  Linear cluster of mutations \n      (OncodriveCLUSTL)", 
               "  3D cluster of mutations \n      (HotMAPS)", 
               "  Cluster in domain \n      (smRegions)",
               "  Functional impact bias \n      (OncodriveFML)", 
               "  Excess of mutations \n      (CBaSE)", 
               "  Trinucleotide-specific bias \n      (MutPanning)", 
               "  Excess of mutations \n      (dNdScv)", 
               
               "+ 3D cluster of mutations \n      (HotMAPS)", 
               "+ Cluster in domain \n      (smRegions)",
               "+ Functional impact bias \n      (OncodriveFML)", 
               "+ Excess of mutations \n      (CBaSE)", 
               "+ Trinucleotide-specific bias \n      (MutPanning)", 
               "+ Excess of mutations \n      (dNdScv)", 
               
               "  Splicing variants \n      (AbSplice)", 
               "  Splicing outliers \n      (FRASER)", 
               "  Expression outliers \n      (OUTRIDER+NB-act)", 
               
               "+ Splicing variants \n      (AbSplice)", 
               "+ Splicing outliers \n      (FRASER)",
               "+ Expression outliers \n      (OUTRIDER+NB-act)"
))

feat_label_dt[, feature_source := 'Genomic']
feat_label_dt[input_feature %in% c(
  "or,ac",
  "fr",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice,fr",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice,fr"
) , feature_source := 'Transcriptomic']

feat_color <- c(
  'Genomic'=RColorBrewer::brewer.pal(12, "Paired")[2],
  'Transcriptomic'=RColorBrewer::brewer.pal(12, "Paired")[6]
)

curation_color <- c('Reported hematologic malignancy driver gene'=RColorBrewer::brewer.pal(8, "Paired")[4], 
                    'Reported cancer driver gene'=RColorBrewer::brewer.pal(8, "Paired")[3], 
                    'Candidate hematologic malignancy driver gene'=RColorBrewer::brewer.pal(8, "Paired")[1])


# get parameters
output_dir <- snakemake@params$htmlOutputPath
project_dir <- snakemake@params$projectPath
intogen_dir <- snakemake@params$intogenDir
vep_dir <- snakemake@params$vep_path

single_group <- snakemake@params$inputDatasets

sample_group <- snakemake@params$outputDatasets
sample_group_by_size <- sample_group
sample_group_by_cell_type <- c('AML', 'MDS', 'MPN', 'MDS_MPN_group', 'CMML', 
                               'CML', 'NK', 'MACS_group', 'Mast_group', 
                               'TL_group', 'CLL', 'Lym_group', 'HZL_group', 'BCP_ALL')

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]
fwrite(gencode, snakemake@output$gencode)
gencode_pr <- gencode[gene_type=='protein_coding', ]

experiment_design <- fread(snakemake@params$experimentDesign)
experiment_design[is.na(experiment_design)] <- ''

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))

mll_cgc_leukemia_gene <- fread(snakemake@params$MLL_CGC_leukemia_gene_list)
fwrite(mll_cgc_leukemia_gene, snakemake@output$MLL_CGC_leukemia_gene)
cgc_leukemia_gene <- fread(snakemake@params$CGC_leukemia_gene_list)
fwrite(cgc_leukemia_gene, snakemake@output$CGC_leukemia_gene)

mll_prc <- fread(snakemake@input$mll_prc)
mll_ap <- fread(snakemake@input$mll_ap)
mll_ap_full <- fread(snakemake@input$mll_ap_full)
fwrite(mll_prc, snakemake@output$mll_prc)
fwrite(mll_ap, snakemake@output$mll_ap)
fwrite(mll_ap_full, snakemake@output$mll_ap_full)

mll_benchmark_prc <- fread(snakemake@input$mll_benchmark_prc)
mll_benchmark_ap <- fread(snakemake@input$mll_benchmark_ap)
mll_benchmark_rp <- fread(snakemake@input$mll_benchmark_rp)

mll_var_filter_prc <- fread(snakemake@input$mll_var_filter_prc)
mll_var_filter_ap <- fread(snakemake@input$mll_var_filter_ap)
mll_var_filter_ap_full <- fread(snakemake@input$mll_var_filter_ap_full)

mll_emb_prc <- fread(snakemake@input$mll_emb_prc)
mll_emb_ap <- fread(snakemake@input$mll_emb_ap)
mll_emb_ap_full <- fread(snakemake@input$mll_emb_ap_full)

mll_method_prc <- fread(snakemake@input$mll_method_prc)
mll_method_ap <- fread(snakemake@input$mll_method_ap)
mll_method_ap_full <- fread(snakemake@input$mll_method_ap_full)

drug_target <- fread(snakemake@params$drug_target)

pharos_target <- fread('/s/project/vale/Resource/pharos_targets_0240205.csv')
colnames(pharos_target) <- gsub(" ", "_", colnames(pharos_target))
pharos_target_pr <- pharos_target[Symbol %in% gencode_pr[, gene_name_orig], ]




### Figure #### 
#### a. ap individual, mll #### 
#' number of RNAseq
samp_anno_exp[DROP_GROUP==single_group, .N]

setnames(mll_ap, 'Training_setup', 'input_feature')
mll_ap <- merge(mll_ap, feat_label_dt, by='input_feature')
mll_ap[, input_feature := factor(input_feature, levels = rev(feat_label_dt[, input_feature]))]

feat_label <- mll_ap[, feat_label]
names(feat_label) <- mll_ap[, input_feature]

# individual tools
individual_tool <- c("fml", "clustl", "hotmaps", "smregions", "cbase", "dndscv", "mutpanning", 
                     "absplice", "or,ac", "fr")
res_ap <- mll_ap[input_feature %in% individual_tool,]
# res_ap[, input_feature := factor(input_feature, levels = ]

p_a_raw <- ggplot(res_ap, aes(input_feature, AP_median, fill=feature_source)) + 
  geom_errorbar(aes(ymin=AP_1q, ymax=AP_9q), color='gray20', width = 0.2, linewidth = 0.3) +
  geom_col(alpha = 0.3) +
  coord_flip()

p_a <- p_a_raw +  
  scale_fill_manual(values = feat_color) +
  scale_x_discrete(labels = feat_label) +
  scale_y_continuous(limits = c(0, 0.27), breaks=c(0, 0.1, 0.2)) +
  xlab("Individual feature") + 
  ylab("Average Precision") +
  theme_vale +
  theme(
    legend.position = 'bottom',
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  guides(fill=guide_legend("Feature source", title.position = "top"))

# p_a




#### b. ap cumulative, mll #### 
setnames(mll_ap_full, 'Training_setup', 'input_feature')
selected_input_feature <- c(
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice,fr"
)
selected_input_feature <- c(
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice,fr",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice,fr"
)

stat.test <- ggpubr::compare_means(formula = AP ~ input_feature, 
                                   data = mll_ap_full[input_feature %in% selected_input_feature], 
                                   method = "wilcox.test")
stat.test$label <- ''

# add on top
on_top <- c(
  "clustl",
  "clustl,hotmaps",
  "clustl,hotmaps,smregions",
  "clustl,hotmaps,smregions,fml",
  "clustl,hotmaps,smregions,fml,cbase",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice,fr",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice,fr")
res_ap <- mll_ap[input_feature %in% on_top,]

p_b_raw <- ggplot(res_ap, aes(input_feature, AP_median)) + 
  geom_errorbar(aes(ymin=AP_1q, ymax=AP_9q), color='gray20', width = 0.2, linewidth = 0.3) +
  geom_col(aes(fill=feature_source), alpha = 0.3) +
  stat_pvalue_manual(stat.test[c(1,5),], label="label", 
                     y.position = 0.25, bracket.shorten = -0.1, color='gray20') +
  annotate("text", x=3.5, y=0.27, label='***', size = 4, angle=90, color='gray20') +
  annotate("text", x=2, y=0.27, label='**', size = 4, angle=90, color='gray20') +
  coord_flip()

# p_b_raw

p_b <- p_b_raw +  
  scale_fill_manual(values = feat_color) +
  scale_x_discrete(labels = feat_label) +
  scale_y_continuous(limits = c(0, 0.27), breaks=c(0, 0.1, 0.2)) +
  xlab("Cumulative feature") + 
  ylab("Average Precision") +
  theme_vale +
  theme(
    legend.position = 'bottom',
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )  +
  guides(fill=guide_legend("Feature source", title.position = "top"))  

# p_b




#### s1. prc individual, mll #### 
setnames(mll_prc, 'Training_setup', 'input_feature')
mll_prc <- merge(mll_prc, mll_ap[, .(input_feature, AP_median)], by='input_feature')
mll_prc <- merge(mll_prc, feat_label_dt, by='input_feature')
mll_prc[, feat_label_ap := paste0(feat_label, " (AP median = ", round(AP_median, 4), ")")]
mll_prc[, input_feature := factor(input_feature, levels = feat_label_dt[, input_feature])]

feat_label_ap <- mll_prc[, unique(feat_label_ap)]
names(feat_label_ap) <- mll_prc[, unique(input_feature)]

# individual tools
res_prc <- mll_prc[input_feature %in% individual_tool, ] # each tool

p_s1_raw <- ggplot(res_prc, aes(recall, precision_median)) + 
  geom_ribbon(aes(ymax=precision_up, ymin=precision_dn, fill=input_feature), alpha=0.3) +
  geom_step(aes(color=input_feature), direction='vh') 

p_s1 <- p_s1_raw + 
  scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired"), 
                     labels = feat_label_ap, name = 'Individual feature') +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Paired"), guide="none") +
  xlab("Recall") + 
  ylab("Precision") +
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6)
  )

# p_s1




#### s2. prc cumulative, mll #### 
# add on top
res_prc <- mll_prc[input_feature %in% on_top, ] # on top

p_s2_raw <- ggplot(res_prc, aes(recall, precision_median)) + 
  geom_ribbon(aes(ymax=precision_up, ymin=precision_dn, fill=input_feature), alpha=0.3) +
  geom_step(aes(color=input_feature), direction='vh') 

p_s2 <- p_s2_raw + 
  scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired"), 
                     labels = feat_label_ap, name = 'Cumulative feature') +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Paired"), guide="none") +
  xlab("Recall") + 
  ylab("Precision") +
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "left",
    legend.margin = margin(6, 6, 6, 6)
  )

# p_s2




#### s3. prc, mll benchmark #### 
res_benchmark_ap <- merge(mll_benchmark_ap, benchmark_label_dt, by='Method')

res_benchmark_ap_sub_1 <- res_benchmark_ap[Method %in% c('vale'), ]
res_benchmark_ap_sub_2 <- res_benchmark_ap[Method %in% c('intogen', 
                                                         'mutsigcv'), ]

res_benchmark_ap_sub <- merge(res_benchmark_ap_sub_1, res_benchmark_ap_sub_2, 
                              by='Sample_group', all = TRUE)
res_benchmark_ap_sub <- merge(res_benchmark_ap_sub, 
                              manuscript_wording[, .(Study_group, Study_group_during_analysis)] %>% unique(), 
                              by.x = 'Sample_group', by.y = 'Study_group_during_analysis', all.x = TRUE)
res_benchmark_ap_sub[Sample_group=='leukemia_14group', Study_group := 'Complete dataset']
res_benchmark_ap_sub[, benchmark_label.y := factor(benchmark_label.y, 
                                                   levels=c('MutSigCV', 'intOGen'))]

p_s3_raw <- ggplot(res_benchmark_ap_sub, aes(x=AP.y, y=AP.x, label=Study_group)) +
  geom_abline(slope=1, intercept=0, linetype='dashed', color='darkgray') +
  geom_point() +
  geom_text_repel(size=2) +
  facet_grid(benchmark_label.x ~ benchmark_label.y) +
  coord_fixed()

# p_s3_raw

p_s3 <- p_s3_raw + 
  ylab('Average precision') + 
  xlab('Average precision') +
  xlim(0, 0.125) +
  ylim(0, 0.125) +
  theme_vale

# p_s3




#### s4. prc var filter, mll #### 
# test if improvement sigficant
setnames(mll_var_filter_ap_full, 'Training_setup', 'input_feature')
stat.test <- ggpubr::compare_means(formula = AP ~ input_feature, 
                                   data = mll_var_filter_ap_full[input_feature %in% c("driver_prediction_202401_0",
                                                                                      "driver_prediction_202401_1",
                                                                                      "driver_prediction_202401_4",
                                                                                      "driver_prediction_202401_2"), ], 
                                   method = "wilcox.test")
stat.test

var_filter_label_dt <- data.table(
  input_feature = c(
    "driver_prediction_202401_0", 
    "driver_prediction_202401_1", 
    "driver_prediction_202401_2", 
    "driver_prediction_202401_3",
    "driver_prediction_202401_4"
  ),
  var_filter_label = c(
    "keep VAF >= 0.1\n",
    "keep VAF >= 0.1 ; selected consequence\n",
    "keep VAF >= 0.15; selected consequence; sequencing depth>=20\n",
    "keep VAF >= 0.1 ; selected consequence; sequencing depth>=20\n",
    "keep VAF >= 0.15; selected consequence\n"
  ))


setnames(mll_var_filter_prc, 'Training_setup', 'input_feature')
setnames(mll_var_filter_ap, 'Training_setup', 'input_feature')
mll_var_filter_prc <- merge(mll_var_filter_prc, mll_var_filter_ap[, .(input_feature, AP_median)], by='input_feature')
mll_var_filter_prc <- merge(mll_var_filter_prc, var_filter_label_dt, by='input_feature')
mll_var_filter_prc[, var_filter_label_ap := paste0(var_filter_label, " (AP median = ", round(AP_median, 4), ")\n")]
mll_var_filter_prc[, input_feature := factor(input_feature, levels = c("driver_prediction_202401_0",
                                                                       "driver_prediction_202401_1",
                                                                       "driver_prediction_202401_4",
                                                                       "driver_prediction_202401_2",
                                                                       "driver_prediction_202401_3"))]

var_filter_label_ap <- mll_var_filter_prc[, unique(var_filter_label_ap)]
names(var_filter_label_ap) <- mll_var_filter_prc[, unique(input_feature)]

# individual tools
res_prc <- mll_var_filter_prc[input_feature %in% c("driver_prediction_202401_0",
                                                   "driver_prediction_202401_1",
                                                   "driver_prediction_202401_4",
                                                   "driver_prediction_202401_2"), ] # each tool

p_s4_raw <- ggplot(res_prc, aes(recall, precision_median)) + 
  geom_ribbon(aes(ymax=precision_up, ymin=precision_dn, fill=input_feature), alpha=0.3) +
  geom_step(aes(color=input_feature), direction='vh') 

p_s4 <- p_s4_raw + 
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), 
                     labels = var_filter_label_ap, name = 'Additional variant filter criteria:') +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), guide="none") +
  xlab("Recall") + 
  ylab("Precision") +
  ggtitle("Variant filter criteria:\nQUALITY==’PASS’; discard MAF >= 0.000,5") +
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = 'right',
    legend.margin = margin(6, 6, 6, 6)
  )

# p_s4




#### s5. prc emb, mll #### 
# test if improvement sigficant
setnames(mll_emb_ap_full, 'Training_setup', 'input_feature')
stat.test <- ggpubr::compare_means(formula = AP ~ input_feature, 
                                   data = mll_emb_ap_full, 
                                   method = "wilcox.test")
stat.test
stat.test[c(1,2,3,5),]

emb_label_dt <- data.table(
  input_feature = c(
    "",
    "coess_cluster",
    "emb_omics",      
    "emb_omics,emb_string", 
    "emb_string"     
  ),
  emb_label = c(
    "7 intOGen tools, OUTRIDER, NB-act, FRASER, AbSplice
    ",
    "7 intOGen tools, OUTRIDER, NB-act, FRASER, AbSplice
    + Co-essential modules
    ",
    "7 intOGen tools, OUTRIDER, NB-act, FRASER, AbSplice
    + Embedding Omics
    ",      
    "7 intOGen tools, OUTRIDER, NB-act, FRASER, AbSplice
    + Embedding Omics + Embedding STRING
    ",     
    "7 intOGen tools, OUTRIDER, NB-act, FRASER, AbSplice
    + Embedding STRING
    "     
  ))


setnames(mll_emb_prc, 'Training_setup', 'input_feature')
setnames(mll_emb_ap, 'Training_setup', 'input_feature')
mll_emb_prc <- merge(mll_emb_prc, mll_emb_ap[, .(input_feature, AP_median)], by='input_feature')
mll_emb_prc <- merge(mll_emb_prc, emb_label_dt, by='input_feature')
mll_emb_prc[, emb_label_ap := paste0(emb_label, " (AP median = ", round(AP_median, 4), ")\n")]
mll_emb_prc[, input_feature := factor(input_feature, levels = rev(emb_label_dt[, input_feature]))]

emb_label_ap <- mll_emb_prc[, unique(emb_label_ap)]
names(emb_label_ap) <- mll_emb_prc[, unique(input_feature)]

# individual tools
res_prc <- mll_emb_prc[input_feature %in% c(
  '',
  'coess_cluster', 
  "emb_omics",      
  "emb_string",     
  "emb_omics,emb_string" 
), ]

p_s5_raw <- ggplot(res_prc, aes(recall, precision_median)) + 
  # geom_ribbon(aes(ymax=precision_up, ymin=precision_dn, fill=input_feature), alpha=0.3) +
  geom_step(aes(color=input_feature), direction='vh')

p_s5 <- p_s5_raw + 
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), 
                     labels = emb_label_ap, name = 'Feature',
                     guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), guide="none") +
  xlab("Recall") + 
  ylab("Precision") +
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = 'right',
    legend.margin = margin(6, 6, 6, 6)
  ) 

# p_s5




#### s8. prc method, mll #### 
# test if improvement sigficant
setnames(mll_method_ap_full, 'Training_setup', 'input_feature')
stat.test <- ggpubr::compare_means(formula = AP ~ input_feature, 
                                   data = mll_method_ap_full, 
                                   method = "wilcox.test")
wilcox.test(mll_method_ap_full[input_feature=='nn', AP],
            mll_method_ap_full[input_feature=='lr', AP])
wilcox.test(mll_method_ap_full[input_feature=='lr', AP],
            mll_method_ap_full[input_feature=='xgb_op', AP])
wilcox.test(mll_method_ap_full[input_feature=='xgb_op', AP],
            mll_method_ap_full[input_feature=='rf', AP])

stat.test

method_label_dt <- data.table(
  input_feature = c(
    "nn",
    "lr",
    "xgb_op",      
    "rf"
  ),
  method_label = c(
    "Neural network",
    "Logistic regression",
    "XGboost",      
    "Random forest"
  ))


setnames(mll_method_prc, 'Training_setup', 'input_feature')
setnames(mll_method_ap, 'Training_setup', 'input_feature')
mll_method_prc <- merge(mll_method_prc, mll_method_ap[, .(input_feature, AP_median)], by='input_feature')
mll_method_prc <- merge(mll_method_prc, method_label_dt, by='input_feature')
mll_method_prc[, method_label_ap := paste0(method_label, " (AP median = ", round(AP_median, 4), ")\n")]
mll_method_prc[, input_feature := factor(input_feature, levels = rev(method_label_dt[, input_feature]))]

method_label_ap <- mll_method_prc[, unique(method_label_ap)]
names(method_label_ap) <- mll_method_prc[, unique(input_feature)]

# individual tools
res_prc <- mll_method_prc

p_s8_raw <- ggplot(res_prc, aes(recall, precision_median)) + 
  # geom_ribbon(aes(ymax=precision_up, ymin=precision_dn, fill=input_feature), alpha=0.3) +
  geom_step(aes(color=input_feature), direction='vh')

p_s8 <- p_s8_raw + 
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), 
                     labels = method_label_ap, name = 'Feature',
                     guide = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), guide="none") +
  xlab("Recall") + 
  ylab("Precision") +
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = 'right',
    legend.margin = margin(6, 6, 6, 6)
  ) 

# p_s8




#### s9. prc individual, mll #### 
mll_benchmark_rp <- merge(mll_benchmark_rp, benchmark_label_dt, by='Method')
mll_benchmark_rp[, rp_label := benchmark_label]
mll_benchmark_rp[, Method := factor(Method,levels=c("intogen", "mutpanning", "dndscv", "fml",
                                                    "cbase", "smregions", "hotmaps", "clustl"))]

manuscript_wording[, Study_group_number := paste0(Study_group, " (n=", Number_of_samples_per_study_group, ")")]

mll_benchmark_rp <- merge(mll_benchmark_rp, 
                          manuscript_wording[, .(Study_group_number, Study_group_during_analysis)] %>% unique(), 
                          by.x = 'Sample_group', by.y = 'Study_group_during_analysis', all.x = TRUE)
mll_benchmark_rp[, Study_group_number := factor(Study_group_number,levels=c(
  "AML (n=730)", "MDS (n=713)", "MatureB group (n=375)"))]

rp_label <- c('Combination (intOGen)', 'MutPanning', 'dNdScv', 'OncodriveFML',
              'CBaSE', 'smRegions', 'HotMAPS', 'OncodriveCLUSTL')
names(rp_label) <- c("intogen", "mutpanning", "dndscv", "fml",
                     "cbase", "smregions", "hotmaps", "clustl")

# individual tools
res_rp <- mll_benchmark_rp[Criteria=='isCGC',] # each tool

p_s9_raw <- ggplot(res_rp, aes(Rank, Proportion, color=Method)) + 
  geom_line() +
  facet_wrap('Study_group_number') 

# p_s9_raw

p_s9 <- p_s9_raw + 
  scale_color_manual(values = c('darkgrey', 'brown', 'green', 'red',
                                'lightgreen', 'lightblue', 'blue', 'purple'), 
                     labels = rp_label, name = 'Method') +
  scale_y_continuous(limits = c(0, 1), breaks=c(0:5)*0.2, labels = scales::percent_format(accuracy = 1)) + 
  scale_x_continuous(limits = c(0, 140), breaks=c(0:7)*20) +
  ylab("Proportion of CGC elements") +
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = 'right',
    legend.margin = margin(6, 6, 6, 6)
  ) 

# p_s9




#### d2. 14 groups, prediction heatmap ####
n_plot <- 100
# read in res_viz
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == single_group &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == '', ]
exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)

res_viz_pred <- fread(exp_viz[, res_post_path])
res_viz_pred[, Sample_group := single_group]
res_viz_pred[, Prediction_post_scaled := Prediction_post/max(Prediction_post)]
res_viz_pred <- res_viz_pred[, .(GeneSymbol, Prediction_post, Prediction_post_scaled, Sample_group, isCGC, isLeukemia)]

selected_genes <- res_viz_pred[1:n_plot, GeneSymbol]
res_viz_pred <- res_viz_pred[GeneSymbol %in% selected_genes, ]

# read in res_viz_groups
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group %in% sample_group_by_size &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == '', ]
exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)

group_label <- manuscript_wording[, unique(Study_group)]
names(group_label) <- manuscript_wording[, unique(Study_group_during_analysis)]

# read in result for visualization
if(exists("res_viz_groups")){remove(res_viz_groups)}

for (i in 1:nrow(exp_viz)) {
  res_temp <- fread(exp_viz[i, res_post_path])
  res_temp[, Sample_group := exp_viz[i, sample_group]]
  res_temp[, Prediction_post_scaled := Prediction_post/max(Prediction_post)]
  
  if(!exists("res_viz_groups")){
    res_viz_groups <- res_temp
  }else{
    res_viz_groups <- rbind(res_viz_groups, res_temp, fill=TRUE)
  }
  
}

res_viz_groups_pred <- res_viz_groups[, .(GeneSymbol, Prediction_post, Prediction_post_scaled, Sample_group, isCGC, isLeukemia)]
res_viz_groups_pred <- res_viz_groups_pred[GeneSymbol %in% selected_genes, ]

res_viz_groups_pred[, Sample_group := factor(Sample_group, levels = rev(sample_group_by_cell_type))]
res_viz_groups_pred[, Sum_prediction_post := sum(Prediction_post), by='GeneSymbol']

res_viz_groups_pred[, curation := 'Candidate hematologic malignancy driver gene']
res_viz_groups_pred[isCGC==TRUE, curation := 'Reported cancer driver gene']
res_viz_groups_pred[GeneSymbol %in% mll_cgc_leukemia_gene[, GeneSymbol], curation := 'Reported hematologic malignancy driver gene']
res_viz_groups_pred[, curation := factor(curation, levels = names(curation_color))]

# heatmap
heatmap_dt <- res_viz_groups_pred[, .(Sample_group, GeneSymbol, Prediction_post_scaled)] 
heatmap_dt <- dcast(heatmap_dt, GeneSymbol ~ Sample_group, value.var = 'Prediction_post_scaled', fill=NA)


num_gene_clusters <- 5
num_sample_clusters <- 2

# We ran the following chunck several times to get a good-looking heatmap. No need to run now
# as we read in the order of final selected heatmap after this chunck
# get row order
# heatmap_mtx <- as.matrix(heatmap_dt[,-1])
# rownames(heatmap_mtx) <- heatmap_dt[, GeneSymbol]
# set.seed(512)
# p_d2_phm <- pheatmap(heatmap_mtx,
#                      row_km = num_gene_clusters,
#                      column_km = num_sample_clusters
# )
# p_d2_phm <- draw(p_d2_phm)

# saveRDS(p_d2_phm, file=snakemake@params$figure_4_heatmap)

heatmap_path <- "/s/project/vale/Resource/figure_4_heatmap.rds"
p_d2_phm <- readRDS(heatmap_path)
#p_d2_phm <- readRDS(snakemake@params$figure_4_heatmap)

r.dend <- row_dend(p_d2_phm)
rcl.list <- row_order(p_d2_phm)


# trying new row orders to make the heatmap more beautiful
rcl.list <- c(rcl.list[1:2], rev(rcl.list[3:5]))
saveRDS(rcl.list, snakemake@output$rcl_list)

# sort heatmap based on dendogram order
heatmap_dt <- heatmap_dt[rev(unlist(rcl.list )),]
gene_cluster_ranges <- rev(as.integer(lapply(rcl.list, function(x) length(x))))

colnames(heatmap_dt)[2:15] <- sapply(colnames(heatmap_dt)[2:15] , function(x){
  manuscript_wording[Study_group_during_analysis==x, unique(Study_group)]
}) %>% unlist()

heatmap_dt$y_min = c(0:(n_plot - 1)) 
heatmap_dt$y_max = c(1 : n_plot ) 


# Specify column orders here. 0 is the most left column and 13 is the most right column. 
# This is done by hand and not by col_order of heatmap.
column_ords <- c(
  'AML' = 0,
  'MDS' = 1,
  'MPN' = 2,
  'MDS/MPN group' = 3,
  'CMML' = 4,
  'Mastocytosis' = 5,
  'CML' = 6,
  'NK' = 7,
  'PCN group' = 8,
  'CLL' = 9,
  'BCP-ALL' = 10,
  'Hairy-cell group' = 11,
  'T-cell group' = 12,
  'MatureB group' = 13
)
# p.s. clustering of sample annotations is also done manually. the first 7 groups are marked as Myeloid
# and the rest as Lymphoid


ml <- melt(heatmap_dt, id.vars = c("GeneSymbol", "y_min", "y_max"))

ml$x_min <- sapply(ml$variable , function(x){
  column_ords[as.character(x)] - 0.5
})
ml$x_max <- ml$x_min + 1


# create gene annotation legend
y_maxes <- cumsum(gene_cluster_ranges) 
left_annot <- data.frame(
  x_min = rep(13.6, num_gene_clusters),
  x_max = rep(14.3, num_gene_clusters),
  y_min = c(0, y_maxes[-length(y_maxes)]),
  y_max = y_maxes,
  Gene_Grouping = factor(paste0("Cluster ", rev(1:num_gene_clusters)))
  
)

drug_annot <- merge(ml[!duplicated(ml$GeneSymbol), ], drug_target, by.x = "GeneSymbol", by.y = "Gene")
drug_annot[, x_min := rep(14.3, nrow(drug_annot))]
drug_annot[, x_max := rep(15, nrow(drug_annot))]
drug_annot[, Approved_drug := ifelse(Target %in% c("Y", "Y/N"), "Available", "Not available")]


# create sample annotation legend 
top_annot <- data.frame(
  x_min = c(-0.5, 6.5),
  x_max = c(6.5, 13.5),
  y_min = c(-0.2, -0.2),
  y_max = c(-1.7, -1.7),
  Cell_origin = factor(c( "Myeloid", "Lymphoid"))
  
)



# define breaks and colors for sample, gene and prediction values
colors <- c('white', RColorBrewer::brewer.pal(9, "Purples"))
sample_annot_colors <- c(RColorBrewer::brewer.pal(9, "BrBG")[[2]], RColorBrewer::brewer.pal(9, "BrBG")[[8]])
gene_annot_colors <- c(RColorBrewer::brewer.pal(8, "Set2")[[6]],  RColorBrewer::brewer.pal(8, "Set2")[[1]], rev(RColorBrewer::brewer.pal(8, "Set2")[2: 4]))
drug_colors <- c(RColorBrewer::brewer.pal(8, "Greens")[[5]], RColorBrewer::brewer.pal(8, "Greys")[[3]])

p_d2_raw <- ggplot() +
  
  scale_x_continuous(breaks=c(-1 : 13), labels=c("", names(column_ords)), guide = guide_axis(angle = 90 )) +
  scale_y_continuous(breaks=(0 : (n_plot )), labels = c(heatmap_dt[,GeneSymbol], ""), position = "right")  +
  
  geom_rect(data=ml, mapping=aes(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill = value), color = 'lightgrey') +
  scale_fill_gradientn(colours = colors) +
  labs(fill = 'Normalized\npredicted\nprobability') + 
  
  ggnewscale::new_scale_fill() +
  geom_rect(data=left_annot, mapping=aes(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill = Gene_Grouping) ) +
  scale_fill_manual(values=rev(gene_annot_colors)) + 
  labs(fill = 'Gene cluster') +
  
  ggnewscale::new_scale_fill() +
  geom_rect(data=top_annot, mapping=aes(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill = Cell_origin) ) +
  scale_fill_manual(values=sample_annot_colors) +
  labs(fill = 'Cell origin') +
  guides(fill = guide_legend(reverse = TRUE)) + 
  
  ggnewscale::new_scale_fill()+
  geom_rect(data=drug_annot, mapping=aes(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill = Approved_drug) ) +
  scale_fill_manual(values=drug_colors) +
  labs(fill = 'Approved drug')

p_d2 <- p_d2_raw + 
  # coord_fixed() +
  theme_vale +
  theme(
    axis.text.x = element_text(margin = margin(-40), size = 8),
    axis.text.y.right = element_text(face = "italic"),
    # axis.text.y.right = element_text(hjust = 0.5, vjust = 0, size = 8),
    # axis.text.y = element_text(margin = margin(-40), size = 8),
    axis.text.y = element_text(vjust = 0, size = 8), 
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.direction = "vertical",
    legend.box = "horizontal",
    legend.position = "bottom",
    legend.margin = margin(1, 1, 40, 0, "points"),
    legend.spacing.y = unit(0.3, 'cm'),
    plot.margin = margin(-1, 0, 0, 1, "cm"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA)
  ) 

# p_d2




#### d1. 1 big group, prediction histogram ####
curation_color_d1 <- c('Reported hematologic\nmalignancy driver gene'=RColorBrewer::brewer.pal(8, "Paired")[4], 
                       'Reported cancer \ndriver gene'=RColorBrewer::brewer.pal(8, "Paired")[3], 
                       'Candidate hematologic\nmalignancy driver gene'=RColorBrewer::brewer.pal(8, "Paired")[1])
x_label <- heatmap_dt[,GeneSymbol]
res_viz_pred[, GeneSymbol := factor(GeneSymbol, levels = x_label)]

res_viz_pred[, curation := 'Candidate hematologic\nmalignancy driver gene']
res_viz_pred[isCGC==TRUE, curation := 'Reported cancer \ndriver gene']
res_viz_pred[GeneSymbol %in% mll_cgc_leukemia_gene[, GeneSymbol], 
             curation := 'Reported hematologic\nmalignancy driver gene']
res_viz_pred[, curation := factor(curation, levels = names(curation_color_d1))]
#' top 100 curation
res_viz_pred[, table(curation)]

p_d1_raw <- ggplot(res_viz_pred, aes(x=GeneSymbol, y=Prediction_post, fill=curation)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=curation_color_d1) +
  coord_flip()

# p_d1_raw

p_d1 <- p_d1_raw + 
  scale_y_continuous(limits = c(0, 1), breaks=c(0, 0.5, 1)) +
  scale_x_discrete(position = "top") +
  theme_vale +
  xlab('') +
  ylab('Predicted probability') +
  labs(fill = 'CGC') +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = c(0.5, -0.145),
    legend.margin = margin(37, 1, 5.5, 5, "points"),
    legend.spacing.y = unit(0.3, 'cm'),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    plot.margin = margin(0.3, 1, 7.9, 0, "cm")
  )  +
  guides(fill=guide_legend('Manual curation',
                           title.position = "top",
                           direction = "vertical"))

# p_d1




#### c. pie #### 
res_100_curation <- res_viz_pred[, table(curation)]

curation_top100 <- data.table(
  Curation = names(curation_color),
  n_gene = as.numeric(res_100_curation)
)

enrichment100 <-
  matrix(c(res_100_curation[1], 
           sum(res_100_curation) - res_100_curation[1], 
           nrow(mll_cgc_leukemia_gene)-res_100_curation[1], 
           gencode[gene_type=='protein_coding', length(unique(gene_id))] - sum(res_100_curation) - (nrow(mll_cgc_leukemia_gene)-res_100_curation[1])),
         nrow = 2,
         dimnames = list(is100 = c("Y", "N"),
                         isDriver = c("Y", "N")))

fisher.test(enrichment100, alternative = "greater")$estimate
fisher.test(enrichment100, alternative = "greater")$p.value

curation_top100[, Curation := factor(Curation, levels = names(curation_color))]
curation_label <- curation_top100[, paste0(Curation, " (n=", n_gene, ')')]

p_c_raw <- ggplot(curation_top100, aes(x="", y=n_gene, fill=Curation)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=curation_color, labels=curation_label) 

p_c <- p_c_raw +
  theme_vale +
  ggtitle('Top 100 predicted genes') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),
        legend.position="bottom",
        legend.direction="vertical",
        legend.spacing.y = unit(0.3, 'cm')
  ) +
  guides(fill=guide_legend("Manual curation"))

# p_c




#### pharos enrichment ####
pharos_target_pr[, table(Target_Development_Level)]
pharos_target_pr[, length(unique(Symbol))]

pharos_target_pr[Target_Development_Level == 'Tclin', .N]

enrichment100 <-
  matrix(c(drug_target[, grep("Y", Target)] %>% length(), 
           pharos_target_pr[Target_Development_Level == 'Tclin', .N] - 18, 
           n_plot - drug_target[, grep("Y", Target)] %>% length(), 
           nrow(gencode_pr) - pharos_target_pr[Target_Development_Level == 'Tclin', .N] - 
             n_plot +drug_target[, grep("Y", Target)] %>% length()),
         nrow = 2,
         dimnames = list(is100 = c("Y", "N"),
                         isDrug = c("Y", "N")))
enrichment100 
sum(enrichment100)

fisher.test(enrichment100, alternative = "greater")




#### check for false positive ####
# https://bitbucket.org/intogen/intogen-plus/src/master/docs/source/postprocessing.rst

library("rjson")
artifacts <- fromJSON(file = '/s/project/mll/intogen/working_directory/local/intogen-plus/datasets/postprocess/artifacts.json')
artifacts[[2]]

# one known artifiact genes
res_viz_pred[, GeneSymbol %in% artifacts[[2]]] %>% table() 
res_viz_pred[GeneSymbol %in% artifacts[[2]], GeneSymbol]
# [1] TTN
# three suspected artifiact genes
res_viz_pred[, GeneSymbol %in% artifacts[[1]]] %>% table()
res_viz_pred[GeneSymbol %in% artifacts[[1]], GeneSymbol]
# [1] IGLL5 TTN   SYNE1

blacklist <- fread('/s/project/mll/intogen/working_directory/local/intogen-plus/datasets/postprocess/black_listed.txt', header = F)
# no black list
res_viz_pred[, GeneSymbol %in% blacklist[, V1]] %>% table()


olfactory <- fread('/s/project/mll/intogen/working_directory/local/intogen-plus/datasets/others/olfactory_receptors.tsv')
# no olfactory
res_viz_pred[, GeneSymbol %in% olfactory[, Symbol]] %>% table()




### arrange plot - heatmap #### 
#' ### raw
#+ plot p4 raw, fig.width=10, fig.height=14
p_abc_raw <- ggarrange(p_a_raw, p_b_raw, p_c_raw, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d_raw <- ggarrange(p_d1_raw, p_d2_raw, nrow = 1, labels = c(' ', ' '), widths = c(1, 0.4))

p4_raw <- ggarrange(p_abc_raw, p_d_raw, nrow = 1, widths = c(1,1.5), labels = c('', 'D'))
p4_raw

#' ### annotated
#+ plot p4, fig.width=10, fig.height=14
p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c(' ', ' '), widths = c(1, 0.4))

p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5), labels = c('', 'D'))
p4


saveRDS(p_c, snakemake@output$p_c)
saveRDS(p_d, snakemake@output$p_d)

# png(paste0(output_dir, "/figure_4_poster.png"), width = 10, height = 16, units = "in", res = 600)
# p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
# p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c(' ', ' '), widths = c(1, 0.4))
# 
# p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5), labels = c('', 'D'))
# p4
# dev.off()


pdf(paste0(output_dir, "/figure_4.pdf"), width = 10, height = 14)
p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c(' ', ' '), widths = c(1, 0.4))

p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5), labels = c('', 'D'))
p4
dev.off()


png(paste0(output_dir, "/figure_4.png"), width = 10, height = 14, units = "in", res = 600)
p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c(' ', ' '), widths = c(1, 0.4))

p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5), labels = c('', 'D'))
p4
dev.off()



### Supplement #####
#' ## Sup
#' ### s1 raw
#+ plot s1 raw, fig.width=6, fig.height=6
p_s1_raw
#' ### s1 annotated
#+ plot s1, fig.width=6, fig.height=6
p_s1

#' ### s2 raw
#+ plot s2 raw, fig.width=6, fig.height=6
p_s2_raw
#' ### s2 annotated
#+ plot s2, fig.width=6, fig.height=6
p_s2

#' ### s3 raw
#+ plot s3 raw, fig.width=8, fig.height=4
p_s3_raw
#' ### s3 annotated
#+ plot s3, fig.width=8, fig.height=4
p_s3

#' ### s4 raw
#+ plot s4 raw, fig.width=10.5, fig.height=6
p_s4_raw
#' ### s4 annotated
#+ plot s4, fig.width=10.5, fig.height=6
p_s4

#' ### s5 raw
#+ plot s5 raw, fig.width=12, fig.height=6
p_s5_raw
#' ### s5 annotated
#+ plot s5, fig.width=12, fig.height=6
p_s5

#' ### s8 raw
#+ plot s8 raw, fig.width=19, fig.height=6
p_s8_raw
#' ### s8 annotated
#+ plot s8, fig.width=9, fig.height=6
p_s8

#' ### s9 raw
#+ plot s9 raw, fig.width=13.5, fig.height=4
p_s9_raw
#' ### s9 annotated
#+ plot s9, fig.width=13.5, fig.height=4
p_s9




### thesis ####
png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "ap_ind.png"), 
    width = 5, height = 6, units = "in", res = 600)
p_a
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "prc_ind.png"), 
    width = 6, height = 6, units = "in", res = 600)
p_s1
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "ap_accu.png"), 
    width = 5, height = 6, units = "in", res = 600)
p_b
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "prc_accu.png"), 
    width = 6, height = 6, units = "in", res = 600)
p_s2
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "model_all_pie.png"), 
    width = 4, height = 4, units = "in", res = 600)
p_c
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "model_benchmark.png"), 
    width = 8, height = 4, units = "in", res = 600)
p_s3
dev.off() 

p_d <- ggarrange(p_d2, p_d1, nrow = 1, widths = c(1, 0.4))

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "model_heatmap.png"), 
    width = 6, height = 14, units = "in", res = 600)
p_d
dev.off() 
