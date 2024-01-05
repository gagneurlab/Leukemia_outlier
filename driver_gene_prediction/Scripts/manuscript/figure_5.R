#'---
#' title: figure_5
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
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'    - MLL_CGC_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/MLL_CGC_leukemia_gene_list.tsv"`'
#'  input:
#'    - mll_prc: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_5/plot_data/mll_prc.tsv"`'
#'    - mll_ap: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_5/plot_data/mll_ap.tsv"`'
#'    - mll_ap_full: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_5/plot_data/mll_ap_full.tsv"`'
#'    - mll_benchmark_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_5/plot_data/mll_benchmark_prc.tsv"`'
#'    - mll_benchmark_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_5/plot_data/mll_benchmark_ap.tsv"`'
#'  output:
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_5.html"`'
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
                             "/processed_data/snakemake/figure_5.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_5.snakemake")
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
                                              "or,ac",
                                              "fr",
                                              
                                              "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-absplice",
                                              "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice",
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
               "  Expression outliers \n      (OUTRIDER+NB-act)", 
               "  Splicing outliers \n      (FRASER)", 
               
               "+ Splicing variants \n      (AbSplice)", 
               "+ Expression outliers \n      (OUTRIDER+NB-act)", 
               "+ Splicing outliers \n      (FRASER)"))

feat_label_dt[, feature_source := 'Genomic']
feat_label_dt[input_feature %in% c(
  "or,ac",
  "fr",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice",
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


experiment_design <- fread(snakemake@params$experimentDesign)
experiment_design[is.na(experiment_design)] <- ''

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))

mll_cgc_leukemia_gene <- fread(snakemake@params$MLL_CGC_leukemia_gene_list)

mll_prc <- fread(snakemake@input$mll_prc)
mll_ap <- fread(snakemake@input$mll_ap)
mll_ap_full <- fread(snakemake@input$mll_ap_full)
mll_benchmark_prc <- fread(snakemake@input$mll_benchmark_prc)
mll_benchmark_ap <- fread(snakemake@input$mll_benchmark_ap)




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
  scale_y_continuous(limits = c(0, 0.2), breaks=c(0, 0.1, 0.2)) +
  xlab("Individual feature") + 
  ylab("Average Precision") +
  theme_vale +
  theme(
    legend.position = 'bottom',
    # axis.text.x = element_text(angle = 45, hjust=1),
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
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice",
  "clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv-or,ac,absplice,fr")
res_ap <- mll_ap[input_feature %in% on_top,]

p_b_raw <- ggplot(res_ap, aes(input_feature, AP_median)) + 
  geom_errorbar(aes(ymin=AP_1q, ymax=AP_9q), color='gray20', width = 0.2, linewidth = 0.3) +
  geom_col(aes(fill=feature_source), alpha = 0.3) +
  stat_pvalue_manual(stat.test[c(1,5),], label="label", 
                     y.position = 0.185, bracket.shorten = -0.1, color='gray20') +
  annotate("text", x=3.5, y=0.2, label='***', size = 4, angle=90, color='gray20') +
  annotate("text", x=2, y=0.2, label='**', size = 4, angle=90, color='gray20') +
  coord_flip()

# p_b_raw

p_b <- p_b_raw +  
  scale_fill_manual(values = feat_color) +
  scale_x_discrete(labels = feat_label) +
  scale_y_continuous(limits = c(0, 0.2), breaks=c(0, 0.1, 0.2)) +
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
  xlim(0, 0.1) +
  ylim(0, 0.1) +
  theme_vale

# p_s3




#### c. pie #### 
curation_top100 <- data.table(
  Curation = names(curation_color),
  n_gene = c(63, 4, 33)
)

enrichment100 <-
  matrix(c(63, 100-63, 
           377-63, gencode[gene_type=='protein_coding', length(unique(gene_id))] - 100 - (377-63)),
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




#### d2. 14 groups, prediction heatmap ####
n_plot <- 100
# read in res_viz
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == single_group &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == '', ]
exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir)

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
exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir)

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

# get row order
# heatmap_mtx <- as.matrix(heatmap_dt[,-1])
# rownames(heatmap_mtx) <- heatmap_dt[, GeneSymbol]
# 

# 
# set.seed(512)
# p_d2_phm <- pheatmap(heatmap_mtx,
#                      row_km = num_gene_clusters,
#                      column_km = num_sample_clusters
#                      )
# p_d2_phm
# 
# p_d2_phm@row_order
# 
# p_d2_phm <- draw(p_d2_phm)

heatmap_path <- "/s/project/vale/Resource/fig5_heatmap.rds"
# save pheatmap object if needed
# saveRDS(p_d2_phm, file=heatmap_path)

p_d2_phm <- readRDS(heatmap_path)


r.dend <- row_dend(p_d2_phm)
rcl.list <- row_order(p_d2_phm)


# trying new row orders to make the heatmap more beautiful
rcl.list <- c(rcl.list[1:2], rev(rcl.list[3:5]))


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
  x_min = rep(13.61, num_gene_clusters),
  x_max = rep(14.48, num_gene_clusters),
  y_min = c(0, y_maxes[-length(y_maxes)]),
  y_max = y_maxes,
  Gene_Grouping = factor(paste0("Cluster ", rev(1:num_gene_clusters)))
  
)

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


p_d2_raw <- ggplot() +
  
  scale_x_continuous(breaks=c(-1 : 13), labels=c("", names(column_ords)), guide = guide_axis(angle = 90 )) +
  scale_y_continuous(breaks=(0 : (n_plot )), labels = c(heatmap_dt[,GeneSymbol], ""), position = "right")  +
  
  geom_rect(data=ml, mapping=aes(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill = value), color = 'lightgrey') +
  scale_fill_gradientn(colours = colors) +
  labs(fill = 'Normalized\npredicted\nprobability') + 
  
  ggnewscale::new_scale_fill()+
  geom_rect(data=left_annot, mapping=aes(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill = Gene_Grouping) ) +
  scale_fill_manual(values=rev(gene_annot_colors)) + 
  labs(fill = 'Gene cluster') + 
  
  ggnewscale::new_scale_fill()+
  geom_rect(data=top_annot, mapping=aes(xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill = Cell_origin) ) +
  scale_fill_manual(values=sample_annot_colors) + 
  labs(fill = 'Cell origin') + 
  guides(fill = guide_legend(reverse = TRUE)) 

# p_d2_raw

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
#+ plot p6 raw, fig.width=10, fig.height=14
p_abc_raw <- ggarrange(p_a_raw, p_b_raw, p_c_raw, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d_raw <- ggarrange(p_d1_raw, p_d2_raw, nrow = 1, labels = c('D', ' '), widths = c(1, 0.4))

p5_raw <- ggarrange(p_abc_raw, p_d_raw, nrow = 1, widths = c(1,1.5))
p5_raw

#' ### annotated
#+ plot p6, fig.width=10, fig.height=14
p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c('D', ' '), widths = c(1, 0.4))

p5 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5))
p5


# png(paste0(output_dir, "/figure_5_poster.png"), width = 10, height = 16, units = "in", res = 600)
# p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
# p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c('D', ' '), widths = c(1, 0.4))
# 
# p5 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5))
# p5
# dev.off()


pdf(paste0(output_dir, "/figure_5.pdf"), width = 10, height = 14)
p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c('D', ' '), widths = c(1, 0.4))

p5 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5))
p5
dev.off()


png(paste0(output_dir, "/figure_5.png"), width = 10, height = 14, units = "in", res = 600)
p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c('D', ' '), widths = c(1, 0.4))

p5 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5))
p5
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
