#'---
#' title: figure_4_coess
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
#'    - MLL_CGC_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/MLL_CGC_leukemia_gene_list.tsv"`'
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
#'    - mll_benchmark_prc: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_benchmark_prc.tsv"`'
#'    - mll_benchmark_ap: '`sm config["projectPath"] + 
#'                      "/manuscript/figure_4/plot_data/mll_benchmark_ap.tsv"`'
#'  output:
#'    - p_c_coess: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_c_coess.Rds"`' 
#'    - p_d_coess: '`sm config["projectPath"] + "/manuscript/figure_4/plot_data/p_d_coess.Rds"`' 
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_4_coess.html"`'
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
                             "/processed_data/snakemake/figure_4_coess.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/figure_4_coess.snakemake")
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

mll_var_filter_prc <- fread(snakemake@input$mll_var_filter_prc)
mll_var_filter_ap <- fread(snakemake@input$mll_var_filter_ap)
mll_var_filter_ap_full <- fread(snakemake@input$mll_var_filter_ap_full)

mll_emb_prc <- fread(snakemake@input$mll_emb_prc)
mll_emb_ap <- fread(snakemake@input$mll_emb_ap)
mll_emb_ap_full <- fread(snakemake@input$mll_emb_ap_full)




### Figure #### 
### prediction only emb #### 
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == single_group &
                               model_method == 'rf' &
                               intogen_input_feature == '' &
                               outlier_input_feature == '' &
                               coess_input_feature == 'coess_cluster', ]
exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)

res_viz_pred <- fread(exp_viz[, res_post_path])
res_viz_pred[, Sample_group := single_group]
res_viz_pred[, Prediction_post_scaled := Prediction_post/max(Prediction_post)]
res_viz_pred <- res_viz_pred[, .(GeneSymbol, Prediction_post, Prediction_post_scaled, Sample_group, isCGC, isLeukemia)]




### prediction post #### 
#### d2. 14 groups, prediction heatmap ####
n_plot <- 100
# read in res_viz
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == single_group &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == 'coess_cluster', ]
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
                               coess_input_feature == 'coess_cluster', ]
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

# get row order
heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[, GeneSymbol]


# set.seed(512)
p_d2_phm <- pheatmap(heatmap_mtx,
                     row_km = num_gene_clusters,
                     column_km = num_sample_clusters
                     )
p_d2_phm

p_d2_phm@row_order

p_d2_phm <- draw(p_d2_phm)

heatmap_path <- "/s/project/vale/Resource/figure_4_coess_heatmap.rds"
# save pheatmap object if needed
saveRDS(p_d2_phm, file=heatmap_path)

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

p_d2




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




### prediction post same gene as no emb #### 
### arrange plot - heatmap #### 
pdf(paste0(output_dir, "/figure_4_coess_top100.pdf"), width = 10, height = 14)
p_abc <- ggarrange(ggplot() + theme(panel.background = element_blank()), 
                   ggplot() + theme(panel.background = element_blank()), 
                   p_c, 
                   ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c('D', ' '), widths = c(1, 0.4))

p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5))
p4
dev.off()


png(paste0(output_dir, "/figure_4_coess_top100.png"), width = 10, height = 14, units = "in", res = 600)
p_abc <- ggarrange(ggplot() + theme(panel.background = element_blank()), 
                   ggplot() + theme(panel.background = element_blank()), 
                   p_c, 
                   ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c('D', ' '), widths = c(1, 0.4))

p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5))
p4
dev.off()




#### d2. 14 groups, prediction heatmap ####
n_plot <- 100
# read in res_viz
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == single_group &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == 'coess_cluster', ]
exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)

res_viz_pred <- fread(exp_viz[, res_post_path])
res_viz_pred[, Sample_group := single_group]
res_viz_pred[, Prediction_post_scaled := Prediction_post/max(Prediction_post)]
res_viz_pred <- res_viz_pred[, .(GeneSymbol, Prediction_post, Prediction_post_scaled, Sample_group, isCGC, isLeukemia)]

heatmap_path <- "/s/project/vale/Resource/figure_4_heatmap.rds"
p_d2_phm <- readRDS(heatmap_path)
selected_genes <- p_d2_phm@ht_list$matrix_3@row_names_param$labels
res_viz_pred <- res_viz_pred[GeneSymbol %in% selected_genes, ]

# read in res_viz_groups
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group %in% sample_group_by_size &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == 'coess_cluster', ]
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

# get row order
heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[, GeneSymbol]

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




#### check overlap ####
heatmap_path <- "/s/project/vale/Resource/figure_4_heatmap.rds"
phm_omics <- readRDS(heatmap_path)
predicted_100genes_omics <- phm_omics@ht_list$matrix_3@row_names_param$labels

heatmap_path <- "/s/project/vale/Resource/figure_4_coess_heatmap.rds"
phm_omics_coess <- readRDS(heatmap_path)
predicted_100genes_omics_coess <- phm_omics_coess@ht_list$matrix_1@row_names_param$labels

# overlap 86 genes
intersect(predicted_100genes_omics, predicted_100genes_omics_coess) %>% length()

# predicted by omics only
setdiff(predicted_100genes_omics, predicted_100genes_omics_coess) 
# [1] "ANKRD36C" "CHD2"     "DIS3"     "EIF3E"    "FAM221A"  "MYC"      "NBPF11"   "PCNT"     "SEL1L3"   "SETBP1"   "SMG1"    
# [12] "SP140"    "VWF"      "ZNF605"  

# predicted after adding emb
setdiff(predicted_100genes_omics_coess, predicted_100genes_omics) 
# [1] "COPS9"   "ETV6"    "IKZF3"   "KIR2DL3" "KIR3DL2" "KRAS"    "NFE2"    "POTED"   "PRDM2"   "RB1"     "SH3KBP1" "TCF12"  
# [13] "WDR31"   "WT1" 




### arrange plot - heatmap #### 
pdf(paste0(output_dir, "/figure_4_coess_same_gene.pdf"), width = 10, height = 14)
p_abc <- ggarrange(ggplot() + theme(panel.background = element_blank()), 
                   ggplot() + theme(panel.background = element_blank()), 
                   ggplot() + theme(panel.background = element_blank()), 
                   ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c(' ', ' '), widths = c(1, 0.4))

p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5))
p4
dev.off()


png(paste0(output_dir, "/figure_4_coess_same_gene.png"), width = 10, height = 14, units = "in", res = 600)
p_abc <- ggarrange(ggplot() + theme(panel.background = element_blank()), 
                   ggplot() + theme(panel.background = element_blank()), 
                   ggplot() + theme(panel.background = element_blank()), 
                   ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c(' ', ' '), widths = c(1, 0.4))

p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5))
p4
dev.off()


saveRDS(p_c, snakemake@output$p_c_coess)
saveRDS(p_d, snakemake@output$p_d_coess)