dir_path <- "/s/project/vale/driver_prediction_202402/manuscript/"

# plot table
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

source("function.R")
source("manuscript_theme.R")

mll_prc <- fread(paste0(dir_path, 'aggregated_data_figure/mll_prc.tsv'))
mll_ap <- fread(paste0(dir_path, 'aggregated_data_figure/mll_ap.tsv'))
mll_ap_full <- fread(paste0(dir_path, 'aggregated_data_figure/mll_ap_full.tsv'))

res_viz_pred <- fread(paste0(dir_path, "supplementary_table/S13_prediction_complete_dataset.csv"))
res_viz_groups_pred <- fread(paste0(dir_path, "supplementary_table/S15_prediction_study_groups.csv"))
cgc_leukemia_gene <- fread(paste0(dir_path, 'aggregated_data_figure/CGC_leukemia_gene.csv'))
mll_cgc_leukemia_gene <- fread(paste0(dir_path, 'aggregated_data_figure/MLL_CGC_leukemia_gene.csv'))
drug_target <- fread(paste0(dir_path, "supplementary_table/S14_drug.csv"))
gencode <- fread(paste0(dir_path, "aggregated_data_figure/gencode.csv"))
rcl.list <- readRDS(paste0(dir_path, "aggregated_data_figure/rcl_list.Rds"))

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

sample_group_by_cell_type <- c('AML', 'MDS', 'MPN', 'MDS/MPN group', 'CMML', 
                               'CML', 'NK', 'PCN group', 'Mastocytosis', 
                               'T-cell group', 'CLL', 'MatureB group', 'Hairy-cell group', 'BCP-ALL')



#### a. ap individual, mll #### 
#' number of RNAseq
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

p_a




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

p_b





#### d2. 14 groups, prediction heatmap ####
n_plot <- 100
# read in res_viz
setnames(res_viz_pred, 
         c('StudyGroup', "Prediction"), 
         c('Sample_group', "Prediction_post"))
res_viz_pred[, Prediction_post_scaled := Prediction_post/max(Prediction_post)]
res_viz_pred[, isLeukemia := GeneID %in% cgc_leukemia_gene[, ENSGid]]
res_viz_pred <- res_viz_pred[, .(GeneSymbol, Prediction_post, Prediction_post_scaled, Sample_group, isCGC, isLeukemia)]
selected_genes <- res_viz_pred[1:n_plot, GeneSymbol]
res_viz_pred <- res_viz_pred[GeneSymbol %in% selected_genes, ]

# read in res_viz_groups
setnames(res_viz_groups_pred, 
         c('StudyGroup', "Prediction"), 
         c('Sample_group', "Prediction_post"))
res_viz_groups_pred[, Prediction_post_scaled := Prediction_post/max(Prediction_post), by='Sample_group']
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

# sort heatmap based on dendogram order
heatmap_dt <- heatmap_dt[rev(unlist(rcl.list )),]
gene_cluster_ranges <- rev(as.integer(lapply(rcl.list, function(x) length(x))))

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

drug_annot <- merge(ml[!duplicated(ml$GeneSymbol), ], drug_target, by = "GeneSymbol")
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

p_d1_raw

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

p_d1




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

p_c




### arrange plot #### 
p_abc <- ggarrange(p_a, p_b, p_c, ncol = 1, labels = c('A', 'B', 'C'), heights = c(1, 1, 0.8))
p_d <- ggarrange(p_d2, p_d1, nrow = 1, labels = c(' ', ' '), widths = c(1, 0.4))

p4 <- ggarrange(p_abc, p_d, nrow = 1, widths = c(1,1.5), labels = c('', 'D'))
p4

