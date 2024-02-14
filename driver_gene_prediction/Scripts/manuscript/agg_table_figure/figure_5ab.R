dir_path <- "/s/project/vale/driver_prediction_202402/manuscript/"

# plot table
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(pheatmap)
  library(RColorBrewer)
  library(ggpubr)
})

source("function.R")
source("manuscript_theme.R")

CGC_cancer_gene <- fread(paste0(dir_path, 'aggregated_data_figure/CGC_cancer_gene_processed.csv'))



#### a. heatmap fill 0 #####
figure_5a_dt <- fread(paste0(dir_path, 'supplementary_table/S16_association.csv'))

heatmap_dt <- figure_5a_dt[Method == 'Activation' &
                             FDR<0.05 &
                             GeneSymbol %in% CGC_cancer_gene[, GeneSymbol], ]
heatmap_dt <- heatmap_dt[, .(DiseaseEntity, GeneSymbol, OddsRatio)] 
heatmap_dt <- dcast(heatmap_dt, DiseaseEntity ~ GeneSymbol, value.var = 'OddsRatio', fill=0)

heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[, DiseaseEntity]
heatmap_mtx[heatmap_mtx>1000] <- 1000

# customize order
sp <- slanter::sheatmap(heatmap_mtx, oclust_rows=FALSE, oclust_cols=FALSE)
gene_order <- sp$gtable$grobs[[4]]$label
DiseaseEntity_order <- sp$gtable$grobs[[5]]$label
heatmap_mtx <- heatmap_mtx[DiseaseEntity_order, gene_order]

# mimic oncoplot order
DiseaseEntity_sum <- apply(heatmap_mtx, 1, function(x){sum(x!=0)}) %>% sort(decreasing = TRUE)
DiseaseEntity_order <- names(DiseaseEntity_sum[c(1:6, 8, 12, 11, 7, 9, 10)])

gene_value <- heatmap_mtx[DiseaseEntity_order[1], ] %>%  sort(decreasing = TRUE)
# the for loop was manually tested
for (i in DiseaseEntity_order[2:8]) {
  print(i)
  temp_value <- heatmap_mtx[i, names(gene_value[gene_value==0])] %>%  sort(decreasing = TRUE)
  gene_value <- c(gene_value[gene_value!=0], temp_value)
}
gene_order <- names(gene_value)

# define color
curation_colors <- c('known' = RColorBrewer::brewer.pal(9, "YlGn")[7],
                     'little' = RColorBrewer::brewer.pal(9, "YlGn")[5],
                     'unknown' = RColorBrewer::brewer.pal(9, "Blues")[5],
                     'NA' = RColorBrewer::brewer.pal(9, "Greys")[3])
curation_labels <- c('known' = "Multiple literature ",
                     'little' = "One literature",
                     'unknown' = "None",
                     'NA' = "NA")

# prepare plot dt
heatmap_dt <- figure_5a_dt[Method == 'Activation' &
                             FDR<0.05 &
                             GeneSymbol %in% CGC_cancer_gene[, GeneSymbol], ]
heatmap_dt <- heatmap_dt[, .(DiseaseEntity, GeneSymbol, OddsRatio)]
heatmap_dt[OddsRatio>1000, OddsRatio := 1000]
heatmap_dt <- dcast(heatmap_dt, DiseaseEntity ~ GeneSymbol, value.var = 'OddsRatio', fill=1)
heatmap_dt <- melt(heatmap_dt, id.vars = 'DiseaseEntity', variable.name = 'GeneSymbol', value.name = 'OddsRatio')
heatmap_dt[, DiseaseEntity := factor(DiseaseEntity, levels = rev(DiseaseEntity_order))]
heatmap_dt[, GeneSymbol := factor(GeneSymbol, levels = gene_order)]

# define breaks and colors
breaks <- c(10^(0:ceiling(log10(max(heatmap_mtx, na.rm = TRUE)))))
colors <- c('white', brewer.pal(length(breaks)-1, "Purples"))
legend_labels <- c(1, 10, 100, '\u22651000')

p_a_raw <- ggplot(heatmap_dt, aes(GeneSymbol, DiseaseEntity, fill = OddsRatio)) +
  geom_tile(linewidth=0.5, color='lightgrey')

p_a <- p_a_raw +
  xlab('CGC cancer gene') +
  ylab('Disease entity') +
  theme_vale +
  scale_fill_gradientn(colors = c('white', brewer.pal(3, "Purples")),
                       trans = "log", name = "Odds ratio",
                       breaks = breaks, labels = legend_labels) +
  coord_fixed() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1, face = "italic"),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.spacing.y = unit(0.5, 'cm')
  )

p_a




#### b. stacked bar plot #####
figure_5b_dt <- fread(paste0(dir_path, 'Catalog_of_transcriptomic_and_genomic_aberrations_of_24_hematologic_malignancy_entities.csv'))
sample_summary_dt <- fread(paste0(dir_path, 'supplementary_table/S1_sample_summary.csv'))

figure_5b_dt <- figure_5b_dt[Method=='NB-act' & GeneSymbol=="LRP1B" &
                               Metric == "Number of samples predicted when filtering for padjust<0.05 and zScore>0", ]
sample_summary_dt[DiseaseEntity %in% figure_5b_dt[, Entity], Number_of_individual]
figure_5b_dt <- merge(figure_5b_dt, sample_summary_dt[, .(DiseaseEntity, Number_of_individual)],
                      by.x='Entity', by.y='DiseaseEntity')

percent_dt <- figure_5b_dt[, .(Entity, Value, Number_of_individual)]
setnames(percent_dt, c("Entity", "Value", "Number_of_individual"), c("Cohort", "N_Cohort_LRP1B", "N_Cohort"))

percent_dt[, ratio := N_Cohort_LRP1B/N_Cohort]
percent_dt[, Cohort := factor(Cohort, levels = c('HCL-V', 'HCL', 'MZL', 'MM', 'Others'))]

percent_dt[, text_label := paste0( label_percent(accuracy=0.1)(ratio), "\n",
                                   "(", N_Cohort_LRP1B, "/", N_Cohort, ")")]
percent_dt[, color_label := factor(Cohort, levels = c('Non-activated', 'HCL-V', 'HCL', 'MZL', 'MM', 'Others'))]

percent_dt_pval <- sapply(percent_dt[, Cohort], function(x){
  contingency_mtx <-
    matrix(c(percent_dt[Cohort==x, N_Cohort_LRP1B], 
             percent_dt[Cohort==x, N_Cohort - N_Cohort_LRP1B], 
             percent_dt[, sum(N_Cohort_LRP1B)] - percent_dt[Cohort==x, N_Cohort_LRP1B], 
             sample_summary_dt[, sum(Number_of_individual)] - 
               percent_dt[Cohort==x, N_Cohort - N_Cohort_LRP1B] - 
               percent_dt[, sum(N_Cohort_LRP1B)]), 
           nrow = 2)
  fisher.test(contingency_mtx, alternative = 'greater')$p.value
})

percent_dt[, pvalue := percent_dt_pval]
percent_dt <- add_significance(percent_dt, "pvalue")

p_b_raw <- ggplot(percent_dt, aes(x=Cohort, y=ratio)) + 
  geom_bar(stat="identity", fill='lightgrey') +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  annotate("text", x=percent_dt[, Cohort], y=1.15, label=percent_dt[, text_label], size = 3) +
  annotate("text", x=percent_dt[, Cohort], y=0.95, label=percent_dt[, pvalue.signif], size = 3) +
  annotate("text", x=0, y=1.15, label='Percentage\n activated', size = 3) +
  coord_cartesian(ylim=c(0,1), xlim=c(1, 4), clip="off")

p_b_raw

p_b <- p_b_raw +
  # scale_fill_manual(values=cohort_color, name='') +
  ylab(expression(paste(bold("Ratio of "), bolditalic("LRP1B"), bold("-activated samples"))))+
  xlab('Disease entity (dataset)') +
  theme_vale +
  theme(
    plot.margin = margin(45, 14, 14, 30, "points"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

p_b




### arrange plot #### 
p5_bottom <- ggarrange(p_b, ggplot() + theme(panel.background = element_blank()), 
                       ncol = 2, nrow = 1, labels = c("B", ""))
p5 <- ggarrange(p_a, p5_bottom,
                ncol = 1, nrow = 2, labels=c("A", ""), heights = c(3.5, 4))

p5

