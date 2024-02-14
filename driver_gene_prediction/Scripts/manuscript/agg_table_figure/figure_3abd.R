dir_path <- "/s/project/vale/driver_prediction_202402/manuscript/"

# plot table
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(rstatix)
})

source("function.R")
source("manuscript_theme.R")

leu_ocg <- fread(paste0(dir_path, 'aggregated_data_figure/leu_ocg.csv'))
leu_tsg <- fread(paste0(dir_path, 'aggregated_data_figure/leu_tsg.csv'))

  

  
#### a. FRASER, leu13, leukemia tsg #### 
figure_3a_dt <- fread(paste0(dir_path, 'aggregated_data_figure/figure_3a.csv'))
exp_gene_fr <- figure_3a_dt[, exp_gene_fr]
fr_gene_0 <- figure_3a_dt[fr_gene_0==TRUE, exp_gene_fr]
fr_gene_1 <- figure_3a_dt[fr_gene_1==TRUE, exp_gene_fr]
fr_gene_2 <- figure_3a_dt[fr_gene_2==TRUE, exp_gene_fr]
fr_gene_5 <- figure_3a_dt[fr_gene_5==TRUE, exp_gene_fr]
fr_top3_gene_0 <- figure_3a_dt[fr_top3_gene_0==TRUE, exp_gene_fr]
fr_top3_gene_1 <- figure_3a_dt[fr_top3_gene_1==TRUE, exp_gene_fr]
fr_top3_gene_2 <- figure_3a_dt[fr_top3_gene_2==TRUE, exp_gene_fr]
fr_top3_gene_5 <- figure_3a_dt[fr_top3_gene_5==TRUE, exp_gene_fr]

#' enrichment overall
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

p_a_raw

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

p_a




#### b. abSplice, leu13, leukemia tsg ####
figure_3b_dt <- fread(paste0(dir_path, 'aggregated_data_figure/figure_3b.csv'))
exp_gene <- figure_3b_dt[, exp_gene]
absplice_gene_0 <- figure_3b_dt[absplice_gene_0==TRUE, exp_gene]
absplice_gene_1 <- figure_3b_dt[absplice_gene_1==TRUE, exp_gene]
absplice_gene_2 <- figure_3b_dt[absplice_gene_2==TRUE, exp_gene]
absplice_gene_5 <- figure_3b_dt[absplice_gene_5==TRUE, exp_gene]
absplice_top3_gene_0 <- figure_3b_dt[absplice_top3_gene_0==TRUE, exp_gene]
absplice_top3_gene_1 <- figure_3b_dt[absplice_top3_gene_1==TRUE, exp_gene]
absplice_top3_gene_2 <- figure_3b_dt[absplice_top3_gene_2==TRUE, exp_gene]
absplice_top3_gene_5 <- figure_3b_dt[absplice_top3_gene_5==TRUE, exp_gene]

# all
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

p_b_raw

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

p_b




#### c. FRASER, abSplice, leu13, leukemia tsg ####
figure_3d_dt <- fread(paste0(dir_path, 'aggregated_data_figure/figure_3d.csv'))
exp_gene <- figure_3d_dt[, exp_gene]
absplice_fr_gene_0 <- figure_3d_dt[absplice_fr_gene_0==TRUE, exp_gene]
absplice_fr_gene_1 <- figure_3d_dt[absplice_fr_gene_1==TRUE, exp_gene]
absplice_fr_gene_2 <- figure_3d_dt[absplice_fr_gene_2==TRUE, exp_gene]
absplice_fr_gene_5 <- figure_3d_dt[absplice_fr_gene_5==TRUE, exp_gene]
absplice_fr_top3_gene_0 <- figure_3d_dt[absplice_fr_top3_gene_0==TRUE, exp_gene]
absplice_fr_top3_gene_1 <- figure_3d_dt[absplice_fr_top3_gene_1==TRUE, exp_gene]
absplice_fr_top3_gene_2 <- figure_3d_dt[absplice_fr_top3_gene_2==TRUE, exp_gene]
absplice_fr_top3_gene_5 <- figure_3d_dt[absplice_fr_top3_gene_5==TRUE, exp_gene]

# all
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

p_d_raw <- ggplot(absplice_fr_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
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

p_d_raw

p_d <- p_d_raw +
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

p_d




### arrange plot #### 
p3_abd <- grid.arrange(p_a, p_b, 
                       ggplot() + theme(panel.background = element_blank()), p_d, nrow = 2)

p3_abd
