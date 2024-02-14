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

  

  
#### a. OUTRIDER, leu13, leukemia tsg #### 
figure_2a_dt <- fread(paste0(dir_path, 'aggregated_data_figure/figure_2a.csv'))
exp_gene <- figure_2a_dt[, exp_gene]
or_dn_gene_0 <- figure_2a_dt[or_dn_gene_0==TRUE, exp_gene]
or_dn_gene_1 <- figure_2a_dt[or_dn_gene_1==TRUE, exp_gene]
or_dn_gene_2 <- figure_2a_dt[or_dn_gene_2==TRUE, exp_gene]
or_dn_gene_5 <- figure_2a_dt[or_dn_gene_5==TRUE, exp_gene]
or_dn_top3_gene_0 <- figure_2a_dt[or_dn_top3_gene_0==TRUE, exp_gene]
or_dn_top3_gene_1 <- figure_2a_dt[or_dn_top3_gene_1==TRUE, exp_gene]
or_dn_top3_gene_2 <- figure_2a_dt[or_dn_top3_gene_2==TRUE, exp_gene]
or_dn_top3_gene_5 <- figure_2a_dt[or_dn_top3_gene_5==TRUE, exp_gene]

# all
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
  odds_ratio_ci_low = c(
    fisher_test(or_dn_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(or_dn_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(or_dn_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(or_dn_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(or_dn_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(or_dn_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(or_dn_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(or_dn_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[2]
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
  odds_ratio_ci_low = c(
    fisher_test(or_dn_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(or_dn_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(or_dn_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[1],
    fisher_test(or_dn_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(or_dn_top3_gene_0, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(or_dn_top3_gene_1, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(or_dn_top3_gene_2, exp_gene, leu_tsg[, ENSGid])$conf.int[2],
    fisher_test(or_dn_top3_gene_5, exp_gene, leu_tsg[, ENSGid])$conf.int[2]
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
or_dn_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=list(rownames(or_dn_en))]
or_dn_en[total < 1000, total_label := as.character(total)]
or_dn_en[category=='0', total_pos := 1]
or_dn_en[category=='1', total_pos := 2]
or_dn_en[category=='2-4', total_pos := 3]
or_dn_en[category=="\u22655", total_pos := 4]
or_dn_en[cutoff=='All', total_pos := total_pos-0.15]
or_dn_en[cutoff=='At most three', total_pos := total_pos+0.15]

or_dn_en[p_val.signif=='ns', odds_ratio_ci_low := NA]
or_dn_en[p_val.signif=='ns', odds_ratio_ci_up := NA]
or_dn_en[, max(odds_ratio)]
or_dn_en[, max(odds_ratio_ci_up, na.rm=TRUE)]

signif_height_a <- 10
total_height_a <- 17

p_a_raw <- ggplot(or_dn_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=odds_ratio_ci_low, ymax=odds_ratio_ci_up),
                width=0.1, position=position_dodge(width=0.7)) +
  geom_text(data = or_dn_en[cutoff=='All'], aes(y=signif_height_a, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = or_dn_en[cutoff=='At most three'], aes(y=signif_height_a, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=or_dn_en[, total_pos], y=total_height_a, label=or_dn_en[, total_label], size = 3) +
  annotate("text", x=0.25, y=total_height_a, label='n =', size = 3) +
  scale_y_log10(
    breaks = c(0.3, 0.5, 1, 3, 5), 
    minor_breaks = c(c(3:9)/10, c(1:5))
  ) +
  coord_cartesian(ylim=c(0.3, 12), xlim=c(1, 4), clip="off")

p_a_raw

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

p_a




#### b. OUTRIDER, leu13, leukemia ocg #### 
figure_2b_dt <- fread(paste0(dir_path, 'aggregated_data_figure/figure_2b.csv'))
exp_gene <- figure_2b_dt[, exp_gene]
or_up_gene_0 <- figure_2b_dt[or_up_gene_0==TRUE, exp_gene]
or_up_gene_1 <- figure_2b_dt[or_up_gene_1==TRUE, exp_gene]
or_up_gene_2 <- figure_2b_dt[or_up_gene_2==TRUE, exp_gene]
or_up_gene_5 <- figure_2b_dt[or_up_gene_5==TRUE, exp_gene]
or_up_top3_gene_0 <- figure_2b_dt[or_up_top3_gene_0==TRUE, exp_gene]
or_up_top3_gene_1 <- figure_2b_dt[or_up_top3_gene_1==TRUE, exp_gene]
or_up_top3_gene_2 <- figure_2b_dt[or_up_top3_gene_2==TRUE, exp_gene]
or_up_top3_gene_5 <- figure_2b_dt[or_up_top3_gene_5==TRUE, exp_gene]


# all 
#' enrichment overall
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
  odds_ratio_ci_low = c(
    fisher_test(or_up_gene_0, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(or_up_gene_1, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(or_up_gene_2, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(or_up_gene_5, exp_gene, leu_ocg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(or_up_gene_0, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(or_up_gene_1, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(or_up_gene_2, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(or_up_gene_5, exp_gene, leu_ocg[, ENSGid])$conf.int[2]
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
  odds_ratio_ci_low = c(
    fisher_test(or_up_top3_gene_0, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(or_up_top3_gene_1, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(or_up_top3_gene_2, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(or_up_top3_gene_5, exp_gene, leu_ocg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(or_up_top3_gene_0, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(or_up_top3_gene_1, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(or_up_top3_gene_2, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(or_up_top3_gene_5, exp_gene, leu_ocg[, ENSGid])$conf.int[2]
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
or_up_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=list(rownames(or_up_en))]
or_up_en[total < 1000, total_label := as.character(total)]
or_up_en[category=='0', total_pos := 1]
or_up_en[category=='1', total_pos := 2]
or_up_en[category=='2-4', total_pos := 3]
or_up_en[category=="\u22655", total_pos := 4]
or_up_en[cutoff=='All', total_pos := total_pos-0.15]
or_up_en[cutoff=='At most three', total_pos := total_pos+0.15]

or_up_en[p_val.signif=='ns', odds_ratio_ci_low := NA]
or_up_en[p_val.signif=='ns', odds_ratio_ci_up := NA]
or_up_en[, max(odds_ratio)]
or_up_en[, max(odds_ratio_ci_up, na.rm=TRUE)]

signif_height_b <- 10
total_height_b <- 17

p_b_raw <- ggplot(or_up_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=odds_ratio_ci_low, ymax=odds_ratio_ci_up),
                width=0.1, position=position_dodge(width=0.7)) +
  geom_text(data = or_up_en[cutoff=='All'], aes(y=signif_height_b, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = or_up_en[cutoff=='At most three'], aes(y=signif_height_b, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=or_up_en[, total_pos], y=total_height_b, label=or_up_en[, total_label], size = 3) +
  annotate("text", x=0.25, y=total_height_b, label='n =', size = 3) +
  scale_y_log10(
    breaks = c(0.3, 0.5, 1, 3, 5), 
    minor_breaks = c(c(3:9)/10, c(1:5))
  ) +
  coord_cartesian(ylim=c(0.3, 12), xlim=c(1, 4), clip="off")

p_b_raw

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

p_b




#### c. Activation, leu13, leukemia ocg ####
figure_2c_dt <- fread(paste0(dir_path, 'aggregated_data_figure/figure_2c.csv'))
exp_gene <- figure_2c_dt[, exp_gene]
ac_gene_0 <- figure_2c_dt[ac_gene_0==TRUE, exp_gene]
ac_gene_1 <- figure_2c_dt[ac_gene_1==TRUE, exp_gene]
ac_gene_2 <- figure_2c_dt[ac_gene_2==TRUE, exp_gene]
ac_gene_5 <- figure_2c_dt[ac_gene_5==TRUE, exp_gene]
ac_top3_gene_0 <- figure_2c_dt[ac_top3_gene_0==TRUE, exp_gene]
ac_top3_gene_1 <- figure_2c_dt[ac_top3_gene_1==TRUE, exp_gene]
ac_top3_gene_2 <- figure_2c_dt[ac_top3_gene_2==TRUE, exp_gene]
ac_top3_gene_5 <- figure_2c_dt[ac_top3_gene_5==TRUE, exp_gene]

# all 
#' enrichment overall
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
  odds_ratio_ci_low = c(
    fisher_test(ac_gene_0, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(ac_gene_1, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(ac_gene_2, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(ac_gene_5, exp_gene, leu_ocg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(ac_gene_0, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(ac_gene_1, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(ac_gene_2, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(ac_gene_5, exp_gene, leu_ocg[, ENSGid])$conf.int[2]
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
  odds_ratio_ci_low = c(
    fisher_test(ac_top3_gene_0, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(ac_top3_gene_1, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(ac_top3_gene_2, exp_gene, leu_ocg[, ENSGid])$conf.int[1],
    fisher_test(ac_top3_gene_5, exp_gene, leu_ocg[, ENSGid])$conf.int[1]
  ),
  odds_ratio_ci_up = c(
    fisher_test(ac_top3_gene_0, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(ac_top3_gene_1, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(ac_top3_gene_2, exp_gene, leu_ocg[, ENSGid])$conf.int[2],
    fisher_test(ac_top3_gene_5, exp_gene, leu_ocg[, ENSGid])$conf.int[2]
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
ac_en[, total_label := paste0(round(total/1000, digits=1), 'k'), by=list(rownames(ac_en))]
ac_en[total < 1000, total_label := as.character(total)]
ac_en[category=='0', total_pos := 1]
ac_en[category=='1', total_pos := 2]
ac_en[category=='2-4', total_pos := 3]
ac_en[category=="\u22655", total_pos := 4]
ac_en[cutoff=='All', total_pos := total_pos-0.15]
ac_en[cutoff=='At most three', total_pos := total_pos+0.15]

ac_en[p_val.signif=='ns', odds_ratio_ci_low := NA]
ac_en[p_val.signif=='ns', odds_ratio_ci_up := NA]
ac_en[, max(odds_ratio)]
ac_en[, max(odds_ratio_ci_up, na.rm=TRUE)]
ac_en[, min(odds_ratio)]
ac_en[, min(odds_ratio_ci_low, na.rm=TRUE)]

signif_height_c <- 18
total_height_c <- 35

p_c_raw <- ggplot(ac_en, aes(x=category, y=odds_ratio, fill = cutoff)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  geom_errorbar(aes(ymin=odds_ratio_ci_low, ymax=odds_ratio_ci_up),
                width=0.1, position=position_dodge(width=0.7)) +
  geom_text(data = ac_en[cutoff=='All'], aes(y=signif_height_c, label = p_val.signif), 
            stat = "identity", nudge_x = -0.18, nudge_y = 0.05, size = 3) +
  geom_text(data = ac_en[cutoff=='At most three'], aes(y=signif_height_c, label = p_val.signif), 
            stat = "identity", nudge_x = 0.18, nudge_y = 0.05, size = 3) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  annotate("text", x=ac_en[, total_pos], y=total_height_c, label=ac_en[, total_label], size = 3) +
  annotate("text", x=0.25, y=total_height_c, label='n =', size = 3) +
  scale_y_log10(
    breaks = c(0.1, 0.3, 0.5, 1, 3, 5), 
    minor_breaks = c(c(1:9)/10, c(1:5))
  ) +
  coord_cartesian(ylim=c(0.1, 20), xlim=c(1, 4), clip="off")

p_c_raw

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

p_c




### arrange plot #### 
p2_abc <- grid.arrange(p_a, p_b, p_c, nrow = 2)

p2_abc

