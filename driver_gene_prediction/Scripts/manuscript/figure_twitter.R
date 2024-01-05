or_en <- data.table(category = 'OUTRIDER\nunderexpression\noutliers',
                    odds_ratio = 4.7679912, 
                    p_val.signif = '***')
fr_en <- data.table(category = 'FRASER\nsplicing\noutliers',
                    odds_ratio = 13.8946646, 
                    p_val.signif = '****')
ab_en <- data.table(category = 'AbSplice\npredicted\nvariants',
                    odds_ratio = 15.3136146, 
                    p_val.signif = '****')

en_dt <- rbind(or_en, fr_en, ab_en)
en_dt[, category := factor(category, levels = c('OUTRIDER\nunderexpression\noutliers',
                                                'FRASER\nsplicing\noutliers',
                                                'AbSplice\npredicted\nvariants'))]

signif_height <- 17
total_height <- 25

p_en_raw <- ggplot(en_dt, aes(x=category, y=odds_ratio)) +
  geom_bar(stat="identity", width=0.7, alpha=0.7) +
  geom_text(data = en_dt, aes(y=signif_height, label = p_val.signif), 
            stat = "identity", nudge_y = 0.05, size = 5) +
  geom_hline(yintercept=1, linetype="dashed", color = "firebrick") +
  scale_y_log10(
    breaks = c(1, 3, 10),
    minor_breaks = c(c(1:9)/10, c(1:10), c(1:10)*10)
  ) +
  coord_cartesian(ylim=c(0.8, 20), xlim=c(1, 3), clip="off")

# p_en_raw

p_en <- p_en_raw +
  xlab('At most three detections/predictions in at least five samples') +
  ylab('Enrichment for hematologic\ntumor supressor genes') +
  theme_vale + 
  theme(
    legend.position = c(0.5, 1.25),
    legend.direction = "vertical",
    plot.margin = margin(14, 14, 14, 20, "points"), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

p_en
