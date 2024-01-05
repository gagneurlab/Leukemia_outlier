theme_vale <- theme_bw() +
  theme(legend.position="right",
        plot.title = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text=element_text(size = 10),
        strip.text=element_text(size=10),
        plot.margin = margin(14, 14, 14, 20, "points")
  )
