snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2_ac_anova.snakemake")

suppressPackageStartupMessages({
  library(OUTRIDER)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]

mutation_cnv <- fread(snakemake@input$mutation_cnv_gene_sample)
ac_lm_dt <- mutation_cnv[tool=='activation', ]
ac_lm_dt[, aberrant := padjust < 0.05]
ac_lm_dt <- merge(ac_lm_dt, gencode[, .(gene_id_unique, gene_id)] %>% unique(),
                  by='gene_id', all.x=TRUE, all.y=FALSE)

ods_filter_out <- readRDS(snakemake@input$ods_filter_out)

ac_aeCount_zscore_dt <- zScore(ods_filter_out) %>% as.data.table()
ac_aeCount_zscore_dt[, gene_id_unique := rownames(ods_filter_out)] 
ac_aeCount_zscore_dt <- melt(ac_aeCount_zscore_dt, id.vars = 'gene_id_unique', 
                             variable.name = "sampleID", value.name = "aeCount_zscore")

ac_aeCount_dt <- counts(ods_filter_out, normalized=TRUE) %>% as.data.table()
ac_aeCount_dt[, gene_id_unique := rownames(ods_filter_out)] 
ac_aeCount_dt <- melt(ac_aeCount_dt, id.vars = 'gene_id_unique', 
                      variable.name = "sampleID", value.name = "aeCount")

ac_sfCount_dt <- counts(ods_filter_out, normalized = FALSE) %>% as.data.table()
size_factor <- sizeFactors(ods_filter_out)
ac_sfCount_dt <- t(t(ac_sfCount_dt)/size_factor) %>% as.data.table()
ac_sfCount_zscore_dt <- t(scale(t(ac_sfCount_dt))) %>% as.data.table()
ac_sfCount_dt[, gene_id_unique := rownames(ods_filter_out)] 
ac_sfCount_zscore_dt[, gene_id_unique := rownames(ods_filter_out)] 
ac_sfCount_dt <- melt(ac_sfCount_dt, id.vars = 'gene_id_unique', 
                      variable.name = "sampleID", value.name = "sfCount")
ac_sfCount_zscore_dt <- melt(ac_sfCount_zscore_dt, id.vars = 'gene_id_unique', 
                             variable.name = "sampleID", value.name = "sfCount_zscore")



ac_lm_dt <- merge(ac_lm_dt, ac_aeCount_zscore_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# ac_lm_dt[, table(is.na(zScore))]

ac_lm_dt <- merge(ac_lm_dt, ac_sfCount_zscore_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# ac_lm_dt[, table(is.na(zScore))]

ac_lm_dt <- merge(ac_lm_dt, ac_aeCount_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# ac_lm_dt[, table(is.na(aeCount))]

ac_lm_dt <- merge(ac_lm_dt, ac_sfCount_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# ac_lm_dt[, table(is.na(sfCount))]




### zscore - copy gain/copy loss ####
ac_lm_sfCount_zscore <- lm(sfCount_zscore ~ promoter_variant + frameshift_variant + 
                             stop_gained + splice_related_variant + structural_variant +
                             copy_number_gain + copy_number_loss,
                           data=ac_lm_dt)
summary(ac_lm_sfCount_zscore)
anova(ac_lm_sfCount_zscore)

ac_lm_aeCount_zscore <- lm(aeCount_zscore ~ promoter_variant + frameshift_variant + 
                             stop_gained + splice_related_variant + structural_variant +
                             copy_number_gain + copy_number_loss,
                           data=ac_lm_dt)
summary(ac_lm_aeCount_zscore)
anova(ac_lm_aeCount_zscore)


### zscore - copy_ratio ####
ac_lm_sfCount_zscore_ratio <- lm(sfCount_zscore ~ promoter_variant + frameshift_variant + 
                                   stop_gained + splice_related_variant + structural_variant +
                                   MEAN_LOG2_COPY_RATIO,
                                 data=ac_lm_dt)
summary(ac_lm_sfCount_zscore_ratio)
anova(ac_lm_sfCount_zscore_ratio)

ac_lm_aeCount_zscore_ratio <- lm(aeCount_zscore ~ promoter_variant + frameshift_variant + 
                                   stop_gained + splice_related_variant + structural_variant +
                                   MEAN_LOG2_COPY_RATIO,
                                 data=ac_lm_dt)
summary(ac_lm_aeCount_zscore_ratio)
anova(ac_lm_aeCount_zscore_ratio)


### aberrant - copy gain/copy loss ####
ac_lm_aberrant <- lm(aberrant ~ promoter_variant + frameshift_variant +
                       stop_gained + splice_related_variant + structural_variant +
                       copy_number_gain + copy_number_loss,
                     data=ac_lm_dt)
summary(ac_lm_aberrant)
anova(ac_lm_aberrant)


### plot ####
anova_plot <- anova(ac_lm_sfCount_zscore)
anova_plot_dt <- data.table(mutation_type = rownames(anova_plot),
                            sum_sq = anova_plot$`Sum Sq`)
anova_plot_dt[, sum_sq_ratio := sum_sq/sum(sum_sq)]
anova_plot_dt[, Category := 'Activation outliers']
anova_plot_sf <- anova_plot_dt[, count_type := 'Size factor corrected count']

anova_plot <- anova(ac_lm_aeCount_zscore)
anova_plot_dt <- data.table(mutation_type = rownames(anova_plot),
                            sum_sq = anova_plot$`Sum Sq`)
anova_plot_dt[, sum_sq_ratio := sum_sq/sum(sum_sq)]
anova_plot_dt[, Category := 'Activation outliers']
anova_plot_ae <- anova_plot_dt[, count_type := 'Autoencoder corrected count']

anova_plot_dt <- rbind(anova_plot_sf, anova_plot_ae)

# fwrite(anova_plot_dt, snakemake@output$ac_anova_res)
fwrite(anova_plot_dt, "/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/ac_anova_res_sergey.csv")
anova_plot_dt <- fread("/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/ac_anova_res_sergey.csv")

anova_plot_dt[, mutation_type := factor(mutation_type, 
                                        levels = c(
                                          "copy_number_loss",   
                                          "copy_number_gain",
                                          "splice_related_variant",
                                          "frameshift_variant",
                                          "stop_gained",
                                          "structural_variant",
                                          "promoter_variant",
                                          "Residuals"
                                        ),
                                        labels = c(
                                          "Copy number loss",   
                                          "Copy number gain",
                                          "VEP splice-related",
                                          "VEP frameshift",
                                          "VEP stop-gained",
                                          "Structural",
                                          "Promoter (TSS±2Kbp)",
                                          "Residuals"
                                        ))]

anova_plot_dt[, count_type := factor(count_type, levels=c('Size factor corrected count', 'Autoencoder corrected count'))]

p_a <- ggplot(anova_plot_dt[mutation_type!="Residuals", ], 
              aes(x=mutation_type, y=sum_sq_ratio)) +
  facet_wrap('count_type', ncol=1) +
  geom_bar(stat="identity", width=0.7) +
  xlab('Mutation type') +
  ylab('Ratio of variance explained') +
  ggtitle('Activation outliers') +
  theme_vale +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

p_a


### plot ratio ####
anova_plot <- anova(ac_lm_sfCount_zscore_ratio)
anova_plot_dt <- data.table(mutation_type = rownames(anova_plot),
                            sum_sq = anova_plot$`Sum Sq`)
anova_plot_dt[, sum_sq_ratio := sum_sq/sum(sum_sq)]
anova_plot_dt[, Category := 'Activation outliers']
anova_plot_sf <- anova_plot_dt[, count_type := 'Size factor corrected count']

anova_plot <- anova(ac_lm_aeCount_zscore_ratio)
anova_plot_dt <- data.table(mutation_type = rownames(anova_plot),
                            sum_sq = anova_plot$`Sum Sq`)
anova_plot_dt[, sum_sq_ratio := sum_sq/sum(sum_sq)]
anova_plot_dt[, Category := 'Activation outliers']
anova_plot_ae <- anova_plot_dt[, count_type := 'Autoencoder corrected count']

anova_plot_dt <- rbind(anova_plot_sf, anova_plot_ae)

# fwrite(anova_plot_dt, snakemake@output$ac_anova_res)
fwrite(anova_plot_dt, "/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/ac_anova_res_sergey_ratio.csv")
anova_plot_dt <- fread("/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/ac_anova_res_sergey_ratio.csv")

anova_plot_dt[, mutation_type := factor(mutation_type, 
                                        levels = c(
                                          "MEAN_LOG2_COPY_RATIO",   
                                          "splice_related_variant",
                                          "frameshift_variant",
                                          "stop_gained",
                                          "structural_variant",
                                          "promoter_variant",
                                          "Residuals"
                                        ),
                                        labels = c(
                                          "Copy ratio",
                                          "VEP splice-related",
                                          "VEP frameshift",
                                          "VEP stop-gained",
                                          "Structural",
                                          "Promoter (TSS±2Kbp)",
                                          "Residuals"
                                        ))]

anova_plot_dt[, count_type := factor(count_type, levels=c('Size factor corrected count', 'Autoencoder corrected count'))]

p_a <- ggplot(anova_plot_dt[mutation_type!="Residuals", ], 
              aes(x=mutation_type, y=sum_sq_ratio)) +
  facet_wrap('count_type', ncol=1) +
  geom_bar(stat="identity", width=0.7) +
  xlab('Mutation type') +
  ylab('Ratio of variance explained') +
  ggtitle('Activation outliers') +
  theme_vale +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

p_a




# # > ### zscore - copy gain/copy loss ####
# > ac_lm_sfCount_zscore <- lm(sfCount_zscore ~ promoter_variant + frameshift_variant + 
#                                +                              stop_gained + splice_related_variant + structural_variant +
#                                +                              copy_number_gain + copy_number_loss,
#                              +                            data=ac_lm_dt)
# > summary(ac_lm_sfCount_zscore)
# 
# Call:
#   lm(formula = sfCount_zscore ~ promoter_variant + frameshift_variant + 
#        stop_gained + splice_related_variant + structural_variant + 
#        copy_number_gain + copy_number_loss, data = ac_lm_dt)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -2.062 -0.289 -0.157 -0.037 61.291 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                -0.0084885  0.0002235 -37.983  < 2e-16 ***
#   promoter_variantTRUE        0.0018367  0.0009396   1.955   0.0506 .  
# frameshift_variantTRUE     -0.0127390  0.0112813  -1.129   0.2588    
# stop_gainedTRUE             0.0152972  0.0116115   1.317   0.1877    
# splice_related_variantTRUE -0.0020091  0.0052804  -0.380   0.7036    
# structural_variantTRUE      0.0209555  0.0019605  10.689  < 2e-16 ***
#   copy_number_gainTRUE        0.2426352  0.0011364 213.510  < 2e-16 ***
#   copy_number_lossTRUE       -0.0061321  0.0012449  -4.926  8.4e-07 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.002 on 22731750 degrees of freedom
# (11315 observations deleted due to missingness)
# Multiple R-squared:  0.002012,	Adjusted R-squared:  0.002012 
# F-statistic:  6548 on 7 and 22731750 DF,  p-value: < 2.2e-16
# 
# > anova(ac_lm_sfCount_zscore)
# Analysis of Variance Table
# 
# Response: sfCount_zscore
# Df   Sum Sq Mean Sq    F value    Pr(>F)    
# promoter_variant              1       12      12    11.7543  0.000607 ***
#   frameshift_variant            1        1       1     0.9733  0.323855    
# stop_gained                   1        2       2     2.0702  0.150200    
# splice_related_variant        1        0       0     0.0139  0.906126    
# structural_variant            1      236     236   235.6902 < 2.2e-16 ***
#   copy_number_gain              1    45711   45711 45563.7699 < 2.2e-16 ***
#   copy_number_loss              1       24      24    24.2649 8.395e-07 ***
#   Residuals              22731750 22805116       1                         
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > 
#   > ac_lm_aeCount_zscore <- lm(aeCount_zscore ~ promoter_variant + frameshift_variant + 
#                                  +                              stop_gained + splice_related_variant + structural_variant +
#                                  +                              copy_number_gain + copy_number_loss,
#                                +                            data=ac_lm_dt)
# > summary(ac_lm_aeCount_zscore)
# 
# Call:
#   lm(formula = aeCount_zscore ~ promoter_variant + frameshift_variant + 
#        stop_gained + splice_related_variant + structural_variant + 
#        copy_number_gain + copy_number_loss, data = ac_lm_dt)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -8.097 -0.661 -0.204  0.501 33.568 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                -0.0068793  0.0002230 -30.849  < 2e-16 ***
#   promoter_variantTRUE        0.0015488  0.0009375   1.652   0.0985 .  
# frameshift_variantTRUE     -0.0231539  0.0112566  -2.057   0.0397 *  
#   stop_gainedTRUE            -0.0121453  0.0115861  -1.048   0.2945    
# splice_related_variantTRUE -0.0132683  0.0052689  -2.518   0.0118 *  
#   structural_variantTRUE      0.0137254  0.0019562   7.016 2.28e-12 ***
#   copy_number_gainTRUE        0.2384681  0.0011339 210.302  < 2e-16 ***
#   copy_number_lossTRUE       -0.0578238  0.0012421 -46.551  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.9994 on 22731750 degrees of freedom
# (11315 observations deleted due to missingness)
# Multiple R-squared:  0.002028,	Adjusted R-squared:  0.002028 
# F-statistic:  6599 on 7 and 22731750 DF,  p-value: < 2.2e-16
# 
# > anova(ac_lm_aeCount_zscore)
# Analysis of Variance Table
# 
# Response: aeCount_zscore
# Df   Sum Sq Mean Sq    F value    Pr(>F)    
# promoter_variant              1        9       9     8.5288  0.003496 ** 
#   frameshift_variant            1        4       4     4.0640  0.043806 *  
#   stop_gained                   1        1       1     0.9035  0.341833    
# splice_related_variant        1        5       5     5.3105  0.021197 *  
#   structural_variant            1      102     102   102.1388 < 2.2e-16 ***
#   copy_number_gain              1    43853   43853 43903.1627 < 2.2e-16 ***
#   copy_number_loss              1     2165    2165  2167.0408 < 2.2e-16 ***
#   Residuals              22731750 22705595       1                         
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# > ### zscore - copy_ratio ####
# > ac_lm_aeCount_zscore_ratio <- lm(aeCount_zscore ~ promoter_variant + frameshift_variant + 
#                                      +                        stop_gained + splice_related_variant + structural_variant +
#                                      +                          MEAN_LOG2_COPY_RATIO,
#                                    +                      data=ac_lm_dt)
# > summary(ac_lm_aeCount_zscore_ratio)
# 
# Call:
#   lm(formula = aeCount_zscore ~ promoter_variant + frameshift_variant + 
#        stop_gained + splice_related_variant + structural_variant + 
#        MEAN_LOG2_COPY_RATIO, data = ac_lm_dt)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -8.104 -0.662 -0.205  0.501 33.562 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.0003144  0.0002170   1.449  0.14743    
# promoter_variantTRUE        0.0027187  0.0009383   2.897  0.00376 ** 
#   frameshift_variantTRUE     -0.0220210  0.0112667  -1.955  0.05064 .  
# stop_gainedTRUE            -0.0110815  0.0115965  -0.956  0.33928    
# splice_related_variantTRUE -0.0125082  0.0052736  -2.372  0.01770 *  
#   structural_variantTRUE      0.0224547  0.0019565  11.477  < 2e-16 ***
#   MEAN_LOG2_COPY_RATIO        0.0338857  0.0004603  73.610  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1 on 22731751 degrees of freedom
# (11315 observations deleted due to missingness)
# Multiple R-squared:  0.0002436,	Adjusted R-squared:  0.0002434 
# F-statistic: 923.2 on 6 and 22731751 DF,  p-value: < 2.2e-16
# 
# > anova(ac_lm_aeCount_zscore_ratio)
# Analysis of Variance Table
# 
# Response: aeCount_zscore
# Df   Sum Sq Mean Sq   F value    Pr(>F)    
# promoter_variant              1        9     8.5    8.5135  0.003525 ** 
#   frameshift_variant            1        4     4.1    4.0568  0.043994 *  
#   stop_gained                   1        1     0.9    0.9019  0.342265    
# splice_related_variant        1        5     5.3    5.3010  0.021313 *  
#   structural_variant            1      102   102.0  101.9565 < 2.2e-16 ***
#   MEAN_LOG2_COPY_RATIO          1     5422  5421.9 5418.4363 < 2.2e-16 ***
#   Residuals              22731751 22746190     1.0                        
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# > ### aberrant - copy gain/copy loss ####
# > ac_lm_aberrant <- lm(aberrant ~ promoter_variant + frameshift_variant + 
#                          +                        stop_gained + splice_related_variant + structural_variant +
#                          +                        copy_number_gain + copy_number_loss,
#                        +                      data=ac_lm_dt)
# > summary(ac_lm_aberrant)
# 
# Call:
#   lm(formula = aberrant ~ promoter_variant + frameshift_variant + 
#        stop_gained + splice_related_variant + structural_variant + 
#        copy_number_gain + copy_number_loss, data = ac_lm_dt)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.00302 -0.00039 -0.00039 -0.00039  0.99961 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                3.895e-04  4.779e-06  81.502   <2e-16 ***
#   promoter_variantTRUE       6.615e-06  2.009e-05   0.329   0.7419    
# frameshift_variantTRUE     2.724e-04  2.413e-04   1.129   0.2589    
# stop_gainedTRUE            3.176e-04  2.483e-04   1.279   0.2008    
# splice_related_variantTRUE 2.119e-04  1.129e-04   1.878   0.0604 .  
# structural_variantTRUE     5.115e-04  4.189e-05  12.213   <2e-16 ***
#   copy_number_gainTRUE       1.528e-03  2.430e-05  62.882   <2e-16 ***
#   copy_number_lossTRUE       2.952e-04  2.661e-05  11.095   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.02142 on 22743065 degrees of freedom
# Multiple R-squared:  0.0001894,	Adjusted R-squared:  0.0001891 
# F-statistic: 615.4 on 7 and 22743065 DF,  p-value: < 2.2e-16
# 
# > anova(ac_lm_aberrant)
# Analysis of Variance Table
# 
# Response: aberrant
# Df  Sum Sq Mean Sq   F value  Pr(>F)    
# promoter_variant              1     0.0 0.00042    0.9117 0.33966    
# frameshift_variant            1     0.0 0.00070    1.5160 0.21823    
# stop_gained                   1     0.0 0.00084    1.8198 0.17734    
# splice_related_variant        1     0.0 0.00185    4.0421 0.04438 *  
#   structural_variant            1     0.1 0.09053  197.2882 < 2e-16 ***
#   copy_number_gain              1     1.8 1.82600 3979.3010 < 2e-16 ***
#   copy_number_loss              1     0.1 0.05649  123.1083 < 2e-16 ***
#   Residuals              22743065 10436.2 0.00046                      
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
