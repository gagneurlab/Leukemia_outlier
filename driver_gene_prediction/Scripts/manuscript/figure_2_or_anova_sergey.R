snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2_or_anova.snakemake")

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
or_lm_dt <- mutation_cnv[tool=='outrider', ]
or_lm_dt[, aberrant := padjust < 0.05]
or_lm_dt <- merge(or_lm_dt, gencode[, .(gene_id_unique, gene_id)] %>% unique(),
                  by='gene_id', all.x=TRUE, all.y=FALSE)

ods <- readRDS(snakemake@input$ods)

or_aeCount_zscore_dt <- zScore(ods) %>% as.data.table()
or_aeCount_zscore_dt[, gene_id_unique := rownames(ods)] 
or_aeCount_zscore_dt <- melt(or_aeCount_zscore_dt, id.vars = 'gene_id_unique', 
                             variable.name = "sampleID", value.name = "aeCount_zscore")

or_aeCount_dt <- counts(ods, normalized=TRUE) %>% as.data.table()
or_aeCount_dt[, gene_id_unique := rownames(ods)] 
or_aeCount_dt <- melt(or_aeCount_dt, id.vars = 'gene_id_unique', 
                      variable.name = "sampleID", value.name = "aeCount")

or_sfCount_dt <- counts(ods, normalized = FALSE) %>% as.data.table()
size_factor <- sizeFactors(ods)
or_sfCount_dt <- t(t(or_sfCount_dt)/size_factor) %>% as.data.table()
or_sfCount_zscore_dt <- t(scale(t(or_sfCount_dt))) %>% as.data.table()
or_sfCount_dt[, gene_id_unique := rownames(ods)] 
or_sfCount_zscore_dt[, gene_id_unique := rownames(ods)] 
or_sfCount_dt <- melt(or_sfCount_dt, id.vars = 'gene_id_unique', 
                      variable.name = "sampleID", value.name = "sfCount")
or_sfCount_zscore_dt <- melt(or_sfCount_zscore_dt, id.vars = 'gene_id_unique', 
                             variable.name = "sampleID", value.name = "sfCount_zscore")


or_lm_dt <- merge(or_lm_dt, or_aeCount_zscore_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# or_lm_dt[, table(is.na(zScore))]

or_lm_dt <- merge(or_lm_dt, or_sfCount_zscore_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# or_lm_dt[, table(is.na(zScore))]

or_lm_dt <- merge(or_lm_dt, or_aeCount_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# or_lm_dt[, table(is.na(aeCount))]

or_lm_dt <- merge(or_lm_dt, or_sfCount_dt, 
                  by=c('gene_id_unique', 'sampleID'), all.x=TRUE, all.y=FALSE)
# or_lm_dt[, table(is.na(sfCount))]


or_lm_dn_dt <- or_lm_dt[zScore < 0, ]
or_lm_up_dt <- or_lm_dt[zScore > 0, ]




### zscore - copy gain/copy loss ####
or_lm_sfCount_zscore <- lm(sfCount_zscore ~ promoter_variant + frameshift_variant + 
                             stop_gained + splice_related_variant + structural_variant +
                             copy_number_gain + copy_number_loss,
                           data=or_lm_dt)
summary(or_lm_sfCount_zscore)
anova(or_lm_sfCount_zscore)

or_lm_aeCount_zscore <- lm(aeCount_zscore ~ promoter_variant + frameshift_variant + 
                             stop_gained + splice_related_variant + structural_variant +
                             copy_number_gain + copy_number_loss,
                           data=or_lm_dt)
summary(or_lm_aeCount_zscore)
anova(or_lm_aeCount_zscore)


### zscore - copy_ratio ####
or_lm_sfCount_zscore_ratio <- lm(sfCount_zscore ~ promoter_variant + frameshift_variant + 
                                   stop_gained + splice_related_variant + structural_variant +
                                   MEAN_LOG2_COPY_RATIO,
                                 data=or_lm_dt)
summary(or_lm_sfCount_zscore_ratio)
anova(or_lm_sfCount_zscore_ratio)

or_lm_aeCount_zscore_ratio <- lm(aeCount_zscore ~ promoter_variant + frameshift_variant + 
                                   stop_gained + splice_related_variant + structural_variant +
                                   MEAN_LOG2_COPY_RATIO,
                                 data=or_lm_dt)
summary(or_lm_aeCount_zscore_ratio)
anova(or_lm_aeCount_zscore_ratio)


### aberrant - copy gain/copy loss ####
or_lm_aberrant <- lm(aberrant ~ promoter_variant + frameshift_variant +
                       stop_gained + splice_related_variant + structural_variant +
                       copy_number_gain + copy_number_loss,
                     data=or_lm_dt)
summary(or_lm_aberrant)
anova(or_lm_aberrant)


### aberrant - downregulated ####
or_lm_dn_aberrant <- lm(aberrant ~ promoter_variant + frameshift_variant + 
                          stop_gained + splice_related_variant + structural_variant +
                          copy_number_gain + copy_number_loss,
                        data=or_lm_dn_dt)
summary(or_lm_dn_aberrant)
anova(or_lm_dn_aberrant)


### aberrant - upregulated ####
or_lm_up_aberrant <- lm(aberrant ~ promoter_variant + frameshift_variant + 
                          stop_gained + splice_related_variant + structural_variant +
                          copy_number_gain + copy_number_loss,
                        data=or_lm_up_dt)
summary(or_lm_dn_aberrant)
anova(or_lm_up_aberrant)


### plot ####
anova_plot <- anova(or_lm_sfCount_zscore)
anova_plot_dt <- data.table(mutation_type = rownames(anova_plot),
                            sum_sq = anova_plot$`Sum Sq`)
anova_plot_dt[, sum_sq_ratio := sum_sq/sum(sum_sq)]
anova_plot_dt[, Category := 'Activation outliers']
anova_plot_sf <- anova_plot_dt[, count_type := 'Size factor corrected count']

anova_plot <- anova(or_lm_aeCount_zscore)
anova_plot_dt <- data.table(mutation_type = rownames(anova_plot),
                            sum_sq = anova_plot$`Sum Sq`)
anova_plot_dt[, sum_sq_ratio := sum_sq/sum(sum_sq)]
anova_plot_dt[, Category := 'Activation outliers']
anova_plot_ae <- anova_plot_dt[, count_type := 'Autoencoder corrected count']

# fwrite(anova_plot_dt, snakemake@output$or_anova_res)
fwrite(anova_plot_dt, "/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/or_anova_res_sergey.csv")
anova_plot_dt <- fread("/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/or_anova_res_sergey.csv")

anova_plot_dt <- rbind(anova_plot_sf, anova_plot_ae)
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
  ggtitle('Expression outliers') +
  theme_vale +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

p_a



### plot ratio ####
anova_plot <- anova(or_lm_sfCount_zscore_ratio)
anova_plot_dt <- data.table(mutation_type = rownames(anova_plot),
                            sum_sq = anova_plot$`Sum Sq`)
anova_plot_dt[, sum_sq_ratio := sum_sq/sum(sum_sq)]
anova_plot_dt[, Category := 'Activation outliers']
anova_plot_sf <- anova_plot_dt[, count_type := 'Size factor corrected count']

anova_plot <- anova(or_lm_aeCount_zscore_ratio)
anova_plot_dt <- data.table(mutation_type = rownames(anova_plot),
                            sum_sq = anova_plot$`Sum Sq`)
anova_plot_dt[, sum_sq_ratio := sum_sq/sum(sum_sq)]
anova_plot_dt[, Category := 'Activation outliers']
anova_plot_ae <- anova_plot_dt[, count_type := 'Autoencoder corrected count']

anova_plot_dt <- rbind(anova_plot_sf, anova_plot_ae)

# fwrite(anova_plot_dt, snakemake@output$or_anova_res)
fwrite(anova_plot_dt, "/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/or_anova_res_sergey_ratio.csv")
anova_plot_dt <- fread("/s/project/vale/driver_prediction_202402/manuscript/figure_2/plot_data/anova/or_anova_res_sergey_ratio.csv")

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
  ggtitle('Expression outliers') +
  theme_vale +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
  )

p_a

