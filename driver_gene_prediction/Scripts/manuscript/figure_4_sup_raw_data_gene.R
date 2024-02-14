source("Scripts/manuscript/function.R")
source("Scripts/manuscript/manuscript_theme.R")

snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2.snakemake")

samp_anno <- fread(snakemake@params$sampAnno)
single_group <- snakemake@params$inputDatasets
samp_anno_14group <- samp_anno[grep(single_group, DROP_GROUP),]

snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2.snakemake")
or_res <- fread(snakemake@input$outriderRes)
or_res[, samp_symbol := paste0(sampleID, "-", sampleID)]
or_res[, geneID_short := strsplit(geneID, "[.]")[[1]][1], by=rownames(or_res)]

ac_res <- fread(snakemake@input$activtionRes)
ac_res[, samp_symbol := paste0(sampleID, "-", sampleID)]
ac_res[, geneID_short := strsplit(geneID, "[.]")[[1]][1], by=rownames(ac_res)]

snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_3.snakemake")
gencode <- fread(snakemake@params$gencode)
for (x in snakemake@input$fraserResJunc) {
  fr_res_temp <- fread(x)
  fr_res_temp <- merge(fr_res_temp, gencode[, .(gene_name, gene_type)], 
                       by.x='hgncSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  fr_res_temp <- fr_res_temp[gene_type=='protein_coding', ]
  
  if (exists("fr_res_junc")) {
    fr_res_junc <- rbind(fr_res_junc, fr_res_temp)
  } else {
    fr_res_junc <- fr_res_temp
  }
}

absplice_res_var <- fread(snakemake@input$abspliceResVar)
absplice_res_var[, samp_symbol := paste0(sampleID, "-", gene_name)]
absplice_res_var <- separate(absplice_res_var, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_res_var <- separate(absplice_res_var, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()

snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2.snakemake")
gencode <- fread(snakemake@params$gencode)
for (x in snakemake@input$vepRes) {
  vep_res_temp <- fread(x)
  vep_res_temp <- merge(vep_res_temp, gencode[, .(gene_name, gene_type)], 
                        by.x='SYMBOL', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  vep_res_temp <- vep_res_temp[gene_type=='protein_coding', ]
  
  if (exists("vep_res")) {
    vep_res <- rbind(vep_res, vep_res_temp, fill=TRUE)
  } else {
    vep_res <- vep_res_temp
  }
}
vep_res <- vep_res %>% separate('#Uploaded_variation', 
                                c("col1", "array_id", "alt", "ref", "col2"), 
                                sep='__') %>% as.data.table()

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

snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/prediction_tab.snakemake")
res_viz_pred_group <- fread(snakemake@output$prediction_diag_tab)
res_viz_pred_group[, Prediction_scaled := Prediction/max(Prediction), by='StudyGroup']

snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/prediction_emb_tab.snakemake")
res_viz_pred_emb_group <- fread(snakemake@output$prediction_emb_diag_tab)
res_viz_pred_emb_group[, Prediction_scaled := Prediction/max(Prediction), by='StudyGroup']



### define gene ####
# gene_to_curate <- 'RUNX1'
# gene_to_curate <- 'PAX5'
# gene_to_curate <- 'IKZF1'
# 
gene_list_to_curate <- c('TET2', 'ASXL1', 'RUNX1')
# gene_list_to_curate <- c('JAK3', 'STAT5B')
# gene_list_to_curate <- c('IKZF1', 'PAX5')




### pred res####
res_viz_pred_group_gene <- res_viz_pred_group[GeneSymbol %in% gene_list_to_curate, 
                                              .(GeneSymbol, StudyGroup, Prediction_scaled)]
res_viz_pred_group_gene[, input := 'Without\nEmbedding\nSTRING']

res_viz_pred_emb_group_gene <- res_viz_pred_emb_group[GeneSymbol %in% gene_list_to_curate, 
                                                      .(GeneSymbol, StudyGroup, Prediction_scaled)]
res_viz_pred_emb_group_gene[, input := 'With\nEmbedding\nSTRING']

res_gene_list <- rbind(res_viz_pred_group_gene, res_viz_pred_emb_group_gene)
res_gene_list[, StudyGroup := factor(StudyGroup, levels = names(column_ords))]
res_gene_list[, input := factor(input, levels = c('Without\nEmbedding\nSTRING',
                                                     'With\nEmbedding\nSTRING'))]

p_s11_raw <- ggplot(res_gene_list, 
                    aes(x=StudyGroup, y=Prediction_scaled)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  # facet_wrap("GeneSymbol", scales = 'free_y', ncol=1)  +
  facet_grid(input ~ GeneSymbol,  scales = 'free_y')

p_s11 <- p_s11_raw +
  # scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2")) +
  # labels = var_filter_label_ap, name = 'Additional variant filter criteria:') +
  # xlab(NULL) +
  ylab("Normalized predicted probability") +
  # ggtitle(gene_to_curate) + 
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = 'right',
    legend.margin = margin(6, 6, 6, 6),
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x=element_blank(),
    axis.text.x=element_blank()
  )

p_s11




### raw input ####
signal_percent_gene_list <- data.table(method=character(), 
                                       Study_group=character(), 
                                       percentage_with_sample=numeric(),
                                       geneSymbol=character())

for (gene_to_curate in gene_list_to_curate) {
  
  signal_dt <- samp_anno_14group[, .(ArrayID, Cohort_group)]
  
  #### or ac ####
  or_res_gene <- or_res[hgncSymbol %in% gene_to_curate, ]
  ac_res_gene <- ac_res[hgncSymbol %in% gene_to_curate, ]
  
  signal_dt <- merge(signal_dt,
                     rbind(or_res_gene[, .N, by='sampleID'], ac_res_gene[, .N, by='sampleID']), 
                     by.x='ArrayID', by.y='sampleID', all.x=TRUE, all.y=FALSE)
  setnames(signal_dt, 'N', 'expression_outlier')
  signal_dt[is.na(expression_outlier), expression_outlier := 0]
  
  signal_dt[, mean(expression_outlier), by='Cohort_group']
  
  
  
  
  #### fr absplice ####
  fr_res_gene <- fr_res_junc[hgncSymbol %in% gene_to_curate, ]
  absplice_res_gene <- absplice_res_var[gene_name %in% gene_to_curate, ]
  
  signal_dt <- merge(signal_dt,fr_res_gene[, .N, by='sampleID'], 
                     by.x='ArrayID', by.y='sampleID', all.x=TRUE, all.y=FALSE)
  setnames(signal_dt, 'N', 'splicing_outlier')
  signal_dt[is.na(splicing_outlier), splicing_outlier := 0]
  signal_dt[, mean(splicing_outlier), by='Cohort_group']
  
  signal_dt <- merge(signal_dt,absplice_res_gene[, .N, by='sampleID'], 
                     by.x='ArrayID', by.y='sampleID', all.x=TRUE, all.y=FALSE)
  setnames(signal_dt, 'N', 'absplice_variant')
  signal_dt[is.na(absplice_variant), absplice_variant := 0]
  signal_dt[, mean(absplice_variant), by='Cohort_group']
  
  
  
  
  #### vep ####
  vep_res_gene <- vep_res[SYMBOL %in% gene_to_curate, ]
  signal_dt <- merge(signal_dt, vep_res_gene[IMPACT == 'HIGH', .N, by='array_id'], 
                     by.x='ArrayID', by.y='array_id', all.x=TRUE, all.y=FALSE)
  setnames(signal_dt, 'N', 'vep_HIGH')
  signal_dt[is.na(vep_HIGH), vep_HIGH := 0]
  
  signal_dt <- merge(signal_dt, vep_res_gene[IMPACT == 'MODERATE', .N, by='array_id'], 
                     by.x='ArrayID', by.y='array_id', all.x=TRUE, all.y=FALSE)
  setnames(signal_dt, 'N', 'vep_MODERATE')
  signal_dt[is.na(vep_MODERATE), vep_MODERATE := 0]
  
  signal_dt <- merge(signal_dt, vep_res_gene[IMPACT == 'LOW', .N, by='array_id'],
                     by.x='ArrayID', by.y='array_id', all.x=TRUE, all.y=FALSE)
  setnames(signal_dt, 'N', 'vep_LOW')
  signal_dt[is.na(vep_LOW), vep_LOW := 0]
  
  
  
  
  #### combine #### 
  # signal_dt[, complete_dataset := sum(expression_outlier, splicing_outlier, absplice_variant, 
  #                                     vep_HIGH, vep_MODERATE, vep_LOW), by=rownames(signal_dt)]
  signal_melt <- melt(signal_dt, id.vars = c('ArrayID', 'Cohort_group'), 
                      variable.name = "method", value.name = "count_per_sample")
  
  snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_4.snakemake")
  manuscript_wording <- fread(snakemake@params$manuscriptWording)
  colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))
  
  signal_melt <- merge(signal_melt, manuscript_wording[, .(Study_group, Study_group_during_analysis)] %>% unique(), 
                       by.x = 'Cohort_group', by.y = 'Study_group_during_analysis', all.x = TRUE)
  
  signal_melt[, Study_group := factor(Study_group, levels = names(column_ords))]
  
  # stat.test <- ggpubr::compare_means(formula = count ~ Cohort_group, 
  #                                    data = signal_melt, 
  #                                    method = "wilcox.test")
  # stat.test
  
  # ggplot(signal_melt, aes(x=Study_group, y=count_per_sample)) +
  #   geom_point(alpha=0.5, size=0.3, position = position_jitter(w = 0.25, h = 0)) +
  #   # geom_boxplot() +
  #   facet_wrap("method", scales = 'free_y') +
  #   # scale_y_log10() +
  #   ggtitle(gene_to_curate) + 
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  signal_melt[, exist := count_per_sample > 0]
  signal_percent_group <- signal_melt[, sum(exist)/.N, by=c("method", "Study_group")]
  signal_percent_all <- signal_melt[, sum(exist)/.N, by="method"]
  signal_percent_all[, Study_group := 'Complete dataset']
  signal_percent <- rbind(signal_percent_group, signal_percent_all)
  setnames(signal_percent, 'V1', 'percentage_with_sample')
  signal_percent[, geneSymbol := gene_to_curate]
  signal_percent_gene_list <- rbind(signal_percent_gene_list, signal_percent)
}


### visualize #### 
p_s10_raw <- ggplot(signal_percent_gene_list[Study_group != 'Complete dataset'], 
                    aes(x=Study_group, y=percentage_with_sample, fill=method)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  facet_wrap("geneSymbol", scale = 'free_y', nrow=1) 

x_label <- c(
  'Expression outlier',
  'Splicing outlier',
  'AbSplice variant',
  'VEP HIGH IMPACT variant',
  'VEP MODERATE IMPACT variant',
  'VEP LOW IMPACT variant'
)
names(x_label) <- c(
  'expression_outlier', 'splicing_outlier', 'absplice_variant',
  'vep_HIGH', 'vep_MODERATE', 'vep_LOW'
)

p_s10 <- p_s10_raw +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Set2"),
                    labels = x_label, name = 'Method') +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("Study group") +
  ylab("Percentage of samples detected") +
  # ggtitle(gene_to_curate) + 
  theme_vale +
  theme(
    legend.box.background = element_rect(color="black", linewidth=0.3),
    legend.position = 'bottom',
    legend.margin = margin(6, 6, 6, 6),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )

p_s10


p_s12 <- ggarrange(p_s11, p_s10, ncol = 1, 
                   labels = c('A', 'B'), heights = c(0.5, 1))
p_s12
