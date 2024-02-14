source("Scripts/manuscript/function.R")
source("Scripts/manuscript/manuscript_theme.R")



#### res pred ####
snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_4.snakemake")

experiment_design <- fread(snakemake@params$experimentDesign)
experiment_design[is.na(experiment_design)] <- ''

project_dir <- snakemake@params$projectPath
intogen_dir <- snakemake@params$intogenDir
vep_dir <- snakemake@params$vep_path
single_group <- snakemake@params$inputDatasets

n_plot <- 200
# read in res_viz
exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == "leukemia_14group" &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == '', ]
exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)
res_viz_pred <- fread(exp_viz[, res_post_path])
predicted_genes_omics <- res_viz_pred[1:200, GeneSymbol]

exp_viz <- experiment_design[label_gene_list == 'MLL_CGC_leukemia_gene' &
                               sample_group == "leukemia_14group" &
                               model_method == 'rf' &
                               intogen_input_feature == 'clustl,hotmaps,smregions,fml,cbase,mutpanning,dndscv' &
                               outlier_input_feature == 'or,ac,absplice,fr' &
                               coess_input_feature == 'emb_string', ]
exp_viz <- add_result_paths(exp_viz, project_dir, intogen_dir, vep_dir)
res_viz_pred <- fread(exp_viz[, res_post_path])
predicted_genes_omics_emb <- res_viz_pred[1:200, GeneSymbol]

# overlap genes
gene_overlap <- intersect(predicted_genes_omics, predicted_genes_omics_emb)

# predicted by omics only
gene_leu_only <- setdiff(predicted_genes_omics, predicted_genes_omics_emb) 

# predicted after adding emb
gene_emb_only <- setdiff(predicted_genes_omics_emb, predicted_genes_omics) 

signal_dt <- rbind(
  data.table(geneSymbol = gene_overlap,
             status = 'overlap'), 
  data.table(geneSymbol = gene_leu_only,
             status = 'leu_only'), 
  data.table(geneSymbol = gene_emb_only,
             status = 'emb_only') 
)


#### or ac ####
snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/figure_2.snakemake")
or_res <- fread(snakemake@input$outriderRes)
or_res[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]
or_res[, geneID_short := strsplit(geneID, "[.]")[[1]][1], by=rownames(or_res)]

ac_res <- fread(snakemake@input$activtionRes)
ac_res[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]
ac_res[, geneID_short := strsplit(geneID, "[.]")[[1]][1], by=rownames(ac_res)]

# signal_dt <- merge(signal_dt,or_res[, .N, by='hgncSymbol'], 
#                    by.x='geneSymbol', by.y='hgncSymbol', all.x=TRUE, all.y=FALSE)
# setnames(signal_dt, 'N', 'or')
# signal_dt[is.na(or), or := 0]
# 
# signal_dt <- merge(signal_dt,ac_res[, .N, by='hgncSymbol'], 
#                    by.x='geneSymbol', by.y='hgncSymbol', all.x=TRUE, all.y=FALSE)
# setnames(signal_dt, 'N', 'ac')
# signal_dt[is.na(ac), ac := 0]

signal_dt <- merge(signal_dt,
                   rbind(or_res[, .N, by='hgncSymbol'], ac_res[, .N, by='hgncSymbol']), 
                   by.x='geneSymbol', by.y='hgncSymbol', all.x=TRUE, all.y=FALSE)
setnames(signal_dt, 'N', 'exp_outlier')
signal_dt[is.na(exp_outlier), exp_outlier := 0]

signal_dt[, mean(exp_outlier), by='status']
# signal_dt[, mean(or), by='status']
# signal_dt[, mean(ac), by='status']




#### fr absplice ####
snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_3.snakemake")
gencode <- fread(snakemake@params$gencode)
for (x in snakemake@input$fraserRes) {
  fr_res_temp <- fread(x)
  fr_res_temp <- merge(fr_res_temp, gencode[, .(gene_name, gene_type)], 
                       by.x='hgncSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  fr_res_temp <- fr_res_temp[gene_type=='protein_coding', ]
  
  if (exists("fr_res")) {
    fr_res <- rbind(fr_res, fr_res_temp, fill=TRUE)
  } else {
    fr_res <- fr_res_temp
  }
}
fr_res[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]

absplice_res <- fread(snakemake@input$abspliceRes)
absplice_res[, samp_symbol := paste0(sampleID, "-", gene_name)]
absplice_res <- separate(absplice_res, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_res <- separate(absplice_res, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()

signal_dt <- merge(signal_dt,fr_res[, .N, by='hgncSymbol'], 
                   by.x='geneSymbol', by.y='hgncSymbol', all.x=TRUE, all.y=FALSE)
setnames(signal_dt, 'N', 'fr')
signal_dt[is.na(fr), fr := 0]

signal_dt <- merge(signal_dt,absplice_res[, .N, by='gene_name'], 
                   by.x='geneSymbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
setnames(signal_dt, 'N', 'absplice')
signal_dt[is.na(absplice), absplice := 0]

signal_dt[, mean(fr), by='status']
signal_dt[, mean(absplice), by='status']




#### vep ####
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

signal_dt <- merge(signal_dt, vep_res[IMPACT == 'HIGH', .N, by='SYMBOL'], 
                   by.x='geneSymbol', by.y='SYMBOL', all.x=TRUE, all.y=FALSE)
setnames(signal_dt, 'N', 'vep_HIGH')
signal_dt[is.na(vep_HIGH), vep_HIGH := 0]

signal_dt <- merge(signal_dt, vep_res[IMPACT == 'MODERATE', .N, by='SYMBOL'], 
                   by.x='geneSymbol', by.y='SYMBOL', all.x=TRUE, all.y=FALSE)
setnames(signal_dt, 'N', 'vep_MODERATE')
signal_dt[is.na(vep_MODERATE), vep_MODERATE := 0]

signal_dt <- merge(signal_dt, vep_res[IMPACT == 'LOW', .N, by='SYMBOL'], 
                   by.x='geneSymbol', by.y='SYMBOL', all.x=TRUE, all.y=FALSE)
setnames(signal_dt, 'N', 'vep_LOW')
signal_dt[is.na(vep_LOW), vep_LOW := 0]

# signal_dt <- merge(signal_dt, vep_res[IMPACT == 'MODIFIER', .N, by='SYMBOL'], 
#                    by.x='geneSymbol', by.y='SYMBOL', all.x=TRUE, all.y=FALSE)
# setnames(signal_dt, 'N', 'vep_MODIFIER')
# signal_dt[is.na(vep_MODIFIER), vep_MODIFIER := 0]
# 
# signal_dt[, mean(vep_MODIFIER), by='status']




#### visualize #### 
signal_melt <- melt(signal_dt, id.vars = c('geneSymbol', 'status'), 
                    variable.name = "method", value.name = "count")
signal_melt[, status := factor(status, levels = c('leu_only', 'overlap', 'emb_only'))]

stat.test <- ggpubr::compare_means(formula = count ~ status, 
                                   data = signal_melt, 
                                   method = "wilcox.test")
stat.test

ggplot(signal_melt, aes(x=status, y=count+1)) +
  geom_boxplot() +
  facet_wrap("method", scales = 'free_y') +
  scale_y_log10()
