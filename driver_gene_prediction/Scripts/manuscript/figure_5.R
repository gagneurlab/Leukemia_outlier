#'---
#' title: figure_5
#' author: Xueqi Cao, Ata Jadid Ahari
#' wb:
#'  py:
#'   - |
#'    annotations = config["geneAnnotation"]
#'    inputDatasets = config["cohort"]["single_group"]
#'    outputDatasets = config["cohort"]["groups"]
#'  params:
#'    - projectPath: '`sm config["projectPath"]`'
#'    - annotations: '`sm annotations`'
#'    - gencode: '`sm config["gencode"]`'
#'    - sampAnno: '`sm config["sampleAnnotation"]`'
#'    - manuscriptWording: '`sm config["manuscriptWording"]`'
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - CGC_cancer_gene_processed: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_processed.tsv"`'
#'    - MLL_CGC_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/MLL_CGC_leukemia_gene_list.tsv"`'
#'    - lrp1b_survival_hclv: '`sm config["lrp1b_survival_hclv"]`'
#'    - lrp1b_survival_hcl_mzl: '`sm config["lrp1b_survival_hcl_mzl"]`'
#'    - activation_enrichment_curation: '`sm config["activation_enrichment_curation"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'  input:
#'    - mll_cnv: '`sm expand(config["projectPath"] + "/manuscript/cnv/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_arriba: '`sm expand(config["projectPath"] + "/manuscript/arriba/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_star_fusion: '`sm expand(config["projectPath"] + "/manuscript/star_fusion/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_manta: '`sm expand(config["projectPath"] + "/manuscript/manta/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_manta_sv: '`sm expand(config["projectPath"] + "/manuscript/manta_sv/{dataset}.tsv",
#'                      dataset=outputDatasets)`'
#'    - vepRes: '`sm expand(config["vep_path"] +
#'              "/MLL_{dataset}.vep.tsv.gz", dataset=outputDatasets)`'
#'    - total_counts_marc: '`sm config["total_counts_marc"]`'
#'    - size_factor_marc: '`sm config["size_factor_marc"]`'
#'    - gene_anno_marc: '`sm config["gene_anno_marc"]`'
#'    - diagFisher: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_5/plot_data/diag_fisher.tsv"`'
#'    - activationOds: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods_filter_out.Rds",
#'                      annotation=annotations, dataset=inputDatasets)`'
#'    - activationRes: '`sm expand(config["outriderDir"] +
#'                      "/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/res_filter_out.tsv",
#'                      annotation=annotations, dataset=inputDatasets)`'
#'  output:
#'    - curation_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/curation_tab.csv"`'
#'    - CGC_cancer_gene: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/CGC_cancer_gene_processed.csv"`'
#'    - counts_marc: '`sm config["projectPath"] + 
#'                     "/manuscript/agg_table_figure/counts_marc.csv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_5.html"`'
#'  type: noindex
#'  resources:
#'    - mem_mb: 16000 
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE 
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_5.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/figure_5.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(mclust)
  library(data.table)
  library(readxl)
  library(tidyr)
  library(scales)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(rstatix)
  library(magrittr)
  library(slanter)
  library(pheatmap)
  library(OUTRIDER)
  library(tidyverse)
  library(survival)
  library(survminer)
  library(openxlsx)
  library(RColorBrewer)
  library(grid)
})
source("Scripts/manuscript/manuscript_theme.R")

options(bitmapType='cairo')

cohort_color <- c("HCL-V"=RColorBrewer::brewer.pal(8, "Paired")[6],
                  "HCL"=RColorBrewer::brewer.pal(8, "Paired")[8],
                  "MZL"=RColorBrewer::brewer.pal(8, "Set2")[6],
                  "Others"=RColorBrewer::brewer.pal(8, "Set2")[8])




# get parameters
set.seed(2023)

output_dir <- snakemake@params$htmlOutputPath

single_group <- snakemake@params$inputDatasets
sample_group <- snakemake@params$outputDatasets
sample_group <- append(sample_group, single_group)

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()
samp_anno_hcl_group <- samp_anno_exp[DROP_GROUP=='HZL_group',]

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))

cohort_dt <- samp_anno[grep(single_group, DROP_GROUP), .N, by = "Cohort_group"]
cohort_dt <- cohort_dt[order(-N),]
cohort_dt <- rbind(data.table(Cohort_group=single_group, 
                              N=samp_anno[grep(single_group, DROP_GROUP), .N]),
                   cohort_dt)

gencode <- fread(snakemake@params$gencode) 

CGC_cancer_gene <- fread(snakemake@params$CGC_cancer_gene_processed)
fwrite(CGC_cancer_gene, snakemake@output$CGC_cancer_gene)
mll_cgc_leukemia_gene <- fread(snakemake@params$MLL_CGC_leukemia_gene_list)

vep_res_hcl_group <- fread(grep('HZL_group', snakemake@input$vepRes, value = TRUE))
vep_res_hcl_group <- vep_res_hcl_group %>% separate('#Uploaded_variation', 
                                                    c("col1", "array_id", "alt", "ref", "col2"), 
                                                    sep='__') %>% as.data.table()

diag_fisher <- fread(snakemake@input$diagFisher)
activation_enrichment_curation <- read_xlsx(snakemake@params$activation_enrichment_curation) %>% as.data.table()

ac_res <- fread(snakemake@input$activationRes)
ods_filter_out <- readRDS(snakemake@input$activationOds)

total_counts_marc <- readRDS(snakemake@input$total_counts_marc)
size_factor_marc <- fread(snakemake@input$size_factor_marc, header = TRUE)
size_factor_marc[, sampleID := regmatches(V1, gregexpr("ohio[^;]+_001|HCL[^;]+_001", V1)), by=rownames(size_factor_marc)]
gene_anno_marc <- fread(snakemake@input$gene_anno_marc)
geneID_marc <- gene_anno_marc[gene_name == 'LRP1B', gene_id]

geneID_LRP1B <- gencode[gene_name == 'LRP1B', gene_id]




### Figure #### 
#### check variants/CNV/SV ####
lrp1b_act_samples <- ac_res[hgncSymbol=='LRP1B', sampleID]

#'### VEP
vep_res_hcl_group[SYMBOL=='LRP1B' & array_id %in% lrp1b_act_samples, ]
vep_res_hcl_group[SYMBOL=='LRP1B', ]
#' one missense variant overlapped with LRP1B activation --> nothing explainable

#'
#'### CNV
cnv_hcl_group <- fread(grep('HZL_group', snakemake@input$mll_cnv, value = TRUE))
cnv_hcl_group[SYMBOL=='LRP1B' & ARRAY_ID %in% lrp1b_act_samples, ]
cnv_hcl_group[SYMBOL=='LRP1B', ]
#' two minor gain of first few exons found with LRP1B activation --> nothing explainable

#'
#'### Fusion
fusion_manta_hcl_group <- fread(grep('HZL_group', snakemake@input$mll_manta, value = TRUE))
fusion_manta_hcl_group[( gene1=='LRP1B' | gene2=='LRP1B' ) & array_id %in% lrp1b_act_samples, ]
#' no fusion found for manta --> nothing explainable

fusion_arriba_hcl_group <- fread(grep('HZL_group', snakemake@input$mll_arriba, value = TRUE))
fusion_arriba_hcl_group[( gene1=='LRP1B' | gene2=='LRP1B' ) & array_id %in% lrp1b_act_samples, ]
#' no fusion found for arriba --> nothing explainable

fusion_star_fusion_hcl_group <- fread(grep('HZL_group', snakemake@input$mll_star_fusion, value = TRUE))
fusion_star_fusion_hcl_group[( gene1=='LRP1B' | gene2=='LRP1B' ) & array_id %in% lrp1b_act_samples, ]
#' no fusion found for star_fusion --> nothing explainable

#'
#'### SV
sv_manta_hcl_group <- fread(grep('HZL_group', snakemake@input$mll_manta_sv, value = TRUE))
sv_manta_hcl_group[( gene1=='LRP1B' | gene2=='LRP1B' ) & array_id %in% lrp1b_act_samples, ]
#' no sv found for star_fusion --> nothing explainable




#### association summary #####
diag_fisher[, isLeukemia := geneSymbol %in% mll_cgc_leukemia_gene[, GeneSymbol] ]

diag_fisher[, table(p_greater_adjust < 0.05)]
diag_fisher[, length(unique(gene_id))]
diag_fisher[, length(unique(Diag))]

diag_fisher[p_greater_adjust < 0.05, length(unique(Diag))]
diag_fisher[p_greater_adjust < 0.05, table(input_res)]

diag_fisher[p_greater_adjust < 0.05, sum(isCGC)/length(isCGC)]
diag_fisher_gene <- diag_fisher[, .(geneSymbol, isCGC)] %>% unique()
diag_fisher_gene[, sum(isCGC)/.N]

diag_fisher[p_greater_adjust < 0.05 & input_res == 'absplice', sum(isLeukemia)/length(isLeukemia)]
diag_fisher[p_greater_adjust < 0.05 & input_res == 'or_dn', sum(isLeukemia)/length(isLeukemia)]
diag_fisher[p_greater_adjust < 0.05 & input_res == 'or_up', sum(isLeukemia)/length(isLeukemia)]
diag_fisher[p_greater_adjust < 0.05 & input_res == 'ac', sum(isLeukemia)/length(isLeukemia)]
diag_fisher[p_greater_adjust < 0.05 & input_res == 'fr', sum(isLeukemia)/length(isLeukemia)]

diag_fisher[p_greater_adjust < 0.05 & input_res == 'absplice', sum(isCGC)/length(isCGC)]
diag_fisher[p_greater_adjust < 0.05 & input_res == 'or_dn', sum(isCGC)/length(isCGC)]
diag_fisher[p_greater_adjust < 0.05 & input_res == 'or_up', sum(isCGC)/length(isCGC)]
diag_fisher[p_greater_adjust < 0.05 & input_res == 'ac', sum(isCGC)/length(isCGC)]
diag_fisher[p_greater_adjust < 0.05 & input_res == 'fr', sum(isCGC)/length(isCGC)]

diag_fisher_absplice <- diag_fisher[p_greater_adjust < 0.05 & input_res == 'absplice', ]

diag_fisher_ac <- diag_fisher[p_greater_adjust < 0.05 & input_res == 'ac', ]




#### a. heatmap fill 0 #####
heatmap_dt <- diag_fisher[input_res == 'ac' &
                            p_greater_adjust<0.05 &
                            geneSymbol %in% CGC_cancer_gene[, GeneSymbol], ]
heatmap_dt <- heatmap_dt[, .(Diag, geneSymbol, oddsr)] 
heatmap_dt <- dcast(heatmap_dt, Diag ~ geneSymbol, value.var = 'oddsr', fill=0)

heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[, Diag]
heatmap_mtx[heatmap_mtx>1000] <- 1000

# name cohort as manuscript
rownames(heatmap_mtx) <- sapply(rownames(heatmap_mtx), function(x){
  manuscript_wording[Cohort_during_analysis==x, Cohort_abbreviation]
}) %>% unlist()

# customize order
sp <- slanter::sheatmap(heatmap_mtx, oclust_rows=FALSE, oclust_cols=FALSE)
gene_order <- sp$gtable$grobs[[4]]$label
diag_order <- sp$gtable$grobs[[5]]$label
heatmap_mtx <- heatmap_mtx[diag_order, gene_order]

# mimic oncoplot order
diag_sum <- apply(heatmap_mtx, 1, function(x){sum(x!=0)}) %>% sort(decreasing = TRUE)
diag_order <- names(diag_sum[c(1:6, 8, 12, 11, 7, 9, 10)])

gene_value <- heatmap_mtx[diag_order[1], ] %>%  sort(decreasing = TRUE)
# the for loop was manually tested
for (i in diag_order[2:8]) {
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
heatmap_dt <- diag_fisher[input_res == 'ac' &
                            p_greater_adjust<0.05 &
                            geneSymbol %in% CGC_cancer_gene[, GeneSymbol], ]
heatmap_dt <- heatmap_dt[, .(Diag, geneSymbol, oddsr)]
curation_tab <- merge(heatmap_dt, activation_enrichment_curation[, .(Diag, geneSymbol, curation, Cohort_abbreviationion)],
                      by=c('Diag', 'geneSymbol'), all.x=TRUE)
keycol <-c("Diag", "geneSymbol")
setorderv(curation_tab, keycol)
curation_tab <- subset(curation_tab, select = -Diag )
curation_tab <- curation_tab[, c(4, 1, 2, 3)] 

setnames(curation_tab, colnames(curation_tab), 
         c("DiseaseEntity",	"GeneSymbol",	"OddsRatioFromFishertest", "Curation"))
fwrite(curation_tab, snakemake@output$curation_tab)

heatmap_dt[oddsr>1000, oddsr := 1000]
heatmap_dt <- dcast(heatmap_dt, Diag ~ geneSymbol, value.var = 'oddsr', fill=1)
heatmap_dt <- melt(heatmap_dt, id.vars = 'Diag', variable.name = 'geneSymbol', value.name = 'oddsr')
heatmap_dt <- merge(heatmap_dt, manuscript_wording[, .(Cohort_during_analysis, Cohort_abbreviation)], 
                    by.x='Diag', by.y='Cohort_during_analysis')
heatmap_dt <- merge(heatmap_dt, activation_enrichment_curation[, .(Diag, geneSymbol, curation)],
                    by=c('Diag', 'geneSymbol'), all.x=TRUE)
heatmap_dt[is.na(curation), curation := 'NA']
heatmap_dt[, curation := factor(curation, levels = names(curation_colors))]
heatmap_dt[, Cohort_abbreviation := factor(Cohort_abbreviation, levels = rev(diag_order))]
heatmap_dt[, geneSymbol := factor(geneSymbol, levels = gene_order)]

# define breaks and colors
breaks <- c(10^(0:ceiling(log10(max(heatmap_mtx, na.rm = TRUE)))))
colors <- c('white', brewer.pal(length(breaks)-1, "Purples"))
legend_labels <- c(1, 10, 100, '\u22651000')

p_a_raw <- ggplot(heatmap_dt, aes(geneSymbol, Cohort_abbreviation, fill = oddsr)) +
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

# p_a




#### b. stacked bar plot #####
samp_anno_lrp1b <- samp_anno_exp[DROP_GROUP==single_group, .(ArrayID, Diag)] %>% unique()

samp_anno_lrp1b[, LRP1B_act := ArrayID %in% lrp1b_act_samples]
samp_anno_lrp1b[is.na(LRP1B_act), LRP1B_act := 'LRP1B_inact']

samp_anno_lrp1b[Diag=='HZL', Cohort:='HCL']
samp_anno_lrp1b[Diag=='vHZL', Cohort:='HCL-V']
samp_anno_lrp1b[Diag=='MZL', Cohort:='MZL']
samp_anno_lrp1b[Diag=='Multiples Myelom', Cohort:='MM']
samp_anno_lrp1b[!Diag %in% c('vHZL', 'HZL', 'MZL', 'Multiples Myelom'), Cohort:='Others']

samp_anno_lrp1b[, N_Cohort_LRP1B := .N, by=c('LRP1B_act', 'Cohort')]
samp_anno_lrp1b[, N_Cohort := .N, by='Cohort']

percent_dt <- samp_anno_lrp1b[, .(LRP1B_act, Cohort, N_Cohort_LRP1B, N_Cohort)] %>% unique()
percent_dt[, ratio := N_Cohort_LRP1B/N_Cohort]
percent_dt[, Cohort := factor(Cohort, levels = c('HCL-V', 'HCL', 'MZL', 'MM', 'Others'))]

percent_dt[, text_label := paste0( label_percent(accuracy=0.1)(ratio), "\n",
                                   "(", N_Cohort_LRP1B, "/", N_Cohort, ")")]
percent_dt <- percent_dt[LRP1B_act==TRUE, ]
percent_dt[, color_label := factor(Cohort, levels = c('Non-activated', 'HCL-V', 'HCL', 'MZL', 'MM', 'Others'))]

percent_dt_pval <- sapply(percent_dt[, Cohort], function(x){
  samp_anno_lrp1b[, fisher_col := (Cohort==x)]
  contingency_df <- samp_anno_lrp1b[, table(.(LRP1B_act, fisher_col))]
  fisher.test(contingency_df, alternative = 'greater')$p.value
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

# p_b_raw

p_b <- p_b_raw +
  scale_fill_manual(values=cohort_color, name='') +
  ylab(expression(paste(bold("Ratio of "), bolditalic("LRP1B"), bold("-activated samples"))))+
  xlab('Disease entity (dataset)') +
  theme_vale +
  theme(
    plot.margin = margin(45, 14, 14, 30, "points"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.2)
  )

# p_b




#### c. pie chart #####
composition_dt <- samp_anno_lrp1b[LRP1B_act==TRUE, .N, by='Diag']
setnames(composition_dt, 'N', 'n_sample')
composition_dt[, ratio := n_sample/sum(n_sample)]
composition_dt <- composition_dt[order(-ratio), ]

composition_dt[, Cohort := Diag]
composition_dt[Diag=='HZL', Cohort:='HCL']
composition_dt[Diag=='vHZL', Cohort:='HCL-V']
composition_dt[Diag=='Multiples Myelom', Cohort:='Others']
composition_dt[, Cohort := factor(Cohort, levels = c('HCL-V', 'HCL', 'MZL', 'Others'))]
composition_dt[, text_label := paste0(Cohort, " (n=", n_sample, ")")]

p_c_raw <- ggplot(composition_dt, aes(x="", y=ratio, fill=Cohort)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=cohort_color, labels=composition_dt[, text_label]) 

p_c <- p_c_raw +
  theme_vale +
  ggtitle(expression(paste(bolditalic("LRP1B"), bold("-activated samples from dataset (n=24)"))))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        axis.text.x=element_blank(),
        legend.position="right",
        legend.direction="vertical"
  ) +
  guides(fill=guide_legend("Disease entity"))

# p_c




#### d. dot/box plot mll #####
counts_mll <- data.table(sampleID = colnames(ods_filter_out),
                         raw_count = OUTRIDER::counts(ods_filter_out, normalized = FALSE)[geneID_LRP1B,],
                         size_factor = OUTRIDER::sizeFactors(ods_filter_out))
counts_mll <- merge(counts_mll, 
                    samp_anno_lrp1b[, .(ArrayID, Diag, Cohort)] %>% unique(), 
                    by.x='sampleID', by.y='ArrayID', all.x=TRUE, all.y=FALSE)
counts_mll[, norm_count := (raw_count + 1)/size_factor]
counts_mll[, Cohort := factor(Cohort, levels = c('HCL-V', 'HCL', 'MZL', 'Others'))]

#' norm count median
counts_mll[Cohort == 'HCL-V', median(norm_count)]
counts_mll[Cohort == 'HCL', median(norm_count)]
counts_mll[Cohort == 'MZL', median(norm_count)]
counts_mll[Cohort == 'Others', median(norm_count)]

stat.test <- ggpubr::compare_means(formula = norm_count ~ Cohort, data = counts_mll, method = "wilcox.test")
# stat.test

p_d_raw <- ggplot(counts_mll, aes(x=Cohort, y=norm_count)) +
  geom_point(aes(colour=Cohort, fill=Cohort), alpha=0.3) + 
  geom_boxplot(fill=NA, width=0.25)

p_d_raw <- ggplot(counts_mll, aes(x=Cohort, y=norm_count)) +
  geom_violin(aes(colour=Cohort, fill=Cohort), alpha=0.3) + 
  geom_boxplot(fill=NA, width=0.25)+
  stat_pvalue_manual(stat.test[c(1,4,6),], label="p.signif", y.position = 5, bracket.shorten = 0.2)+
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
  ) 

# p_d_raw

p_d <- p_d_raw +
  xlab('Dataset') +
  ylab(expression(paste(bolditalic("LRP1B"), bold(" normalized counts"))))+
  scale_color_manual(values=cohort_color, guide = 'none') +
  scale_fill_manual(values=cohort_color) +
  theme_vale +
  guides(fill=guide_legend('Disease entity', override.aes=list(color=NA)))

# p_d




#### fit gaussian mixture #####
counts_marc <- data.table(sampleID = colnames(total_counts_marc),
                          raw_count = assay(total_counts_marc)[geneID_marc,])
counts_marc <- merge(counts_marc, size_factor_marc[, .(sampleID, `size_factor (deseq)`, sample)])
setnames(counts_marc, 'size_factor (deseq)', 'size_factor')
setnames(counts_marc, 'sample', 'Diag')
counts_marc[, size_factor := sub(',', '.', size_factor) %>% as.numeric()]
counts_marc[, norm_count := (raw_count + 1)/size_factor]
counts_marc[Diag=='cl_HCL', Cohort:='HCL']
counts_marc[Diag=='HCL_var', Cohort:='HCL-V']
counts_marc[, sampleID := NULL]
fwrite(counts_marc, snakemake@output$counts_marc)

counts_marc_wilcox <- copy(counts_marc)
counts_marc[, Cohort := factor(Cohort, levels = c('HCL-V', 'HCL', 'MZL', 'Others'))]

# ggplot(counts_marc, aes(sample = log(norm_count))) +
#   geom_qq() + facet_wrap(~Cohort)

# x <- counts_marc[, norm_count]
x <- counts_marc[, log(norm_count)]

# fit <- Mclust(x, G=2, model="V") # Gaussian mixtures with unequal variance
fit <- Mclust(x, G=2, model="E") # Gaussian mixtures with equal variance
summary(fit)
predicted_class <- predict(fit)$classification

counts_marc[predicted_class == 1, fitted_activation := FALSE]
counts_marc[predicted_class == 2, fitted_activation := TRUE]
counts_marc[, sum(fitted_activation)/.N, by='Cohort']
counts_marc[, sum(fitted_activation), by='Cohort']
counts_marc[, .N, by='Cohort']

ggplot(counts_marc, aes(x=Cohort, y=norm_count, color=fitted_activation)) +
  geom_point() +
  scale_y_log10()




#### e. dot/box plot marc #####
counts_marc[, sum(fitted_activation)/.N, by='Cohort'][, V1] %>% as.character()

ratio_vec <- counts_marc[, paste0( label_percent(accuracy=0.1)( sum(fitted_activation)/.N), "\n",
                                   "(", sum(fitted_activation), "/", .N, ")"), by='Cohort'][, V1]

p_e_raw <- ggplot(counts_marc, aes(x=Cohort, y=norm_count + 1)) +
  geom_jitter(aes(colour=fitted_activation), width=0.1) + 
  scale_y_log10(
    breaks = c(1, 10, 100, 1000), 
    minor_breaks = rep(1:9, 3)*(10^rep(0:3, each = 9))
  ) +
  annotate("text", x=c('HCL', 'HCL-V'), y=3000, label=ratio_vec, size = 3) +
  annotate("text", x=0.2, y=3000, label='Percentage\n activated', size = 3) +
  coord_cartesian(ylim=c(1,1000), xlim=c(1, 2), clip="off")

# p_e_raw

p_e <- p_e_raw +
  xlab('Disease entity (validation dataset)') +
  ylab(expression(paste(bolditalic("LRP1B"), 
                        bold(" normalized counts + 1"))))+
  scale_color_manual(values=c('darkgrey', 'firebrick')) +
  theme_vale +
  guides(color=guide_legend(expression(atop(bolditalic("LRP1B"),
                                             bold("activated"))), 
                            reverse = TRUE) ) + 
  
  theme(
    plot.margin = margin(45, 14, 14, 30, "points"),
    panel.border = element_rect(colour = "black", 
                                fill=NA, linewidth=0.2)
  )

# p_e




# #### f. coverage #####
# hg38: 
# exon12-14: chr2:141013556-141020102
# exon12-14+-500: chr2:141013056-141020602
# hg19: 
# exon12-14: chr2:141771125-141777671
# exon12-14+-500: chr2:141770625-141778171

# # region is LRP1B exon 11 - 16
# gr_chr <- "chr2"
# gr_start <- 141750000
# gr_end <- 141810000
# gr_strand <- "-"
# exon_range <- c(11:16)
# 
# # # region is LRP1B exon 12 - 14
# # gr_chr <- "chr2"
# # gr_start <- 141762066
# # gr_end <- 141779012
# # gr_strand <- "-"
# # exon_range <- c(12:15)
# 
# # define region
# gr_region <- GRanges(seqnames = gr_chr,
#                      ranges = IRanges(gr_start, width = gr_end-gr_start+1),
#                      strand = gr_strand)
# 
# # get exons
# exonic_regions <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
# exon_region <- exonic_regions[seqnames(exonic_regions) == gr_chr & 
#                                 start(exonic_regions) > gr_start & 
#                                 end(exonic_regions) < gr_end & 
#                                 strand(exonic_regions) == gr_strand]
# exon_pos <- sapply(seq_along(exon_region), function(i){
#   c(start(exon_region[i]):end(exon_region[i]))
# }) %>% unlist()
# 
# # get transcript 
# region_dt <- exon_region %>% as.data.table() 
# region_dt[, transcript_name:='LRP1B']
# 
# region_dt <- shorten_gaps(
#   region_dt, 
#   to_intron(region_dt, "transcript_name"), 
#   group_var = "transcript_name",
#   target_gap_width = 150L
# ) %>% as.data.table()
# region_dt[, start_mirror:=-start]
# region_dt[, end_mirror:=-end]
# region_dt[, number_pos := (start_mirror+end_mirror)/2, by=list(rownames(region_dt))]
# 
# exon_ref_dt <- region_dt[type=='exon',]
# exon_ref_dt[, transcript_name := 'Reference']
# exon_ref_dt[, exon_number := rev(exon_range)]
# exon_ref_dt <- exon_ref_dt[exon_id %in% c(43120:43123), ]
# intron_ref_dt <- region_dt[type=='intron',]
# intron_ref_dt[, transcript_name := 'Reference']
# 
# exon_asm_dt <- exon_ref_dt[exon_id!='43123', ]
# exon_asm_dt[, transcript_name := 'Assembled']
# intron_asm_dt <- intron_ref_dt[exon_id %in% c(43120:43122), ]
# intron_asm_dt[, transcript_name := 'Assembled']
# 
# exon_region_dt <- rbind(exon_asm_dt, exon_ref_dt)
# exon_region_dt[, transcript_name := factor(transcript_name, levels=c('Reference', 'Assembled'))]
# intron_region_dt <- rbind(intron_asm_dt, intron_ref_dt)
# intron_region_dt[, transcript_name := factor(transcript_name, levels=c('Reference', 'Assembled'))]
# 
# # plot transcript
# p_transcript_raw <- ggplot(exon_region_dt, 
#                            aes(xstart=start_mirror, xend=end_mirror, y=transcript_name)) +
#   geom_intron(data = intron_region_dt, aes(strand=strand), color='blue') +
#   geom_range(fill='blue') +
#   geom_text(aes(x=number_pos, label=exon_number), color='white')
# 
# # p_transcript_raw
# 
# p_transcript <- p_transcript_raw + 
#   theme(
#     axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     axis.title.y=element_blank(),
#     axis.text.y=element_text(hjust = 0),
#     axis.ticks.y=element_blank(),
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank()
#   )
# 
# # p_transcript
# 
# 
# # get coverage and subset for exon
# p_coverage <- autoplot("/s/project/mll/rawdata/alignments/MLL_57538-M005.alignments.bam", which = gr_region)
# coverage_dt <- ggplot_build(p_coverage@ggplot)$data[[1]] %>% as.data.table()
# coverage_region <- coverage_dt[x %in% exon_pos, ]
# setnames(coverage_region, c('x', 'y'), c('pos_ori', 'coverage'))
# coverage_region <- coverage_region[, .(pos_ori, coverage)]
# 
# # get truncated gap of each pos
# transcript_region <- copy(exon_region_dt[transcript_name=='Reference', ])
# exon_region <- exon_region[2:5,]
# transcript_region[, start_ori := start(exon_region)]
# transcript_region[, end_ori := end(exon_region)]
# transcript_region[, truncated_gaps := start_ori - start]
# transcript_region <- transcript_region[, .(start_ori, end_ori, start, truncated_gaps)]
# 
# transcript_region <- do.call("rbind", replicate(length(exon_pos), transcript_region, simplify = FALSE))
# transcript_region[, pos_ori:=rep(exon_pos, each=length(exon_region))]
# subset_vec <- transcript_region[, pos_ori %in% c(start_ori:end_ori), by=rownames(transcript_region)]$V1
# transcript_region <- transcript_region[subset_vec, ]
# 
# coverage_region <- merge(coverage_region, transcript_region, by='pos_ori')
# coverage_region[, pos := pos_ori-truncated_gaps]
# coverage_region[, pos_mirror := -pos]
# 
# # plot coverage
# p_coverage_raw <- ggplot(coverage_region, aes(x=pos_mirror,y=coverage))+
#   geom_bar(stat = "identity")
# 
# # p_coverage_raw
# 
# p_coverage <- p_coverage_raw + 
#   theme(
#     axis.title.x=element_blank(),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     axis.title.y=element_blank(),
#     axis.text.y=element_blank(),
#     axis.ticks.y=element_blank(),
#     panel.grid.major = element_blank(), 
#     panel.grid.minor = element_blank(),
#     panel.background = element_blank(),
#     plot.margin = margin(0,1,0,2.5, "cm")
#   )
# 
# # p_coverage
# 
# p_e_raw <- ggarrange(p_coverage_raw, p_transcript_raw, ncol=1)
# # p_e_raw
# 
# p_e <- ggarrange(p_coverage, p_transcript, ncol=1)
# # p_e




### survival ####
lrp1b_survival_hcl <- read.xlsx(snakemake@params$lrp1b_survival_hcl_mzl, sheet = 'HZL') %>% as.data.table()
colnames(lrp1b_survival_hcl) <- gsub("[.]", "_", colnames(lrp1b_survival_hcl))

lrp1b_survival_mzl <- read.xlsx(snakemake@params$lrp1b_survival_hcl_mzl, sheet = 'MZL') %>% as.data.table()
colnames(lrp1b_survival_mzl) <- gsub("[.]", "_", colnames(lrp1b_survival_mzl))

lrp1b_survival_hclv <- read.xlsx(snakemake@params$lrp1b_survival_hclv)  %>% as.data.table()
colnames(lrp1b_survival_hclv) <- gsub("[.]", "_", colnames(lrp1b_survival_hclv))

lrp1b_survival <- rbind(lrp1b_survival_hclv, lrp1b_survival_hcl, lrp1b_survival_mzl, fill=TRUE)

counts_survival <- merge(counts_mll, lrp1b_survival, by.x='sampleID', by.y='Array_ID')
counts_survival[, LRP1B_high := sampleID %in% lrp1b_act_samples]
counts_other <- counts_mll[!sampleID %in% counts_survival[, sampleID], ]
counts_other[, LRP1B_high := sampleID %in% lrp1b_act_samples]

hist(counts_other[, norm_count], 100)

cor.test(counts_survival[, raw_count], counts_survival[, OS_days], method = 'spearman')
cor.test(counts_survival[, norm_count], counts_survival[, OS_days], method = 'spearman')

ggplot(counts_survival, aes(raw_count, OS_days, color=Cohort)) +
  geom_point() +
  xlim(0, NA) +
  ylim(0, NA) +
  xlab('LRP1B raw count') +
  ylab('Overall survival (days)') +
  stat_cor(method = "spearman") +
  theme_vale

ggplot(counts_survival, aes(norm_count, OS_days)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap('Cohort') +
  xlim(0, NA) +
  ylim(0, NA) +
  xlab('LRP1B size-factor normalized count') + 
  ylab('Overall survival (days)') + 
  stat_cor(method = "spearman") +
  theme_vale

p_s7_raw <- ggplot(counts_survival, aes(norm_count, OS_days, color=LRP1B_high)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlim(0, NA) +
  ylim(0, NA) 

p_s7 <- p_s7_raw + 
  guides(color=guide_legend("LRP1B activated")) +
  xlab('LRP1B size-factor normalized count') + 
  ylab('Overall survival (days)') + 
  stat_cor(method = "spearman") +
  theme_vale

p_s8_raw <- ggplot(counts_survival[LRP1B_high==TRUE, ], 
       aes(norm_count, OS_days)) +
  facet_wrap('Cohort') +
  geom_point() +
  geom_smooth(method = 'lm') +
  xlim(0, NA) +
  ylim(0, 20000) 

p_s8 <- p_s8_raw + 
  xlab('LRP1B size-factor normalized count') + 
  ylab('Overall survival (days)') + 
  stat_cor(method = "spearman") +
  theme_vale




### arrange plot ####
#' ### raw
#+ plot p5 raw, fig.width=10, fig.height=7.5
p5_raw_bottom <- ggarrange(p_b_raw, p_e_raw, 
                           ncol = 2, nrow = 1, labels = c("B", "C"))
p5_raw <- ggarrange(p_a_raw, p5_raw_bottom,
                    ncol = 1, nrow = 2, labels=c("A", ""), heights = c(3.5, 4))

p5_raw

#' ### annotated
#+ plot p5, fig.width=10, fig.height=7.5
p5_bottom <- ggarrange(p_b, p_e, 
                       ncol = 2, nrow = 1, labels = c("B", "C"))
p5 <- ggarrange(p_a, p5_bottom,
                ncol = 1, nrow = 2, labels=c("A", ""), heights = c(3.5, 4))

p5


pdf(paste0(output_dir, "/figure_5.pdf"), width = 10, height = 7.5)
p5_bottom <- ggarrange(p_b, p_e, 
                       ncol = 2, nrow = 1, labels = c("B", "C"))
p5 <- ggarrange(p_a, p5_bottom,
                ncol = 1, nrow = 2, labels=c("A", ""), heights = c(3.5, 4))

p5
dev.off()


png(paste0(output_dir, "/figure_5.png"), width = 10, height = 7.5, units = "in", res = 600)
p5_bottom <- ggarrange(p_b, p_e, 
                       ncol = 2, nrow = 1, labels = c("B", "C"), widths = c(1, 1.2))
p5 <- ggarrange(p_a, p5_bottom,
                ncol = 1, nrow = 2, labels=c("A", ""), heights = c(3.5, 4))

p5
dev.off()




#' ## Sup
### Supplement #####
#### s1. survival curve #####
# Read combined KM data and format
lrp1b_survival_hclv <- read.xlsx(snakemake@params$lrp1b_survival_hclv) 
# udpate lrp1b high def
lrp1b_survival_hclv$LRP1B_high <- as.numeric(lrp1b_survival_hclv$Array.ID %in% lrp1b_act_samples) 

lrp1b_survival_hclv <- lrp1b_survival_hclv %>%
  mutate(OS_years = OS_days/365,
         cohort = factor(LRP1B_high, levels = c("0", "1")))

# Fit survival models
survival_fit <- survfit(Surv(OS_years, Zensor_OS) ~ cohort, data = lrp1b_survival_hclv)

# Plot  
#' ### s1
p_s1 <- ggsurvplot(survival_fit,
                  pval = TRUE,
                  conf.int = FALSE,
                  censor= TRUE,
                  censor.shape = "I",
                  title = NULL,
                  xlab = "Time (years)",
                  ylab = "Survival probability",
                  legend = "right",
                  legend.title = "Categories",
                  legend.labs = c("LRP1B non-activated", "LRP1B activated"),
                  palette = c("grey", "red"),
                  axes.offset = FALSE, 
                  break.x.by = 1,
                  break.y.by = 0.1,
                  xlim = c(0, 10.5),
                  ylim = c(0, 1),
                  risk.table = TRUE,
                  surv.median.line = "h"
)

#+ plot s1, fig.width=8, fig.height=6
p_s1




#### s2. ac enrichment heatmap #####
#' ### s2 ac
heatmap_dt <- diag_fisher[input_res == 'ac', ]

heatmap_dt <- heatmap_dt[, .(Diag, geneSymbol, oddsr)] 
heatmap_dt <- dcast(heatmap_dt, geneSymbol ~ Diag, value.var = 'oddsr', fill=0)
heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[,geneSymbol]
heatmap_mtx[heatmap_mtx==Inf] <- max(heatmap_mtx[heatmap_mtx!=Inf])

# name cohort as manuscript
colnames(heatmap_mtx) <- sapply(colnames(heatmap_mtx), function(x){
  manuscript_wording[Cohort_during_analysis==x, Cohort_abbreviation]
}) %>% unlist()

# define breaks and colors
breaks <- c(0, 10^(0:ceiling(log10(max(heatmap_mtx, na.rm = TRUE)))))
colors <- c('white', brewer.pal(length(breaks)-1, "Purples"))

#+ plot s2, fig.width=8, fig.height=10
p_s2 <- pheatmap(heatmap_mtx,
                 color=colors, breaks=breaks, legend_breaks=breaks, legend=T,
                 # cluster_rows=F, cluster_cols=F, na_col = "grey",
                 angle_col = 90
)
grid.text('Enrichment of activated outliers', x = 0.5, y = 0.99)




#### s3. or_up enrichment heatmap #####
#' ### s3 or_up
heatmap_dt <- diag_fisher[input_res == 'or_up', ]

heatmap_dt <- heatmap_dt[, .(Diag, geneSymbol, oddsr)] 
heatmap_dt <- dcast(heatmap_dt, geneSymbol ~ Diag, value.var = 'oddsr', fill=0)
heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[,geneSymbol]
heatmap_mtx[heatmap_mtx==Inf] <- max(heatmap_mtx[heatmap_mtx!=Inf])

# name cohort as manuscript
colnames(heatmap_mtx) <- sapply(colnames(heatmap_mtx), function(x){
  manuscript_wording[Cohort_during_analysis==x, Cohort_abbreviation]
}) %>% unlist()

# define breaks and colors
breaks <- c(0, 10^(0:ceiling(log10(max(heatmap_mtx, na.rm = TRUE)))))
colors <- c('white', brewer.pal(length(breaks)-1, "Purples"))

#+ plot s3, fig.width=8, fig.height=10
p_s3 <- pheatmap(heatmap_mtx,
                 color=colors, breaks=breaks, legend_breaks=breaks, legend=T,
                 # cluster_rows=F, cluster_cols=F, na_col = "grey",
                 angle_col = 90
)$gtable
grid.text('Enrichment of over-expressed outliers', x = 0.5, y = 0.99)




#### s4. or_dn enrichment heatmap #####
#' ### s4 or_dn
heatmap_dt <- diag_fisher[input_res == 'or_dn', ]

heatmap_dt <- heatmap_dt[, .(Diag, geneSymbol, oddsr)] 
heatmap_dt <- dcast(heatmap_dt, geneSymbol ~ Diag, value.var = 'oddsr', fill=0)
heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[,geneSymbol]
heatmap_mtx[heatmap_mtx==Inf] <- max(heatmap_mtx[heatmap_mtx!=Inf])

# name cohort as manuscript
colnames(heatmap_mtx) <- sapply(colnames(heatmap_mtx), function(x){
  manuscript_wording[Cohort_during_analysis==x, Cohort_abbreviation]
}) %>% unlist()

# define breaks and colors
breaks <- c(0, 10^(0:ceiling(log10(max(heatmap_mtx, na.rm = TRUE)))))
colors <- c('white', brewer.pal(length(breaks)-1, "Purples"))

#+ plot s4, fig.width=8, fig.height=10
p_s4 <- pheatmap(heatmap_mtx,
                 color=colors, breaks=breaks, legend_breaks=breaks, legend=T,
                 # cluster_rows=F, cluster_cols=F, na_col = "grey",
                 angle_col = 90
)$gtable
grid.text('Enrichment of under-expressed outliers', x = 0.5, y = 0.99)




#### s5. fr enrichment heatmap #####
#' ### s5 fr
heatmap_dt <- diag_fisher[input_res == 'fr', ]

heatmap_dt <- heatmap_dt[, .(Diag, geneSymbol, oddsr)] 
heatmap_dt <- dcast(heatmap_dt, geneSymbol ~ Diag, value.var = 'oddsr', fill=0)
heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[,geneSymbol]
heatmap_mtx[heatmap_mtx==Inf] <- max(heatmap_mtx[heatmap_mtx!=Inf])

# name cohort as manuscript
colnames(heatmap_mtx) <- sapply(colnames(heatmap_mtx), function(x){
  manuscript_wording[Cohort_during_analysis==x, Cohort_abbreviation]
}) %>% unlist()

# define breaks and colors
breaks <- c(0, 10^(0:ceiling(log10(max(heatmap_mtx, na.rm = TRUE)))))
colors <- c('white', brewer.pal(length(breaks)-1, "Purples"))

#+ plot s5, fig.width=8, fig.height=10
p_s5 <- pheatmap(heatmap_mtx,
                 color=colors, breaks=breaks, legend_breaks=breaks, legend=T,
                 # cluster_rows=F, cluster_cols=F, na_col = "grey",
                 angle_col = 90
)$gtable
grid.text('Enrichment of splicing outliers', x = 0.5, y = 0.99)




#### s6. absplice enrichment heatmap #####
#' ### s6 absplice
heatmap_dt <- diag_fisher[input_res == 'absplice', ]

heatmap_dt <- heatmap_dt[, .(Diag, geneSymbol, oddsr)] 
heatmap_dt <- dcast(heatmap_dt, geneSymbol ~ Diag, value.var = 'oddsr', fill=0)
heatmap_mtx <- as.matrix(heatmap_dt[,-1])
rownames(heatmap_mtx) <- heatmap_dt[,geneSymbol]
heatmap_mtx[heatmap_mtx==Inf] <- max(heatmap_mtx[heatmap_mtx!=Inf])

# name cohort as manuscript
colnames(heatmap_mtx) <- sapply(colnames(heatmap_mtx), function(x){
  manuscript_wording[Cohort_during_analysis==x, Cohort_abbreviation]
}) %>% unlist()

# define breaks and colors
breaks <- c(0, 10^(0:ceiling(log10(max(heatmap_mtx, na.rm = TRUE)))))
colors <- c('white', brewer.pal(length(breaks)-1, "Purples"))

#+ plot s6, fig.width=8, fig.height=10
p_s6 <- pheatmap(heatmap_mtx,
                 color=colors, breaks=breaks, legend_breaks=breaks, legend=T,
                 # cluster_rows=F, cluster_cols=F, na_col = "grey",
                 angle_col = 90
)$gtable
grid.text('Enrichment of splicing variants', x = 0.5, y = 0.99)

#' ### s7 raw
#+ plot s7 raw, fig.width=9, fig.height=6
p_s7_raw
#' ### s7 annotated
#+ plot s7, fig.width=9, fig.height=6
p_s7

#' ### s8 raw
#+ plot s8 raw, fig.width=9, fig.height=6
p_s8_raw
#' ### s8 annotated
#+ plot s8, fig.width=9, fig.height=6
p_s8




### thesis ####
png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "lrp1b_heatmap.png"), 
    width = 8, height = 3, units = "in", res = 600)
p_a
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "lrp1b_mll.png"), 
    width = 4, height = 4, units = "in", res = 600)
p_b
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "lrp1b_marc.png"), 
    width = 4, height = 4, units = "in", res = 600)
p_e
dev.off() 

png(paste0("/s/project/vale/driver_prediction_202304/Output/html/thesis/", 
           "lrp1b_survival_hclv.png"), 
    width = 8, height = 6, units = "in", res = 600)
p_s1
dev.off() 