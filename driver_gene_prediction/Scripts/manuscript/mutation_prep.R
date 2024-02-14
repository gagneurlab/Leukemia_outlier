#'---
#' title: mutation_prep
#' author: Xueqi Cao
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
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - predictedConsequence: '`sm config["predictedConsequence"]`'
#'    - AML_variants: '`sm config["AML_variants"]`'
#'    - AML_enformer: '`sm config["AML_enformer"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'  input:
#'    - mll_cnv: '`sm expand(config["projectPath"] + "/manuscript/cnv/{dataset}.tsv",
#'                   dataset=outputDatasets)`'
#'    - mll_cnv_full: '`sm expand(config["projectPath"] + "/manuscript/cnv_full/{dataset}.tsv",
#'                    dataset=outputDatasets)`'
#'    - ods: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds",
#'             annotation=annotations, dataset=inputDatasets)`'
#'    - ods_filter_out: '`sm expand(config["outriderDir"] +"/processed_results/aberrant_expression/{annotation}/outrider/{dataset}/ods_filter_out.Rds",
#'                        annotation=annotations, dataset=inputDatasets)`'
#'    - fraserRes: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv",
#'                    annotation=annotations, dataset=outputDatasets)`'
#'    - fraserResJunc: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv",
#'                    annotation=annotations, dataset=outputDatasets)`'
#'    - fds: '`sm expand(config["fraserDir"] +
#'            "/processed_results/aberrant_splicing/datasets/savedObjects/{dataset}--{annotation}/fds-object.RDS",
#'             annotation=annotations, dataset=outputDatasets)`'
#'    - abspliceRes: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice.tsv"`'
#'    - abspliceResVar: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice_var.tsv"`'
#'  output:
#'    - mutation_gene_sample: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/mutation_gene_sample.csv"`'
#'    - mutation_fr: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/mutation_fr.csv"`'
#'    - prop_mutation_or: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/prop_mutation_or.csv"`'
#'    - prop_mutation_fr: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_2/plot_data/anova/prop_mutation_fr.csv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/mutation_prep.html"`'
#'  type: noindex
#'  threads: 1  
#'  resources:
#'    - mem_mb: 256000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/mutation_prep.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202402/processed_data/snakemake/mutation_prep.snakemake")

.libPaths("~/R/4.1/FRASER2")
suppressPackageStartupMessages({
  library(OUTRIDER)
  library(FRASER)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id_unique := gene_id]
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]
gencode <- gencode[gene_type == 'protein_coding', ]

samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()
samp_anno_14groups <- samp_anno_exp[DROP_GROUP %in% snakemake@params$inputDatasets, ] %>% unique()
single_group_id <- samp_anno_14groups[, ArrayID]

mutation_dt <- fread('/s/project/mll/sergey/effect_prediction/promoter_mutations/res/mutations_accumulated.tsv.gz')

setnames(mutation_dt, c('geneName', 'array_id'), c('gene_id', 'sampleID'))
mutation_dt <- merge(mutation_dt, gencode[, .(gene_id, gene_name_orig)] %>% unique(),
                     by='gene_id', all.x=TRUE, all.y=FALSE)
# mutation_dt[, aberrant := padjust < 0.05]
mutation_dt[, samp_symbol := paste0(sampleID, "-", gene_name_orig)]
# mutation_dt[, samp_ensgid := paste0(sampleID, "-", gene_id)]
  
fwrite(mutation_dt, snakemake@output$mutation_gene_sample)

# cnv/all mutation among outliers
mutation_aberrant <- mutation_dt[padjust < 0.05, ]
mutation_aberrant[, mutation_type_count := sum(promoter_variant, frameshift_variant, stop_gained, splice_related_variant, 
                                               structural_variant, copy_number_gain, copy_number_loss), 
                  by=rownames(mutation_aberrant)]
mutation_aberrant[, mutation_exist := mutation_type_count >=1]
mutation_aberrant[, cnv_exist := (copy_number_gain + copy_number_loss)>=1]
mutation_aberrant[, sum(cnv_exist)/sum(mutation_exist)]




### or ac ####
mutation_or_dn <- mutation_dt[tool=='outrider' & l2fc < 0 & padjust < 0.05, ]
mutation_or_dn[, mutation_type_count := sum(promoter_variant, frameshift_variant, stop_gained, splice_related_variant, 
                                            structural_variant, copy_number_gain, copy_number_loss), 
               by=rownames(mutation_or_dn)]
mutation_or_dn[, table(mutation_type_count)]
mutation_or_dn[, multiple_type := mutation_type_count >=2]
mutation_or_dn[multiple_type == TRUE, promoter_variant := FALSE]
mutation_or_dn[multiple_type == TRUE, frameshift_variant := FALSE]
mutation_or_dn[multiple_type == TRUE, stop_gained := FALSE]
mutation_or_dn[multiple_type == TRUE, splice_related_variant := FALSE]
mutation_or_dn[multiple_type == TRUE, structural_variant := FALSE]
mutation_or_dn[multiple_type == TRUE, copy_number_gain := FALSE]
mutation_or_dn[multiple_type == TRUE, copy_number_loss := FALSE]

mutation_or_dn_melt <- melt(mutation_or_dn, id.vars=c('gene_id', 'sampleID', 'padjust', 'l2fc', 'tool'),
                            measure.vars = c('promoter_variant', 'frameshift_variant', 'stop_gained', 'splice_related_variant', 
                                             'structural_variant', 'copy_number_gain', 'copy_number_loss', 'multiple_type'),
                            variable.name = 'mutation_type', value.name='mutation_exist')
mutation_or_dn_melt <- mutation_or_dn_melt[mutation_exist == TRUE,]
or_dn_overlap <- mutation_or_dn_melt[, .N, by='mutation_type']
setnames(or_dn_overlap, "N", "count")
or_dn_overlap <- rbind(or_dn_overlap, 
                       data.table(mutation_type = 'none', 
                                  count=(nrow(mutation_or_dn)-or_dn_overlap[, sum(count)]))
)
or_dn_overlap[, Category := 'Underexpression outlier']
or_dn_overlap[, total := sum(count)]

mutation_or_up <- mutation_dt[tool=='outrider' & l2fc > 0 & padjust < 0.05, ]
mutation_or_up[, mutation_type_count := sum(promoter_variant, frameshift_variant, stop_gained, splice_related_variant, 
                                            structural_variant, copy_number_gain, copy_number_loss), 
               by=rownames(mutation_or_up)]
mutation_or_up[, table(mutation_type_count)]
mutation_or_up[, multiple_type := mutation_type_count >=2]
mutation_or_up[multiple_type == TRUE, promoter_variant := FALSE]
mutation_or_up[multiple_type == TRUE, frameshift_variant := FALSE]
mutation_or_up[multiple_type == TRUE, stop_gained := FALSE]
mutation_or_up[multiple_type == TRUE, splice_related_variant := FALSE]
mutation_or_up[multiple_type == TRUE, structural_variant := FALSE]
mutation_or_up[multiple_type == TRUE, copy_number_gain := FALSE]
mutation_or_up[multiple_type == TRUE, copy_number_loss := FALSE]

mutation_or_up_melt <- melt(mutation_or_up, id.vars=c('gene_id', 'sampleID', 'padjust', 'l2fc', 'tool'),
                            measure.vars = c('promoter_variant', 'frameshift_variant', 'stop_gained', 'splice_related_variant', 
                                             'structural_variant', 'copy_number_gain', 'copy_number_loss', 'multiple_type'),
                            variable.name = 'mutation_type', value.name='mutation_exist')
mutation_or_up_melt <- mutation_or_up_melt[mutation_exist == TRUE,]
or_up_overlap <- mutation_or_up_melt[, .N, by='mutation_type']
setnames(or_up_overlap, "N", "count")
or_up_overlap <- rbind(or_up_overlap, 
                       data.table(mutation_type = 'none', 
                                  count=(nrow(mutation_or_up)-or_up_overlap[, sum(count)]))
)
or_up_overlap[, Category := 'Overexpression outlier']
or_up_overlap[, total := sum(count)]

mutation_ac <- mutation_dt[tool=='activation' & l2fc > 0 & padjust < 0.05, ]
mutation_ac[, mutation_type_count := sum(promoter_variant, frameshift_variant, stop_gained, splice_related_variant, 
                                         structural_variant, copy_number_gain, copy_number_loss), 
            by=rownames(mutation_ac)]
mutation_ac[, table(mutation_type_count)]
mutation_ac[, multiple_type := mutation_type_count >=2]
mutation_ac[multiple_type == TRUE, promoter_variant := FALSE]
mutation_ac[multiple_type == TRUE, frameshift_variant := FALSE]
mutation_ac[multiple_type == TRUE, stop_gained := FALSE]
mutation_ac[multiple_type == TRUE, splice_related_variant := FALSE]
mutation_ac[multiple_type == TRUE, structural_variant := FALSE]
mutation_ac[multiple_type == TRUE, copy_number_gain := FALSE]
mutation_ac[multiple_type == TRUE, copy_number_loss := FALSE]

mutation_ac_melt <- melt(mutation_ac, id.vars=c('gene_id', 'sampleID', 'padjust', 'l2fc', 'tool'),
                         measure.vars = c('promoter_variant', 'frameshift_variant', 'stop_gained', 'splice_related_variant', 
                                          'structural_variant', 'copy_number_gain', 'copy_number_loss', 'multiple_type'),
                         variable.name = 'mutation_type', value.name='mutation_exist')
mutation_ac_melt <- mutation_ac_melt[mutation_exist == TRUE,]
ac_overlap <- mutation_ac_melt[, .N, by='mutation_type']
setnames(ac_overlap, "N", "count")
ac_overlap <- rbind(ac_overlap, 
                    data.table(mutation_type = 'none', 
                               count=(nrow(mutation_ac)-ac_overlap[, sum(count)]))
)
ac_overlap[, Category := 'Activation outlier']
ac_overlap[, total := sum(count)]

mutation_none <- mutation_dt[padjust > 0.05, ]
mutation_none[, mutation_type_count := sum(promoter_variant, frameshift_variant, stop_gained, splice_related_variant, 
                                           structural_variant, copy_number_gain, copy_number_loss), 
              by=rownames(mutation_none)]
mutation_none[, table(mutation_type_count)]
mutation_none[, multiple_type := mutation_type_count >=2]
mutation_none[multiple_type == TRUE, promoter_variant := FALSE]
mutation_none[multiple_type == TRUE, frameshift_variant := FALSE]
mutation_none[multiple_type == TRUE, stop_gained := FALSE]
mutation_none[multiple_type == TRUE, splice_related_variant := FALSE]
mutation_none[multiple_type == TRUE, structural_variant := FALSE]
mutation_none[multiple_type == TRUE, copy_number_gain := FALSE]
mutation_none[multiple_type == TRUE, copy_number_loss := FALSE]

mutation_none_melt <- melt(mutation_none, id.vars=c('gene_id', 'sampleID', 'padjust', 'l2fc', 'tool'),
                           measure.vars = c('promoter_variant', 'frameshift_variant', 'stop_gained', 'splice_related_variant', 
                                            'structural_variant', 'copy_number_gain', 'copy_number_loss', 'multiple_type'),
                           variable.name = 'mutation_type', value.name='mutation_exist')
mutation_none_melt <- mutation_none_melt[mutation_exist == TRUE,]
none_overlap <- mutation_none_melt[, .N, by='mutation_type']
setnames(none_overlap, "N", "count")
none_overlap <- rbind(none_overlap, 
                      data.table(mutation_type = 'none', 
                                 count=(nrow(mutation_none)-none_overlap[, sum(count)]))
)
none_overlap[, Category := 'None outlier']
none_overlap[, total := sum(count)]

prop_mutation_or <- rbind(or_dn_overlap, or_up_overlap, ac_overlap, none_overlap)
prop_mutation_or[, Percentage := count/total]

fwrite(prop_mutation_or, snakemake@output$prop_mutation_or)



### fr ####
absplice_res <- fread(snakemake@input$abspliceRes)
absplice_res[, samp_symbol := paste0(sampleID, "-", gene_name)]
absplice_res <- separate(absplice_res, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_res <- separate(absplice_res, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()
absplice_res_sig <- absplice_res[AbSplice_DNA >= 0.2, ]

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

exp_gene_fr <- c()
for (fds_path in snakemake@input$fds) {
  print(fds_path)
  fds_temp <- loadFraserDataSet(file = fds_path)
  exp_gene_fr_temp <- rowData(fds_temp)$hgnc_symbol
  exp_gene_fr <- c(exp_gene_fr, exp_gene_fr_temp)
}
exp_gene_fr <- exp_gene_fr[!is.na(exp_gene_fr)] %>% unique()
exp_gene_fr <- sapply(exp_gene_fr, function(x){strsplit(x, ";")[[1]]}) %>% unlist()
names(exp_gene_fr) <- NULL
exp_gene_fr <- unique(exp_gene_fr)


mutation_fr <- mutation_dt[gene_name_orig %in% exp_gene_fr, ]
mutation_fr_sig <- mutation_fr[samp_symbol %in% fr_res[, samp_symbol], ]
mutation_fr_sig[, absplice_variant := samp_symbol %in% absplice_res_sig[, samp_symbol]]
mutation_fr_sig[, absplice_vep := (absplice_variant & splice_related_variant)]
mutation_fr_sig[, absplice_only := (absplice_variant & !absplice_vep)]
mutation_fr_sig[, vep_only := (splice_related_variant & !absplice_vep)]
mutation_fr_sig[, mutation_type_count := sum(promoter_variant, frameshift_variant, stop_gained, 
                                             absplice_vep, absplice_only, vep_only,
                                             structural_variant, copy_number_gain, copy_number_loss), 
                by=rownames(mutation_fr_sig)]

mutation_fr_sig[, table(mutation_type_count)]
mutation_fr_sig[, multiple_type := mutation_type_count >=2]
mutation_fr_sig[multiple_type == TRUE, promoter_variant := FALSE]
mutation_fr_sig[multiple_type == TRUE, frameshift_variant := FALSE]
mutation_fr_sig[multiple_type == TRUE, stop_gained := FALSE]
mutation_fr_sig[multiple_type == TRUE, absplice_vep := FALSE]
mutation_fr_sig[multiple_type == TRUE, absplice_only := FALSE]
mutation_fr_sig[multiple_type == TRUE, vep_only := FALSE]
mutation_fr_sig[multiple_type == TRUE, structural_variant := FALSE]
mutation_fr_sig[multiple_type == TRUE, copy_number_gain := FALSE]
mutation_fr_sig[multiple_type == TRUE, copy_number_loss := FALSE]

mutation_fr_sig_melt <- melt(mutation_fr_sig, id.vars=c('gene_id', 'sampleID', 'padjust', 'l2fc', 'tool'),
                             measure.vars = c('promoter_variant', 'frameshift_variant', 'stop_gained', 
                                              'absplice_vep', 'absplice_only', 'vep_only', 
                                              'structural_variant', 'copy_number_gain', 'copy_number_loss',
                                              'multiple_type'),
                             variable.name = 'mutation_type', value.name='mutation_exist')
mutation_fr_sig_melt <- mutation_fr_sig_melt[mutation_exist == TRUE,]
fr_sig_overlap <- mutation_fr_sig_melt[, .N, by='mutation_type']
setnames(fr_sig_overlap, "N", "count")
fr_sig_overlap <- rbind(fr_sig_overlap, 
                        data.table(mutation_type = 'none', 
                                   count=(nrow(fr_res)-fr_sig_overlap[, sum(count)]))
)
fr_sig_overlap[, Category := 'Splicing outlier']
fr_sig_overlap[, total := sum(count)]

mutation_fr_non_sig <- mutation_fr[!samp_symbol %in% fr_res[, samp_symbol], ]
mutation_fr_non_sig[, absplice_variant := samp_symbol %in% absplice_res_sig[, samp_symbol]]
mutation_fr_non_sig[, absplice_vep := absplice_variant & splice_related_variant]
mutation_fr_non_sig[, absplice_only := absplice_variant & !absplice_vep]
mutation_fr_non_sig[, vep_only := splice_related_variant & !absplice_vep]
mutation_fr_non_sig[, mutation_type_count := sum(promoter_variant, frameshift_variant, stop_gained, 
                                                 absplice_vep, absplice_only, vep_only,
                                                 structural_variant, copy_number_gain, copy_number_loss), 
                    by=rownames(mutation_fr_non_sig)]

mutation_fr_non_sig[, table(mutation_type_count)]
mutation_fr_non_sig[, multiple_type := mutation_type_count >=2]
mutation_fr_non_sig[multiple_type == TRUE, promoter_variant := FALSE]
mutation_fr_non_sig[multiple_type == TRUE, frameshift_variant := FALSE]
mutation_fr_non_sig[multiple_type == TRUE, stop_gained := FALSE]
mutation_fr_non_sig[multiple_type == TRUE, absplice_vep := FALSE]
mutation_fr_non_sig[multiple_type == TRUE, absplice_only := FALSE]
mutation_fr_non_sig[multiple_type == TRUE, vep_only := FALSE]
mutation_fr_non_sig[multiple_type == TRUE, structural_variant := FALSE]
mutation_fr_non_sig[multiple_type == TRUE, copy_number_gain := FALSE]
mutation_fr_non_sig[multiple_type == TRUE, copy_number_loss := FALSE]

mutation_fr_non_sig_melt <- melt(mutation_fr_non_sig, id.vars=c('gene_id', 'sampleID', 'padjust', 'l2fc', 'tool'),
                                 measure.vars = c('promoter_variant', 'frameshift_variant', 'stop_gained', 
                                                  'absplice_vep', 'absplice_only', 'vep_only', 
                                                  'structural_variant', 'copy_number_gain', 'copy_number_loss',
                                                  'multiple_type'),
                                 variable.name = 'mutation_type', value.name='mutation_exist')
mutation_fr_non_sig_melt <- mutation_fr_non_sig_melt[mutation_exist == TRUE,]
fr_non_sig_overlap <- mutation_fr_non_sig_melt[, .N, by='mutation_type']
setnames(fr_non_sig_overlap, "N", "count")
fr_non_sig_overlap <- rbind(fr_non_sig_overlap, 
                            data.table(mutation_type = 'none', 
                                       count=(nrow(mutation_fr_non_sig)-fr_sig_overlap[, sum(count)]))
)
fr_non_sig_overlap[, Category := 'None outlier']
fr_non_sig_overlap[, total := sum(count)]

prop_mutation_fr <- rbind(fr_sig_overlap, fr_non_sig_overlap)
prop_mutation_fr[, Percentage := count/total]

fwrite(prop_mutation_fr, snakemake@output$prop_mutation_fr)


fr_res <- merge(fr_res, fr_res[, .(hgncSymbol, rank(padjust)), by='sampleID'], 
                by = c('hgncSymbol', 'sampleID'))
setnames(fr_res, 'V2', 'padjustRankInSample')
fr_top3 <- fr_res[padjustRankInSample <= 3, ]

mutation_fr[, aberrant := samp_symbol %in% fr_res[, samp_symbol]]
mutation_fr[, top3 := samp_symbol %in% fr_top3[, samp_symbol]]
mutation_fr[, absplice_variant := samp_symbol %in% absplice_res_sig[, samp_symbol]]

fwrite(mutation_fr, snakemake@output$mutation_fr)