#'---
#' title: figure_3_prep_vep
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
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - predictedConsequence: '`sm config["predictedConsequence"]`'
#'    - mll_panel_genes: '`sm config["mll_panel_genes"]`'
#'    - AML_variants: '`sm config["AML_variants"]`'
#'    - htmlOutputPath: '`sm config["htmlOutputPath"] + "/manuscript"`'
#'    - CGC_leukemia_OCG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_OCG_list.tsv"`'
#'    - CGC_leukemia_TSG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_TSG_list.tsv"`'
#'    - anonymization_table_mll: '`sm config["anonymization_table_mll"]`'
#'    - abspliceDirNoHeader: '`sm config["abspliceDirNoHeader"]`'
#'  input:
#'    - vepRes: '`sm expand(config["vep_path"] +
#'              "/MLL_{dataset}.vep.tsv.gz", dataset=outputDatasets)`'
#'    - fraserResJunc: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results_per_junction.tsv",
#'                    annotation=annotations, dataset=outputDatasets)`'
#'    - abspliceResVar: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice_var.tsv"`'
#'  output:
#'    - venn_list: '`sm config["projectPath"] + 
#'                     "/manuscript/figure_3/plot_data/venn_list.Rds"`'
#'    - vep_splice: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_splice.tsv"`'
#'    - vep_res_abSplice: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_res_abSplice"`'
#'    - vep_full: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/vep_full.tsv"`'
#'  type: script
#'  resources:
#'    - mem_mb: 256000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_3.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_3.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})




# get parameters
set.seed(2023)

gencode <- fread(snakemake@params$gencode)

absplice_res_var <- fread(snakemake@input$abspliceResVar)
absplice_res_var[, samp_symbol := paste0(sampleID, "-", gene_name)]
absplice_res_var <- separate(absplice_res_var, col='variant', into=c('chr', 'pos', 'ref_alt'), sep=":", remove = FALSE)
absplice_res_var <- separate(absplice_res_var, col='ref_alt', into=c('ref', 'alt'), sep=">") %>% as.data.table()


# read in fr junc
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
fr_res_junc[, samp_symbol := paste0(sampleID, "-", hgncSymbol)]
fr_res_junc[, junction_id := paste0("Intron ", seqnames, ": ", 
                                    format(start, nsmall=1, big.mark=","), "-", 
                                    format(end, nsmall=1, big.mark=",")), 
            by = rownames(fr_res_junc)]


# read in vep 
for (x in snakemake@input$vepRes) {
  vep_res_temp <- fread(x)
  vep_res_temp <- vep_res_temp %>% unique()
  vep_res_temp <- merge(vep_res_temp, gencode[, .(gene_name, gene_type)], 
                        by.x='SYMBOL', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  vep_res_temp <- vep_res_temp[gene_type=='protein_coding', ]
  
  if (exists("vep_res")) {
    vep_res <- rbind(vep_res, vep_res_temp)
  } else {
    vep_res <- vep_res_temp
  }
}
vep_res <- vep_res %>% separate('#Uploaded_variation', 
                                c("col1", "array_id", "alt", "ref", "pos"), 
                                sep='__', remove = FALSE) %>% as.data.table()
vep_res[, samp_symbol := paste0(array_id, "-", SYMBOL)]
vep_splice <- vep_res[grep('splice', Consequence), ]
vep_splice[, table(Consequence)]

fwrite(vep_splice, snakemake@output$vep_splice)
fwrite(vep_res, snakemake@output$vep_full)


vcfs <- list.files(path=snakemake@params$abspliceDirNoHeader, full.names = TRUE)

for (vcf in vcfs) {
  
  vep_res_temp <- fread(vcf)
  vep_res_temp <-vep_res_temp[grep("splice", CSQT), ]
  sample_id <- colnames(vep_res_temp)[10]
  setnames(vep_res_temp, c(sample_id), c("sample"))
  vep_res_temp <- vep_res_temp[,sample_id := sample_id]
  
  if (exists("vep_res_abSplice")) {
    vep_res_abSplice <- rbind(vep_res_abSplice, vep_res_temp)
  } else {
    vep_res_abSplice <- vep_res_temp
  }
}

fwrite(vep_res_abSplice, snakemake@output$vep_res_abSplice)


#### overlap venn ####
x <- list( 
  'Splice-affecting variants' = absplice_res_var[AbSplice_DNA >= 0.2, samp_symbol], 
  'Splicing outliers' = fr_res_junc[, samp_symbol],
  'VEP splice variants' = vep_splice[, samp_symbol]
)

saveRDS(x, snakemake@output$venn_list)


