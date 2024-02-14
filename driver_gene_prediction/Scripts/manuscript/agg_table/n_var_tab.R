#'---
#' title: n_var_tab
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
#'    - manuscriptWording: '`sm config["manuscriptWording"]`'
#'    - anonymization_table_mll: '`sm config["anonymization_table_mll"]`' 
#'    - CGC_cancer_gene_processed: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_processed.tsv"`'
#'    - abspliceDirNoHeader: '`sm config["abspliceDirNoHeader"]`'
#'  output:
#'    - n_var_samp_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/n_var_samp_tab.csv"`'
#'    - n_var_gene_tab: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/n_var_gene_tab.csv"`'
#'    - n_var_samp_list: '`sm config["projectPath"] + 
#'                     "/manuscript/sup_table/n_var_samp_list.Rds"`'
#'  type: script
#'  threads: 96  
#'  resources:
#'    - mem_mb: 96000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/n_var_tab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202401/processed_data/snakemake/n_var_tab.snakemake")
print("Snakemake saved") 

suppressPackageStartupMessages({
  library(data.table) 
  library(BiocParallel)
  library(magrittr)
  library(tidyr)
})




#+ load 
single_group <- snakemake@params$inputDatasets

gencode <- fread(snakemake@params$gencode)
gencode[, gene_id := strsplit(gene_id, '[.]')[[1]][1], by=rownames(gencode)]
gencode_pr <- gencode[gene_type=='protein_coding', ]

manuscript_wording <- fread(snakemake@params$manuscriptWording)
colnames(manuscript_wording) <- gsub(" ", "_", colnames(manuscript_wording))
diag_order <- manuscript_wording[order(-Number_of_samples_per_cohort), unique(Cohort_abbreviation)]

anonymization_table_mll <- fread(snakemake@params$anonymization_table_mll)


#+ parallel setup
register(MulticoreParam(snakemake@threads))
message(date(), ": Running with ", bpworkers(), " workers ...")


#+ number of input variant
absplice_var_input_dir <- snakemake@params$abspliceDirNoHeader
samp_anno <- fread(snakemake@params$sampAnno)
samp_anno_exp <- separate_rows(samp_anno, 'DROP_GROUP', sep=',') %>% as.data.table()
samp_anno_14group <- samp_anno_exp[DROP_GROUP==single_group, ]
samp_anno_14group[, absplice_input_file := gsub('somatic.vcf.gz', 'somatic.abSplice.txt', VCF_file), by=rownames(samp_anno_14group)]
samp_anno_14group[, absplice_input_path := paste0(absplice_var_input_dir, "/", absplice_input_file), by=rownames(samp_anno_14group)]
samp_anno_14group[, table(file.exists(absplice_input_path))]
samp_anno_14group[, Diag_sum := .N, by='Diag']

n_var_samp_list <- bplapply(samp_anno_14group[, ArrayID],
                            function(i){ 
                              abs_temp <- fread(samp_anno_14group[ArrayID == i, absplice_input_path])
                              abs_temp_explode <- separate_rows(abs_temp, 'CSQT', sep=',') %>% as.data.table()
                              abs_temp_explode <- abs_temp_explode %>% separate('CSQT', 
                                                                                c("col1", "SYMBOL", "ENST_id", "Consequence"), 
                                                                                sep='[|]', remove = FALSE) %>% as.data.table()
                              
                              # fix outdated symbol
                              abs_temp_explode[SYMBOL=='AC018720.10', SYMBOL := 'SPDYE5']
                              abs_temp_explode[SYMBOL=='RP11-1055B8.7', SYMBOL := 'BAHCC1']
                              abs_temp_explode[SYMBOL=='NPHP3', SYMBOL := 'NPHP3-ACAD11']
                              
                              abs_temp_explode_pr <- merge(abs_temp_explode, gencode_pr, by.x = 'SYMBOL', by.y = 'gene_name', all.x = FALSE, all.y = FALSE)
                              
                              n_var_samp <- nrow(abs_temp_explode_pr)
                              gene_samp <- abs_temp_explode_pr[, .N, by=c('gene_id','SYMBOL')]
                              gene_samp[, sampleID := i]
                              return(list(n_var_samp, gene_samp))
                            })

saveRDS(n_var_samp_list, snakemake@output$n_var_samp_list)

print(paste0('Number of abSplice input vcf:', length(samp_anno_14group[, absplice_input_path])))
print(paste0('Number of abSplice input variants:', sum(sapply(n_var_samp_list,"[[",1))))

n_var_gene_dt <- lapply(n_var_samp_list,function(x) x[[2]]) %>% rbindlist()
n_var_gene_dt <- merge(n_var_gene_dt, samp_anno_14group[, .(ArrayID, Diag)], by.x='sampleID', by.y='ArrayID')
n_var_gene_dt <- merge(n_var_gene_dt, 
                       manuscript_wording[, .(Cohort_abbreviation, Cohort_during_analysis)],
                       by.x='Diag', by.y='Cohort_during_analysis')
n_var_gene_dt[, n_var_gene := sum(N), by='SYMBOL']

n_var_gene_entity_dt <- n_var_gene_dt[, sum(N), by=c('gene_id', 'SYMBOL', 'Cohort_abbreviation')]
n_var_gene_entity_dt <- dcast(n_var_gene_entity_dt, gene_id + SYMBOL ~ Cohort_abbreviation, value.var='V1')
n_var_gene_entity_dt[is.na(n_var_gene_entity_dt)] <- 0 

n_var_gene_entity_dt <- merge(n_var_gene_entity_dt, n_var_gene_dt[, .(SYMBOL, n_var_gene)] %>% unique(), by='SYMBOL')
n_var_gene_entity_dt <- n_var_gene_entity_dt[order(-n_var_gene),]

setnames(n_var_gene_entity_dt, c("gene_id", "SYMBOL", "n_var_gene"), c('GeneID', 'GeneSymbol', "Total"))
setcolorder(n_var_gene_entity_dt, c("GeneID", "GeneSymbol", diag_order, "Total"))

fwrite(n_var_gene_entity_dt, snakemake@output$n_var_gene_tab)




n_var_samp_dt <- samp_anno_14group[, .(ArrayID, Diag, Diag_sum)]
n_var_samp_dt[, Number_of_variant := sapply(n_var_samp_list,"[[",1)]

hist(n_var_samp_dt[, Number_of_variant], 100)
hist(n_var_samp_dt[, log10(Number_of_variant)], 100)

n_var_samp_dt <- merge(n_var_samp_dt, 
                       manuscript_wording[, .(Cohort_abbreviation, Cohort_during_analysis)],
                       by.x='Diag', by.y='Cohort_during_analysis')
n_var_samp_dt[, Diag := NULL]
setnames(n_var_samp_dt, 'Cohort_abbreviation', 'DiseaseEntity')

n_var_samp_dt <- merge(anonymization_table_mll, n_var_samp_dt, by='ArrayID')
n_var_samp_dt[, ArrayID := NULL]

n_var_samp_dt <- n_var_samp_dt[order(-Diag_sum, -Number_of_variant),]
n_var_samp_dt[, Diag_sum := NULL]

setcolorder(n_var_samp_dt, c("AnonamizedID", "DiseaseEntity", "Number_of_variant"))

fwrite(n_var_samp_dt, snakemake@output$n_var_samp_tab)