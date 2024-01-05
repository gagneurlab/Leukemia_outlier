#'---
#' title: figure_3 preparation absplice
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
#'    - CGC_cancer_gene_processed: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_processed.tsv"`'
#'    - abspliceDirNoHeader: '`sm config["abspliceDirNoHeader"]`'
#'  input:
#'    - abspliceRes: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/res_absplice.tsv"`'
#'  output:
#'    - enrichmentTab: '`sm config["projectPath"] + 
#'                    "/manuscript/figure_3/plot_data/enrichments_absplice.csv"`'
#'    - wBhtml: '`sm config["htmlOutputPath"] + "/manuscript/figure_3_prep_absplice.html"`'
#'  type: noindex
#'  threads: 96  
#'  resources:
#'    - mem_mb: 96000
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_3_prep_ab.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/figure_3_prep_ab.snakemake")
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

cgc_cancer_gene <- fread(snakemake@params$CGC_cancer_gene_processed)
cgc_pr <- cgc_cancer_gene[gene_type=='protein_coding', ]

absplice_res <- fread(snakemake@input$abspliceRes)

n_rep <- 99


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

n_var_samp_list <- bplapply(samp_anno_14group[, absplice_input_path],
         function(i){ 
           abs_temp <- fread(i)
           abs_temp_explode <- separate_rows(abs_temp, 'CSQT', sep=',') %>% as.data.table()
           abs_temp_explode <- abs_temp_explode %>% separate('CSQT', 
                                           c("col1", "SYMBOL", "ENST_id", "Consequence"), 
                                           sep='[|]', remove = FALSE) %>% as.data.table()
           
           abs_temp_explode_pr <- merge(abs_temp_explode, gencode_pr, by.x = 'SYMBOL', by.y = 'gene_name', all.x = FALSE, all.y = FALSE)
           
           n_var_samp <- nrow(abs_temp_explode_pr)
           gene_samp <- unique(abs_temp_explode_pr[, SYMBOL])
           return(list(n_var_samp, gene_samp))
         })

print(paste0('Number of abSplice input vcf:', length(samp_anno_14group[, absplice_input_path])))
print(paste0('Number of abSplice input variants:', sum(sapply(n_var_samp_list,"[[",1))))
print(paste0('Number of abSplice input genes:', length(unique(unlist(sapply(n_var_samp_list,"[[",2))))))

#+ read in
absplice_res_02 <- absplice_res[AbSplice_DNA>=0.2, ]
absplice_res_02[, aberrantByGene := length(sampleID), by='gene_id']

#+ categorize
absplice_sum <- absplice_res_02[, .(gene_id, gene_name, aberrantByGene)] %>% unique()
absplice_sum <- merge(absplice_sum[, gene_name := NULL], 
                      gencode_pr[, .(gene_id, gene_name)], 
                      by='gene_id', all.x=TRUE, all.y=TRUE)

absplice_sum[is.na(aberrantByGene), aberrantByGene := 0]
absplice_sum[, rank_aberrantByGene := rank(-aberrantByGene, ties.method = 'min')]
absplice_sum <- absplice_sum[order(rank_aberrantByGene), ]

absplice_sum[, num_samp_abs_per_gene := as.character(aberrantByGene)]
absplice_sum[aberrantByGene >=8, num_samp_abs_per_gene := '>=8']

#+ simulate random
absplice_cat_simu_ls <- bplapply(1:n_rep,
                                 function(i){ 
                                   cgc_simu <- copy(cgc_pr)
                                   cgc_simu <- cgc_simu[, ENSGid := sample(gencode_pr[, gene_id], nrow(cgc_pr))]
                                   
                                   absplice_sum[, isCGC := gene_id %in% cgc_simu[, ENSGid]]
                                   absplice_sum[, isOCG := gene_id %in% cgc_simu[grep('oncogene', RoleinCancer), ENSGid]]
                                   absplice_sum[, isTSG := gene_id %in% cgc_simu[grep('TSG', RoleinCancer), ENSGid]]
                                   absplice_sum[, isCGCleu := gene_id %in% cgc_simu[grep('L', TissueType), ENSGid]]
                                   absplice_sum[, isOCGleu := gene_id %in% cgc_simu[intersect(grep("L", TissueType), grep("oncogene", RoleinCancer)), ENSGid]]
                                   absplice_sum[, isTSGleu := gene_id %in% cgc_simu[intersect(grep("L", TissueType), grep("TSG", RoleinCancer)), ENSGid]]
                                   
                                   absplice_cat_simu <- data.table(
                                     num_samp_abs_per_gene=unique(absplice_sum[, num_samp_abs_per_gene])
                                   )
                                   
                                   sapply(absplice_cat_simu[, num_samp_abs_per_gene], function(x){
                                     absplice_cat_simu[num_samp_abs_per_gene==x, CGC := absplice_sum[num_samp_abs_per_gene==x, sum(isCGC)]] 
                                     absplice_cat_simu[num_samp_abs_per_gene==x, OCG := absplice_sum[num_samp_abs_per_gene==x, sum(isOCG)]] 
                                     absplice_cat_simu[num_samp_abs_per_gene==x, TSG := absplice_sum[num_samp_abs_per_gene==x, sum(isTSG)]] 
                                     absplice_cat_simu[num_samp_abs_per_gene==x, CGCleu := absplice_sum[num_samp_abs_per_gene==x, sum(isCGCleu)]] 
                                     absplice_cat_simu[num_samp_abs_per_gene==x, OCGleu := absplice_sum[num_samp_abs_per_gene==x, sum(isOCGleu)]] 
                                     absplice_cat_simu[num_samp_abs_per_gene==x, TSGleu := absplice_sum[num_samp_abs_per_gene==x, sum(isTSGleu)]] 
                                   })
                                   
                                   return(absplice_cat_simu)
                                 })
absplice_cat_simu <- rbindlist(absplice_cat_simu_ls)
absplice_cat_simu <- melt(absplice_cat_simu, id.vars = 'num_samp_abs_per_gene',
                          variable.name = 'gene_list', value.name = 'count')

absplice_cat_simu[, `:=` (random = median(count),
                          random95 = quantile(count, 0.95),
                          random05 = quantile(count, 0.05)), by=c('num_samp_abs_per_gene', 'gene_list')]
absplice_cat_simu[, count:= NULL] 
en_dt <- absplice_cat_simu %>% unique()


#+ calculate real
absplice_sum[, isCGC := gene_id %in% cgc_pr[, ENSGid]]
absplice_sum[, isOCG := gene_id %in% cgc_pr[grep('oncogene', RoleinCancer), ENSGid]]
absplice_sum[, isTSG := gene_id %in% cgc_pr[grep('TSG', RoleinCancer), ENSGid]]
absplice_sum[, isCGCleu := gene_id %in% cgc_pr[grep('L', TissueType), ENSGid]]
absplice_sum[, isOCGleu := gene_id %in% cgc_pr[intersect(grep("L", TissueType), grep("oncogene", RoleinCancer)), ENSGid]]
absplice_sum[, isTSGleu := gene_id %in% cgc_pr[intersect(grep("L", TissueType), grep("TSG", RoleinCancer)), ENSGid]]

sapply(en_dt[, num_samp_abs_per_gene], function(x){
  en_dt[num_samp_abs_per_gene==x & gene_list=='CGC', 
        real := absplice_sum[num_samp_abs_per_gene==x, sum(isCGC)]] 
  en_dt[num_samp_abs_per_gene==x & gene_list=='OCG', 
        real := absplice_sum[num_samp_abs_per_gene==x, sum(isOCG)]] 
  en_dt[num_samp_abs_per_gene==x & gene_list=='TSG', 
        real := absplice_sum[num_samp_abs_per_gene==x, sum(isTSG)]] 
  en_dt[num_samp_abs_per_gene==x & gene_list=='CGCleu', 
        real := absplice_sum[num_samp_abs_per_gene==x, sum(isCGCleu)]] 
  en_dt[num_samp_abs_per_gene==x & gene_list=='OCGleu', 
        real := absplice_sum[num_samp_abs_per_gene==x, sum(isOCGleu)]] 
  en_dt[num_samp_abs_per_gene==x & gene_list=='TSGleu', 
        real := absplice_sum[num_samp_abs_per_gene==x, sum(isTSGleu)]] 
})

en_dt[, enrichment := real/random]
en_dt[, en95 := real/random95]
en_dt[, en05 := real/random05]


#+ save enrichments
fwrite(en_dt, snakemake@output$enrichmentTab)

