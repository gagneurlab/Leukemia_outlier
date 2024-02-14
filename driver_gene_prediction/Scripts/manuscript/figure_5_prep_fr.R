#'---
#' title: figure_5 preparation OUTRIDER
#' author: xueqicao
#' wb:
#'  py:
#'   - |
#'    annotations = config["geneAnnotation"]
#'    inputDatasets = config["cohort"]["groups"]
#'    outputDatasets = config["cohort"]["groups"]
#'  params:
#'    - projectPath: '`sm config["projectPath"]`'
#'    - annotations: '`sm annotations`'
#'    - gencode: '`sm config["gencode"]`'
#'    - sampAnno: '`sm config["sampleAnnotation"]`'
#'    - inputDatasets: '`sm inputDatasets`'
#'    - outputDatasets: '`sm outputDatasets`'
#'    - CGC_cancer_gene_processed: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_processed.tsv"`'
#'  input:
#'    - fraserRes: '`sm expand(config["fraserDir"] +
#'                    "/processed_results/aberrant_splicing/results/{annotation}/fraser/{dataset}/results.tsv",
#'                    annotation=annotations, dataset=groupsDatasets)`'
#'  output:
#'    - fr_res: '`sm config["projectPath"] + 
#'                      "/processed_data/fr_res.tsv"`'
#'  type: script
#'  threads: 84  
#'  resources:
#'    - mem_mb: 400000
#'---

#+ echo=FALSE
saveRDS(snakemake, file.path(snakemake@params$projectPath,
                             "/processed_data/snakemake/figure_5_prep_fr.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202303_2/processed_data/snakemake/figure_5_prep_fr.snakemake")
print("Snakemake saved") 




library(data.table)

groups <- snakemake@params$inputDatasets
fr_res_all <-  data.table(hgncSymbol=character(), sampleID=character())

for (x in groups) {
  fr_res <- fread(snakemake@input$fraserRes[groups==x])
  fr_res <- fr_res[, .(hgncSymbol, sampleID)]
  fr_res_all <- rbind(fr_res_all, fr_res)
}


#+ save fr_res_all
fwrite(fr_res_all, snakemake@output$fr_res, sep='\t')
