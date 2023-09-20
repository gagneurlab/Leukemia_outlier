#'---
#' title: generate_driver_gene_list
#' author: xueqicao
#' wb:
#'  params:
#'    - projectPath: '`sm config["projectPath"]`'
#'  input:
#'    - gencode: '`sm config["gencode"]`'
#'    - CGC_cancer_gene: '`sm config["CGC_cancer_gene"]`'
#'    - IntOGen_cancer_gene: '`sm config["IntOGen_cancer_gene"]`'
#'    - MLL_leukemia_gene: '`sm config["MLL_leukemia_gene"]`'
#'  output:
#'    - CGC_cancer_gene_processed: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_processed.tsv"`'
#'    - CGC_AML_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_AML_gene_list.tsv"`'
#'    - CGC_AML_OCG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_AML_OCG_list.tsv"`'
#'    - CGC_AML_TSG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_AML_TSG_list.tsv"`'
#'    - CGC_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_gene_list.tsv"`'
#'    - CGC_leukemia_OCG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_OCG_list.tsv"`'
#'    - CGC_leukemia_TSG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_leukemia_TSG_list.tsv"`'
#'    - CGC_cancer_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_gene_list.tsv"`'
#'    - CGC_cancer_OCG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_OCG_list.tsv"`'
#'    - CGC_cancer_TSG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/CGC_cancer_TSG_list.tsv"`'
#'    - IntOGen_cancer_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/IntOGen_cancer_gene_list.tsv"`'
#'    - IntOGen_cancer_OCG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/IntOGen_cancer_OCG_list.tsv"`'
#'    - IntOGen_cancer_TSG_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/IntOGen_cancer_TSG_list.tsv"`'
#'    - MLL_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/MLL_leukemia_gene_list.tsv"`'
#'    - MLL_CGC_leukemia_gene_list: '`sm config["projectPath"] + "/processed_data/driver_gene_list/MLL_CGC_leukemia_gene_list.tsv"`'
#' output:   
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

#+ echo=F
saveRDS(snakemake, file.path(snakemake@params$projectPath, 
                             "/processed_data/snakemake/generate_label_list.snakemake"))
# snakemake <- readRDS("/s/project/vale/driver_prediction_202304/processed_data/snakemake/generate_label_list.snakemake")




#+ read in data
suppressPackageStartupMessages({
  library(data.table)
  library(magrittr)
  library(tidyr)
  library(readxl)
})

gencode <- fread(snakemake@input$gencode)
gencode[, gene_id_unique := gene_id]
ensg_id <- sapply(gencode[, gene_id], function(x){strsplit(x, "[.]")[[1]][1]})
gencode[, gene_id := ensg_id]



# process CGC
CGC_cancer_gene <- fread(snakemake@input$CGC_cancer_gene)
colnames(CGC_cancer_gene) <- gsub(" ","", colnames(CGC_cancer_gene))

# manual fix
CGC_cancer_gene[GeneSymbol=='CARS', GeneSymbol:='CARS1']
CGC_cancer_gene[GeneSymbol=='GAS7', GeneSymbol:='AC005747.1']
CGC_cancer_gene[GeneSymbol=='H3F3A', GeneSymbol:='H3-3A']
CGC_cancer_gene[GeneSymbol=='H3F3B', GeneSymbol:='H3-3B']
CGC_cancer_gene[GeneSymbol=='SEPT5', GeneSymbol:='SEPTIN5']
CGC_cancer_gene[GeneSymbol=='SEPT6', GeneSymbol:='SEPTIN6']
CGC_cancer_gene[GeneSymbol=='SEPT9', GeneSymbol:='SEPTIN9']

# the ENSGid can not be found for gene "IGH" "IGK" "IGL" "TRA" "TRB" "TRD"
# merge on GeneSymbol: !the GeneSymbol is correct but Synonyms contains mistakes!
CGC_cancer_gene <- merge(CGC_cancer_gene, gencode[, .(gene_name, gene_id, gene_type)], by.x='GeneSymbol', by.y='gene_name')
setnames(CGC_cancer_gene, 'gene_id', 'ENSGid')

fwrite(CGC_cancer_gene,
       file = snakemake@output$CGC_cancer_gene_processed,
       sep = "\t", quote = F)

# get protein coding gene for the model
CGC_cancer_gene <- CGC_cancer_gene[ENSGid %in% gencode[gene_type=='protein_coding', gene_id],]



# process IntOGen
IntOGen_cancer_gene <- fread(snakemake@input$IntOGen_cancer_gene)

IntOGen_cancer_gene[SYMBOL == "CARS", SYMBOL:= "CARS1"]
IntOGen_cancer_gene[SYMBOL == "FAM46C", SYMBOL:= "TENT5C"]
IntOGen_cancer_gene[SYMBOL == "H3F3A", SYMBOL:= "H3-3A"]
IntOGen_cancer_gene[SYMBOL == "HIST1H3B", SYMBOL:= "H3C2"]
IntOGen_cancer_gene[SYMBOL == "HIST1H4I", SYMBOL:= "H4C9"]
IntOGen_cancer_gene[SYMBOL == "SEPT9", SYMBOL:= "SEPTIN9"]

IntOGen_cancer_gene <- merge(IntOGen_cancer_gene, gencode[, .(gene_id, gene_name)],
                             by.x = "SYMBOL", by.y = "gene_name", all.x = TRUE)
setnames(IntOGen_cancer_gene, "gene_id", "ENSGid")
setnames(IntOGen_cancer_gene, "SYMBOL", "GeneSymbol")

IntOGen_cancer_gene <- IntOGen_cancer_gene[!is.na(ENSGid)]
IntOGen_cancer_gene <- IntOGen_cancer_gene[ENSGid %in% gencode[gene_type=='protein_coding', gene_id],]



# process MLL genes
mll_leukemia_gene <- read_xlsx(snakemake@input$MLL_leukemia_gene) %>% as.data.table()

mll_leukemia_gene[, gene_name := Gene]
mll_leukemia_gene[Gene == "FAM46C", gene_name:= "TENT5C"]
mll_leukemia_gene[Gene == "GPR98", gene_name:= "ADGRV1"]
mll_leukemia_gene[Gene == "WHSC1", gene_name:= "NSD2"]

mll_leukemia_gene <- merge(mll_leukemia_gene, gencode[, .(gene_id, gene_name)],
                           by = "gene_name", all.x = TRUE)
setnames(mll_leukemia_gene, "gene_id", "ENSGid")
setnames(mll_leukemia_gene, "gene_name", "GeneSymbol")

mll_cgc_leukemia_gene <- rbind(CGC_cancer_gene[grep("L", TissueType), .(GeneSymbol, ENSGid)] ,
                               mll_leukemia_gene[, .(GeneSymbol, ENSGid)]) %>% unique()




##### CGC cancer gene list #####
#+ CGC cancer gene list
fwrite(CGC_cancer_gene[, .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_cancer_gene_list,
       sep = "\t", quote = F)

#+ CGC cancer oncogene list
fwrite(CGC_cancer_gene[grep("oncogene", RoleinCancer), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_cancer_OCG_list,
       sep = "\t", quote = F)

#+ CGC cancer TSG list
fwrite(CGC_cancer_gene[grep("TSG", RoleinCancer), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_cancer_TSG_list,
       sep = "\t", quote = F)




##### CGC leukemia gene list #####
#+ CGC leukemia gene list
fwrite(CGC_cancer_gene[grep("L", TissueType), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_leukemia_gene_list,
       sep = "\t", quote = F)

#+ CGC leukemia oncogene list
fwrite(CGC_cancer_gene[intersect(grep("L", TissueType), grep("oncogene", RoleinCancer)), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_leukemia_OCG_list,
       sep = "\t", quote = F)

#+ CGC leukemia TSG list
fwrite(CGC_cancer_gene[intersect(grep("L", TissueType), grep("TSG", RoleinCancer)), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_leukemia_TSG_list,
       sep = "\t", quote = F)




##### CGC AML gene list #####
#+ CGC AML gene list
fwrite(CGC_cancer_gene[grep("AML", `TumourTypes(Somatic)`), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_AML_gene_list,
       sep = "\t", quote = F)

#+ CGC AML oncogene list
fwrite(CGC_cancer_gene[intersect(grep("AML", `TumourTypes(Somatic)`), grep("oncogene", RoleinCancer)), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_AML_OCG_list,
       sep = "\t", quote = F)

#+ CGC AML TSG list
fwrite(CGC_cancer_gene[intersect(grep("AML", `TumourTypes(Somatic)`), grep("TSG", RoleinCancer)), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$CGC_AML_TSG_list,
       sep = "\t", quote = F) 




##### IntOGen cancer gene list #####
#+ IntOGen cancer gene list
fwrite(IntOGen_cancer_gene[, .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$IntOGen_cancer_gene_list,
       sep = "\t", quote = F)

#+ IntOGen cancer oncogene list
fwrite(IntOGen_cancer_gene[grep("Act", ROLE), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$IntOGen_cancer_OCG_list,
       sep = "\t", quote = F)

#+ IntOGen cancer TSG list
fwrite(IntOGen_cancer_gene[grep("LoF", ROLE), .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$IntOGen_cancer_TSG_list,
       sep = "\t", quote = F)




##### MLL leukemia gene list #####
#+ MLL leukemia gene list
fwrite(mll_leukemia_gene[, .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$MLL_leukemia_gene_list,
       sep = "\t", quote = F)

#+ MLL CGC leukemia gene list
fwrite(mll_cgc_leukemia_gene[, .(GeneSymbol, ENSGid)] %>% unique(),
       file = snakemake@output$MLL_CGC_leukemia_gene_list,
       sep = "\t", quote = F)
