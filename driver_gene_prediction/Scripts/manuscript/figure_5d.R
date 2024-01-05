
print('loading packages')
suppressPackageStartupMessages({
  library(data.table)
  library(ggbio)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(tidyverse)
  library(ggforce)
  library(gridExtra)
  library(ggrastr)
  library(ggpubr)
})

source("Scripts/manuscript/manuscript_theme.R")





ac_res <- fread("/s/project/mll/drop_2023feb/processed_results/aberrant_expression/v33b/outrider/leukemia_14group/res_filter_out.tsv")
samp_anno_mll <- fread("/s/project/mll/sample_info/sample_annotation_for_drop.tsv")
samp_anno_marc <- fread("/s/project/vale/drop_marc_lrp1b/sample_annotation.tsv")
lrp1b_act_samples <- ac_res[hgncSymbol=='LRP1B', sampleID]
mll_anonymized_ids <- fread("/s/project/vale/Resource/anonymization_table_mll.csv")
marc_anonymized_ids <- fread("/s/project/vale/Resource/res_marc_seifert_lrp1b/anonymization_table_marc.csv")




# hg19 transcript



gr_chr <- "chr2"
gr_start <- 140986990
gr_end <- 142890599
gr_strand <- "-"
exon_range <- c(1:92)




exonic_regions <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
exon_region <- exonic_regions[seqnames(exonic_regions) == gr_chr &
                                start(exonic_regions) > gr_start &
                                end(exonic_regions) < gr_end &
                                strand(exonic_regions) == gr_strand]
exon_pos <- sapply(seq_along(exon_region), function(i){
  c(start(exon_region[i]):end(exon_region[i]))
}) %>% unlist()

# get transcript
region_dt <- exon_region %>% as.data.table()
region_dt[, transcript_name:='LRP1B']
region_dt <- as.data.table(region_dt)

region_dt <- region_dt[exon_id!='43123', ]
# first exon in igv: chr2: 142888217 - 142888585



transcript_plot <- ggplot() +
  #scale_x_continuous(name='Genomic coordinate (chr2), "-" strand') + 
  scale_y_continuous() +
  geom_rect(data=region_dt, mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1), color="black", alpha=0.5) +
  geom_hline(yintercept = 0.5) + 
  xlim(min(region_dt$start) , max(region_dt$end)) +
  # ggtitle("Exons of LRP1b Gene") +
  xlab("") +
  ylab("") +
  geom_text(aes(x=142130000, y = 0.58, label="<")) + 
  geom_text(aes(x=142400000, y = 0.58, label="<")) + 
  geom_text(aes(x=142700000, y = 0.58, label="<")) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = margin(0.6, 0.63, 0, 2.13, "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) 





#### marc coverage ####
genome_38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
lrp1b_38 <- genes(genome_38)[which(genes(genome_38)$gene_id ==  53353),]

lrp1b_38 <- GRanges(seqnames = "2",ranges = IRanges(start(lrp1b_38), width = end(lrp1b_38)-start(lrp1b_38)+1), strand = "-")
lrp1b_38_chr <- lrp1b_38
seqlevels(lrp1b_38_chr) <- "chr2" # some samples have "chr2" as their seqnames, not "2"
pos_ori <- start(lrp1b_38): end(lrp1b_38)
coverage_dt_marc <- data.table(pos_ori)







arrayID <- "ohio5_TAGTGCCA_L001_R1_001"
bam_sample <- samp_anno_marc[ArrayID == arrayID, RNA_BAM_FILE]
sample_name <- sapply(str_split(bam_sample, "/"), tail, 1)
anonymized_ID <- marc_anonymized_ids[ArrayID == arrayID, AnonamizedID]

p_coverage <- autoplot(bam_sample, which = lrp1b_38)

  
sample_coverage_dt <- ggplot_build(p_coverage@ggplot)$data[[1]] %>% as.data.table()
  
sample_coverage_dt <- sample_coverage_dt[, .(x, y)]
setnames(sample_coverage_dt, c('x', 'y'), c('pos_ori', anonymized_ID))
coverage_dt_marc <- merge(coverage_dt_marc, sample_coverage_dt, all.x=TRUE, all.y=TRUE)
  
  

# plotting whole gene
melted_dt_marc <- melt(coverage_dt_marc, "pos_ori")

melted_dt_marc$value[is.na(melted_dt_marc$value)] <- 0

exon_13_hg38 <- 141012000
marc_coverage <- ggplot(melted_dt_marc, aes(x = pos_ori, y = value)) + 
  geom_bar(stat="identity", color = "grey")  + 
  geom_vline(xintercept = exon_13_hg38, linetype = "dashed")  +
  
  facet_wrap("variable", ncol=1, scales = "free_y", strip.position = "left") + 
  ylab("") + xlab("") +   
  scale_x_continuous(breaks = NULL) +  
  #ylim(0, 2000) + 
  theme_vale + 
  theme(
        axis.text.y = element_text(size = 10), 
        strip.text.y = element_text(size = 8) ,
        legend.position="none", 
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        plot.margin = margin(-0.3, 0, 0, 0.5, "cm"),) 

marc_coverage <- rasterise(marc_coverage)
#marc_coverage

#### mll coverage ####
genome_19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
lrp1b_19 <- genes(genome_19)[which(genes(genome_19)$gene_id ==  53353),]
seqlevels(lrp1b_19) <- "chr2"

lrp1b_act_samples_mll <- samp_anno_mll[samp_anno_mll$ArrayID %in% lrp1b_act_samples][1:24]
pos_ori <- start(lrp1b_19): end(lrp1b_19)
coverage_dt <- data.table(pos_ori)

# Adding two mll samples


arrayID <- "MLL_25594"
anonymized_ID<- mll_anonymized_ids[ArrayID == arrayID, AnonamizedID]
bam_sample <-  samp_anno_mll[ArrayID == arrayID, RNA_BAM_FILE]
sample_name <- sapply(str_split(bam_sample, "/"), tail, 1)
p_coverage <- autoplot(bam_sample, which = lrp1b_19)
sample_coverage_dt <- ggplot_build(p_coverage@ggplot)$data[[1]] %>% as.data.table()
sample_coverage_dt <- sample_coverage_dt[, .(x, y)]
setnames(sample_coverage_dt, c('x', 'y'), c('pos_ori', anonymized_ID))
coverage_dt <- merge(coverage_dt, sample_coverage_dt, all.x=TRUE, all.y=TRUE)




arrayID <- "MLL_22363"
anonymized_ID<- mll_anonymized_ids[ArrayID == arrayID, AnonamizedID]
bam_sample <-  samp_anno_mll[ArrayID == arrayID, RNA_BAM_FILE]
sample_name <- sapply(str_split(bam_sample, "/"), tail, 1)


p_coverage <- autoplot(bam_sample, which = lrp1b_19)
sample_coverage_dt <- ggplot_build(p_coverage@ggplot)$data[[1]] %>% as.data.table()
sample_coverage_dt <- sample_coverage_dt[, .(x, y)]
setnames(sample_coverage_dt, c('x', 'y'), c('pos_ori', anonymized_ID))
coverage_dt <- merge(coverage_dt, sample_coverage_dt, all.x=TRUE, all.y=TRUE)



# plotting whole gene

melted_dt <- melt(coverage_dt, "pos_ori")

melted_dt$value[is.na(melted_dt$value)] <- 0

#exon_13_hg19: 14177
exon_13 <- 141773365

# reorder the plots, first mll truncated, then mll full and then marc truncated
re_orderd <- transform(melted_dt,
                       variable=factor(variable, rev(levels(melted_dt$variable))))

coverage <- ggplot(re_orderd, aes(x = pos_ori, y = value)) + 
  geom_bar(stat="identity", color = "grey") + 
  facet_wrap("variable", ncol=1, scales = "free_y", strip.position = "left") + 
  scale_x_continuous(breaks = NULL) + 
  ylab("") + xlab("") +   
  geom_vline(xintercept = exon_13, linetype = "dashed") + 
  theme_vale + 
  theme(
        axis.text.y = element_text(size = 10), 
        strip.text.y = element_text(size = 8) ,
        legend.position="none", 
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

coverage <- rasterise(coverage)

p_coverage <- ggarrange(coverage, marc_coverage, ncol=1,  heights = c(3.1, 1), align = "v")

p_coverage_transcript <- ggarrange(p_coverage, transcript_plot, ncol=1,  heights = c(4, 1))

                                                        
device=cairo_pdf




png("~/figure5d.png", 
    width = 15, height = 3.5, units = "in", res = 600,
    type="cairo")


p_coverage_transcript
dev.off()

output_file <- "~/figure5d.pdf"

pdf(output_file, width = 15, height = 3.5)
print(p_coverage_transcript)
dev.off()




