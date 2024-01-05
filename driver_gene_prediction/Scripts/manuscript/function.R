add_result_paths <- function(exp_viz, project_dir, intogen_dir){
  exp_viz[, res_pre_path := paste0(project_dir, '/processed_results/pre/result_', experiment_no, '.tsv')]
  exp_viz[, coeff_path := paste0(project_dir, '/processed_results/coeff_/result_', experiment_no, '.tsv')]
  exp_viz[, res_post_path := paste0(project_dir, '/processed_results/post/result_', experiment_no, '.tsv')]
  exp_viz[, rp_path := paste0(project_dir, '/processed_results/post/rp_', experiment_no, '.tsv')]
  exp_viz[, vep_path := paste0(intogen_dir, '/steps/vep/MLL_WGS_MLL_', toupper(sample_group), '.tsv.gz')]
  exp_viz[, dndscv_path := paste0(intogen_dir, '/steps/dndscv/MLL_WGS_MLL_', toupper(sample_group), '.dndscv.tsv.gz')]
}

add_input_feature <- function(exp_viz){
  exp_viz[, input_feature := paste(c(intogen_input_feature, outlier_input_feature, coess_input_feature)[
    c(intogen_input_feature, outlier_input_feature, coess_input_feature)!=''], 
    collapse = ','), by=list(rownames(exp_viz))]
}

calculate_precision_recall <- function(score, label, decr=TRUE, total=sum(label)){
  dt <- data.table(score=score, label=label)
  dt <- dt[order(score, decreasing=decr, na.last=TRUE)]
  dt[,rank      := 1:.N]
  dt[,TP        := cumsum(label)]
  dt[,c('TP', 'rank'):=list(max(TP), max(rank)), by=score]
  dt[,precision := TP/rank]
  dt[,recall    := TP/total]
  dt[,recall_diff := diff(c(0, recall))]
  dt[,au_prc := sum(precision * recall_diff)]
  dt[,recall_diff:=NULL]
  dt
}
# Chromosome lengths for hg19/GRCh37
chromlens <- tibble(
  chrom = as.character(c(1:22, "X", "Y")),
  len = as.numeric(c(
    249250621,
    243199373,
    198022430,
    191154276,
    180915260,
    171115067,
    159138663,
    146364022,
    141213431,
    135534747,
    135006516,
    133851895,
    115169878,
    107349540,
    102531392,
    90354753,
    81195210,
    78077248,
    59128983,
    63025520,
    48129895,
    51304566,
    155270560,
    59373566
  ))
)

add_chromlens <- chromlens %>%
  mutate(cumpos = cumsum(len),
         addpos = (cumpos - len) / 10^6)

names(add_chromlens)[names(add_chromlens) == "chrom"] <- "CONTIG"



plotCumulativeChromOUTRIDERzScore <- function(ods, sample, value, chrom_to_plot, color_chr = c("black", "firebrick")){
  stopifnot("Sample not in ods" = sample %in% colnames(ods))
  stopifnot("Value should be either pvalue, zscore or l2fc" =
              value %in% c('pvalue', 'pValue', 'pv', 
                           'padjust', 'pAdjust', 'pa',
                           'zscore', 'zScore', 
                           'l2fc', 'L2FC', 'log2fc', 
                           'counts'))
  # get granges
  if(all(c('start', 'seqnames') %in% colnames(rowData(ods)))){
    # if generated through DROP, ods should already contain seqnames and position
    gr <- GRanges(seqnames = rowData(ods)$seqnames, IRanges(start = rowData(ods)$start, end = rowData(ods)$end))
    # Add names to display on plot
    names(gr) <- row.names(ods)
  } else if(is(rowRanges(ods), "GRangesList")){
    # gr <- unlist(endoapply(rowRanges(ods), range))
    gr <- unlist(endoapply(rowRanges(ods), function(rr) rr[1,])) # faster than range
  } else if(is(rowRanges(ods), "GRanges")){
    gr <- GRanges(seqnames = seqnames(rowRanges(ods)), 
                  IRanges(start = start(rowRanges(ods)), end = end(rowRanges(ods))))
  }
  # seqlevelsStyle(gr) <- 'NCBI'
  # Add values to granges
  if(value %in% c('pvalue', 'pValue', 'pv')){
    gr$value <- -log10(assays(ods)$pValue[, sample])
    value <- '-log10(pvalue)'
  }
  if(value %in% c('padjust', 'pAdjust', 'pa')){
    gr$value <- -log10(assays(ods)$padjust[, sample])
    value <- '-log10(padjust)'
  }
  if(value %in% c('zscore', 'zScore'))
    gr$value <- assays(ods)$zScore[, sample]
  if(value %in% c('l2fc', 'L2FC', 'log2fc'))
    gr$value <- assays(ods)$l2fc[, sample]
  if(value == 'counts')
    gr$value <- log10(assays(ods)$counts[, sample]+1)
  gr$aberrant <- aberrant(ods)[,sample]
  # Sort granges for plot
  seqlevels(gr) <- gsub('chr','',seqlevels(gr))
  gr <- sortSeqlevels(gr)
  gr <- gr[gr@seqnames == chrom_to_plot]
  gr <- sort(gr)
  gr <- as.data.table(gr)
  max_value <- ceiling( max(abs(gr$value)))
  add_chromlens <- add_chromlens[add_chromlens$CONTIG == '4',]
  plotdata <- gr %>%
    filter(!is.na(seqnames), 
           seqnames != "MT") %>%
    left_join(add_chromlens, by = c("seqnames" = "CONTIG")) %>%
    mutate(START_mbp = start / 1e06,
           START_mbp_cum = START_mbp + addpos) #%>%
  #sample_frac(0.01) # sampling 1% of points is detailed enough
  
  
  p <- ggplot() +
    geom_point(data = plotdata,
               aes(x = START_mbp_cum, y = value, color = aberrant),
               alpha = 0.5, size = 0.5) +
    labs(y = value) + 
    ylim(-max_value, max_value) + 
    scale_color_manual(values = color_chr)
  
  return(p)
}



plotCumulativeChromCNV <- function(cnv, chrom_to_plot, color_chr = c("black", "firebrick")){
  add_chromlens <- add_chromlens[add_chromlens$CONTIG == chrom_to_plot,]
  cnv_chrom <- cnv[cnv$CONTIG == chrom_to_plot,]
  plotdata <- cnv_chrom %>%
    filter(!is.na(CONTIG), 
           CONTIG != "MT") %>%
    left_join(add_chromlens, by = c("CONTIG" = "CONTIG")) %>%
    mutate(START_mbp = START / 1e06,
           START_mbp_cum = START_mbp + addpos) #%>%
  #sample_frac(0.01) # sampling 1% of points is detailed enough
  
  # Plot: Coverage + chrom rectangles + chrom labels
  p <- ggplot() +
    geom_point(data = plotdata,
               aes(x = START_mbp_cum, y = LOG2_COPY_RATIO),
               alpha = 0.5, size = 0.5) 
  
  return(p)
}

plotChromOUTRIDERzScore <- function(ods, sample, value, chrom_to_plot, color_chr = c("black", "firebrick")){
  stopifnot("Sample not in ods" = sample %in% colnames(ods))
  stopifnot("Value should be either pvalue, zscore or l2fc" =
              value %in% c('pvalue', 'pValue', 'pv', 
                           'padjust', 'pAdjust', 'pa',
                           'zscore', 'zScore', 
                           'l2fc', 'L2FC', 'log2fc', 
                           'counts'))
  # get granges
  if(all(c('start', 'seqnames') %in% colnames(rowData(ods)))){
    # if generated through DROP, ods should already contain seqnames and position
    gr <- GRanges(seqnames = rowData(ods)$seqnames, IRanges(start = rowData(ods)$start, end = rowData(ods)$end))
    # Add names to display on plot
    names(gr) <- row.names(ods)
  } else if(is(rowRanges(ods), "GRangesList")){
    # gr <- unlist(endoapply(rowRanges(ods), range))
    gr <- unlist(endoapply(rowRanges(ods), function(rr) rr[1,])) # faster than range
  } else if(is(rowRanges(ods), "GRanges")){
    gr <- GRanges(seqnames = seqnames(rowRanges(ods)), 
                  IRanges(start = start(rowRanges(ods)), end = end(rowRanges(ods))),
                  symbol = names(rowRanges(ods)))
  }
  # seqlevelsStyle(gr) <- 'NCBI'
  # Add values to granges
  if(value %in% c('pvalue', 'pValue', 'pv')){
    gr$value <- -log10(assays(ods)$pValue[, sample])
    value <- '-log10(pvalue)'
  }
  if(value %in% c('padjust', 'pAdjust', 'pa')){
    gr$value <- -log10(assays(ods)$padjust[, sample])
    value <- '-log10(padjust)'
  }
  if(value %in% c('zscore', 'zScore'))
    gr$value <- assays(ods)$zScore[, sample]
  if(value %in% c('l2fc', 'L2FC', 'log2fc'))
    gr$value <- assays(ods)$l2fc[, sample]
  if(value == 'counts')
    gr$value <- log10(assays(ods)$counts[, sample]+1)
  gr$aberrant <- aberrant(ods)[,sample]
  # Sort granges for plot
  seqlevels(gr) <- gsub('chr','',seqlevels(gr))
  gr <- sortSeqlevels(gr)
  gr <- gr[gr@seqnames == chrom_to_plot]
  gr <- sort(gr)
  gr <- as.data.table(gr)
  max_value <- ceiling( max(abs(gr$value)))
  # add_chromlens <- add_chromlens[add_chromlens$CONTIG == '4',]
  plotdata <- gr %>%
     filter(!is.na(seqnames), 
            seqnames != "MT") %>%
  
     mutate(START_mbp = start / 1e06)
  vline_pos <- plotdata[aberrant==TRUE, START_mbp]
  # #sample_frac(0.01) # sampling 1% of points is detailed enough

  p <- ggplot(data = plotdata, aes(x = START_mbp, y = value, color = aberrant)) +
    geom_vline(xintercept=vline_pos, alpha=0.5, color='firebrick') + 
    geom_point(alpha = 0.5, size = 0.5) +
    labs(y = value) + 
    ylim(-max_value, max_value) + 
    scale_color_manual(values = color_chr)
  
  return(list(p, plotdata))
}



plotCumulativeChromCNV <- function(cnv, chrom_to_plot, color_chr = c("black", "firebrick")){
  add_chromlens <- add_chromlens[add_chromlens$CONTIG == chrom_to_plot,]
  cnv_chrom <- cnv[cnv$CONTIG == chrom_to_plot,]
  plotdata <- cnv_chrom %>%
    filter(!is.na(CONTIG), 
           CONTIG != "MT") %>%
    left_join(add_chromlens, by = c("CONTIG" = "CONTIG")) %>%
    mutate(START_mbp = START / 1e06,
           START_mbp_cum = START_mbp + addpos) #%>%
  #sample_frac(0.01) # sampling 1% of points is detailed enough
  
  # Plot: Coverage + chrom rectangles + chrom labels
  p <- ggplot() +
    geom_point(data = plotdata,
               aes(x = START_mbp_cum, y = LOG2_COPY_RATIO),
               alpha = 0.5, size = 0.5) 
  
  return(p)
}


plotChromCNV <- function(cnv, chrom_to_plot, color_chr = c("black", "firebrick")){
  cnv_chrom <- cnv[cnv$CONTIG == chrom_to_plot,]
  plotdata <- cnv_chrom %>%
    filter(!is.na(CONTIG), 
           CONTIG != "MT") %>%
    mutate(START_mbp = START / 1e06)
  #sample_frac(0.01) # sampling 1% of points is detailed enough
  
  # Plot: Coverage + chrom rectangles + chrom labels
  p <- ggplot() +
    geom_point(data = plotdata,
               aes(x = START_mbp, y = COPY_RATIO),
               alpha = 0.5, size = 0.5, color='darkgrey') 
  
  return(p)
}

fisher_test <- function(test_gene, all_gene, test_list) {
  non_test_gene <- all_gene[!all_gene %in% test_gene]
  
  contingency_df <- data.frame(
    "test_gene" = c(sum(test_gene %in% test_list), 
                    sum(!test_gene %in% test_list)),
    "non_test_gene" = c(sum(non_test_gene %in% test_list), 
                        sum(!non_test_gene %in% test_list)),
    row.names = c("leu", "non_leu"), stringsAsFactors = FALSE
  )
  
  fisher.test(contingency_df)
}
