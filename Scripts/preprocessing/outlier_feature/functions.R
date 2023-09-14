##### function #####
nr_score_is_topK_by_zScore <- function(k, ranks, zScore, zScore_direction){
  
  switch(zScore_direction,
         up = nr_is_top_socre <- rowSums(ranks <= k & zScore > 0),
         dn = nr_is_top_socre <- rowSums(ranks <= k & zScore < 0)
  )
  
  dt <- data.table(geneID = rownames(ranks),
                   nrIsTopScore = nr_is_top_socre,
                   k = k)
}

nr_score_is_pvalCutoff_by_zScore <- function(cutoff, padj, gene_dt, zScore, zScore_direction){
  
  switch(zScore_direction,
         up = nr_sig <- rowSums(padj <= cutoff & zScore > 0),
         dn = nr_sig <- rowSums(padj <= cutoff & zScore < 0)
  )
  
  data.table(geneID = rownames(zScore),
             nrSignificant = nr_sig,
             cutoff = cutoff)
}

nr_score_is_zScore_by_padj <- function(cutoff, gene_dt, zScore_mtx, padj_mtx){
  
  zScore_up_vec <- rowSums((zScore_mtx > cutoff) & (padj_mtx < 0.05))
  zScore_dn_vec <- rowSums((zScore_mtx < -cutoff) & (padj_mtx < 0.05))
  
  zScore_dt <- data.table(geneID = sapply(rownames(zScore_mtx), function(x){strsplit(x, "[.]")[[1]][1]}),
                          zScore_up = zScore_up_vec,
                          zScore_dn = zScore_dn_vec)
  setnames(zScore_dt, c('zScore_up', 'zScore_dn'), paste("zScore", cutoff, c('up', 'dn'), sep = "_"))
  
  zScore_dt <- left_join(gene_dt, zScore_dt, by = "geneID")
  zScore_dt_wide <- melt(zScore_dt, id.vars=c('geneID', 'geneSymbol'))
  zScore_dt_wide[is.na(value), value := 0]
}

getRanks <- function(rc_scores){
  apply(rc_scores, 1,
        function(scores, ties.method){
          rank(-scores, ties.method = ties.method)
        }, ties.method= "min")
}

# create gene_dt
createGeneTable <- function(ods, gencode){
  gene_dt <- data.table(gene_id = rownames(rowData(ods)))
  gene_dt <- left_join(gene_dt, gencode[, .(gene_id, gene_name)] %>% unique(),
                       by = "gene_id") %>% as.data.table()
  gene_dt[, geneID := sapply(gene_dt[, gene_id],
                             function(x){strsplit(x, "[.]")[[1]][1]})]
  setnames(gene_dt, "gene_name", "geneSymbol")
  gene_dt[, gene_id := NULL]
  
  return(gene_dt)
}

# create top K with zScore directioin
createTopKbyZScore <- function(ods, k, gene_dt, filepath){
  pval <- assay(ods, 'pValue') %>% t()
  zscore <- zScore(ods)
  
  ranks <- getRanks(-pval)
  rownames(ranks) <- rowData(ods)$geneID
  
  res_list_zScore_up <- lapply(k, nr_score_is_topK_by_zScore, ranks=ranks, zScore=zscore, zScore_direction="up")
  res_list_zScore_dn <- lapply(k, nr_score_is_topK_by_zScore, ranks=ranks, zScore=zscore, zScore_direction="dn")
  
  dt_long <- rbindlist(res_list_zScore_up)
  dt_wide <- dcast(dt_long, geneID ~ paste0("k=", k, "_up"), value.var = "nrIsTopScore")
  dt_wide <- gene_dt[dt_wide, on="geneID"]
  dt_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  dt_wide_up <- copy(dt_wide)
  
  dt_long <- rbindlist(res_list_zScore_dn)
  dt_wide <- dcast(dt_long, geneID ~ paste0("k=", k, "_dn"), value.var = "nrIsTopScore")
  dt_wide <- gene_dt[dt_wide, on="geneID"]
  dt_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  dt_wide_dn <- copy(dt_wide)
  
  dt_wide <- merge(dt_wide_up, dt_wide_dn, by = c("geneSymbol", "geneID"))
  setnames(dt_wide, 'geneID', 'gene_id')
  dt_wide[, geneSymbol := NULL]
  
  prefix <- sub("^([^-]*-[^-]*).*", "\\1", basename(filepath))
  colnames(dt_wide)[-1] <- paste(prefix, colnames(dt_wide)[-1], sep="-")
  
  fwrite(dt_wide,
         file = filepath,
         sep = "\t", quote = F)
  
  return(dt_wide)
}

# create pval cutoff
createPadjCutoffbyZScore <- function(ods, cutoffs, gene_dt, filepath){
  padj <- assay(ods, 'padjust') 
  zscore <- zScore(ods) 
  
  padj_cutoffs_up <- lapply(cutoffs, nr_score_is_pvalCutoff_by_zScore, padj = padj, gene_dt = gene_dt, zScore=zscore, zScore_direction="up")
  padj_cutoffs_dn <- lapply(cutoffs, nr_score_is_pvalCutoff_by_zScore, padj = padj, gene_dt = gene_dt, zScore=zscore, zScore_direction="dn")
  
  padj_res <- rbindlist(padj_cutoffs_up)
  padj_res_wide <- dcast(padj_res, geneID ~ paste0("cutoff=", cutoff, "_up"),
                         value.var = "nrSignificant")
  padj_res_wide <- gene_dt[padj_res_wide, on="geneID"]
  padj_res_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  padj_res_wide_up <- copy(padj_res_wide)
  
  padj_res <- rbindlist(padj_cutoffs_dn)
  padj_res_wide <- dcast(padj_res, geneID ~ paste0("cutoff=", cutoff, "_dn"),
                         value.var = "nrSignificant")
  padj_res_wide <- gene_dt[padj_res_wide, on="geneID"]
  padj_res_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  padj_res_wide_dn <- copy(padj_res_wide)
  
  padj_res_wide <- merge(padj_res_wide_up, padj_res_wide_dn, by = c("geneSymbol", "geneID"))
  setnames(padj_res_wide, 'geneID', 'gene_id')
  padj_res_wide[, geneSymbol := NULL]
  
  prefix <- sub("^([^-]*-[^-]*).*", "\\1", basename(filepath))
  colnames(padj_res_wide)[-1] <- paste(prefix, colnames(padj_res_wide)[-1], sep="-")
  
  fwrite(padj_res_wide,
         file = filepath,
         sep = "\t", quote = F)
  
  return(padj_res_wide)
}

# create zScore count
createZscoreSig <- function(ods, cutoffs, gene_dt, filepath){
  zScore_mtx <- zScore(ods)
  padj_mtx <- padj(ods)
  
  zScore_dt_list <- lapply(cutoffs, nr_score_is_zScore_by_padj, gene_dt=gene_dt, zScore_mtx=zScore_mtx, padj_mtx = padj_mtx)
  zScore_dt_long <- rbindlist(zScore_dt_list)
  zScore_dt <- dcast(zScore_dt_long, geneID+geneSymbol ~ variable, value.var = "value")
  
  setnames(zScore_dt, 'geneID', 'gene_id')
  zScore_dt[, geneSymbol := NULL]
  
  prefix <- sub("^([^-]*-[^-]*).*", "\\1", basename(filepath))
  colnames(zScore_dt)[-1] <- paste(prefix, colnames(zScore_dt)[-1], sep="-")
  
  #+ write table
  fwrite(zScore_dt, 
         file = filepath,
         sep = "\t", quote = F)
  
  return(zScore_dt)
}

get_pval_gene <- function(fds, gencode){
  pval_junc <- pVals(fds) %>% as.data.table()
  gr <- rowRanges(fds, type=type)
  pval_junc[is.na(pval_junc)] <- 1
  pval_junc[, hgnc_symbol := gr$hgnc_symbol]
  pval_junc <- merge(pval_junc, gencode[, .(gene_id, gene_name)], 
                     by.x='hgnc_symbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  pval_junc[, hgnc_symbol := NULL]  
  pval_junc <- pval_junc[!is.na(gene_id), ]
  pval_junc <- separate_rows(pval_junc, 'gene_id', sep = ";") %>% as.data.table()
  pval_junc_melt <- melt(pval_junc, id.vars='gene_id')
  pval_dcast <- dcast(pval_junc_melt, gene_id ~ variable, fun=min)
  row_symbol <- pval_dcast[, gene_id]
  pval <- pval_dcast[, gene_id := NULL] %>% as.matrix()
  rownames(pval) <- row_symbol
  return(pval)
}

get_padj_gene <- function(fds, gencode){
  padj_junc <- padjVals(fds) %>% as.data.table()
  gr <- rowRanges(fds, type=type)
  padj_junc[is.na(padj_junc)] <- 1
  padj_junc[, hgnc_symbol := gr$hgnc_symbol]
  padj_junc <- merge(padj_junc, gencode[, .(gene_id, gene_name)], 
                     by.x='hgnc_symbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  padj_junc[, hgnc_symbol := NULL]  
  padj_junc <- padj_junc[!is.na(gene_id), ]
  padj_junc <- separate_rows(padj_junc, 'gene_id', sep = ";") %>% as.data.table()
  padj_junc_melt <- melt(padj_junc, id.vars='gene_id')
  padj_dcast <- dcast(padj_junc_melt, gene_id ~ variable, fun=min)
  row_symbol <- padj_dcast[, gene_id]
  padj <- padj_dcast[, gene_id := NULL] %>% as.matrix()
  rownames(padj) <- row_symbol
  return(padj)
}

get_deltaJac_gene <- function(fds, gencode){
  deltaJac_junc <- deltaPsiValue(fds, type) %>% as.data.table()
  gr <- rowRanges(fds, type=type)
  deltaJac_junc[is.na(deltaJac_junc)] <- 0
  deltaJac_junc[, hgnc_symbol := gr$hgnc_symbol]
  deltaJac_junc <- merge(deltaJac_junc, gencode[, .(gene_id, gene_name)], 
                         by.x='hgnc_symbol', by.y='gene_name', all.x=TRUE, all.y=FALSE)
  deltaJac_junc[, hgnc_symbol := NULL] 
  deltaJac_junc <- deltaJac_junc[!is.na(gene_id), ]
  deltaJac_junc <- separate_rows(deltaJac_junc, 'gene_id', sep = ";") %>% as.data.table()
  deltaJac_junc_melt <- melt(deltaJac_junc, id.vars='gene_id')
  max_abs <- function(x){x[which.max(abs(x))]}
  deltaJac_junc_melt <- deltaJac_junc_melt[, max_abs(value), by=c('gene_id', 'variable')]
  deltaJac_dcast <- dcast(deltaJac_junc_melt, gene_id ~ variable, value.var='V1')
  row_symbol <- deltaJac_dcast[, gene_id]
  deltaJac <- deltaJac_dcast[, gene_id := NULL] %>% as.matrix()
  rownames(deltaJac) <- row_symbol
  return(deltaJac)
}


createGeneTableFds <- function(fds, gencode){
  gr <- rowRanges(fds, type=type)
  
  gene_dt <- data.table(gene_name = unique(gr$hgnc_symbol))
  gene_dt <- separate_rows(gene_dt,1,sep = ";") %>% as.data.table()
  gene_dt <- left_join(gene_dt, gencode[, .(gene_id, gene_name)] %>% unique(),
                       by = "gene_name") %>% as.data.table()
  gene_dt <- gene_dt[!is.na(gene_id),]
  gene_dt <- unique(gene_dt)
  gene_dt[, geneID := sapply(gene_dt[, gene_id],
                             function(x){strsplit(x, "[.]")[[1]][1]})]
  setnames(gene_dt, "gene_name", "geneSymbol")
  gene_dt[, gene_id := NULL]
  
  return(gene_dt)
}


createTopKbyDeltaJac <- function(fds, k, gene_dt, gencode, filepath){
  pval <- get_pval_gene(fds, gencode) %>% t()
  deltaJac <- get_deltaJac_gene(fds, gencode)
  
  ranks <- getRanks(-pval)
  rownames(ranks) <- rownames(deltaJac)
  
  res_list_deltaJac_up <- lapply(k, nr_score_is_topK_by_zScore, ranks=ranks, zScore=deltaJac, zScore_direction="up")
  res_list_deltaJac_dn <- lapply(k, nr_score_is_topK_by_zScore, ranks=ranks, zScore=deltaJac, zScore_direction="dn")
  
  dt_long <- rbindlist(res_list_deltaJac_up)
  dt_wide <- dcast(dt_long, geneID ~ paste0("k=", k, "_up"), value.var = "nrIsTopScore")
  dt_wide <- gene_dt[dt_wide, on="geneID"]
  dt_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  dt_wide_up <- copy(dt_wide)
  
  dt_long <- rbindlist(res_list_deltaJac_dn)
  dt_wide <- dcast(dt_long, geneID ~ paste0("k=", k, "_dn"), value.var = "nrIsTopScore")
  dt_wide <- gene_dt[dt_wide, on="geneID"]
  dt_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  dt_wide_dn <- copy(dt_wide)
  
  dt_wide <- merge(dt_wide_up, dt_wide_dn, by = c("geneSymbol", "geneID"))
  setnames(dt_wide, 'geneID', 'gene_id')
  dt_wide[, geneSymbol := NULL]
  
  prefix <- sub("^([^-]*-[^-]*).*", "\\1", basename(filepath))
  colnames(dt_wide)[-1] <- paste(prefix, colnames(dt_wide)[-1], sep="-")
  
  fwrite(dt_wide,
         file = filepath,
         sep = "\t", quote = F)
  
  return(dt_wide)
}

createPadjCutoffbyDeltaJac <- function(fds, cutoffs, gene_dt, gencode, filepath){
  padj <- get_padj_gene(fds, gencode)
  deltaJac <- get_deltaJac_gene(fds, gencode)
  
  padj_cutoffs_up <- lapply(cutoffs, nr_score_is_pvalCutoff_by_zScore, padj = padj, gene_dt = gene_dt, zScore=deltaJac, zScore_direction="up")
  padj_cutoffs_dn <- lapply(cutoffs, nr_score_is_pvalCutoff_by_zScore, padj = padj, gene_dt = gene_dt, zScore=deltaJac, zScore_direction="dn")
  
  padj_res <- rbindlist(padj_cutoffs_up)
  padj_res_wide <- dcast(padj_res, geneID ~ paste0("cutoff=", cutoff, "_up"),
                         value.var = "nrSignificant")
  padj_res_wide <- gene_dt[padj_res_wide, on="geneID"]
  padj_res_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  padj_res_wide_up <- copy(padj_res_wide)
  
  padj_res <- rbindlist(padj_cutoffs_dn)
  padj_res_wide <- dcast(padj_res, geneID ~ paste0("cutoff=", cutoff, "_dn"),
                         value.var = "nrSignificant")
  padj_res_wide <- gene_dt[padj_res_wide, on="geneID"]
  padj_res_wide[,geneID:=sapply(strsplit(geneID, ".", fixed=TRUE), '[', 1)]
  padj_res_wide_dn <- copy(padj_res_wide)
  
  padj_res_wide <- merge(padj_res_wide_up, padj_res_wide_dn, by = c("geneSymbol", "geneID"))
  setnames(padj_res_wide, 'geneID', 'gene_id')
  padj_res_wide[, geneSymbol := NULL]
  
  prefix <- sub("^([^-]*-[^-]*).*", "\\1", basename(filepath))
  colnames(padj_res_wide)[-1] <- paste(prefix, colnames(padj_res_wide)[-1], sep="-")
  
  fwrite(padj_res_wide,
         file = filepath,
         sep = "\t", quote = F)
  
  return(padj_res_wide)
}

createDeltaJacSig <- function(fds, cutoffs, gene_dt, gencode, filepath){
  padj_mtx <- get_padj_gene(fds, gencode)
  deltaJac_mtx <- get_deltaJac_gene(fds, gencode)
  
  zScore_dt_list <- lapply(cutoffs, nr_score_is_zScore_by_padj, gene_dt=gene_dt, zScore_mtx=deltaJac_mtx, padj_mtx = padj_mtx)
  zScore_dt_long <- rbindlist(zScore_dt_list)
  zScore_dt <- dcast(zScore_dt_long, geneID+geneSymbol ~ variable, value.var = "value")
  
  setnames(zScore_dt, 'geneID', 'gene_id')
  zScore_dt[, geneSymbol := NULL]
  
  prefix <- sub("^([^-]*-[^-]*).*", "\\1", basename(filepath))
  colnames(zScore_dt)[-1] <- paste(prefix, colnames(zScore_dt)[-1], sep="-")
  
  #+ write table
  fwrite(zScore_dt, 
         file = filepath,
         sep = "\t", quote = F)
  
  return(zScore_dt)
}