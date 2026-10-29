#!/usr/bin/env Rscript
## convert_glyco_AD_coloc.R -- turn Ru's glycoQTL x AD ColocBoost summary tables
## (Oct 16 patch, *_AD_summary_pv1e4.tsv: one row per colocalization set) into the
## ColocBoost export layout the build reads (one row per variant).
##
## Rules:
## - Only sets that include an AD outcome are written (sets of glyco/pQTL alone are not AD colocalization).
## - Outcome names are kept as given (DLPFC_Klein_gpQTL_*_gp_*, DLPFC_Bennett_pQTL_*); region_ID is the set's gene.
## - Gene x context pairs already present in an existing ColocBoost export are dropped, so nothing is counted twice.
## - coef is not in the input and is written as NA.
##
## Rscript convert_glyco_AD_coloc.R <out.bed> <existing_export.bed> <in1.tsv> [in2.tsv ...]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) >= 3)
out <- a[1]; prev <- fread(a[2]); ins <- a[-(1:2)]
ctx <- function(e) sub('_(ENSG[0-9]+|gp_[0-9]+)$', '', e)
seen <- unique(prev[, .(e = trimws(unlist(strsplit(event_ID, ';')))), by = region_ID][!grepl('^AD_', e), paste(region_ID, ctx(e))])
rows <- list(); nd <- 0L
for (f in ins) {
  x <- fread(f); tag <- sub('_summary_pv1e4.tsv$', '', basename(f))
  for (r in seq_len(nrow(x))) {
    ev <- trimws(strsplit(x$colocalized_outcomes[r], ';')[[1]]); ev <- ev[ev != '']
    ad <- ev[grepl('^AD_', ev)]; if (!length(ad)) next
    q <- ev[!grepl('^AD_', ev)]; g <- x$gene[r]
    keep <- !(paste(g, ctx(q)) %in% seen); nd <- nd + sum(!keep); q <- q[keep]
    if (!length(q)) next
    v <- trimws(strsplit(x$colocalized_variables[r], ';')[[1]]); w <- as.numeric(trimws(strsplit(as.character(x$colocalized_variables_vcp[r]), ';')[[1]]))
    stopifnot(length(v) == length(w))
    e <- c(q, ad)
    rows[[length(rows) + 1]] <- data.table(variant_ID = v, region_ID = g, event_ID = paste(e, collapse = '; '),
      cos_ID = paste0(g, ':', tag, ':', x$cos_id[r]), vcp = w, cos_npc = x$cos_npc[r], min_npc_outcome = x$min_npc_outcome[r],
      coef = paste(rep('NA', length(e)), collapse = ';'))
  }
}
y <- rbindlist(rows); sp <- tstrsplit(y$variant_ID, ':', fixed = TRUE)
y[, `:=`(`#chr` = sp[[1]], start = as.integer(sp[[2]]), end = as.integer(sp[[2]]), a1 = sp[[4]], a2 = sp[[3]])]
setcolorder(y, c('#chr','start','end','a1','a2','variant_ID','region_ID','event_ID','cos_ID','vcp','cos_npc','min_npc_outcome','coef'))
fwrite(y, out, sep = '\t', quote = FALSE)
message(sprintf('[glycoAD] %d sets, %d rows, %d genes, %d duplicate gene x context pairs dropped', uniqueN(y$cos_ID), nrow(y), uniqueN(y$region_ID), nd))
print(unique(y[, .(region_ID, event_ID)]))
