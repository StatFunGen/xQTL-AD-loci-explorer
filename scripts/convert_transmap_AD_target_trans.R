#!/usr/bin/env Rscript
## convert_transmap_AD_target_trans.R -- write the trans side of the transmap
## AD-target colocalization table (Xuewei, transmap_data) in the per-context
## top-loci layout the build reads, so trans genes, programs and biomarkers that
## colocalize with AD reach the trans columns.
##
## Rules:
## - One output file per outcome context: ROSMAP eQTL tissues (genes and local or
##   global programs) as <tissue>_DeJager_eQTL, biomarkers as <set>_transmap.
## - One row per variant x trans outcome; cs_coverage_0.95 is the set number,
##   PIP the variant's vcp, z left empty (not in the input).
## - cis outcomes are not written here (convert_transmap_AD_target.R handles them).
##
## Rscript convert_transmap_AD_target_trans.R <in.tsv> <out_dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); x <- fread(a[1]); out <- a[2]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
rows <- vector('list', nrow(x))
for (r in seq_len(nrow(x))) {
  v <- trimws(strsplit(x$colocalized_variables[r], ';')[[1]]); w <- as.numeric(trimws(strsplit(x$colocalized_variables_vcp[r], ';')[[1]]))
  stopifnot(length(v) == length(w))
  ev <- trimws(strsplit(x$colocalized_outcomes[r], ';')[[1]]); ev <- ev[grepl('_trans_', ev)]
  if (!length(ev)) next
  rows[[r]] <- CJ(i = seq_along(v), event_ID = ev, sorted = FALSE)[, `:=`(variant_ID = v[i], PIP = w[i], cs_coverage_0.95 = r, purity = x$purity[r], region_ID = x$region[r])][, i := NULL]
}
y <- rbindlist(rows)
set <- sub('_trans_.*$', '', y$event_ID)
y[, context := ifelse(grepl('_trans_biomarker_', event_ID), paste0(set, '_transmap'), paste0(set, '_DeJager_eQTL'))]
sp <- tstrsplit(y$variant_ID, ':', fixed = TRUE)
y[, `:=`(`#chr` = sp[[1]], start = as.integer(sp[[2]]), end = as.integer(sp[[2]]), a1 = sp[[4]], a2 = sp[[3]], z = NA_real_)]
for (cx in sort(unique(y$context))) {
  f <- file.path(out, paste0(cx, '.transmap_ADtarget.exported.toploci.bed.gz'))
  fwrite(y[context == cx, .(`#chr`, start, end, a1, a2, variant_ID, region_ID, event_ID, cs_coverage_0.95, purity, PIP, z)][order(`#chr`, start)], f, sep = '\t')
  message(sprintf('[tmAD trans] %s: %d rows, %d outcomes, %d sets', cx, sum(y$context == cx), uniqueN(y[context == cx]$event_ID), uniqueN(y[context == cx]$cs_coverage_0.95)))
}
