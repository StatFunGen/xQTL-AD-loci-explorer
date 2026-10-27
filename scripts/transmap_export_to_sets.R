#!/usr/bin/env Rscript
## transmap_export_to_sets.R -- collapse a TransMap x AD ColocBoost export (one row per
## variant) into one row per colocalization set, in the column layout of
## transmap_AD_target_results.tsv, so convert_transmap_AD_target.R (cis genes) and
## convert_transmap_AD_target_trans.R (trans genes/programs/biomarkers) can read it unchanged.
## Rscript transmap_export_to_sets.R <export.bed> <out.tsv>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); x <- fread(a[1])
y <- x[, .(focal_outcome = grep('^AD_', trimws(strsplit(event_ID[1], ';')[[1]]), value = TRUE)[1],
  colocalized_outcomes = event_ID[1], purity = NA_real_, cos_npc = cos_npc[1], min_npc_outcome = min_npc_outcome[1],
  colocalized_variables = paste(variant_ID, collapse = '; '), colocalized_variables_vcp = paste(vcp, collapse = '; '),
  region = region_ID[1]), by = .(cos_id = cos_ID)]
fwrite(y, a[2], sep = '\t')
message(sprintf('[tm sets] %d variant rows -> %d sets, %d regions', nrow(x), nrow(y), uniqueN(y$region)))
