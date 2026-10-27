#!/usr/bin/env Rscript
## build_record.R -- after a build, save what it used and what it made, and what changed
## since a previous build. Writes <out>/_build_record/:
##   inputs.csv            the input list this build used (copy)
##   outputs_md5.csv       checksum and size of every output file
##   summary.csv           headline counts (loci, genes per tier)
##   changes_since_previous.txt   inputs, outputs and counts that differ from <prev_out>
## Rscript build_record.R <out> <inputs.csv> [prev_out]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); out <- a[1]; inp <- a[2]; prev <- if (length(a) > 2) a[3] else ''
rd <- file.path(out, '_build_record'); dir.create(rd, showWarnings = FALSE)
invisible(file.copy(inp, file.path(rd, "inputs.csv"), overwrite = TRUE))
fs <- list.files(out, full.names = TRUE); fs <- fs[!dir.exists(fs)]
o <- data.table(file = basename(fs), size_bytes = file.size(fs), md5 = unname(tools::md5sum(fs)))
fwrite(o, file.path(rd, 'outputs_md5.csv'))
s <- list()
lf <- file.path(out, 'AD_loci_unified_cs95orColocs.csv.gz')
if (file.exists(lf)) { l <- fread(lf); s[['AD loci before p<1e-5 filter']] <- uniqueN(l$ADlocus) }
vf <- file.path(out, 'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz')
if (file.exists(vf)) { v <- fread(vf); s[['AD loci (final, p<1e-5)']] <- uniqueN(v$ADlocus); s[['AD loci variants (final)']] <- uniqueN(v$variant_ID) }
tf <- file.path(out, 'gene_tier_assignment.csv')
if (file.exists(tf)) { t <- fread(tf); s[['genes (any tier)']] <- uniqueN(t$gene_ID); for (k in sort(unique(t$tier))) s[[paste('genes', k)]] <- uniqueN(t[tier == k]$gene_ID) }
S <- data.table(quantity = names(s), value = unlist(s)); fwrite(S, file.path(rd, 'summary.csv'))
print(S)
if (nzchar(prev) && dir.exists(file.path(prev, '_build_record'))) {
  pr <- file.path(prev, '_build_record'); L <- c(sprintf('Changes from %s to %s', basename(prev), basename(out)), '')
  pi <- fread(file.path(pr, 'inputs.csv'), colClasses = 'character'); ci <- fread(inp, colClasses = 'character')
  k <- function(d) paste(d$Path, d$md5)
  L <- c(L, 'Inputs added or changed:', paste0('  ', ci$Path[!k(ci) %in% k(pi)]), 'Inputs removed or changed (previous version):', paste0('  ', pi$Path[!k(pi) %in% k(ci)]), '')
  ps <- fread(file.path(pr, 'summary.csv')); m <- merge(ps, S, by = 'quantity', all = TRUE, suffixes = c('_prev', '_now'))
  L <- c(L, 'Counts (previous -> now):', sprintf('  %s: %s -> %s', m$quantity, m$value_prev, m$value_now), '')
  po <- fread(file.path(pr, 'outputs_md5.csv')); mo <- merge(po, o, by = 'file', all = TRUE, suffixes = c('_prev', '_now'))
  L <- c(L, 'Output files that differ:', paste0('  ', mo[is.na(md5_prev) | is.na(md5_now) | md5_prev != md5_now]$file))
  writeLines(L, file.path(rd, 'changes_since_previous.txt')); cat(L, sep = '\n')
} else cat('[build_record] no previous build record to compare with\n')
