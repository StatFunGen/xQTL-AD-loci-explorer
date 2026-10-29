#!/usr/bin/env Rscript
## compare_builds.R -- compare the scientific tables of two build output folders.
## Tables are compared after sorting rows and columns; numbers within a tolerance.
## Workbooks are not byte-compared (they carry metadata such as timestamps).
## Rscript compare_builds.R <out_a> <out_b> [tolerance=1e-8]
## Exits with status 1 when any table differs or is missing.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); A <- a[1]; B <- a[2]; tol <- if (length(a) > 2) as.numeric(a[3]) else 1e-8
find <- function(d, stem) { f <- list.files(d, pattern = paste0('^', stem), full.names = TRUE); f[1] }
stems <- c('AD_loci_unified_cs95orColocs_Pval1e5_variant_level', 'AD_loci_unified_cs95orColocs(_[0-9]{8})?\\.csv', 'gene_tier_assignment',
           'res_AD_variants_xQTL', 'xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5(_[0-9]{8})?\\.csv')
norm <- function(x) { x <- copy(x); setcolorder(x, sort(names(x)))
  for (c in names(x)) if (is.character(x[[c]])) set(x, j = c, value = fifelse(x[[c]] == '', NA_character_, x[[c]]))
  setorderv(x, names(x), na.last = TRUE); x }
ok_all <- TRUE
for (s in stems) {
  fa <- find(A, s); fb <- find(B, s)
  if (is.na(fa) || is.na(fb)) { cat(sprintf('%-60s missing in %s\n', s, if (is.na(fa)) 'A' else 'B')); ok_all <- FALSE; next }
  x <- norm(fread(fa)); y <- norm(fread(fb))
  same_cols <- identical(names(x), names(y))
  eq <- same_cols && nrow(x) == nrow(y) && isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))
  cat(sprintf('%-60s rows %8d vs %8d  cols %s  %s\n', basename(fa), nrow(x), nrow(y), if (same_cols) 'same' else 'DIFFER', if (eq) 'IDENTICAL' else 'DIFFERENT'))
  if (!eq) { ok_all <- FALSE
    if (!same_cols) cat('   only in A:', setdiff(names(x), names(y)), ' only in B:', setdiff(names(y), names(x)), '\n')
    else if (nrow(x) == nrow(y)) { d <- all.equal(x, y, tolerance = tol, check.attributes = FALSE); cat('  ', head(d, 5), sep = '\n   ') } }
}
for (d in c(A, B)) { t <- find(d, 'gene_tier_assignment'); if (!is.na(t)) { z <- fread(t); cat(basename(d), ': genes', uniqueN(z$gene_ID), ' by tier:', paste(names(table(z$tier)), table(z$tier), sep = '=', collapse = ' '), '\n') } }
cat(if (ok_all) 'RESULT: all compared tables identical\n' else 'RESULT: differences found (see above)\n')
if (!ok_all) quit(status = 1)
