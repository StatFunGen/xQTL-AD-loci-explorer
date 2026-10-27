#!/usr/bin/env Rscript
## dedup_against_export.R -- remove from the tpQTL x AD ColocBoost export every ROSMAP
## eQTL outcome (gene + tissue) that is already in the main ROSMAP AD-xQTL export, so the
## same colocalization is not counted twice. tpQTL (factor) outcomes and eQTL outcomes not
## in the main export are kept. The matching coef entry is removed too; rows left with
## only AD outcomes are dropped.
## The two files name the same eQTL differently, matched here on gene + tissue:
##   tpQTL export:  <T>_DeJager_eQTL_<gene>          (T = AC, DLPFC, PCC, Ast, Exc, ...)
##   main export:   ROSMAP_<T>_<gene> (AC/DLPFC/PCC)  or  DeJager_<T>_<gene> (cell types)
## Rscript dedup_against_export.R <in.bed> <out.bed> <main_ROSMAP_export.bed>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) == 3)
x <- fread(a[1]); m <- fread(a[3])
me <- unique(trimws(unlist(strsplit(m$event_ID, ';'))))
seen <- sub('^(ROSMAP|DeJager)_([A-Za-z]+)_(ENSG[0-9]+)$', '\\2 \\3', grep('^(ROSMAP_(AC|DLPFC|PCC)|DeJager_[A-Za-z]+)_ENSG[0-9]+$', me, value = TRUE))
key <- function(e) sub('^([A-Za-z]+)_DeJager_eQTL_(ENSG[0-9]+)$', '\\1 \\2', e)
nrm <- 0L
for (i in seq_len(nrow(x))) {
  e <- trimws(strsplit(x$event_ID[i], ';')[[1]]); cf <- trimws(strsplit(as.character(x$coef[i]), ';')[[1]])
  drop <- grepl('_DeJager_eQTL_ENSG', e) & key(e) %in% seen; nrm <- nrm + sum(drop)
  x$event_ID[i] <- paste(e[!drop], collapse = '; '); if (length(cf) == length(e)) x$coef[i] <- paste(cf[!drop], collapse = ';')
}
n0 <- nrow(x); x <- x[sapply(strsplit(event_ID, ';'), function(e) any(!grepl('^AD_', trimws(e))))]
fwrite(x, a[2], sep = '\t', quote = FALSE)
message(sprintf('[dedup] %d outcome entries removed; %d -> %d rows', nrm, n0, nrow(x)))
