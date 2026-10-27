#!/usr/bin/env Rscript
## drop_events.R -- remove named outcomes from a ColocBoost export (one row per variant)
## without touching anything else: each listed event is taken out of event_ID and its
## coef entry is removed; edited rows left with only AD outcomes (or none) are dropped.
## Rows that were not edited are written back unchanged.
## Used for agreed manual fixes, e.g. TMEM106B glyco signals gp_8491/gp_8492 that the
## Oct 16 glycoQTL rerun no longer has (Gao/Ru, Slack Oct 26).
## Rscript drop_events.R <in.bed> <out.bed> <region_ID> <event1> [event2 ...]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) >= 4)
x <- fread(a[1]); reg <- a[3]; drop <- a[-(1:3)]
hit <- which(x$region_ID == reg & sapply(strsplit(x$event_ID, ';'), function(e) any(trimws(e) %in% drop)))
for (i in hit) {
  e <- trimws(strsplit(x$event_ID[i], ';')[[1]]); cf <- trimws(strsplit(as.character(x$coef[i]), ';')[[1]])
  k <- !(e %in% drop); x$event_ID[i] <- paste(e[k], collapse = '; ')
  if (length(cf) == length(e)) x$coef[i] <- paste(cf[k], collapse = ';')
}
left <- sapply(strsplit(x$event_ID, ';'), function(e) sum(!grepl('^AD_', trimws(e))))
n0 <- nrow(x); keep <- rep(TRUE, n0); keep[hit] <- left[hit] > 0; x <- x[keep]
fwrite(x, a[2], sep = '\t', quote = FALSE)
message(sprintf('[drop_events] %s: %d rows edited, %d rows dropped (no xQTL outcome left), %d -> %d rows', reg, length(hit), n0 - nrow(x), n0, nrow(x)))
