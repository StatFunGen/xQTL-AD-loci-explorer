#!/usr/bin/env Rscript
## drop_cos.R -- remove whole colocalization sets (CoS) from a ColocBoost export (one
## row per variant): every CoS on <chr> whose outcomes include <gene> as a cis gene.
## Used for the agreed TMEM106B fix in Xuewei's TransMap x AD export (Slack, Oct 27).
## Rscript drop_cos.R <in.bed> <out.bed> <chr> <gene_id>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) == 4)
b <- fread(a[1])
d <- unique(b[`#chr` == a[3] & grepl(paste0('cis_gene_group[0-9]+_[A-Za-z]+_', a[4]), event_ID)]$cos_ID)
y <- b[!cos_ID %in% d]
fwrite(y, a[2], sep = '\t', quote = FALSE)
message(sprintf('[drop_cos] removed %d CoS (%d rows); %d -> %d rows', length(d), nrow(b) - nrow(y), nrow(b), nrow(y)))
