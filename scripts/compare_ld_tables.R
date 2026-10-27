## Compare two LD tables (old, new) on the variants present in both, for pairs within 1 Mb.
## Usage: Rscript scripts/compare_ld_tables.R <old.csv.gz> <new.csv.gz>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
rd <- function(f) { d <- unique(fread(f, select = c('variant_id','variant2','r')), by = c('variant_id','variant2'))
  d[, dist := abs(as.numeric(sub('^[^:]+:([0-9]+):.*','\\1',variant_id)) - as.numeric(sub('^[^:]+:([0-9]+):.*','\\1',variant2)))][] }
o <- rd(a[1]); n <- rd(a[2])
cat('old pairs', nrow(o), 'beyond 1Mb', o[dist > 1e6, .N], '| new pairs', nrow(n), '\n')
vs <- intersect(unique(o$variant_id), unique(n$variant_id))
cat('variants: old', uniqueN(o$variant_id), 'new', uniqueN(n$variant_id), 'both', length(vs), '\n')
k <- merge(o[variant_id %in% vs & dist <= 1e6], n[variant_id %in% vs], by = c('variant_id','variant2'), all = TRUE, suffixes = c('.old','.new'))
cat('pairs on shared variants: both', k[!is.na(r.old) & !is.na(r.new), .N], 'only old', k[is.na(r.new), .N], 'only new', k[is.na(r.old), .N],
    'max |r diff|', k[, max(abs(r.old - r.new), na.rm = TRUE)], '\n')
if (k[is.na(r.new) | is.na(r.old), .N]) print(head(k[is.na(r.new) | is.na(r.old)], 10))
