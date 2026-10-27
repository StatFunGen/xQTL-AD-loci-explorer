## Pairwise LD table used by the AD loci build for set extension (staging gwas_variants_cor0.5.csv.gz).
## For every variant in <variants.txt>, lists all reference-panel variants within 1 Mb with |r| > 0.5,
## using the ADSP R4 EUR LD sketch through pecotmr::loadLdMatrix. Each variant is looked up in a 3 Mb window
## that has it in the middle 1 Mb, so every partner within 1 Mb is in the window.
## Usage: Rscript scripts/compute_gwas_ld_table.R <variants.txt> <sketch ld_meta_file.tsv> <out.csv.gz> [workers]
suppressPackageStartupMessages({library(pecotmr); library(data.table); library(parallel)})
a <- commandArgs(TRUE); if (length(a) < 3) stop('usage: compute_gwas_ld_table.R <variants.txt> <ld_meta_file.tsv> <out.csv.gz> [workers]')
nw <- if (length(a) >= 4) as.integer(a[4]) else 1L
v <- data.table(variant_id = unique(readLines(a[1])))
v[, chr := sub(':.*', '', variant_id)][, pos := as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', variant_id))]
v <- v[!is.na(pos)][, bin := floor((pos - 1) / 1e6)]
W <- unique(v[, .(chr, bin)])
cat('[ld] variants', nrow(v), 'windows', nrow(W), 'workers', nw, '\n')
one <- function(i) {
  w <- W[i]; s <- max(1, w$bin * 1e6 - 1e6 + 1); e <- w$bin * 1e6 + 2e6
  reg <- sprintf('%s:%d-%d', w$chr, as.integer(s), as.integer(e))
  m <- tryCatch(loadLdMatrix(a[2], region = reg)@correlation, error = function(err) NULL)
  q <- v[chr == w$chr & bin == w$bin]$variant_id
  if (is.null(m)) return(list(pairs = NULL, missing = q, failed = reg))
  pn <- as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', colnames(m)))
  have <- q[q %in% rownames(m)]
  out <- rbindlist(lapply(have, function(x) {
    px <- as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', x)); r <- m[x, ]
    k <- which(abs(r) > 0.5 & abs(pn - px) <= 1e6)
    data.table(variant_id = x, variant2 = colnames(m)[k], r = unname(r[k]))
  }))
  if (nrow(out)) out[, region := reg]
  list(pairs = out, missing = setdiff(q, have), failed = NULL)
}
res <- mclapply(seq_len(nrow(W)), one, mc.cores = nw, mc.preschedule = FALSE)
bad <- which(sapply(res, function(x) inherits(x, 'try-error') || is.null(x)))
if (length(bad)) stop('worker errors in windows: ', paste(head(bad, 10), collapse = ','))
T <- rbindlist(lapply(res, `[[`, 'pairs'))
T[, `:=`(variant_ID = variant_id, variant_ID2 = variant2)]
fwrite(T, a[3])
miss <- unlist(lapply(res, `[[`, 'missing')); fail <- unlist(lapply(res, `[[`, 'failed'))
writeLines(miss, sub('\\.csv\\.gz$', '.not_in_panel.txt', a[3]))
cat('[ld] pairs', nrow(T), 'variants with pairs', uniqueN(T$variant_id), 'not in panel', length(miss), 'failed windows', length(fail), '\n')
if (length(fail)) { cat('[ld] FAILED windows:', fail, '\n'); quit(status = 1) }
