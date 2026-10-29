## Pairwise LD table used by the AD loci build for set extension (staging gwas_variants_cor0.5.csv.gz).
## For every variant in <variants.txt>, lists all reference-panel variants within 1 Mb with |r| > 0.5,
## using the ADSP R4 EUR LD sketch through pecotmr::loadLdMatrix. Each variant is looked up in a 3 Mb window
## that has it in the middle 1 Mb, so every partner within 1 Mb is in the window.
## Usage: Rscript scripts/compute_gwas_ld_table.R <variants.txt> <sketch ld_meta_file.tsv> <out.csv.gz> [workers] [max distance, default 1e6] ['listed': partners only from variants.txt]
suppressPackageStartupMessages({library(pecotmr); library(data.table); library(parallel)})
a <- commandArgs(TRUE); if (length(a) < 3) stop('usage: compute_gwas_ld_table.R <variants.txt> <ld_meta_file.tsv> <out.csv.gz> [workers]')
nw <- if (length(a) >= 4) as.integer(a[4]) else 1L
D <- if (length(a) >= 5) as.numeric(a[5]) else 1e6        ## max distance between the two variants of a pair
only_listed <- length(a) >= 6 && a[6] == 'listed'          ## keep only partners that are in <variants.txt>
v <- data.table(variant_id = unique(readLines(a[1])))
v[, chr := sub(':.*', '', variant_id)][, pos := as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', variant_id))]
v <- v[!is.na(pos)][, bin := floor((pos - 1) / 1e6)]
W <- unique(v[, .(chr, bin)])
cat('[ld] variants', nrow(v), 'windows', nrow(W), 'workers', nw, '\n')
## 'listed' mode (fast): partners come only from <variants.txt>, so load each 1 Mb bin once with returnGenotype=TRUE,
## keep the listed columns, and compute r between listed variants per chromosome (same r as @correlation; checked to 1e-14).
if (only_listed) {
  B <- unique(v[, .(chr, bin)])
  ld1 <- function(i) {
    w <- B[i]; reg <- sprintf('%s:%d-%d', w$chr, as.integer(w$bin * 1e6 + 1), as.integer(w$bin * 1e6 + 1e6))
    q <- v[chr == w$chr & bin == w$bin]$variant_id
    g <- tryCatch(getGenotypes(loadLdMatrix(a[2], region = reg, returnGenotype = TRUE)), error = function(err) NULL)
    if (is.null(g)) return(list(G = NULL, missing = q, failed = reg))
    have <- q[q %in% colnames(g)]
    list(G = g[, have, drop = FALSE], missing = setdiff(q, have), failed = NULL)
  }
  res <- mclapply(seq_len(nrow(B)), ld1, mc.cores = nw, mc.preschedule = FALSE)
  bad <- which(sapply(res, function(x) inherits(x, 'try-error') || is.null(x)))
  if (length(bad)) stop('worker errors in bins: ', paste(head(bad, 10), collapse = ','))
  T <- rbindlist(lapply(unique(B$chr), function(ch) {
    gl <- Filter(Negate(is.null), lapply(res[B$chr == ch], `[[`, 'G')); gl <- gl[sapply(gl, ncol) > 0]
    if (!length(gl)) return(NULL)
    if (!all(sapply(gl, function(x) identical(rownames(x), rownames(gl[[1]]))))) stop('sample order differs across bins on ', ch)
    G <- do.call(cbind, gl); R <- suppressWarnings(cor(G)); p <- as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', colnames(G)))
    k <- which(abs(R) > 0.5 & abs(outer(p, p, '-')) <= D, arr.ind = TRUE)
    data.table(variant_id = colnames(G)[k[, 1]], variant2 = colnames(G)[k[, 2]], r = R[k])
  }))
  if (!nrow(T)) T <- data.table(variant_id = character(), variant2 = character(), r = numeric())
  T[, region := v$chr[match(variant_id, v$variant_id)]]
  T[, `:=`(variant_ID = variant_id, variant_ID2 = variant2)]
  fwrite(T, a[3])
  miss <- unlist(lapply(res, `[[`, 'missing')); fail <- unlist(lapply(res, `[[`, 'failed'))
  writeLines(miss, sub('\\.csv\\.gz$', '.not_in_panel.txt', a[3]))
  cat('[ld] pairs', nrow(T), 'variants with pairs', uniqueN(T$variant_id), 'not in panel', length(miss), 'failed bins', length(fail), '\n')
  if (length(fail)) { cat('[ld] FAILED bins:', fail, '\n'); quit(status = 1) }
  quit(status = 0)
}
one <- function(i) {
  w <- W[i]; s <- max(1, w$bin * 1e6 - D + 1); e <- w$bin * 1e6 + 1e6 + D
  reg <- sprintf('%s:%d-%d', w$chr, as.integer(s), as.integer(e))
  m <- tryCatch(loadLdMatrix(a[2], region = reg)@correlation, error = function(err) NULL)
  q <- v[chr == w$chr & bin == w$bin]$variant_id
  if (is.null(m)) return(list(pairs = NULL, missing = q, failed = reg))
  pn <- as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', colnames(m)))
  have <- q[q %in% rownames(m)]
  out <- rbindlist(lapply(have, function(x) {
    px <- as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', x)); r <- m[x, ]
    k <- which(abs(r) > 0.5 & abs(pn - px) <= D & (!only_listed | colnames(m) %in% v$variant_id))
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
