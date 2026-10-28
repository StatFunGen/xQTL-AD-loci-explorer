## LD among the AD variants, computed when the build runs from the ADSP R4 EUR LD sketch with pecotmr (no precomputed LD table).
## Each 1 Mb bin that holds AD variants is loaded once (loadLdMatrix(returnGenotype = TRUE)) and only the AD variant columns are
## kept; r is then computed among the AD variants of each extension region. Same r as loadLdMatrix(...)@correlation
## (checked to 1e-14), without building full-panel matrices. Returns every pair with |r| > min_abs_r within a region.
ld_pairs_sketch <- function(ids, meta, regions, workers = 1L, min_abs_r = 0.5) {
  suppressPackageStartupMessages({library(pecotmr); library(data.table); library(parallel)})
  v <- data.table(variant_id = unique(ids))
  v[, chr := sub(':.*', '', variant_id)][, pos := as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', variant_id))]
  v <- v[!is.na(pos)]
  rg <- data.table(region = unique(regions))[, c('chr', 's', 'e') := tstrsplit(region, '_', fixed = TRUE)][, `:=`(s = as.numeric(s), e = as.numeric(e))]
  v[rg, region := i.region, on = .(chr, pos >= s, pos <= e)]
  if (anyNA(v$region)) stop('[LD] ', sum(is.na(v$region)), ' AD variants fall outside the extension regions')
  v[, bin := floor((pos - 1) / 1e6)]
  B <- unique(v[, .(chr, bin)]); t0 <- Sys.time()
  cat('[LD] sketch', meta, ':', nrow(v), 'AD variants,', nrow(rg), 'regions,', nrow(B), '1 Mb bins,', workers, 'workers\n')
  res <- mclapply(seq_len(nrow(B)), function(i) {
    w <- B[i]; reg <- sprintf('%s:%d-%d', w$chr, as.integer(w$bin * 1e6 + 1), as.integer(w$bin * 1e6 + 1e6))
    q <- v[chr == w$chr & bin == w$bin]$variant_id
    g <- tryCatch(getGenotypes(loadLdMatrix(meta, region = reg, returnGenotype = TRUE)), error = function(err) NULL)
    if (is.null(g)) return(list(G = NULL, missing = q, failed = reg))
    have <- q[q %in% colnames(g)]
    list(G = g[, have, drop = FALSE], missing = setdiff(q, have), failed = NULL)
  }, mc.cores = workers, mc.preschedule = FALSE)
  if (any(vapply(res, function(x) inherits(x, 'try-error') || is.null(x), TRUE))) stop('[LD] worker error while loading the sketch')
  fail <- unlist(lapply(res, `[[`, 'failed')); if (length(fail)) stop('[LD] sketch load failed for ', paste(fail, collapse = ' '))
  B[, key := paste(chr, bin)]; v[, key := paste(chr, bin)]
  out <- rbindlist(lapply(rg$region, function(r) {
    gl <- lapply(res[B$key %in% v[region == r, unique(key)]], `[[`, 'G'); gl <- gl[vapply(gl, ncol, 1L) > 0]
    if (!length(gl)) return(NULL)
    if (!all(vapply(gl, function(x) identical(rownames(x), rownames(gl[[1]])), TRUE))) stop('[LD] sample order differs across sketch bins in ', r)
    G <- do.call(cbind, gl); R <- suppressWarnings(cor(G)); k <- which(abs(R) > min_abs_r, arr.ind = TRUE)
    data.table(variant_id = colnames(G)[k[, 1]], variant2 = colnames(G)[k[, 2]], r = R[k], region = r)
  }))
  if (!nrow(out)) out <- data.table(variant_id = character(), variant2 = character(), r = numeric(), region = character())
  out[, `:=`(variant_ID = variant_id, variant_ID2 = variant2)]
  miss <- unlist(lapply(res, `[[`, 'missing'))
  cat('[LD] pairs |r| >', min_abs_r, ':', nrow(out), '; AD variants not in the sketch:', length(miss), '; minutes:', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), '\n')
  attr(out, 'not_in_panel') <- miss
  out
}
