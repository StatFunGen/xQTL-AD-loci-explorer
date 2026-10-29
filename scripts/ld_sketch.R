## LD among the AD variants, computed when the build runs from the ADSP R4 EUR LD sketch with pecotmr (no precomputed LD table).
## Only the AD variants are read from each chromosome's pgen, with pecotmr's own genotype reader (handle -> variant
## index -> dosage -> stochastic-scale inversion, the steps of loadGenotypeRegion() but with the variant index applied
## before reading). r is then computed among the AD variants of each extension region. Same r as
## loadLdMatrix(...)@correlation (checked to 1e-15 on chr22). Returns every pair with |r| > min_abs_r and distance <= max_dist
## within a region.
.sketch_genotypes <- function(meta, chr, ids) {
  ns <- asNamespace('pecotmr'); m <- data.table::fread(meta)
  p <- m[[4]][sub('^chr', '', as.character(m[[1]])) == sub('^chr', '', chr)]
  if (length(p) != 1) stop('[LD] sketch metadata has ', length(p), ' entries for ', chr)
  pfx <- file.path(dirname(meta), p)
  h <- ns$.loadGenoHandle(pfx); si <- ns$getSnpInfo(h); id <- ns$.snpInfoToVariantInfo(si)$id
  idx <- which(id %in% ids)
  if (anyDuplicated(id[idx])) stop('[LD] duplicated variant IDs in the sketch on ', chr)
  if (!length(idx)) return(list(G = NULL, missing = ids))
  res <- list(X = ns$.dosageMatrix(h, idx, meanImpute = FALSE),
              variant_info = ns$.loadGenoAttachAfreq(h, ns$.snpInfoToVariantInfo(dplyr::slice(si, idx))))
  res <- ns$.loadGenoInvertStochastic(res, pfx, NULL, NULL)
  list(G = res$X, missing = setdiff(ids, colnames(res$X)))
}
ld_pairs_sketch <- function(ids, meta, regions, workers = 1L, min_abs_r = 0.5, max_dist = Inf) {
  suppressPackageStartupMessages({library(pecotmr); library(data.table); library(parallel)})
  v <- data.table(variant_id = unique(ids))
  v[, chr := sub(':.*', '', variant_id)][, pos := as.numeric(sub('^[^:]+:([0-9]+):.*', '\\1', variant_id))]
  v <- v[!is.na(pos)]
  rg <- data.table(region = unique(regions))[, c('chr', 's', 'e') := tstrsplit(region, '_', fixed = TRUE)][, `:=`(s = as.numeric(s), e = as.numeric(e))]
  v[rg, region := i.region, on = .(chr, pos >= s, pos <= e)]
  if (anyNA(v$region)) stop('[LD] ', sum(is.na(v$region)), ' AD variants fall outside the extension regions')
  chrs <- unique(v$chr); t0 <- Sys.time()
  cat('[LD] sketch', meta, ':', nrow(v), 'AD variants,', nrow(rg), 'regions,', length(chrs), 'chromosomes, max distance', max_dist, ',', workers, 'workers\n')
  res <- mclapply(chrs, function(ch) .sketch_genotypes(meta, ch, v[chr == ch]$variant_id), mc.cores = workers, mc.preschedule = FALSE)
  bad <- vapply(res, function(x) inherits(x, 'try-error') || is.null(x), TRUE)
  if (any(bad)) stop('[LD] sketch read failed for ', paste(chrs[bad], collapse = ' '), ': ', paste(unlist(res[bad]), collapse = ' '))
  names(res) <- chrs
  out <- rbindlist(lapply(rg$region, function(r) {
    G <- res[[rg[region == r]$chr]]$G; k <- intersect(v[region == r]$variant_id, colnames(G))
    if (!length(k)) return(NULL)
    R <- suppressWarnings(cor(G[, k, drop = FALSE])); p <- v[match(k, variant_id)]$pos
    x <- which(abs(R) > min_abs_r & abs(outer(p, p, '-')) <= max_dist, arr.ind = TRUE)
    data.table(variant_id = k[x[, 1]], variant2 = k[x[, 2]], r = R[x], region = r)
  }))
  if (!nrow(out)) out <- data.table(variant_id = character(), variant2 = character(), r = numeric(), region = character())
  out[, `:=`(variant_ID = variant_id, variant_ID2 = variant2)]
  miss <- unlist(lapply(res, `[[`, 'missing'))
  cat('[LD] pairs |r| >', min_abs_r, ':', nrow(out), '; AD variants not in the sketch:', length(miss), '; minutes:', round(as.numeric(difftime(Sys.time(), t0, units = 'mins')), 1), '\n')
  attr(out, 'not_in_panel') <- miss
  out
}
