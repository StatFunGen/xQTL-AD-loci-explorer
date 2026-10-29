# Fast LD extension used by build_AD_locus_table.R (same output as the original loop).
# Call only after loading the existing coordinate helpers (seqid/start/end).
extend_ad_sets <- function(res_gwf, res_cosad, regions, variants_corsf,
                           workers = as.integer(Sys.getenv('NSLOTS', '1')),
                           reuse_pairs = TRUE, timing_dir = NULL) {
  stopifnot(requireNamespace('data.table', quietly = TRUE))
  if (length(workers) != 1L || is.na(workers) || workers < 1L) stop('workers must be positive')
  slots <- suppressWarnings(as.integer(Sys.getenv('NSLOTS', as.character(workers))))
  if (is.na(slots) || slots < 1L) stop('NSLOTS must be positive')
  workers <- min(workers, slots, max(1L, length(regions)))
  if (workers > 1L && .Platform$OS.type == 'windows') stop('Parallel region execution requires Unix')
  utc <- function() format(Sys.time(), '%Y-%m-%dT%H:%M:%OS6Z', tz = 'UTC')
  elapsed <- function() unname(proc.time()[['elapsed']])
  overall_start <- elapsed(); partition_start <- elapsed()
  previous_threads <- data.table::getDTthreads()
  data.table::setDTthreads(1L)
  on.exit(data.table::setDTthreads(previous_threads), add = TRUE)

  ## GWAS sets are assigned to LD blocks by variant position, like the colocalization sets below. Matching the
  ## fine-mapping block name to the LD block name skipped every GWAS set whose fine-mapping block differs from the
  ## LD reference block (most blocks since the Sep 2026 GWAS export), so those sets were never extended.
  gchr <- seqid(res_gwf$variant_ID); gpos <- pos(res_gwf$variant_ID)
  initial <- lapply(regions, function(region_id) {
    x <- data.table::melt(res_gwf[gchr == seqid(region_id) & gpos >= start(region_id) & gpos <= end(region_id)],
      measure.vars = c('cs_coverage_0.95', 'cs_coverage_0.7', 'cs_coverage_0.5'),
      variable.name = 'coverage', value.name = 'cs_num')[cs_num != 0]
    y <- res_cosad[chr == seqid(region_id) & pos >= start(region_id) & pos <= end(region_id)]
    data.table::rbindlist(list(x, y), fill = TRUE)
  })
  membership <- data.table::rbindlist(lapply(seq_along(initial), function(k) {
    if (!nrow(initial[[k]])) return(NULL)
    data.table::data.table(locuscontext_id = unique(initial[[k]]$locuscontext_id), region_id = regions[k])
  }))
  shared_sets <- if (nrow(membership)) membership[, .(regions = data.table::uniqueN(region_id)), by = locuscontext_id][regions > 1L] else data.table::data.table()
  # Some real set IDs span LD regions. The legacy algorithm nevertheless reads
  # only gwas_new[region == r] on every pass: these are region-scoped copies,
  # not shared mutable state. Keep every cache and membership vector inside its
  # region, and report the shared labels rather than merging their states.
  if (nrow(shared_sets)) message(nrow(shared_sets), ' set IDs span regions; retaining legacy region-scoped state')
  # Sets only acquire variants from other sets in their own region. Therefore this
  # initial union contains every variant that any later pair lookup can request.
  correlations <- variants_corsf[, .(variant_ID, variant_ID2, r)]
  data.table::setindexv(correlations, 'variant_ID')
  ld <- lapply(initial, function(x) {
    universe <- unique(x$variant_ID)
    z <- correlations[.(universe), on = 'variant_ID', nomatch = 0L][
      .(universe), on = 'variant_ID2', nomatch = 0L]
    data.table::setkeyv(z, c('variant_ID', 'variant_ID2'))
    z
  })
  groups <- lapply(initial, function(x) {
    data.table::rbindlist(lapply(unique(x$locuscontext_id), function(id) unique(x[locuscontext_id == id, .(locuscontext_id, source)])))
  })
  caches <- lapply(groups, function(g) list(values = matrix(FALSE, nrow(g), length(unique(g$locuscontext_id))),
                                           valid = matrix(FALSE, nrow(g), length(unique(g$locuscontext_id))),
                                           a = NULL, b = NULL))
  partition_seconds <- elapsed() - partition_start
  states <- initial; pass <- 0L; pass_times <- list(); region_times <- list()
  work <- vapply(seq_along(regions), function(k) as.numeric(nrow(groups[[k]])) * data.table::uniqueN(initial[[k]]$locuscontext_id), numeric(1))
  heavy <- if (workers > 1L) which(work >= 0.1 * sum(work)) else integer(0)
  if (length(heavy)) message('extension: ', length(heavy), ' heavy region(s) run with set-level parallelism: ', paste(regions[heavy], collapse = ', '))

  region_step <- function(k, pass_id, inner = 1L) {
    data.table::setDTthreads(1L)
    begun <- utc(); clock <- elapsed()
    x <- states[[k]]; g <- groups[[k]]; cache <- caches[[k]]
    loci <- unique(initial[[k]]$locuscontext_id)
    if (!length(loci)) {
      # Legacy empty regions alternate: rbindlist(NULL)[, region := r] creates
      # a placeholder on odd passes; grouping that NA set produces no row on
      # the next pass. Preserve this lifecycle, including a zero-row pass 4.
      empty_state <- data.table::data.table(locuscontext_id = NA_character_, source = NA_character_,
          variant_ID = NA_character_, original_cs = NA, region = regions[k])
      if (pass_id %% 2L == 0L) empty_state <- empty_state[0L]
      return(list(state = empty_state, cache = cache,
        timing = data.table::data.table(pass = pass_id, region = regions[k], start = begun,
          end = utc(), seconds = elapsed() - clock, pair_checks = 0L, cache_hits = 0L, added = 0L)))
    }
    b <- lapply(loci, function(id) x[locuscontext_id == id]$variant_ID)
    a_rows <- lapply(seq_len(nrow(g)), function(j) {
      id <- g$locuscontext_id[j]; src <- g$source[j]
      which(x$locuscontext_id == id & if (is.na(src)) is.na(x$source) else !is.na(x$source) & x$source == src)
    })
    a <- lapply(a_rows, function(idx) x$variant_ID[idx])
    # Pass 1 may have duplicates from melt(). Comparing the full vectors also
    # detects their collapse in pass 2, even when unique membership is unchanged.
    changed_a <- vapply(seq_along(a), function(j) !identical(a[[j]], cache$a[[j]]), logical(1))
    changed_b <- vapply(seq_along(b), function(j) !identical(b[[j]], cache$b[[j]]), logical(1))
    # Each set j reads the same previous-pass state and fixed LD and writes only its
    # own output and its own cache row, so sets can run in parallel (used for the few
    # regions that dominate run time, e.g. the chr17 MAPT region). Results are
    # reassembled in set order, so the output is identical to the serial loop.
    step_j <- function(j) {
      id <- g$locuscontext_id[j]; vars_a <- a[[j]]
      flags <- if (pass_id == 1L) vars_a else x$variant_ID[a_rows[[j]]][x$original_cs[a_rows[[j]]]]
      accept <- rep(FALSE, length(loci)); vals <- cache$values[j, ]; vld <- cache$valid[j, ]; ch <- 0L; hi <- 0L
      for (h in seq_along(loci)) {
        if (loci[h] == id) next
        if (!reuse_pairs || !vld[h] || changed_a[j] || changed_b[h]) {
          pairs <- ld[[k]][.(unique(vars_a)), on = 'variant_ID', nomatch = 0L][
            .(unique(b[[h]])), on = 'variant_ID2', nomatch = 0L]
          ok <- nrow(pairs) > 0L && isTRUE(any(abs(pairs$r) > 0.8) &&
            nrow(pairs) == length(unique(vars_a)) * length(unique(b[[h]])))  ## unique counts (review Q4)
          vals[h] <- ok; vld[h] <- TRUE; ch <- ch + 1L
        } else hi <- hi + 1L
        accept[h] <- vals[h]
      }
      # Preserve the original table order of B's rows, including duplicates.
      other <- x[x$locuscontext_id != id]
      to_add <- other$variant_ID[accept[match(other$locuscontext_id, loci)]]
      merged <- union(vars_a, to_add)
      list(out = data.table::data.table(locuscontext_id = id, source = g$source[j],
        variant_ID = merged, original_cs = merged %in% flags),
        vals = vals, vld = vld, checks = ch, hits = hi, added = length(setdiff(to_add, vars_a)))
    }
    res <- if (inner > 1L && nrow(g) > 1L) parallel::mclapply(seq_len(nrow(g)), step_j,
      mc.cores = inner, mc.preschedule = FALSE, mc.set.seed = FALSE) else lapply(seq_len(nrow(g)), step_j)
    # a killed worker returns NULL rather than an error, so require a full result from every set
    if (length(res) != nrow(g) || any(vapply(res, function(z) is.null(z) || inherits(z, 'try-error') || !is.list(z) || is.null(z$out), logical(1))))
      stop('Set worker failed in region ', regions[k])
    out <- lapply(res, '[[', 'out')
    for (j in seq_along(res)) { cache$values[j, ] <- res[[j]]$vals; cache$valid[j, ] <- res[[j]]$vld }
    checks <- sum(vapply(res, function(r) r$checks, integer(1)))
    hits <- sum(vapply(res, function(r) r$hits, integer(1)))
    added <- sum(vapply(res, function(r) r$added, integer(1)))
    cache$a <- a; cache$b <- b
    result <- data.table::rbindlist(out, fill = TRUE)[, region := regions[k]]
    list(state = result, cache = cache, timing = data.table::data.table(pass = pass_id,
      region = regions[k], start = begun, end = utc(), seconds = elapsed() - clock,
      pair_checks = checks, cache_hits = hits, added = added))
  }
  repeat {
    pass <- pass + 1L; begun <- utc(); clock <- elapsed()
    # states and caches are read-only for this entire synchronous pass. Publish
    # the next state only after all workers return, retaining region order.
    values <- vector('list', length(regions))
    # the heaviest regions run one at a time with their sets spread over all workers;
    # the rest run as before, one region per worker
    for (kk in heavy) values[[kk]] <- region_step(kk, pass_id = pass, inner = workers)
    light <- setdiff(seq_along(regions), heavy)
    values[light] <- if (workers == 1L) lapply(light, region_step, pass_id = pass) else
      parallel::mclapply(light, region_step, pass_id = pass,
        mc.cores = workers, mc.preschedule = TRUE, mc.set.seed = FALSE)
    if (any(vapply(values, function(z) is.null(z) || inherits(z, 'try-error') || !is.list(z) || is.null(z$state), logical(1)))) stop('Region worker failed')
    states <- lapply(values, '[[', 'state'); caches <- lapply(values, '[[', 'cache')
    rt <- data.table::rbindlist(lapply(values, '[[', 'timing'))
    region_times[[pass]] <- rt
    pass_times[[pass]] <- data.table::data.table(pass = pass, start = begun, end = utc(),
      seconds = elapsed() - clock, added = sum(rt$added), pair_checks = sum(rt$pair_checks),
      cache_hits = sum(rt$cache_hits),
      empty_placeholder_rows = sum(vapply(states, function(z) sum(is.na(z$variant_ID) & is.na(z$locuscontext_id)), numeric(1))))
    if (!is.null(timing_dir)) {
      dir.create(timing_dir, recursive = TRUE, showWarnings = FALSE)
      data.table::fwrite(data.table::rbindlist(pass_times), file.path(timing_dir, 'pass-times.csv'))
      data.table::fwrite(data.table::rbindlist(region_times), file.path(timing_dir, 'region-times.csv'))
    }
    # The legacy algorithm always executes pass 2, even if pass 1 added nothing.
    if (pass >= 2L && sum(rt$added) == 0L) break
  }
  list(result = data.table::rbindlist(states, fill = TRUE), passes = pass,
       pass_times = data.table::rbindlist(pass_times), region_times = data.table::rbindlist(region_times),
       partition_seconds = partition_seconds, total_seconds = elapsed() - overall_start,
       workers = workers, set_region_membership = membership, shared_set_regions = shared_sets)
}
