#!/usr/bin/env Rscript
## build_record.R -- the one record of a build. After a build, writes <out>/_build_record/:
##   run.csv                       run id, code commit (and whether the checkout had local
##                                 changes), input list and its sha256, config sha256, R and
##                                 package versions, host, user, time
##   inputs.csv                    the input list this build used (copy)
##   config/ + config/MD5SUMS      the config files the build read (copy)
##   environment.txt               R session information
##   outputs_md5.csv               checksum and size of every output file
##   summary.csv                   headline counts (loci, genes per tier)
##   changes_since_previous.txt    inputs, outputs and counts that differ from <prev_out>
## The build itself adds registry_runtime.csv (the registry with run-time paths), and the
## SCC job adds check_inputs.txt and validation.txt.
## Rscript build_record.R <out> <inputs.csv> [prev_out]
## Reads AD_LOCI_CONFIG (config folder) and AD_LOCI_COMMIT / AD_LOCI_DIRTY (set by the job).
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); out <- a[1]; inp <- a[2]; prev <- if (length(a) > 2) a[3] else ''
if (!dir.exists(out) || !file.exists(inp)) stop('build_record: missing output folder or input list')
rd <- file.path(out, '_build_record'); dir.create(rd, showWarnings = FALSE); if (!dir.exists(rd)) stop('build_record: cannot create ', rd)
sha <- function(f) { r <- tryCatch(system2('sha256sum', shQuote(f), stdout = TRUE), error = function(e) NA_character_); sub(' .*', '', r[1]) }

if (!file.copy(inp, file.path(rd, 'inputs.csv'), overwrite = TRUE)) stop('build_record: could not copy the input list')
cfg <- Sys.getenv('AD_LOCI_CONFIG'); cd <- file.path(rd, 'config')
if (nzchar(cfg) && dir.exists(cfg)) {
  dir.create(cd, showWarnings = FALSE); cf <- list.files(cfg, full.names = TRUE); cf <- cf[!dir.exists(cf)]
  if (!all(file.copy(cf, cd, overwrite = TRUE))) stop('build_record: could not copy the config')
  writeLines(paste(unname(tools::md5sum(file.path(cd, basename(cf)))), basename(cf), sep = '  '), file.path(cd, 'MD5SUMS'))
}
writeLines(capture.output(sessionInfo()), file.path(rd, 'environment.txt'))
pk <- c('data.table', 'stringr', 'dplyr', 'tidyr', 'openxlsx', 'parallel')
pv <- vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), character(1))
run <- data.table(field = c('run_id', 'code_commit', 'code_has_local_changes', 'input_list', 'input_list_sha256',
                            'config_dir', 'config_sha256', 'R_version', paste0('pkg_', pk), 'host', 'user', 'recorded_at'),
  value = c(basename(normalizePath(out)), Sys.getenv('AD_LOCI_COMMIT', NA), Sys.getenv('AD_LOCI_DIRTY', NA),
            normalizePath(inp), sha(inp), cfg, if (file.exists(file.path(cd, 'MD5SUMS'))) sha(file.path(cd, 'MD5SUMS')) else NA,
            R.version.string, pv, Sys.info()[['nodename']], Sys.info()[['user']], format(Sys.time(), '%Y-%m-%dT%H:%M:%S%z')))
fwrite(run, file.path(rd, 'run.csv'))

fs <- list.files(out, full.names = TRUE); fs <- fs[!dir.exists(fs)]
o <- data.table(file = basename(fs), size_bytes = file.size(fs), md5 = unname(tools::md5sum(fs)))
fwrite(o, file.path(rd, 'outputs_md5.csv'))
s <- list()
find1 <- function(stem) { f <- list.files(out, pattern = paste0('^', stem), full.names = TRUE); if (length(f)) f[1] else '' }
lf <- find1('AD_loci_unified_cs95orColocs(_[0-9]{8})?\\.csv'); if (nzchar(lf)) { l <- fread(lf); s[['AD loci before p<1e-5 filter']] <- uniqueN(l$ADlocus) }
vf <- find1('AD_loci_unified_cs95orColocs_Pval1e5_variant_level'); if (nzchar(vf)) { v <- fread(vf); s[['AD loci (final, p<1e-5)']] <- uniqueN(v$ADlocus); s[['AD loci variants (final)']] <- uniqueN(v$variant_ID) }
tf <- find1('gene_tier_assignment'); if (nzchar(tf)) { t <- fread(tf); s[['genes (any tier)']] <- uniqueN(t$gene_ID); for (k in sort(unique(t$tier))) s[[paste('genes', k)]] <- uniqueN(t[tier == k]$gene_ID) }
S <- data.table(quantity = names(s), value = unlist(s)); fwrite(S, file.path(rd, 'summary.csv'))
print(S)
if (nzchar(prev) && dir.exists(file.path(prev, '_build_record'))) {
  pr <- file.path(prev, '_build_record'); L <- c(sprintf('Changes from %s to %s', basename(prev), basename(out)), '')
  pi <- fread(file.path(pr, 'inputs.csv'), colClasses = 'character'); ci <- fread(inp, colClasses = 'character')
  k <- function(d) paste(d$Path, d$md5)
  L <- c(L, 'Inputs added or changed:', paste0('  ', ci$Path[!k(ci) %in% k(pi)]), 'Inputs removed or changed (previous version):', paste0('  ', pi$Path[!k(pi) %in% k(ci)]), '')
  if (file.exists(file.path(pr, 'summary.csv'))) { ps <- fread(file.path(pr, 'summary.csv')); m <- merge(ps, S, by = 'quantity', all = TRUE, suffixes = c('_prev', '_now'))
    L <- c(L, 'Counts (previous -> now):', sprintf('  %s: %s -> %s', m$quantity, m$value_prev, m$value_now), '') }
  if (file.exists(file.path(pr, 'outputs_md5.csv'))) { po <- fread(file.path(pr, 'outputs_md5.csv')); mo <- merge(po, o, by = 'file', all = TRUE, suffixes = c('_prev', '_now'))
    L <- c(L, 'Output files that differ:', paste0('  ', mo[is.na(md5_prev) | is.na(md5_now) | md5_prev != md5_now]$file)) }
  writeLines(L, file.path(rd, 'changes_since_previous.txt')); cat(L, sep = '\n')
} else cat('[build_record] no previous build record to compare with\n')
