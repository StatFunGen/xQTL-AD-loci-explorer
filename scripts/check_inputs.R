#!/usr/bin/env Rscript
## check_inputs.R -- verify, before a build reads anything, that every file in the
## build's input list is present and byte-identical to the recorded checksum, and
## (when the config is given) that every file the build will read is in the list.
## Folders (no single checksum) must exist. Any missing, changed, unlisted or
## unchecksummed file stops the build.
## Rscript check_inputs.R <inputs.csv> <project_root> [config_dir] [build_script]
##   config_dir    the config folder the build will read (its metadata_analysis.csv)
##   build_script  build_AD_locus_table.R, for the fixed input paths written in it
## Staging rows (staging/... or staging_gene_prio/...) are checked where the build
## reads them: $AD_LOCI_STAGING/<file> when AD_LOCI_STAGING is set.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) %in% 2:4)
x <- fread(a[1], colClasses = 'character'); root <- a[2]
cfg <- if (length(a) >= 3) a[3] else ''; bscript <- if (length(a) >= 4) a[4] else ''
## the shared input list masks one collaborator's user name; the real name is kept outside the repo
col <- Sys.getenv('AD_LOCI_COLLABORATOR', unset = tryCatch(readLines('~/.ad_loci_collaborator', 1), error = function(e) ''))
unmask <- function(p) gsub('<collaborator>', col, p, fixed = TRUE)
stg <- Sys.getenv('AD_LOCI_STAGING')
is_stg <- grepl('^staging(_gene_prio)?/', x$Path)
x[, p := file.path(root, unmask(Path))]
if (nzchar(stg)) x[is_stg, p := file.path(stg, basename(Path))]
x[, kind := fifelse(dir.exists(p), 'folder', fifelse(file.exists(p), 'file', 'missing'))]
f <- x[kind == 'file']
f[, now := unname(tools::md5sum(p))]
nomd5 <- f[is.na(md5) | md5 == '']
bad <- f[!(is.na(md5) | md5 == '') & (is.na(now) | now != md5)]
miss <- x[kind == 'missing']

## files the build will read that the input list does not name
unlisted <- character(0); gone <- character(0)
if (nzchar(cfg)) {
  m <- fread(file.path(cfg, 'metadata_analysis.csv'), colClasses = 'character')
  ## keep only the rows the build's loaders read (build_AD_locus_table.R):
  ## - interaction rows are not read: the build uses the precomputed staging summaries
  ## - LR rows are read only for sQTL
  ## - ColocBoost rows are read only for AD x xQTL and AD meta colocalization
  m <- m[!Method %in% c('msex interaction', 'APOE interaction') &
         !(Method == 'LR' & `Data Type` != 'sQTL') &
         !(Method == 'ColocBoost' & !Modality %in% c('AD_xQTL_colocalization', 'AD_meta_colocalization'))]
  need <- unique(unmask(m$Path[!is.na(m$Path) & m$Path != '']))
  if (nzchar(bscript)) {
    src <- readLines(bscript); src <- src[!grepl('^\\s*#', src)]
    ## the build replaces the registry's QR path with a fixed one
    qr <- regmatches(src, regexpr("Method=='QR',Path:='[^']+'", src))
    if (length(qr)) { need <- need[!need %in% m$Path[m$Method == 'QR']]; need <- c(need, sub(".*Path:='([^']+)'", '\\1', qr)) }
    ## fixed single-file inputs assigned at the top of the build (name=ROOTP('...'))
    fx <- regmatches(src, regexpr("^[A-Za-z0-9_.]+ *= *ROOTP\\('[^']+'\\)", src))
    fx <- sub(".*ROOTP\\('([^']+)'\\).*", '\\1', fx)
    need <- c(need, fx[!grepl('/$', fx)])
    ## per-chromosome files the build lists from a folder
    cv <- fx[grepl('cv2f', fx) & grepl('/$', fx)]
    if (!length(cv)) stop('check_inputs: could not find the cv2f folder in the build script; update this check')
    for (d in cv) need <- c(need, file.path(sub('/$', '', d), list.files(file.path(root, d), pattern = '\\.cv2f$')))
  }
  ## every dependency must be in the input list, whether or not it exists on disk
  ondisk <- file.exists(file.path(root, need)) | dir.exists(file.path(root, need))
  listed <- unmask(x$Path)
  inlist <- need %in% listed | vapply(need, function(n) any(startsWith(n, paste0(sub('/$', '', listed[x$kind == 'folder']), '/'))), logical(1))
  unlisted <- need[!inlist]
  gone <- need[!inlist & !ondisk]   # neither listed nor on disk (listed but missing files are reported under MISSING)
  ## the two staging tables the build reads
  if (nzchar(stg)) for (s in c('res_APOE_interaction_summ.csv.gz', 'res_msex_interaction_summ.csv.gz'))
    if (!any(is_stg & basename(x$Path) == s)) unlisted <- c(unlisted, file.path('staging', s))
}

cat(sprintf('[check_inputs] %d rows: %d files checked, %d match, %d changed, %d without md5, %d folders, %d missing, %d read by the build but not listed\n',
  nrow(x), nrow(f), nrow(f) - nrow(bad) - nrow(nomd5), nrow(bad), nrow(nomd5), sum(x$kind == 'folder'), nrow(miss), length(unlisted)))
if (nrow(bad)) cat('CHANGED (recorded -> now):\n', paste0('  ', bad$Path, '  ', bad$md5, ' -> ', bad$now, collapse = '\n'), '\n')
if (nrow(nomd5)) cat('NO RECORDED MD5:\n', paste0('  ', nomd5$Path, collapse = '\n'), '\n')
if (nrow(miss)) cat('MISSING:\n', paste0('  ', miss$p, collapse = '\n'), '\n')
if (length(unlisted)) cat('NOT IN THE INPUT LIST:\n', paste0('  ', unlisted, ifelse(unlisted %in% gone, '  (also missing on disk)', ''), collapse = '\n'), '\n')
if (nrow(bad) || nrow(nomd5) || nrow(miss) || length(unlisted)) {
  cat('[check_inputs] STOP: inputs do not match the input list. Update the list (with the data owner\'s confirmation) or restore the file.\n'); quit(status = 1) }
cat('[check_inputs] OK\n')
