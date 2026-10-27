#!/usr/bin/env Rscript
## check_inputs.R -- verify, before a build reads anything, that every file in the
## build's input list is present and byte-identical to the recorded checksum.
## Folders (no single checksum) must exist. Any missing or changed file stops the build.
## Rscript check_inputs.R <inputs.csv> <project_root>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) == 2)
x <- fread(a[1], colClasses = 'character'); root <- a[2]
## the shared input list masks one collaborator's user name; the real name is kept outside the repo
col <- Sys.getenv('AD_LOCI_COLLABORATOR', unset = tryCatch(readLines('~/.ad_loci_collaborator', 1), error = function(e) ''))
x[, p := file.path(root, gsub('<collaborator>', col, Path, fixed = TRUE))]
x[grepl('^staging/', Path), p := file.path(Sys.getenv('AD_LOCI_OUT'), sub('^staging/', '', Path))]
x[, kind := fifelse(dir.exists(p), 'folder', fifelse(file.exists(p), 'file', 'missing'))]
f <- x[kind == 'file' & md5 != '']
f[, now := unname(tools::md5sum(p))]
bad <- f[now != md5]
miss <- x[kind == 'missing']
cat(sprintf('[check_inputs] %d rows: %d files checked, %d match, %d changed, %d folders, %d missing, %d files without md5\n',
  nrow(x), nrow(f), nrow(f) - nrow(bad), nrow(bad), sum(x$kind == 'folder'), nrow(miss), sum(x$kind == 'file' & x$md5 == '')))
if (nrow(bad)) cat('CHANGED (recorded -> now):\n', paste0('  ', bad$Path, '  ', bad$md5, ' -> ', bad$now, collapse = '\n'), '\n')
if (nrow(miss)) cat('MISSING:\n', paste0('  ', miss$p, collapse = '\n'), '\n')
if (nrow(bad) || nrow(miss)) { cat('[check_inputs] STOP: inputs do not match the input list. Update the list (with the data owner\'s confirmation) or restore the file.\n'); quit(status = 1) }
cat('[check_inputs] OK\n')
