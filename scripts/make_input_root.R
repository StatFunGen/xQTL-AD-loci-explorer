#!/usr/bin/env Rscript
## make_input_root.R -- build a project folder of links laid out exactly as the build expects,
## pointing at wherever each input actually lives on this cluster. The build then runs unchanged
## with AD_LOCI_ROOT=<root>, and check_inputs.R confirms every file matches the input list.
##   Rscript scripts/make_input_root.R <inputs.csv> <root> [column] [prefix]
##   column  the input-list column holding each file's location on this cluster
##           (default columbia_path). Use Path with a prefix to link a copy that keeps the
##           project layout, e.g. on BU SCC: Path /restricted/projectnb/xqtl/jaempawi/xqtl/AD_loci_xQTL
## Staging tables are linked under <root>/staging_gene_prio/, where the job script looks for them.
## Only links are made; nothing is copied, and existing links in <root> are replaced.
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); stopifnot(length(a) >= 2)
x <- fread(a[1], colClasses = 'character'); root <- a[2]
col <- if (length(a) >= 3) a[3] else 'columbia_path'; pre <- if (length(a) >= 4) a[4] else ''
if (!col %in% names(x)) stop('no column ', col, ' in ', a[1])
## the shared input list masks one collaborator's user name; the real name is kept outside the repo
cl <- Sys.getenv('AD_LOCI_COLLABORATOR', unset = tryCatch(readLines('~/.ad_loci_collaborator', 1), error = function(e) ''))
unmask <- function(p) gsub('<collaborator>', cl, p, fixed = TRUE)
stg <- grepl('^staging(_gene_prio)?/', x$Path)
x[, dest := file.path(root, fifelse(stg, file.path('staging_gene_prio', basename(Path)), unmask(Path)))]
x[, src := sub(' \\(.*\\)$', '', unmask(get(col)))]  # drop a trailing note such as '(not yet verified)'
if (nzchar(pre)) x[, src := file.path(pre, fifelse(stg, file.path('staging_gene_prio', basename(Path)), src))]
x[, src := sub('/$', '', src)]; x[, dest := sub('/$', '', dest)]
x[, kind := fifelse(dir.exists(src), 'folder', fifelse(file.exists(src), 'file', 'missing'))]
## a file inside a linked folder is reached through that folder; never create paths inside a link
fold <- x[kind == 'folder']$dest
x[, inside := vapply(dest, function(d) { o <- fold[fold != d]; length(o) > 0 && any(startsWith(d, paste0(o, '/'))) }, logical(1))]
dir.create(root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, 'repro', 'tmp'), recursive = TRUE, showWarnings = FALSE)  # the build's scratch for large reads
n <- 0L
for (i in which(x$kind != 'missing' & !x$inside)) {
  d <- x$dest[i]; dir.create(dirname(d), recursive = TRUE, showWarnings = FALSE)
  if (nzchar(Sys.readlink(d))) file.remove(d) else if (file.exists(d)) stop('not a link, refusing to replace: ', d)
  if (!file.symlink(x$src[i], d)) stop('could not link ', d)
  n <- n + 1L
}
cat(sprintf('[make_input_root] %d inputs: %d linked, %d reached through a linked folder, %d missing\n',
  nrow(x), n, sum(x$inside & x$kind != 'missing'), sum(x$kind == 'missing')))
if (any(x$kind == 'missing')) { cat('MISSING on this cluster:\n', paste0('  ', x[kind == 'missing']$src, collapse = '\n'), '\n'); quit(status = 1) }
cat('[make_input_root] next: Rscript scripts/check_inputs.R <inputs.csv> ', root, ' config scripts/build_AD_locus_table.R  (with AD_LOCI_STAGING=', file.path(root, 'staging_gene_prio'), ')\n', sep = '')
