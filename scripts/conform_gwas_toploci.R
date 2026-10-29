## Convert an AD GWAS fine-mapping export (<study>.gwas.top_loci.bed.gz) into the
## column layout the AD loci build reads (same rows, same values):
##   #chr (no 'chr' prefix), start = pos - 1, end = pos, region_ID and gene_ID = the
##   fine-mapping block, event_ID = the study name, empty fields written as NA.
## Usage: Rscript scripts/conform_gwas_toploci.R <in.top_loci.bed.gz> <out.top_loci.bed.gz>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
if (length(a) != 2) stop('usage: conform_gwas_toploci.R <in> <out>')
study <- sub('\\.gwas\\.top_loci\\.bed\\.gz$', '', basename(a[1]))
d <- fread(a[1], sep = '\t', colClasses = 'character', na.strings = NULL)
d[, (names(d)) := lapply(.SD, function(x) fifelse(x == '', NA_character_, x))]
d[, `#chr` := sub('^chr', '', chr)]
d[, start := as.character(as.numeric(pos) - 1)]
d[, end := pos]
d[, region_ID := blockId]
d[, gene_ID := blockId]
d[, event_ID := study]
cols <- c('#chr','start','end','a1','a2','variant_ID','region_ID','event_ID','chr','pos','af','betahat','sebetahat','z','gene_ID',
          'cs_coverage_0.95','cs_coverage_0.7','cs_coverage_0.5','cs_95_purity','PIP','logBF','conditional_effect','lfsr','blockId','context','within_cs_pip','block_id')
miss <- setdiff(cols, names(d)); if (length(miss)) stop('missing columns: ', paste(miss, collapse = ', '))
fwrite(d[, ..cols], a[2], sep = '\t', quote = FALSE, na = 'NA', scipen = 50)
cat('[conform_gwas_toploci]', study, nrow(d), 'rows ->', a[2], '\n')
