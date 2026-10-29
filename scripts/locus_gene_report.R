## Multi-function of AD loci and context specificity (Gao, Oct 28): from the tiered xQTL overlap table of a build,
## per AD locus count the genes and cell types with tiered support, and per gene count the loci it is supported at.
## A gene's tier at a locus is its strongest tier among its rows at that locus (trans evidence has no tier).
## Cell type = first word of context_short (Ast, End, Exc, Inh, Mic, OPC, Oli, Peri; bulk = bulk brain tissue;
## bMic / bMono / bMac = bulk microglia / monocyte / macrophage). Myeloid = Mic, bMic, bMono, bMac.
## Writes <out>/locus_gene_report.txt and <out>/locus_gene_support.csv.gz.
## Usage: Rscript scripts/locus_gene_report.R <build out dir>
suppressPackageStartupMessages(library(data.table))
out <- commandArgs(TRUE)[1]; if (is.na(out)) stop('usage: locus_gene_report.R <build out dir>')
x <- fread(file.path(out, 'xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5_noTrans.csv.gz'),
           select = c('ADlocusID', 'gene_ID', 'gene_name', 'context_short', 'confidence_lvl', 'APOE_region'))
x <- x[confidence_lvl %in% paste0('T', 1:5)][, tier := as.integer(sub('T', '', confidence_lvl))]
x[, cell := sub(' .*', '', context_short)][cell == '', cell := NA_character_]
myeloid <- c('Mic', 'bMic', 'bMono', 'bMac')
lines <- character(); say <- function(...) { s <- paste0(...); lines <<- c(lines, s); cat(s, '\n') }
pct <- function(a, b) sprintf('%d of %d (%.1f%%)', a, b, 100 * a / b)
report <- function(d, label) {
  g <- d[, .(tier = min(tier)), by = .(ADlocusID, gene_ID)]                 # gene tier at each locus
  L4 <- unique(g[tier <= 4]$ADlocusID); n <- length(L4)
  say('== ', label)
  say('AD loci with T1-T4 gene support: ', n)
  ng4 <- g[tier <= 4, .N, by = ADlocusID]; ng5 <- g[ADlocusID %in% L4 & tier <= 5, .N, by = ADlocusID]
  c4 <- d[tier <= 4 & !is.na(cell), .(n = uniqueN(cell)), by = ADlocusID]; c5 <- d[ADlocusID %in% L4 & tier <= 5 & !is.na(cell), .(n = uniqueN(cell)), by = ADlocusID]
  say('  multiple T1-T4 genes: ', pct(sum(ng4$N > 1), n), ' | genes per locus 1/2/3/4+: ',
      paste(sum(ng4$N == 1), sum(ng4$N == 2), sum(ng4$N == 3), sum(ng4$N >= 4), sep = '/'))
  say('  multiple cell types (T1-T4): ', pct(sum(c4$n > 1), n))
  c4s <- d[tier <= 4 & !is.na(cell) & cell != 'bulk', .(n = uniqueN(cell)), by = ADlocusID]
  c5s <- d[ADlocusID %in% L4 & tier <= 5 & !is.na(cell) & cell != 'bulk', .(n = uniqueN(cell)), by = ADlocusID]
  say('  multiple cell types, not counting bulk brain tissue: T1-T4 ', pct(sum(c4s$n > 1), n), ', including T5 ', pct(sum(c5s$n > 1), n))
  say('  including T5 at these loci: multiple genes ', pct(sum(ng5$N > 1), n), ', multiple cell types ', pct(sum(c5$n > 1), n))
  say('  myeloid (Mic, bMic, bMono, bMac) T1-T4 support: ', pct(uniqueN(d[tier <= 4 & cell %in% myeloid]$ADlocusID), n))
  ct <- d[tier <= 4 & !is.na(cell), .(loci = uniqueN(ADlocusID)), by = cell][order(-loci)]
  say('  loci with T1-T4 support by cell type: ', paste0(ct$cell, ' ', ct$loci, collapse = ', '))
  one <- d[tier <= 4 & !is.na(cell), .(n = uniqueN(cell), c = cell[1]), by = ADlocusID][n == 1, .N, by = c][order(-N)]
  say('  loci whose T1-T4 support is from one cell type only: ', paste0(one$c, ' ', one$N, collapse = ', '))
  for (t in c(4, 5)) { k <- g[tier <= t, .(loci = uniqueN(ADlocusID)), by = gene_ID]
    say('  genes with T1-T', t, ' support at more than one locus: ', pct(sum(k$loci > 1), nrow(k))) }
  g
}
g <- report(x, 'all loci')
invisible(report(x[!APOE_region %in% TRUE], 'excluding the APOE region (chr19:43.9-45.9 Mb)'))
w <- g[, .(T1_T4_genes = sum(tier <= 4), T1_T5_genes = .N), by = ADlocusID]
w <- merge(w, x[tier <= 4 & !is.na(cell), .(T1_T4_cell_types = paste(sort(unique(cell)), collapse = ';')), by = ADlocusID], by = 'ADlocusID', all.x = TRUE)
fwrite(w[order(ADlocusID)], file.path(out, 'locus_gene_support.csv.gz'))
say('rows without a cell type (not counted in cell-type numbers): ', x[is.na(cell), .N], ' of ', nrow(x))
writeLines(lines, file.path(out, 'locus_gene_report.txt'))
