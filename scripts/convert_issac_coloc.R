## Convert ISSAC single-nucleus sQTL x AD ColocBoost exports into the gene-keyed
## export format used by the other AD-xQTL ColocBoost files: region_ID is the gene,
## each event is "<cell type>_<gene>", and the AD event is listed last.
## A splice site is assigned to a gene by the sn-sQTL association table's own
## annotation (spliced_genes) when the site is listed there; otherwise by the gene
## on the same strand whose body contains the site, when there is exactly one.
## Usage: Rscript convert_issac_coloc.R <gene_gtf> <sn_sQTL_assoc.txt.gz> <out_dir> <export.bed> ...
suppressMessages({library(data.table); library(stringr)})
a <- commandArgs(TRUE)
gtf <- a[1]; assoc <- a[2]; od <- a[3]; ins <- a[-(1:3)]
stopifnot(file.exists(gtf), file.exists(assoc), dir.exists(od), length(ins) > 0)

g <- fread(cmd = paste("grep -v '^#'", shQuote(gtf)), header = FALSE, sep = '\t', quote = '')[V3 == 'gene']
g[, gene_ID := sub('\\..*$', '', sub('.*gene_id "([^"]+)".*', '\\1', V9))]
g[, gene_name := sub('.*gene_name "([^"]+)".*', '\\1', V9)]
g[, gene_type := sub('.*gene_type "([^"]+)".*', '\\1', V9)]
g <- g[, .(chr = V1, start = V4, end = V5, strand = V7, gene_ID, gene_name, gene_type)]
setkey(g, chr, strand, start, end)

s <- fread(assoc, col.names = c('splice_site','variant_ID','pval','beta','se','FDR','if_significant',
                                'spliced_genes','gene_type','splice_type','cell_type'))
s <- unique(s[!is.na(spliced_genes) & spliced_genes != '', .(site = splice_site, gene_name = sub(',.*$', '', spliced_genes))])
s <- merge(s, unique(g[, .(gene_name, gene_ID)])[, if (.N == 1) .SD, by = gene_name], by = 'gene_name')
s <- s[, if (uniqueN(gene_ID) == 1) .SD[1], by = site]

x <- rbindlist(lapply(ins, function(f) fread(f)[, src := basename(f)]))
setnames(x, '#chr', 'chr_', skip_absent = TRUE)
x[, row := .I]
L <- x[, .(ev = trimws(strsplit(event_ID, ';', fixed = TRUE)[[1]]),
           cf = trimws(strsplit(as.character(coef), ';', fixed = TRUE)[[1]])), by = row]
L[, is_ad := startsWith(ev, 'AD_')]
ad <- L[(is_ad), .(ad_ev = ev[1], ad_cf = cf[1]), by = row]
q <- L[!(is_ad)]
q[, ctx := sub('_chr.*$', '', ev)]
q[, site := sub('^.*_(chr[0-9XY]+:[+-]:[0-9]+)$', '\\1', ev)]
q[, `:=`(chr = sub(':.*', '', site), strand = sub('^chr[0-9XY]+:([+-]):.*', '\\1', site),
         pos = as.integer(sub('.*:', '', site)))]
## (1) annotation from the sn-sQTL association table
q <- merge(q, s[, .(site, gene_assoc = gene_ID)], by = 'site', all.x = TRUE)
## (2) gene body on the same strand
us <- unique(q[, .(site, chr, strand, start = pos, end = pos)])
setkey(us, chr, strand, start, end)
ov <- foverlaps(us, g, nomatch = NULL)[, .(n = uniqueN(gene_ID), gene_body = gene_ID[1]), by = site]
q <- merge(q, ov, by = 'site', all.x = TRUE)
q[, gene := fifelse(!is.na(gene_assoc), gene_assoc, fifelse(!is.na(n) & n == 1, gene_body, NA_character_))]
q[, how := fifelse(!is.na(gene_assoc), 'sn-sQTL annotation',
           fifelse(is.na(n), 'no gene at site', fifelse(n == 1, 'single gene body', 'several genes (dropped)')))]
cat('== splice sites by how they were assigned to a gene\n')
print(unique(q[, .(site, how)])[, .N, by = how][order(-N)])

q <- q[!is.na(gene)]
q[, cfn := suppressWarnings(as.numeric(cf))]
## one event per cell type and gene: keep the splice site with the largest |coef|
q <- q[order(row, gene, ctx, -abs(cfn))][, .SD[1], by = .(row, gene, ctx)]
o <- q[, .(ev2 = paste(paste0(ctx, '_', gene), collapse = '; '), cf2 = paste(cf, collapse = ';')), by = .(row, gene)]
o <- merge(o, ad, by = 'row')
o <- merge(o, x, by = 'row')
o[, `:=`(region_ID = gene, event_ID = paste(ev2, ad_ev, sep = '; '),
         coef = paste(cf2, ad_cf, sep = ';'), cos_ID = paste0(gene, ':', cos_ID))]
cols <- c('chr_','start','end','a1','a2','variant_ID','region_ID','event_ID','cos_ID','vcp','cos_npc','min_npc_outcome','coef','gwas_1e5')
cat('== rows in -> rows out, by file\n')
for (f in unique(x$src)) {
  y <- o[src == f, ..cols]; setnames(y, 'chr_', '#chr')
  fo <- file.path(od, sub('\\.bed$', '.gene.bed', f))
  fwrite(y, fo, sep = '\t', quote = FALSE)
  cat(sprintf('%s: %d -> %d rows, %d genes, %d CoS -> %s\n', f, x[src == f, .N], nrow(y),
              uniqueN(y$region_ID), uniqueN(x[src == f, cos_ID]), basename(fo)))
}
cat('CONVERT-DONE\n')
