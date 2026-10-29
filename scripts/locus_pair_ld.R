## |r| diagnostic between AD loci (review comment on build:1588, Ru Feng): for every pair of loci with at least one
## variant pair in the LD table (|r| > 0.5), report how many cross-locus variant pairs are listed, the smallest and
## largest listed |r|, and the exact minimum when all n1 x n2 pairs are listed (otherwise '<0.5': some pair is below 0.5).
## max_found > 0.8 would mean two loci that the merge rule could have joined.
## Usage: Rscript scripts/locus_pair_ld.R <build out dir> [LD table, default <out>/gwas_variants_cor0.5.csv.gz]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); out <- a[1]; ldf <- if (length(a) > 1) a[2] else file.path(out, 'gwas_variants_cor0.5.csv.gz')
vl <- unique(fread(file.path(out, 'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'), select = c('variant_ID', 'ADlocus')))
if (vl[, uniqueN(ADlocus), by = variant_ID][V1 > 1, .N]) stop('variants mapped to more than one locus')
sz <- vl[, .(n = .N), by = ADlocus]
ld <- fread(ldf, select = c('variant_ID', 'variant_ID2', 'r'))
cat('[locus_pair_ld] locus variants missing from the LD table:', vl[!variant_ID %in% ld$variant_ID, .N], 'of', nrow(vl), '\n')
lc <- unique(rbind(ld[, .(v1 = variant_ID, v2 = variant_ID2, a = abs(r))], ld[, .(v1 = variant_ID2, v2 = variant_ID, a = abs(r))]), by = c('v1', 'v2'))
lc <- merge(lc, vl[, .(v1 = variant_ID, L1 = ADlocus)], by = 'v1')
lc <- merge(lc, vl[, .(v2 = variant_ID, L2 = ADlocus)], by = 'v2')[L1 < L2]
p <- lc[, .(n_found = .N, min_found = min(a), max_found = max(a)), by = .(ADlocus_1 = L1, ADlocus_2 = L2)]
p <- merge(merge(p, sz[, .(ADlocus_1 = ADlocus, n1 = n)], by = 'ADlocus_1'), sz[, .(ADlocus_2 = ADlocus, n2 = n)], by = 'ADlocus_2')
p[, n_pairs := as.numeric(n1) * n2][, min_absr := fifelse(n_found == n_pairs, min_found, NA_real_)]
p[, label := fifelse(is.na(min_absr), '<0.5', sprintf('%.3f', min_absr))]
setcolorder(p, c('ADlocus_1', 'ADlocus_2', 'n_found', 'min_found', 'max_found', 'n1', 'n2', 'n_pairs', 'min_absr', 'label'))
fwrite(p[order(ADlocus_1, ADlocus_2)], file.path(out, 'locus_pair_min_abs_ld.csv.gz'))
cat('[locus_pair_ld] locus pairs with any |r| > 0.5:', nrow(p), '| all pairs > 0.5 (exact min):', p[!is.na(min_absr), .N],
    '| any pair |r| > 0.8:', p[max_found > 0.8, .N], '\n')
print(p[max_found > 0.8][order(-max_found)])
