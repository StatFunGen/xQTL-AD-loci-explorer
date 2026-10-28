## TransMap trans-xQTL-only colocalization (disease agnostic, Xuewei Cao) overlapped with the AD loci, for discovery.
## Reads the export as one file, keeps the TransMap CoS whose variants sit at an AD locus variant (same chr:pos,
## the overlap rule used for every xQTL method), and labels each outcome: TransMap-eGene (trans gene, ENSG),
## TransMap-eProgram (trans_local_program ... hotspot_module) or TransMap-biomarker (glyco, protein, metabolite).
## It counts as trans xQTL support for an AD locus (with the TransMap x AD colocalization), but does not define loci
## or set gene tiers (Gao, Oct 28: TransMap is a discovery tool, it just does not count into the gene tiers).
## Usage: Rscript scripts/annotate_transmap_xqtl_only.R <trans_xQTL_only_colocalization_summary_table.bed> <build out dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE); if (length(a) != 2) stop('usage: annotate_transmap_xqtl_only.R <file> <out dir>')
out <- a[2]
v <- fread(file.path(out, 'AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz'), select = c('ADlocusID', 'variant_ID'))
v[, `:=`(chr = sub(':.*', '', variant_ID), pos = as.integer(sub('^[^:]+:([0-9]+):.*', '\\1', variant_ID)))]
t <- fread(a[1], select = c('#chr', 'end', 'variant_ID', 'cos_ID', 'event_ID', 'vcp'))
setnames(t, c('#chr', 'end', 'variant_ID'), c('chr', 'pos', 'transmap_variant_ID'))
m <- merge(unique(v), t, by = c('chr', 'pos'))
m[, rid := .I]
o <- m[, .(outcome = strsplit(event_ID, '; ', fixed = TRUE)[[1]]), by = rid]
m <- merge(m[, -'event_ID'], o, by = 'rid')[, rid := NULL]
m[, type := fcase(grepl('_trans_gene_', outcome) & grepl('ENSG', outcome), 'TransMap-eGene',
                  grepl('_trans_local_program_', outcome) & grepl('hotspot_module', outcome), 'TransMap-eProgram',
                  grepl('trans_biomarker', outcome), 'TransMap-biomarker', default = 'other')]
m[, context := fcase(grepl('^adjusted_glyco_', outcome), 'glyco', grepl('^protein_', outcome), 'protein',
                     grepl('^metabolomics_', outcome), 'metabolite', default = sub('_trans_.*', '', outcome))]
if (m[type == 'other', .N]) stop('unrecognised TransMap outcomes: ', paste(head(unique(m[type == 'other']$outcome)), collapse = ', '))
fwrite(m[order(ADlocusID, variant_ID, cos_ID, outcome)], file.path(out, 'transmap_xQTL_only_overlap_ADloci.csv.gz'))
x <- fread(file.path(out, 'xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz'), select = c('ADlocusID', 'Method'))
cisM <- c('AD_xQTL_colocalization', 'Coloc', 'fSuSiE_finemapping', 'multi_context_finemapping', 'multi_gene_finemapping', 'single_context_finemapping')
s <- data.table(ADlocusID = unique(v$ADlocusID))
s[, cis_xQTL_support := ADlocusID %in% x[Method %in% cisM]$ADlocusID]
s[, TransMap_ADtarget := ADlocusID %in% x[Method == 'transmap']$ADlocusID]
k <- m[, .(TransMap_xQTLonly_variants = uniqueN(variant_ID), eGene_CoS = uniqueN(cos_ID[type == 'TransMap-eGene']),
           eProgram_CoS = uniqueN(cos_ID[type == 'TransMap-eProgram']), biomarker_CoS = uniqueN(cos_ID[type == 'TransMap-biomarker']),
           contexts = paste(sort(unique(context)), collapse = ';')), by = ADlocusID]
s <- merge(s, k, by = 'ADlocusID', all.x = TRUE)
for (c in c('TransMap_xQTLonly_variants', 'eGene_CoS', 'eProgram_CoS', 'biomarker_CoS')) s[is.na(get(c)), (c) := 0L]
## xQTL support per locus: cis (fine-mapping and colocalization) and trans (TransMap x AD colocalization or TransMap xQTL-only overlap)
s[, trans_xQTL_support := TransMap_ADtarget | TransMap_xQTLonly_variants > 0]
s[, xQTL_support := fcase(cis_xQTL_support & trans_xQTL_support, 'cis and trans', cis_xQTL_support, 'cis only', trans_xQTL_support, 'trans only', default = 'none')]
fwrite(s[order(ADlocusID)], file.path(out, 'transmap_xQTL_only_locus_summary.csv'))
cat('[xqtl_support] loci', nrow(s), '| with xQTL support', s[xQTL_support != 'none', .N], '| cis only', s[xQTL_support == 'cis only', .N],
    '| trans only', s[xQTL_support == 'trans only', .N], '| cis and trans', s[xQTL_support == 'cis and trans', .N], '\n')
cat('[transmap_xqtl_only] loci', nrow(s), '| with TransMap xQTL-only overlap', s[TransMap_xQTLonly_variants > 0, .N],
    '| no cis-xQTL support but TransMap xQTL-only', s[!cis_xQTL_support & TransMap_xQTLonly_variants > 0, .N], '\n')
print(s[!cis_xQTL_support & TransMap_xQTLonly_variants > 0])
