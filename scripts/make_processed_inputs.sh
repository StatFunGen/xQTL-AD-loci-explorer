#!/bin/bash
## make_processed_inputs.sh -- regenerate every post-processed build input from its
## original files. Each step names its parents; outputs go to OUTDIR and can be
## compared byte for byte with the files the build used (see the input list).
## bash scripts/make_processed_inputs.sh <project_root> <outdir> [inputs.csv]
## With an input list, every regenerated file is compared with the md5 recorded for it
## there, and the script fails on any difference or on a file the list does not name.
set -euo pipefail
ROOT=$1; OUT=$2; INPUTS=${3:+$(realpath "$3")}; S=$(cd "$(dirname "$0")" && pwd); RS=${AD_LOCI_RSCRIPT:-$HOME/.pixi/bin/Rscript}; [ -x "$RS" ] || RS=$(command -v Rscript)
COL=${AD_LOCI_COLLABORATOR:-$(cat ~/.ad_loci_collaborator)}
mkdir -p $OUT/xuewei_slack $OUT/issac_coloc_gene $OUT/k5_transfer $OUT/ru_glyco_patch_20261016 $OUT/tmp
cd $ROOT
# ROSMAP AD x xQTL: remove the two TMEM106B glyco signals not in Ru's Oct 16 rerun
$RS $S/drop_events.R xuewei_slack/AD_xQTL_ROSMAP_colocboost_export_release_20260422.bed $OUT/xuewei_slack/AD_xQTL_ROSMAP_colocboost_export_release_20260422.TMEM106B_fix.bed ENSG00000106460 DLPFC_Klein_gpQTL_adjusted_gp_8491 DLPFC_Klein_gpQTL_unadjusted_gp_8492
# tpQTL x AD: remove ROSMAP eQTL outcomes already in the main ROSMAP export
$RS $S/dedup_against_export.R AD_UofCsQTL_colocboost_export_release.bed $OUT/AD_UofCsQTL_colocboost_export_release.ROSMAP_eQTL_dups_removed.bed xuewei_slack/AD_xQTL_ROSMAP_colocboost_export_release_20260422.bed
# TransMap x AD: drop the TMEM106B-locus CoS, collapse to sets, split trans by context, keep AD + cis
$RS $S/drop_cos.R xuewei_slack/AD_xQTL_transmap_export_release.bed $OUT/tmp/AD_xQTL_transmap_export_release.noTMEM106B.bed chr7 ENSG00000106460
$RS $S/transmap_export_to_sets.R $OUT/tmp/AD_xQTL_transmap_export_release.noTMEM106B.bed $OUT/tmp/AD_xQTL_transmap_export_release.noTMEM106B.sets.tsv
$RS $S/convert_transmap_AD_target_trans.R $OUT/tmp/AD_xQTL_transmap_export_release.noTMEM106B.sets.tsv $OUT/trans_toploci_transmap_ADtarget_20261027
$RS $S/convert_transmap_AD_target.R $OUT/tmp/AD_xQTL_transmap_export_release.noTMEM106B.sets.tsv $OUT/AD_transmapADtarget_colocboost_export_20261027.bed $OUT/xuewei_slack/AD_xQTL_ROSMAP_colocboost_export_release_20260422.TMEM106B_fix.bed
# ISSAC x AD: splice sites -> genes (GENCODE v32 genes + sn-sQTL association table)
zcat resource/references/gencode.v32.basic.annotation.gtf.gz | awk -F'\t' '$3=="gene"' > $OUT/tmp/gencode.v32.genes.gtf
$RS $S/convert_issac_coloc.R $OUT/tmp/gencode.v32.genes.gtf interactive_analysis/$COL/data/sn-sQTL/major_cell_FDR_0_05.txt.gz $OUT/issac_coloc_gene AD_*_issac_colocboost_export_release.bed
# glycoQTL x AD (Ru, Oct 16): sets with AD, one row per variant
$RS $S/convert_glyco_AD_coloc.R $OUT/ru_glyco_patch_20261016/AD_glycoQTL_colocboost_export_ru20261016.bed xuewei_slack/AD_xQTL_ROSMAP_colocboost_export_release_20260422.bed ru_glyco_patch_20261016/ROSMAP_pQTL_glycoQTL_adjusted_AD_summary_pv1e4.tsv ru_glyco_patch_20261016/ROSMAP_pQTL_glycoQTL_unadjusted_AD_summary_pv1e4.tsv
# tpQTL fine-mapping (Ru K5): split the genome-wide table by region, rows unchanged
for c in AC DLPFC PCC; do zcat k5_transfer/ROSMAP_sQTL_factor.GenomeWide.exported.toploci.bed.gz | awk -F'\t' -v c=ROSMAP_${c}_sQTL_factor 'NR==1 || index($8,c)==1' | gzip -n > $OUT/k5_transfer/ROSMAP_${c}_sQTL_factor.exported.toploci.bed.gz; done
if [ -n "$INPUTS" ]; then
  cat > $OUT/tmp/verify.R <<'REOF'
## each regenerated file must equal the file the build read: same bytes, or for .gz the
## same decompressed content (gzip headers differ between gzip versions)
suppressPackageStartupMessages(library(data.table)); a <- commandArgs(TRUE)
x <- fread(a[1], colClasses = 'character'); out <- normalizePath(a[2]); root <- a[3]
fs <- list.files(out, recursive = TRUE); fs <- fs[!startsWith(fs, 'tmp/')]
m <- data.table(file = fs, now = unname(tools::md5sum(file.path(out, fs))))
m[, recorded := x$md5[match(file, x$Path)]]
## md5 of the decompressed content; NA if the file does not decompress cleanly
zmd5 <- function(p) { r <- suppressWarnings(system2('bash', c('-c', shQuote(paste('set -o pipefail; gzip -dc', shQuote(p), '| md5sum'))), stdout = TRUE))
  if (!is.null(attr(r, 'status')) || !length(r)) NA_character_ else sub(' .*', '', r[1]) }
m[, result := fifelse(is.na(recorded), 'not in input list', fifelse(now == recorded, 'same bytes', 'DIFFERENT'))]
## a gzip file may differ only in its header: compare content with the project copy, but only
## after confirming that copy is the recorded version
for (i in which(m$result == 'DIFFERENT' & grepl('\\.gz$', m$file))) {
  ref <- file.path(root, m$file[i])
  if (!file.exists(ref) || unname(tools::md5sum(ref)) != m$recorded[i]) { m$result[i] <- 'DIFFERENT (project copy is not the recorded version)'; next }
  a1 <- zmd5(file.path(out, m$file[i])); a2 <- zmd5(ref)
  if (!is.na(a1) && !is.na(a2) && a1 == a2) m$result[i] <- 'same content (gzip header differs)' }
print(m[, .(file, result)], nrows = 200)
if (any(startsWith(m$result, 'DIFFERENT') | m$result == 'not in input list')) quit(status = 1)
cat('MAKE-VERIFY: all', nrow(m), 'regenerated files match the files the build read\n')
REOF
  "$RS" $OUT/tmp/verify.R "$INPUTS" "$OUT" "$ROOT" || { echo "MAKE-VERIFY failed"; exit 1; }
fi
echo MAKE-DONE
