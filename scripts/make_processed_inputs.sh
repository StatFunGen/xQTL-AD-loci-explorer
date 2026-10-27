#!/bin/bash
## make_processed_inputs.sh -- regenerate every post-processed build input from its
## original files. Each step names its parents; outputs go to OUTDIR and can be
## compared byte for byte with the files the build used (see the input list).
## bash scripts/make_processed_inputs.sh <project_root> <outdir>
set -euo pipefail
ROOT=$1; OUT=$2; S=$(cd "$(dirname "$0")" && pwd); RS=${AD_LOCI_RSCRIPT:-$HOME/.pixi/bin/Rscript}; [ -x "$RS" ] || RS=$(command -v Rscript)
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
$RS $S/convert_transmap_AD_target.R $OUT/tmp/AD_xQTL_transmap_export_release.noTMEM106B.sets.tsv $OUT/AD_transmapADtarget_colocboost_export_20261027.bed xuewei_slack/AD_xQTL_ROSMAP_colocboost_export_release_20260422.TMEM106B_fix.bed
# ISSAC x AD: splice sites -> genes (GENCODE v32 genes + sn-sQTL association table)
zcat resource/references/gencode.v32.basic.annotation.gtf.gz | awk -F'\t' '$3=="gene"' > $OUT/tmp/gencode.v32.genes.gtf
$RS $S/convert_issac_coloc.R $OUT/tmp/gencode.v32.genes.gtf interactive_analysis/$COL/data/sn-sQTL/major_cell_FDR_0_05.txt.gz $OUT/issac_coloc_gene AD_*_issac_colocboost_export_release.bed
# glycoQTL x AD (Ru, Oct 16): sets with AD, one row per variant
$RS $S/convert_glyco_AD_coloc.R $OUT/ru_glyco_patch_20261016/AD_glycoQTL_colocboost_export_ru20261016.bed xuewei_slack/AD_xQTL_ROSMAP_colocboost_export_release_20260422.bed ru_glyco_patch_20261016/ROSMAP_pQTL_glycoQTL_adjusted_AD_summary_pv1e4.tsv ru_glyco_patch_20261016/ROSMAP_pQTL_glycoQTL_unadjusted_AD_summary_pv1e4.tsv
# tpQTL fine-mapping (Ru K5): split the genome-wide table by region, rows unchanged
for c in AC DLPFC PCC; do zcat k5_transfer/ROSMAP_sQTL_factor.GenomeWide.exported.toploci.bed.gz | awk -F'\t' -v c=ROSMAP_${c}_sQTL_factor 'NR==1 || index($8,c)==1' | gzip -n > $OUT/k5_transfer/ROSMAP_${c}_sQTL_factor.exported.toploci.bed.gz; done
echo MAKE-DONE
