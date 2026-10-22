# Released tables - 196-loci build

Locus-level results for the 196 candidate Alzheimer's disease loci in this
release, produced by build_AD_locus_table.R. See the repository README for how
to regenerate them.

Every table here comes from the same build. Verified locus and row counts:

| file | rows | loci | contents |
|---|---|---|---|
| unified_AD_loci_xQTL_summary_20261022.xlsx | - | 196 | unified summary workbook: locus, gene and variant evidence across fine-mapping, colocalization, TWAS/MR and cTWAS, with evidence tiers T1-T6 |
| AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz | 20,572 | 196 | one row per variant-locus: variants entering a locus through a 95% credible set or a colocalization and reaching p < 1e-5, with per-variant GWAS summary columns |
| ADlocus_name_map_195_to_196.csv | 199 | 196 | maps each locus to its name in the 195-locus build: same, boundary changed, merged, or new |

## Not included here

xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz is the
long-form table with one row per variant-locus-method-context-gene. At 160 MB it exceeds
what a git repository can hold, as does the noTrans variant at 136 MB.
Both are available on request: open an issue on this repository or email
jaempawi@bu.edu.

## Provenance

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.

GWAS credible sets in this build were fine-mapped with SuSiE-RSS EB-mix. Colocalization evidence adds single-nucleus sQTL (sn-sQTL) colocalization with AD from ColocBoost in six brain cell types; each splice event is assigned to a gene from the sn-sQTL annotation, or otherwise to the single same-strand gene whose body contains it, and events with no single gene are left out. Evidence tiers follow the T1-T6 definitions in the repository README. Trans evidence is reported for reference and does not set a tier. See `../README.md` for what differs between the builds.
