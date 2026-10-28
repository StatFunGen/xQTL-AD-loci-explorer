# Released tables - 20261028 build

Locus-level results for the candidate Alzheimer's disease loci in the final
20261028 build (input list `build_inputs/20261028.2.inputs.csv`), produced by
build_AD_locus_table.R. See the repository README for how to regenerate them.

Every table here comes from the same build.

| file | contents |
|---|---|
| unified_AD_loci_xQTL_summary_20261028.xlsx | unified summary workbook: locus, gene and variant evidence across fine-mapping, colocalization, TWAS/MR and cTWAS, with evidence tiers T1-T6 |
| AD_loci_unified_cs95orColocs_Pval1e5_variant_level.csv.gz | one row per variant-locus: variants entering a locus through a 95% credible set or a colocalization, with per-variant GWAS summary columns. The file name is kept from earlier builds; this build applies no p-value filter |
| ADlocus_name_map_20261022_to_20261028.csv | maps each locus to its name in the 20261022 build by shared variants: same, boundary changed, merged, or new |

## Not included here

xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz is the
long-form table with one row per variant-locus-method-context-gene. It exceeds
what a git repository can hold. It is on Synapse (AD_xQTL_integration_summary,
syn69865823), or available on request: open an issue on this repository.

## Provenance

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.

GWAS credible sets come from the final SuSiE-RSS EB-mix run of the eight AD GWAS
studies (R4_L10). LD among AD variants is computed during the build from the
ADSP R4 EUR LD sketch, in 2 Mb windows. Evidence tiers follow the T1-T6
definitions in the repository README. TransMap colocalization counts as trans
xQTL support for a locus and does not set a gene tier. fSuSiE colocalization of
AD with bulk methylation and histone QTL is included and counts as CoS support. See `../README.md` for
what differs between the builds.
