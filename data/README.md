# Released tables

Four builds are published here. They are kept in separate directories because the
same workbook filename has been used across builds, and the `202605` / `202609`
suffix is the build version of the summary script, not the number of loci.

| directory | build | workbook | rows | loci |
|---|---|---|---|---|
| 20260514/ | earlier release | unified_AD_loci_xQTL_summary_202605.xlsx | 4,114 | 183 |
| 20260917/ | earlier release | unified_AD_loci_xQTL_summary_202609.xlsx | 4,195 | 188 |
| 20261002/ | earlier release | unified_AD_loci_xQTL_summary_20261002.xlsx | 5,016 | 195 |
| 20261022/ | current release, shown in the Explorer | unified_AD_loci_xQTL_summary_20261022.xlsx | 5,874 | 196 |

Counts are measured from the `Gene Locus table` sheet of each workbook and agree
on both `locus index` and `ADlocus`.

## What differs between 20260514 and 20260917

Both builds draw on the same eight AD GWAS studies. They differ in how those
studies were fine-mapped, and in one added evidence source.

| | 20260514 (202605) | 20260917 (202609) |
|---|---|---|
| GWAS studies | the same eight | the same eight, in the current format |
| GWAS fine-mapping | SuSiE-RSS version 1 | SuSiE-RSS EB-mix |
| colocalization | as in that release | adds transmap colocalization |

The credible sets were recomputed with a different method rather than simply
extended, so the two locus sets are not nested and neither contains the other.
Comparing the two workbooks on variant IDs: 143 loci carry the same identifier,
8 loci in the 20260917 build share no variant with any locus in the 20260514
build, and 8 loci in the 20260514 build share no variant with any locus in the
20260917 build; the remainder shift their boundaries. Of the 1,182 variants
carrying a maximum inclusion score in both, 869 are identical, 120 are higher in
the 20260917 build and 193 are lower.

That some scores are *lower* is expected here and is worth understanding. Within
a single build the inclusion score is a maximum across methods and sources, so
adding a source can only raise it. It can fall only when the underlying
probabilities are themselves recomputed, which is exactly what changing the
fine-mapping method does. A drop between these two builds is therefore a
property of the method change, not evidence that something was removed.

Locus numbers are ordinals within one build, so the same number does not denote
the same locus across the two. Always compare on variant IDs.

## What differs between 20260917 and 20261002

Both builds use the same eight AD GWAS studies, fine-mapped the same way. The 20261002 build changes the xQTL evidence:

| | 20260917 (202609) | 20261002 (20261002) |
|---|---|---|
| colocalization | ColocBoost exports and transmap colocalization | adds tpQTL and CUIMC1 eQTL colocalization with AD, and transmap AD-target colocalization |
| trans evidence | all trans fine-mapping results | the selected single-context trans fine-mapping set plus trans-xQTL-only and AD-target colocalization; genome-wide trans-eGene, trans-hotspot, trans pQTL and trans gpQTL fine-mapping are not included |

Loci are defined by a 95% credible set or an AD colocalization, so the added AD colocalizations define 7 new loci, and 4 existing loci gain a wider boundary. All 188 earlier loci are kept. Locus names are ordinals and are renumbered; `20261002/ADlocus_name_map_20260917_to_20261002.csv` maps each locus to its earlier name. Gene tiers use cis evidence only.


## What differs between 20261002 and 20261022

Both builds use the same eight AD GWAS studies, fine-mapped the same way. The 20261022 build adds evidence and changes how genes are tiered:

| | 20261002 (20261002) | 20261022 (20261022) |
|---|---|---|
| colocalization | ColocBoost exports, tpQTL and CUIMC1 eQTL colocalization, transmap AD-target colocalization | adds single-nucleus sQTL colocalization with AD in six brain cell types |
| evidence tiers | earlier tier rules | T1-T6 as defined in the repository README: T1-T5 need an AD GWAS 95% credible set or a 95% colocalization; trans evidence does not set a tier |

The added colocalizations define 4 new loci, merge 5 earlier loci into 2, and change the boundary of 9 others. All 195 earlier loci are kept or merged. Locus names are ordinals and are renumbered; `20261022/ADlocus_name_map_20261002_to_20261022.csv` maps each locus to its earlier name.

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.
