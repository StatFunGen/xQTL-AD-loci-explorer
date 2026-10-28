# xQTL AD Loci Explorer

Code and released tables behind the **AD Loci Explorer**, a browsable view of the
ADSP FunGen-xQTL release across 196 candidate Alzheimer's disease loci.

Live app: https://xqtl.niagads.org/AD-loci-explorer/ (direct link: https://jenny-empawi.shinyapps.io/xQTL-AD-loci-explorer/)

## Layout

```
scripts/
  build_AD_locus_table.R   locus-level build: evidence integration and tier assignment
  gene_prio_utils.R        helper functions, including the T1-T6 tier rules
  convert_issac_coloc.R    keys single-nucleus sQTL (ISSAC) x AD colocalization exports to genes
  add_evidence_source.R    register a new evidence source (see 'Adding data')
  check_inputs.R           checks every file in a build input list exists and matches its md5
  build_record.R           writes the build record, <out>/_build_record/ (see Running it)
  compare_builds.R         compares two builds table by table
  make_processed_inputs.sh  regenerates the processed inputs from the raw exports
  convert_*.R, drop_*.R, dedup_against_export.R  the processing steps make_processed_inputs.sh runs
  validate_outputs.R       check a release against the expectations of this build
  run_pipeline.qsub        end-to-end job script for the SCC scheduler
config/                    metadata tables the build reads
data/                      released tables, one directory per build date
app/                       Shiny application, runnable from a clone
DEPENDENCIES.md            R version and package versions
```

## Requirements

R 4.5.3 with the packages listed in DEPENDENCIES.md:

```r
install.packages(c("data.table","stringr","openxlsx","tidyverse","ggplot2",
                   "scales","circlize","shiny","bslib","DT","plotly",
                   "shinycssloaders"))
remotes::install_github("StatFunGen/pecotmr")   # not on CRAN
```

## Running it

There are two ways in, depending on whether you want to reproduce the published
release or add your own data to it.

### Option 1: reproduce the release from scratch

Rebuilds the locus table from the registered inputs. This needs access to the
tables listed in `config/metadata_analysis.csv`, which are controlled-access; see
Inputs for how to obtain them.

```bash
export AD_LOCI_ROOT=/path/to/AD_loci_xQTL      # where the input trees live
export AD_LOCI_STAGING=/path/to/staging        # three precomputed tables, see Inputs
export AD_LOCI_OUT=$AD_LOCI_ROOT/out_$(date +%Y%m%d)    # optional
export AD_LOCI_CONFIG=/path/to/private/config          # optional, see Inputs

Rscript scripts/build_AD_locus_table.R                 # -> $AD_LOCI_OUT
Rscript scripts/validate_outputs.R  $AD_LOCI_OUT       # optional count check: AD_LOCI_EXPECT_LOCI, AD_LOCI_EXPECT_GENES
Rscript app/build_shiny_data.R $AD_LOCI_OUT app/data.csv
```

On the SCC the same steps run as one job. From the root of a clean checkout of the commit you want to build:

```bash
qsub -v AD_LOCI_INPUTS=data/build_inputs/YYYYMMDD.inputs.csv,AD_LOCI_OUT=/path/to/out_YYYYMMDD scripts/run_pipeline.qsub
```

The job:

1. refuses to write into an output folder that already has files (override with `AD_LOCI_ALLOW_EXISTING_OUT=1`);
2. runs `check_inputs.R`, which stops the job if any file the build will read is missing, has changed, or is not in the input list;
3. builds and validates the tables;
4. writes the build record (below);
5. writes `_COMPLETE` only if the build, the validation and the build record all succeeded.

Before anything runs, the job copies the input list and the config (`config/`, or `AD_LOCI_CONFIG`) into the build record; the input check, the build and the record all use those copies, so editing the originals during a run cannot change it. The three precomputed tables are read in place from `$AD_LOCI_ROOT/staging_gene_prio` (or `AD_LOCI_STAGING`) and are never written to. Scratch files go to `<out>/_scratch`.

Each build has one record, `<out>/_build_record/`:

| file | contents |
|---|---|
| `run.csv` | run id, code commit and whether the checkout had local changes, input list and its sha256, config sha256, R and package versions, host, user, time |
| `inputs.csv` | the input list the build used (copied before the build) |
| `config/` | the config files the build read (copied before the build), with `MD5SUMS` |
| `registry_runtime.csv` | the registry as the build filled it in (run-time paths of intermediate tables) |
| `check_inputs.txt`, `validation.txt` | output of the input check and of the validation |
| `outputs_md5.csv`, `summary.csv` | checksum of every output file; headline counts |
| `changes_since_previous.txt` | inputs, counts and outputs that differ from the build in `AD_LOCI_PREV`, when set |
| `environment.txt` | R session information |

The build uses `NSLOTS` threads. The per-file overlap step runs in parallel; set `AD_LOCI_OVERLAP_WORKERS` to change the number of workers. Set `AD_LOCI_RSCRIPT` to choose which Rscript runs the steps (default: `~/.pixi/bin/Rscript`, otherwise the one on `PATH`).

The build checks every required input before starting and stops with a list of
whatever is missing, rather than failing part-way through a long run.

### Option 2: add your own evidence to the published build

Register new results as an evidence source and rebuild with them included.
`scripts/add_evidence_source.R` converts a ColocBoost-style colocalization
summary and adds its rows to the registry; see Adding data for the input format,
and for what a source in any other format needs.

```bash
Rscript scripts/add_evidence_source.R --name mysource --bed coloc_summary.bed
Rscript scripts/build_AD_locus_table.R
```

Adding a source is a configuration change, not a code change. Be aware that the
rebuild still reads the whole registry, so this route needs the same input
access as option 1; it saves you from modifying the pipeline, not from having
the data.

If all you want is the Explorer, neither option is needed: `app/data.csv` ships
with the repository. See Running the app.

## Versions

Three things are versioned separately:

| what | how |
|---|---|
| code, config, processing recipe | git commit (each build records its commit in `_build_record/run.csv`) |
| accepted inputs | the dated input list in `data/build_inputs/` (path, md5, source, parents and recipe); files on project storage are not modified once listed |
| published releases | git tag `release-YYYYMMDD` on the commit that added the release's tables to `data/`; corrections get a new tag |

Work happens on short-lived branches that are merged and then deleted.

## Inputs

`config/metadata_analysis.csv` is the registry of every input the build reads.
One row per exported analysis table:

| column | meaning |
|---|---|
| context | dataset label, joined to config/contexts_metadata.csv |
| Data Type | eQTL, sQTL, pQTL, ... |
| Cohort | ROSMAP, MSBB, KNIGHT, ... |
| Modality | assay or analysis modality |
| Method | how the table was produced; selects the loader in the build |
| Path | location of the table, relative to `$AD_LOCI_ROOT` |
| variant_level_method | whether the method contributes variant-level evidence |

Paths under `$AD_LOCI_ROOT` refer to the FunGen-xQTL release, Synapse
`syn68872650`. The layout on Synapse does not necessarily mirror these relative
paths; email jaempawi@bu.edu if you need help locating a particular table. Some
rows carry `<user>` and `<collaborator>` placeholders where
tables were exported from per-user analysis directories; replace these with the
directory names in your own copy.

Three precomputed tables are too large to distribute with the code and are read
from `$AD_LOCI_STAGING`:

| file | size | source |
|---|---|---|
| gwas_variants_cor0.5.csv.gz | 295 MB | Synapse `syn75082260` |
| res_APOE_interaction_summ.csv.gz | 45 MB | on request |
| res_msex_interaction_summ.csv.gz | 1.6 MB | on request |

The two interaction summaries are not redistributed publicly. To request them,
open an issue on this repository or email jaempawi@bu.edu.

The correlation table drives the LD-based credible-set extension, so the
released variant sets cannot be reproduced without it.

### Build input lists

Each build has an input list in `data/build_inputs/`, named by build date (`YYYYMMDD.inputs.csv`). One row per file the build reads, with its path, md5, size, who provided it and where it came from. Files that were processed before use (for example, with signals removed or duplicates dropped) also list their `parents` (path and md5) and the `recipe` that produced them; `scripts/make_processed_inputs.sh` reruns every recipe and reproduces those files exactly. Rows with `role` = `reference` are lookup files rather than evidence.

#### Processed inputs

A few delivered files are changed before the build reads them. The original file is kept and listed as the parent; every step below is run by `scripts/make_processed_inputs.sh`.

| input | what is done | why | script |
|---|---|---|---|
| ROSMAP AD x xQTL colocalization | remove the two TMEM106B glycoprotein (DLPFC_Klein gpQTL) signals | they are not in Ru's Oct 16 glycoQTL rerun | `drop_events.R` |
| tpQTL x AD colocalization | remove ROSMAP eQTL outcomes that are already in the main ROSMAP export | avoid counting the same colocalization twice | `dedup_against_export.R` |
| TransMap x AD colocalization | remove the CoS at the TMEM106B locus, group rows into colocalization sets, then write the trans part as one file per context and the AD + cis part as one ColocBoost-style file | the build reads trans evidence one context at a time, and reads cis colocalization with the other AD x xQTL files; trans is shown as annotation only and does not set tiers | `drop_cos.R`, `transmap_export_to_sets.R`, `convert_transmap_AD_target_trans.R`, `convert_transmap_AD_target.R` |
| single-nucleus sQTL (ISSAC) x AD colocalization | map each splice site to its gene (GENCODE v32 genes and the sn-sQTL association table) | the build links colocalization to genes, and these exports carry splice sites | `convert_issac_coloc.R` |
| glycoQTL x AD colocalization | reshape to one row per variant per colocalization set | match the format of the other AD x xQTL files | `convert_glyco_AD_coloc.R` |
| tpQTL fine-mapping (genome-wide table) | split by context (AC, DLPFC, PCC); rows unchanged | the build reads fine-mapping one context at a time | in `make_processed_inputs.sh` |

To check a set of inputs before a run:

```bash
Rscript scripts/check_inputs.R data/build_inputs/YYYYMMDD.inputs.csv $AD_LOCI_ROOT
```

To see what changed between two builds:

```bash
Rscript scripts/compare_builds.R /path/to/out_A /path/to/out_B
```

## Adding data

To add a new evidence source you register it in `config/metadata_analysis.csv`.
The build auto-detects any Method there that points at per-context toploci files
and has no dedicated loader, so registration is the only step.

For a ColocBoost-style colocalization summary, `scripts/add_evidence_source.R` does the
conversion and registration for you:

```bash
Rscript scripts/add_evidence_source.R --name transmap \
        --bed trans_xQTL_only_colocalization_summary_table.bed
Rscript scripts/build_AD_locus_table.R
```

It is safe to re-run: conversion is skipped when outputs already exist, and rows
for the same `--name` are replaced rather than duplicated.

### Input format

The `--bed` file is tab-delimited with a header, one row per variant, and event
lists delimited by `;`:

| # | column | notes |
|---|---|---|
| 1 | #chr | |
| 2 | start | |
| 3 | end | |
| 4 | a1 | |
| 5 | a2 | |
| 6 | variant_ID | chr:pos:ref:alt |
| 7 | region_ID | |
| 8 | event_ID | semicolon-delimited list; context is the part before `_cis_` or `_trans_` |
| 9 | cos_ID | must contain `cosN`; rows with `cos0` are skipped |
| 10 | vcp | becomes PIP in the emitted file |
| 11 | cos_npc | becomes purity |
| 14 | zscore | semicolon-delimited, positionally matched to event_ID |

It emits one toploci file per context with the columns the build expects:
`#chr start end a1 a2 variant_ID region_ID event_ID cs_coverage_0.95 purity PIP z`.

A source in any other format needs its own loader in `scripts/build_AD_locus_table.R`;
follow one of the existing `Method ==` blocks.

## Running on another cluster

To change the code and test the change against a baseline build, see [CONTRIBUTING.md](CONTRIBUTING.md).

The build reads every input through one project folder (`AD_LOCI_ROOT`). On a cluster
where the inputs sit in different places, make that folder out of links, using the
location column of the input list (`columbia_path` for the Columbia cluster):

```bash
pixi install                                  # R and packages at the pinned versions (pixi.toml, pixi.lock)
Rscript scripts/make_input_root.R data/build_inputs/YYYYMMDD.inputs.csv /path/to/input_root
AD_LOCI_STAGING=/path/to/input_root/staging_gene_prio \
  Rscript scripts/check_inputs.R data/build_inputs/YYYYMMDD.inputs.csv /path/to/input_root config scripts/build_AD_locus_table.R
sbatch --export=ALL,AD_LOCI_ROOT=/path/to/input_root,AD_LOCI_INPUTS=data/build_inputs/YYYYMMDD.inputs.csv,AD_LOCI_OUT=/path/to/out_YYYYMMDD,AD_LOCI_RSCRIPT=$PWD/.pixi/envs/default/bin/Rscript scripts/run_pipeline.slurm
```

`make_input_root.R` only makes links; it stops if any input is missing on that cluster.
The input check then confirms every file is byte-identical to the input list, and
`scripts/compare_builds.R <release> <new build>` shows whether the result matches a
release exactly. `run_pipeline.slurm` runs the same steps as `run_pipeline.qsub` and
needs about 64 GB of memory.

## Methodology

This section describes what the build does, step by step, and the thresholds it
applies. Two terms recur. A *credible set* is the small group of variants that
fine-mapping says most likely contains the causal one; 95% coverage (cs95) means
the set is built to have a 95% chance of containing it, and cs70 and cs50 are
weaker versions of the same idea. *Colocalization* asks whether an AD signal and
a molecular signal in the same place are driven by the same variant rather than
by two neighbouring ones.

### Step 1: the AD locus set

The locus set is imported, not derived here: GWAS fine-mapping is run upstream
and is not repeated by this pipeline. The credible sets behind the current
build come from eight AD GWAS studies fine-mapped with SuSiE-RSS EB-mix; the
earlier 20260514 release used SuSiE-RSS version 1 on the same studies, which is
why the two are not nested (see `data/README.md`).
At startup the build scans the registry
for the analysis methods that point at per-context top-loci files and reports
what it found, so the source list is discovered from the configuration rather
than hardcoded. The job log records it as `[toploci] auto-detected sources`.

### Step 2: harmonization

Nothing is merged until the identifiers agree.

*Contexts.* Dataset and assay labels are mapped to one harmonized context name
through `config/contexts_metadata.csv`. Single-nucleus tables that carry a cell
type instead of a context are mapped by cell type.

*Variants.* Identifiers arrive in several spellings. Underscores are converted
to colons and a `chr` prefix is added where it is missing, but the merge itself
is done on chromosome and position rather than on the identifier string, and the
AD locus table's spelling is then adopted for the merged row. This is what
absorbs allele-order differences: where two sources name the same position with
the alleles written the other way round, the merge still succeeds, and the build
counts the disagreeing spellings instead of dropping those variants (213 of
5,318 in this build).

*Events and genes.* Event identifiers have trailing `_chr`, `_ENSG` and `_gp`
fragments stripped. Genes are matched on the Ensembl ID with the version suffix
removed, so `ENSG00000141510.16` and `ENSG00000141510` are treated as the same
gene.

*GWAS studies.* Source and study names are harmonized so that one GWAS is not
counted twice under two labels.

### Step 3: evidence families and their significance rules

| evidence | what makes it count |
|---|---|
| fine-mapping | membership in a credible set; the coverage level is recorded on each row as cs95, cs70 or cs50, and the strongest level available is kept per locus |
| colocalization | a variant-level colocalization probability (VCP), or the SNP-level PPH4 where that is what the method reports |
| TWAS | gene-level p < 2.5e-6, either in the selected best method or in at least half of the methods run, evaluated per gene, block, context and GWAS source |
| MR | TWAS-significant, and cPIP > 0.5, and at least two credible sets, and I2 < 0.5 |
| cTWAS | a non-empty credible set and SuSiE PIP > 0.75 |

A row that names a method but lacks that method's statistic is dropped, so a
colocalization row with no VCP and a fine-mapping row with no PIP never enter
the table.

### Step 4: one number per variant

Each surviving row is reduced to a single `variant_inclusion_probability` by
precedence rather than by averaging: the colocalization VCP when present,
otherwise the fine-mapping PIP, otherwise the SNP-level PPH4. The maximum across
every method and source is then taken per variant, and the method that produced
that maximum is recorded beside it. Because it is a maximum, adding a source can
only raise a variant's score and never lower it; a score that falls between two
builds therefore means the inputs or the filtering changed, not merely that
something was added. The GWAS methods, sources and effect directions behind each
variant are collapsed into parallel `|`-separated fields on the same row.

### Step 5: assembling and filtering loci

A locus identifier is built from its chromosome and the first and last positions
of its member variants, and the loci are then numbered in order. That number is
an ordinal within a single build, so the same number does not denote the same
locus in another build and cross-build comparisons have to be made on variant
IDs.

Three filters decide what survives. Only sets that qualify on their own form loci: AD GWAS
credible sets at 95% coverage and colocalization sets. cs70 and cs50 sets are
dropped as whole sets before sets that share a variant are merged into loci, so
a weaker set cannot join a locus by overlapping a qualifying one. A locus is kept only if its best GWAS p-value
is below 1e-5. Each surviving locus is then represented by up to five variants
whose GWAS PIP or VCP exceeds 0.1; where no variant clears 0.1, the single best
one is shown, ranked on the GWAS probability first and the xQTL probability
second. GWAS significance is reported in bands: p < 1e-5, p < 1e-6, p < 5e-8.

Two flags are carried per locus: whether it is supported only by proxy-based
GWAS, and whether it falls in the APOE region (chr19:43,905,790-45,905,791),
which is reported separately because the extreme linkage disequilibrium there
makes fine-mapping and colocalization hard to read.

### Step 6: tiers, outputs and provenance

Tiers are assigned per variant, gene and context group by the rules in the next
section, and each gene is reported at the strongest tier it reaches anywhere.
The build writes the locus-level summary, the variant-level membership table and
the unified workbook to `$AD_LOCI_OUT`; what the run used (code commit, inputs,
config, environment) goes into the build record, `_build_record/`. `scripts/validate_outputs.R` then
checks that release: the expected tables are present, and the locus and gene
counts match `AD_LOCI_EXPECT_LOCI` / `AD_LOCI_EXPECT_GENES` when those are set.

The output directory is deliberately created empty. The GWAS credible-set
harmonization recomputes from its inputs rather than reusing a previous
release's artifacts, so a rerun cannot silently inherit stale intermediate
files.

## Confidence tiers

Tiers T1–T5 need localized cis-xQTL support; genes implicated only by trans-xQTL
are reported for reference. AD GWAS credible sets (CS) come from SuSiE-rss.
xQTL credible sets come from single-context SuSiE, regional fSuSiE, and
multi-context or multi-gene mvSuSiE. A 95% CS overlap is at least one variant
shared by an AD GWAS 95% CS and an xQTL 95% CS. CoS support is AD and a cis
molecular trait of the gene in the same ColocBoost event with a 95% CoS.
Gene-level TWAS, MR and cTWAS evidence follows the regulome-wide association
companion study. Rules are checked in order.

| tier | evidence |
|---|---|
| T1 | 95% CS overlap and CoS support, with MR or cTWAS |
| T2 | 95% CS overlap or CoS support, with MR or cTWAS |
| T3 | 95% CS overlap or CoS support, with TWAS |
| T4 | CoS support, and an xQTL 95% CS for the gene in a colocalized context (same context, any dataset) that shares a variant with the CoS; no separate AD GWAS CS needed |
| T5 | other localized evidence: CoS alone, or an AD GWAS 95% CS overlapping an xQTL CS at 95% or 70% |
| T6 | gene-level evidence only: TWAS, MR, or cis quantile/interaction QTL at variants in an AD GWAS 95% CS |

`scripts/gene_prio_utils.R` assigns `top_confidence` per row during the build.
Evidence outside these rules gives no tier: 50% credible sets, AD credible sets
below 95%, single-nucleus sQTL association without colocalization, and trans
evidence. Genes with TWAS, MR or cTWAS results but no row in the xQTL overlap
table are given T6 from the gene-level tables produced earlier in the same run.

A gene is reported at its strongest tier. The app takes its tiers from the
release, so the table and the Explorer cannot diverge.

## Running the app

```r
shiny::runApp("app")
```

`app/data.csv` ships with the repository, so the Explorer runs without the
pipeline. Note that `app/build_shiny_data.R` carries a number of gene-level
columns forward from the existing `app/data.csv`, so that file is an input as
well as an output.

## Data

Derived from ADSP/NIAGADS study data. Use of the underlying controlled-access
resources is governed by their own data use terms.
