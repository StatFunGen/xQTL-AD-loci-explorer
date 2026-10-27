# Changing the pipeline and testing the change

This is the loop for changing the build code and showing what the change does: set up once,
run a baseline from the unchanged code, then run your change and compare the two builds.
The example uses the Columbia cluster (Slurm); on BU SCC use `scripts/run_pipeline.qsub`
instead (see the README).

## 1. Set up once

```bash
git clone https://github.com/StatFunGen/xQTL-AD-loci-explorer
cd xQTL-AD-loci-explorer
curl -fsSL https://pixi.sh/install.sh | bash   # only if pixi is not installed yet
export PATH=$HOME/.pixi/bin:$PATH
pixi install                                   # R and packages at the pinned versions
RS=$PWD/.pixi/envs/default/bin/Rscript
LIST=data/build_inputs/20261027.2.inputs.csv   # the input list of the build you start from
ROOT=/your/path/ad_loci_root                   # a new, empty folder
```

One input path contains a masked user name (`<collaborator>`). Ask the maintainer for the
real name and put it in `~/.ad_loci_collaborator` (one line). It is kept out of the repository.

Where the inputs are on the Columbia cluster (each file's exact path is in the `columbia_path`
column of the input list, so nothing has to be looked up by hand):

| location | what |
|---|---|
| `/mnt/lustre/lab/gwang/mmcloud_2026/` | the original analysis results and references, as released |
| `/mnt/lustre/home/jae2200/AD_loci_inputs/` | the processed inputs (made by `scripts/make_processed_inputs.sh`) and the three staging tables |
| `/mnt/lustre/lab/gwang/users/rf2872/` | two glycoQTL files kept by their owner |

Make the project folder. It contains only links to each input's location on the cluster
(the `columbia_path` column of the input list); nothing is copied:

```bash
$RS scripts/make_input_root.R $LIST $ROOT
AD_LOCI_STAGING=$ROOT/staging_gene_prio $RS scripts/check_inputs.R $LIST $ROOT config scripts/build_AD_locus_table.R
```

The check must end with `[check_inputs] OK`: every input is then byte-identical to the one
the release used. If it reports a changed or missing file, stop and tell the maintainer.

## 2. Run a baseline from the unchanged code

```bash
sbatch --export=ALL,AD_LOCI_ROOT=$ROOT,AD_LOCI_INPUTS=$LIST,AD_LOCI_OUT=/your/path/out_baseline,AD_LOCI_RSCRIPT=$RS scripts/run_pipeline.slurm
```

The job needs about 64 GB of memory and takes about 1.5 hours. It writes
`out_baseline/_COMPLETE` only if the input check, the build, the output validation and the
build record all succeed. Its log is `adloci_e2e_<job id>.log` in the submit folder.

To confirm the baseline reproduces the release, compare it with the release outputs:

```bash
$RS scripts/compare_builds.R /path/to/release_outputs /your/path/out_baseline
```

## 3. Make a change and test it

```bash
git checkout -b my-change
# edit scripts/..., then
git commit -am "What the change does and why"
sbatch --export=ALL,AD_LOCI_ROOT=$ROOT,AD_LOCI_INPUTS=$LIST,AD_LOCI_OUT=/your/path/out_my_change,AD_LOCI_RSCRIPT=$RS scripts/run_pipeline.slurm
$RS scripts/compare_builds.R /your/path/out_baseline /your/path/out_my_change
```

- Commit before submitting. The build record stores the commit and flags uncommitted edits
  (`code_has_local_changes` in `_build_record/run.csv`), so an uncommitted run cannot be traced.
- Use a new output folder for every run; the job refuses a folder that already has files.
- `compare_builds.R` lists every output file that differs and exits with status 1 if any does.
  No difference means the change only reorganized code. Otherwise
  `out_my_change/_build_record/changes_since_previous.txt` (submit with
  `AD_LOCI_PREV=/your/path/out_baseline` in `--export`) summarizes the changed counts.

## 4. Hand it over

Push the branch and open a pull request. Say in the description what changed in the results
(from `compare_builds.R` and `changes_since_previous.txt`). The maintainer reruns the branch on
BU SCC and checks that it gives the same outputs as your run.

## If you change the inputs

Inputs are fixed by the input list, not by the code. To add or replace an input, add a new dated
list in `data/build_inputs/` with the file's location and md5 (the README explains the columns
and how processed inputs are made), rather than editing the list a release was built from.
