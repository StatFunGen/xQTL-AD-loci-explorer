#!/bin/bash
## Slurm version of run_pipeline.qsub (e.g. the Columbia cluster). It sets the same
## variables the SGE job gets and then runs run_pipeline.qsub, so both clusters run the
## identical steps: snapshot, input check, build, validation, build record.
##   Rscript scripts/make_input_root.R data/build_inputs/YYYYMMDD.inputs.csv /path/to/input_root
##   sbatch --export=ALL,AD_LOCI_ROOT=/path/to/input_root,AD_LOCI_INPUTS=data/build_inputs/YYYYMMDD.inputs.csv,AD_LOCI_OUT=/path/to/out_YYYYMMDD scripts/run_pipeline.slurm
## Submit from a clean checkout of the commit you want to build, with the pixi environment
## of this repository on PATH (or AD_LOCI_RSCRIPT pointing at its Rscript).
#SBATCH -J adloci_e2e
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o adloci_e2e_%j.log
: "${AD_LOCI_ROOT:?set AD_LOCI_ROOT to the input folder made by scripts/make_input_root.R}"
export NSLOTS=${SLURM_CPUS_PER_TASK:-1}
export AD_LOCI_REPO=${AD_LOCI_REPO:-$SLURM_SUBMIT_DIR}
export AD_LOCI_RSCRIPT=${AD_LOCI_RSCRIPT:-$(command -v Rscript)}
exec bash "$AD_LOCI_REPO/scripts/run_pipeline.qsub"
