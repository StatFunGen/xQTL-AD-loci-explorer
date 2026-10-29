#!/usr/bin/env Rscript
## Check a release directory against the expectations of the 188-loci build.
## usage: Rscript validate_outputs.R <release_dir>
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript validate_outputs.R <release_dir>", call. = FALSE)
release <- args[1]
if (!dir.exists(release)) stop("no such release directory: ", release, call. = FALSE)

fail <- character()
note <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fail <<- c(fail, msg)
}

cat("release:", release, "\n\n")

## ---- expected outputs ------------------------------------------------------
expected <- c("res_AD_variants_xQTL.csv.gz",
              "res_all_single_gwas_finemapping_cs50orgreater.csv.gz",
              "res_coloc_AD_xQTL_unified.csv.gz",
              "res_AD_XWAS.csv.gz",
              "res_AD_cTWAS_pip075.csv.gz")
for (f in expected) note(file.exists(file.path(release, f)), paste("present:", f))

xl <- list.files(release, pattern = "^unified_AD_loci_xQTL_summary.*\\.xlsx$")
note(length(xl) > 0, "present: unified_AD_loci_xQTL_summary*.xlsx")

## ---- content checks --------------------------------------------------------
vf <- file.path(release, "res_AD_variants_xQTL.csv.gz")
if (file.exists(vf)) {
  v <- fread(vf, select = c("locus_index", "gene_name", "top_confidence"))
  n_loci <- uniqueN(v[!is.na(locus_index)]$locus_index)
  ## expected counts are given per build (e.g. from the release being reproduced), never hard-coded;
## without them the counts are reported and the build record compares them with the previous build
.el <- suppressWarnings(as.integer(Sys.getenv("AD_LOCI_EXPECT_LOCI", "")))
if (!is.na(.el)) note(n_loci == .el, sprintf("locus count is %d (found %d)", .el, n_loci)) else cat(sprintf("[info] locus count %d\n", n_loci))

  t <- v[!is.na(top_confidence) & top_confidence != "" &
         !is.na(gene_name) & gene_name != ""]
  g <- t[, .(tier = min(top_confidence)), by = .(gene = gene_name)]
  .eg <- suppressWarnings(as.integer(Sys.getenv("AD_LOCI_EXPECT_GENES", "")))
if (!is.na(.eg)) note(nrow(g) == .eg, sprintf("tiered genes is %d (found %d)", .eg, nrow(g))) else cat(sprintf("[info] tiered genes %d\n", nrow(g)))

  tb <- table(g$tier)
  cat("\n  tier distribution:", paste(sprintf("%s=%d", names(tb), tb), collapse = " "), "\n")
  note("T6" %in% names(tb) && tb[["T6"]] > 0,
       "T6 is populated (genes with gene-level evidence only)")
  cat("\n")
}

## ---- every AD-xQTL colocalization input reaches the coloc table -----------
root <- Sys.getenv("AD_LOCI_ROOT"); cfg <- Sys.getenv("AD_LOCI_CONFIG")
cf <- file.path(release, "res_coloc_AD_xQTL_unified.csv.gz")
if (nzchar(root) && nzchar(cfg) && file.exists(cf)) {
  mtd <- fread(file.path(cfg, "metadata_analysis.csv"))
  got <- unique(fread(cf, select = "gene_ID")$gene_ID)
  for (f in mtd[Method == "ColocBoost" & Modality == "AD_xQTL_colocalization"]$Path) {
    g <- unique(fread(file.path(root, f), select = "region_ID")$region_ID)
    miss <- setdiff(g, got)
    note(!length(miss), sprintf("coloc input %s: %d of %d genes in the coloc table",
                                basename(f), length(g) - length(miss), length(g)))
  }
} else cat(" [skip] coloc input check: set AD_LOCI_ROOT and AD_LOCI_CONFIG\n")

## ---- tiers agree with the evidence they are built from -------------------
of <- file.path(release, "xQTL_all_methods_overlap_with_AD_loci_unified_cs95orColocs_Pval1e5.csv.gz")
if (file.exists(of)) {
  x <- fread(of, select = c("ADlocus", "gene_ID", "Method", "context", "context_short",
                            "credibleset", intersect(c("susie_coverage", "is.cs95.var"), names(fread(of, nrows = 0))), "MR_signif", "cTWAS_signif",
                            "TWAS_signif", "confidence_lvl"))
  x <- x[!grepl("^AD", context) & !grepl("trans", Method) & !grepl("^tmAD_", gene_ID)]
  fm <- c("single_context_finemapping", "multi_context_finemapping",
          "multi_gene_finemapping") ## epigenomic (fSuSiE) sets do not count at gene level
  x[, ov95 := ((Method %in% fm & grepl("^cs95", credibleset))) &
              (if ("is.cs95.var" %in% names(x)) is.cs95.var %in% TRUE else susie_coverage %in% "cs95")]
  x[, causal := MR_signif %in% TRUE | cTWAS_signif %in% TRUE]
  x[, twas := TWAS_signif %in% TRUE]
  ## gene has a 95% CS overlap and MR/cTWAS in the same context -> at most T2
  e <- x[, .(ov = any(ov95), c = any(causal), t = any(twas)), by = .(gene_ID, context_short)]
  need2 <- unique(e[ov & c]$gene_ID); need3 <- unique(e[ov & t]$gene_ID)
  tier <- x[!is.na(confidence_lvl) & confidence_lvl != "",
            .(tier = min(as.character(confidence_lvl))), by = gene_ID]
  bad2 <- setdiff(need2, tier[tier %in% c("T1", "T2")]$gene_ID)
  bad3 <- setdiff(need3, tier[tier %in% c("T1", "T2", "T3")]$gene_ID)
  note(!length(bad2), sprintf("overlap + MR/cTWAS in the same context gives T1/T2 (%d genes do not)", length(bad2)))
  note(!length(bad3), sprintf("overlap + TWAS in the same context gives T1-T3 (%d genes do not)", length(bad3)))
  if (length(bad2) || length(bad3)) cat("   e.g.", head(c(bad2, bad3), 10), "\n")
}

## ---- no locus overlap file is older than the table it was built from ------
ov <- list.files(release, pattern = "_overlapADloci\\.csv\\.gz$", full.names = TRUE)
src <- sub("_overlapADloci\\.csv\\.gz$", ".csv.gz", ov)
stale <- ov[file.exists(src) & file.mtime(src) > file.mtime(ov)]
note(!length(stale), sprintf("no stale overlap files (%d of %d older than their source)", length(stale), length(ov)))
if (length(stale)) cat("   stale:", basename(stale), sep = "\n    ")

if (length(fail)) {
  cat(sprintf("FAILED %d check(s)\n", length(fail)))
  quit(status = 1)
}
cat("all checks passed\n")
