"""Gene-level TWAS / MR / cTWAS flags per AD locus, read from the released summary workbook.

The app's data.csv carries these flags per variant row, and for most genes the
significant contexts sit on rows the app does not show, so a gene such as BIN1
read "not recorded" although the release marks it TWAS and MR significant.
The workbook's "xQTLs summary per gene" block is the released, gene-level answer.

Usage: python scripts/gene_evidence_from_release.py data/<build>/unified_AD_loci_xQTL_summary_<build>.xlsx app/data_aux/gene_evidence.csv
"""
import sys
import pandas as pd

src, out = sys.argv[1], sys.argv[2]
v = pd.read_excel(src, sheet_name="Gene Locus table", engine="calamine", header=None)
hdr = v.iloc[2].tolist()
col = {name: hdr.index(name) for name in
       ["xQTL target gene", "TWAS significant", "MR significant", "cTWAS significant"]}
d = v.iloc[3:]
d = pd.DataFrame({
    "ADlocus": d[1],
    "gene": d[col["xQTL target gene"]],
    "twas": d[col["TWAS significant"]] == True,
    "mr": d[col["MR significant"]] == True,
    "ctwas": d[col["cTWAS significant"]] == True,
})
d = d[d["gene"].notna() & d["ADlocus"].notna()]
d = d.groupby(["ADlocus", "gene"], as_index=False).any()
d.to_csv(out, index=False)
print(f"{len(d)} locus-gene pairs; TWAS {d.twas.sum()}, MR {d.mr.sum()}, cTWAS {d.ctwas.sum()}")
