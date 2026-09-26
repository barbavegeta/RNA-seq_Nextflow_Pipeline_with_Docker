#!/usr/bin/env python3
"""Check the airway DESeq2 results against well-known glucocorticoid-response genes.

Dexamethasone strongly induces these genes in airway smooth muscle cells
(Himes et al. 2014, PLoS One; the Bioconductor 'airway' workflow). If the
pipeline is correct they should all be up-regulated (log2FC > 1, padj < 0.05).
Also writes a labelled volcano plot and the top-20 table with gene symbols.

    python airway/check_known_genes.py --results results_airway/deseq2_results.tsv \
        --gtf data/airway/ref/Homo_sapiens.GRCh38.112.gtf --out results_airway
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

KNOWN_UP = ["FKBP5", "TSC22D3", "ZBTB16", "PER1", "KLF15", "DUSP1", "CRISPLD2", "SPARCL1"]


def gene_names(gtf):
    ids, names = [], []
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8])
            gn = re.search(r'gene_name "([^"]+)"', f[8])
            if gid:
                ids.append(gid.group(1))
                names.append(gn.group(1) if gn else gid.group(1))
    return pd.Series(names, index=ids, name="gene_name")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    res = pd.read_csv(a.results, sep="\t")
    res = res.join(gene_names(a.gtf), on="gene_id")
    res.to_csv(out / "deseq2_results_with_symbols.tsv", sep="\t", index=False)

    sig = res[(res.padj < 0.05) & (res.log2FoldChange.abs() > 1)]
    print(f"{len(res):,} genes tested; {len(sig):,} with padj < 0.05 and |log2FC| > 1 "
          f"({(sig.log2FoldChange > 0).sum():,} up, {(sig.log2FoldChange < 0).sum():,} down with dex)\n")

    known = res[res.gene_name.isin(KNOWN_UP)].set_index("gene_name").reindex(KNOWN_UP)
    known["pass"] = (known.log2FoldChange > 1) & (known.padj < 0.05)
    print(known[["log2FoldChange", "padj", "pass"]].to_string(float_format=lambda v: f"{v:.3g}"))
    n_pass = int(known["pass"].sum())
    print(f"\n{n_pass}/{len(KNOWN_UP)} known dexamethasone-response genes recovered")
    if n_pass < len(KNOWN_UP) - 1:
        print("WARNING: fewer than expected. Check the contrast direction (--deseq2_reference control) "
              "and the design (--deseq2_formula '~ cell + condition').")
    known.to_csv(out / "known_gene_check.tsv", sep="\t")

    top = res.dropna(subset=["padj"]).sort_values("padj").head(20)
    top[["gene_id", "gene_name", "baseMean", "log2FoldChange", "padj"]].to_csv(out / "top20_genes.tsv", sep="\t", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = res.dropna(subset=["padj"])
    y = -np.log10(d.padj.clip(lower=1e-300))
    up = (d.padj < 0.05) & (d.log2FoldChange > 1)
    down = (d.padj < 0.05) & (d.log2FoldChange < -1)
    fig, ax = plt.subplots(figsize=(7, 5), dpi=150)
    ax.scatter(d.log2FoldChange[~(up | down)], y[~(up | down)], s=4, c="#c3c2b7", edgecolors="none", label="not significant")
    ax.scatter(d.log2FoldChange[up], y[up], s=6, c="#eb6834", edgecolors="none", label=f"up with dex ({up.sum():,})")
    ax.scatter(d.log2FoldChange[down], y[down], s=6, c="#2a78d6", edgecolors="none", label=f"down with dex ({down.sum():,})")
    for g, r in d[d.gene_name.isin(KNOWN_UP)].set_index("gene_name").iterrows():
        ax.annotate(g, (r.log2FoldChange, -np.log10(max(r.padj, 1e-300))), xytext=(4, 3),
                    textcoords="offset points", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("log2 fold change (dexamethasone vs control)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Airway smooth muscle cells: dexamethasone response (GSE52778)", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "volcano_airway.png")
    print(f"\nwrote {out / 'volcano_airway.png'}, {out / 'top20_genes.tsv'}, {out / 'known_gene_check.tsv'}")


if __name__ == "__main__":
    main()
