#!/usr/bin/env python3
"""Sum Salmon transcript-level estimated counts to gene level and write a DESeq2-ready matrix.

Gene IDs come from the Ensembl transcript FASTA headers (``gene:ENSG...``), with
the version suffix removed so they match the GTF ``gene_id``. Summing NumReads per
gene and rounding is what tximport does with countsFromAbundance = "no".

    salmon_gene_counts.py --transcripts transcripts.fa.gz --out counts.tsv sampleA/ sampleB/ ...
"""
import argparse
import gzip
import os
import re
import sys
from collections import defaultdict


def tx2gene(fasta):
    opener = gzip.open if fasta.endswith(".gz") else open
    mapping = {}
    with opener(fasta, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                tx = line[1:].split()[0]
                m = re.search(r"\bgene:(\S+)", line)
                gene = m.group(1) if m else tx
                mapping[tx] = gene.split(".")[0]
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--transcripts", required=True)
    ap.add_argument("--out", default="counts.tsv")
    ap.add_argument("quant_dirs", nargs="+", help="Salmon output folders, one per sample (folder name = sample)")
    a = ap.parse_args()

    t2g = tx2gene(a.transcripts)
    samples = sorted(os.path.basename(os.path.normpath(d)) for d in a.quant_dirs)
    dirs = {os.path.basename(os.path.normpath(d)): d for d in a.quant_dirs}
    counts = defaultdict(lambda: [0.0] * len(samples))
    unmapped = 0
    for j, s in enumerate(samples):
        with open(os.path.join(dirs[s], "quant.sf")) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_name, i_reads = header.index("Name"), header.index("NumReads")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                gene = t2g.get(f[i_name])
                if gene is None:
                    unmapped += 1
                    gene = f[i_name]
                counts[gene][j] += float(f[i_reads])
    with open(a.out, "w") as out:
        out.write("gene_id\t" + "\t".join(samples) + "\n")
        for gene in sorted(counts):
            out.write(gene + "\t" + "\t".join(str(int(round(v))) for v in counts[gene]) + "\n")
    print(f"{len(counts):,} genes x {len(samples)} samples written to {a.out}"
          + (f" ({unmapped} transcript rows had no gene in the FASTA header)" if unmapped else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
