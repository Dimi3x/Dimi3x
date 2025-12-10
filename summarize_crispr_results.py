#!/usr/bin/env python3
"""
Lightweight visualization helper for CRISPR-GPU miner outputs.

Generates summary plots:
- Cluster size histogram (coarse clusters).
- Novelty quadrant: known (sequence) vs structure (FoldSeek) hits per cluster representative.
- Optional UMAP of ESM3 embeddings if available.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict


def read_clusters(path: str) -> Counter:
    sizes = Counter()
    with open(path, "r", encoding="ascii") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            rep, member = row[0], row[1]
            sizes[rep] += 1 if member else 0
            if member != rep:
                sizes[rep] += 0  # member already counted
    return sizes


def reps_with_hits(path: str) -> set[str]:
    reps: set[str] = set()
    if not path or not os.path.isfile(path):
        return reps
    with open(path, "r", encoding="ascii") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            reps.add(row[0])
    return reps


def load_embeddings(path: str):
    import numpy as np

    names = []
    mats = []
    for fname in os.listdir(path):
        if not fname.endswith(".npy"):
            continue
        names.append(os.path.splitext(fname)[0])
        mats.append(np.load(os.path.join(path, fname)))
    if not mats:
        return None, None
    return names, np.vstack(mats)


def plot_cluster_sizes(sizes: Counter, out_png: str):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping cluster size plot", file=sys.stderr)
        return
    vals = list(sizes.values())
    plt.figure(figsize=(6, 4))
    plt.hist(vals, bins=50, color="#4e79a7")
    plt.xlabel("Cluster size (coarse reps)")
    plt.ylabel("Count")
    plt.title("Coarse cluster size distribution")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def plot_novelty_quadrant(reps: set[str], seq_hits: set[str], struct_hits: set[str], out_png: str):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping novelty quadrant", file=sys.stderr)
        return
    quad = defaultdict(int)
    for rep in reps:
        s = rep in seq_hits
        t = rep in struct_hits
        if s and t:
            quad["Seq+ / Struct+"] += 1
        elif s and not t:
            quad["Seq+ / Struct-"] += 1
        elif not s and t:
            quad["Seq- / Struct+"] += 1
        else:
            quad["Seq- / Struct-"] += 1
    labels = ["Seq+ / Struct+", "Seq+ / Struct-", "Seq- / Struct+", "Seq- / Struct-"]
    counts = [quad[l] for l in labels]
    plt.figure(figsize=(5, 4))
    plt.bar(labels, counts, color=["#59a14f", "#edc948", "#76b7b2", "#e15759"])
    plt.ylabel("Rep count")
    plt.xticks(rotation=20)
    plt.title("Novelty quadrant (coarse reps)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def plot_umap(names, embeddings, out_png):
    try:
        import umap
        import numpy as np
        import matplotlib.pyplot as plt
    except ImportError:
        print("umap-learn or matplotlib not installed; skipping UMAP plot", file=sys.stderr)
        return
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, metric="cosine")
    coords = reducer.fit_transform(embeddings)
    plt.figure(figsize=(6, 5))
    plt.scatter(coords[:, 0], coords[:, 1], s=10, alpha=0.7)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("ESM3 embeddings (coarse reps)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Summarize CRISPR-GPU miner outputs with quick plots.")
    parser.add_argument("--clusters-loose", required=True, help="Path to clusters_loose.tsv (rep/member pairs).")
    parser.add_argument("--known-hits", help="known_hits.tsv from mmseqs search.")
    parser.add_argument("--foldseek", help="foldseek_hits.tsv or foldseek_self.tsv.")
    parser.add_argument("--esm3-dir", help="Directory with ESM3 embeddings (.npy per rep).")
    parser.add_argument("--outdir", required=True, help="Directory to write plots.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sizes = read_clusters(args.clusters_loose)
    plot_cluster_sizes(sizes, os.path.join(args.outdir, "cluster_sizes.png"))

    reps = set(sizes.keys())
    seq_hits = reps_with_hits(args.known_hits)
    struct_hits = reps_with_hits(args.foldseek)
    plot_novelty_quadrant(reps, seq_hits, struct_hits, os.path.join(args.outdir, "novelty_quadrant.png"))

    if args.esm3_dir:
        names, emb = load_embeddings(args.esm3_dir)
        if names is not None:
            plot_umap(names, emb, os.path.join(args.outdir, "esm3_umap.png"))


if __name__ == "__main__":
    main()
