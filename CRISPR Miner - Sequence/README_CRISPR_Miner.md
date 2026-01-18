CRISPR-GPU Miner (standalone)
=============================

Purpose
-------
- Standalone Bash + Python pipeline to mine novel CRISPR/Cas operons from assembled MAGs.
- Uses MMseqs2 GPU clustering for sensitivity; Python handles CRISPR array detection, ORF calls, operon assembly, and ranking.
- Outputs full operons: CRISPR array + flanking proteins + adjacency metadata for downstream FoldSeek/ESM3 or database building.

Pipeline overview
-----------------
- Prodigal meta-mode ORF calling per MAG to improve protein accuracy (required in the pipeline; miner can fall back only if run standalone).
- Parse assembled MAG FASTA files.
- Detect CRISPR arrays using repeat/spacer heuristics (min 3 repeats, repeat 23-47 bp, spacers 20-60 bp); CRONUS (required) for noncanonical/imperfect arrays and CRISPRCasFinder (required) for canonical arrays are both ingested by the pipeline.
- If running the miner standalone without Prodigal: it can call ORFs on both strands with the built-in lightweight caller.
- Assemble operons: arrays plus ORFs within a configurable flank window (default 10 kbp).
- Score/rank operons by array quality, proximity and density of ORFs, length sanity, and repeat/spacer consistency.
- Two-tier MMseqs2 GPU clustering for sensitivity and deduplication:
  - Tight pass (e.g., 0.9 id/0.8 cov) to collapse near-duplicates.
  - Coarse pass (e.g., 0.5 id/0.33 cov) on tight representatives for family-level grouping.
- Generate cluster representatives (tight and coarse) for downstream modeling/search.
- Optional GPU search against a provided CRISPR/Cas reference DB to tag known families and highlight novel clusters.
- Optional downstream passes:
  - ESM3 embeddings on coarse reps to prioritize novel candidates (if `--esm3-cmd` or `esm3_embed` is available).
  - AlphaFold3 on coarse reps to obtain structures (if `--alphafold3-cmd` is provided), followed by FoldSeek if installed.
  - Phylogenetic tree on coarse reps (MAFFT + FastTree required; built automatically).

Key files
---------
- `crispr_operon_miner.py`: Python core for array ingestion, operon scoring, JSON/TSV/FASTA outputs.
- `crispr_gpu_pipeline.sh`: Full pipeline driver (Prodigal, CRONUS, CRISPRCasFinder, MMseqs2 GPU, tree, optional extras).
- `run_pipeline_steps.sh`: Minimal step-by-step runner for core stages.
- `summarize_crispr_results.py`: Optional plotting helper.
- `PIPELINE_DIAGRAM.md`: ASCII + Mermaid pipeline diagram.

External tools (not bundled)
----------------------------
- Prodigal (required).
- CRONUS (required).
- CRISPRCasFinder (required).
- MMseqs2 GPU build (CUDA-enabled).
- MAFFT + FastTree/FastTreeMP.
- Optional: ESM3, AlphaFold3, FoldSeek.
Install these separately and ensure they’re on PATH or referenced explicitly.

Ranking rationale
-----------------
- Array quality (repeat length in canonical range, ≥3 spacers, spacer length variability) dominates the score.
- Proximity and density of ORFs around the array are boosted; unrealistic protein lengths penalized.
- Optional MMseqs2 hit density to known DB (if provided) is recorded so novel clusters stand out.
- This ranking helps triage candidates; keep it but always review high-scoring outliers.

Usage (quick start)
-------------------
```bash
# Requires: python3, prodigal, CRONUS, CRISPRCasFinder, mmseqs (GPU build), mafft, fasttree/FastTree
./crispr_gpu_pipeline.sh /path/to/mags /path/to/output \
  --known-db /path/to/crispr_db  # optional MMseqs2 DB
# optional extras (recommended if available):
#   --alphafold3-cmd /path/to/af3_runner
#   --foldseek-target /path/to/structural_db
#   --mafft-cmd /path/to/mafft
#   --fasttree-cmd /path/to/fasttree
#   --esm3-cmd /path/to/esm3_embed
#   --casfinder-cmd /path/to/CRISPRCasFinder.pl
#   --cronus-cmd /path/to/cronus
```
- Outputs:
  - `operons.json`: structured operon calls with scores.
  - `operons.tsv`: tabular summary for quick filtering.
  - `operon_proteins.faa`: proteins near arrays (for clustering/search).
  - `clusters_tight.tsv` and `cluster_reps_tight.faa`: near-duplicate collapsed clusters/reps.
  - `clusters_loose.tsv` and `cluster_reps_loose.faa`: coarse (family-level) clusters/reps.
  - `known_hits.tsv`: optional GPU search annotations if a DB is given.
  - `tmp/prodigal/*.faa|gff`: Prodigal outputs (kept for reuse).
  - `cluster_reps_loose.faa`: coarse cluster representative sequences (for modeling/search).
  - `CRISPRs_REPORT.combined.tsv`: aggregated CRISPRCasFinder output when provided/available.
  - `cronus_combined.tsv`: aggregated CRONUS output (required; produced by pipeline).
  - `esm3_embeddings/`: ESM3 embedding outputs when `--esm3-cmd` (or `esm3_embed` on PATH) is provided.
  - `af3_models/`: AlphaFold3 models when `--alphafold3-cmd` is provided.
  - `foldseek_hits.tsv` or `foldseek_self.tsv`: FoldSeek results on AlphaFold3 models (target DB or self-search).
  - `cluster_reps_loose.aln.faa` and `cluster_reps_loose.tree.nwk`: alignment and phylogenetic tree of coarse reps.

Quick visualization (optional)
------------------------------
- Generate summary plots from outputs:
  ```bash
  python3 summarize_crispr_results.py \
    --clusters-loose /path/to/output/clusters_loose.tsv \
    --known-hits /path/to/output/known_hits.tsv \   # optional
    --foldseek /path/to/output/foldseek_hits.tsv \  # optional
    --esm3-dir /path/to/output/esm3_embeddings \    # optional
    --outdir /path/to/output/plots
  ```
- Produces:
  - `cluster_sizes.png`: coarse cluster size distribution.
  - `novelty_quadrant.png`: counts of reps with sequence/structure hits vs none.
  - `esm3_umap.png`: UMAP of ESM3 embeddings (if embeddings present and umap+matplotlib installed).

Flags worth noting
------------------
- `--known-db`: path to MMseqs2 DB for known CRISPR/Cas tagging (optional).
- `--alphafold3-cmd`: path to your AlphaFold3 runner; if provided, models cluster reps automatically.
- `--foldseek-target`: structural DB for FoldSeek search on AlphaFold3 models; if omitted and FoldSeek is installed, self-search is run.
- `--mafft-cmd` / `--fasttree-cmd`: paths to alignment/tree tools; must resolve (or be on PATH) because tree building is mandatory.
- `--esm3-cmd`: path to ESM3 embedding runner; if provided (or `esm3_embed` on PATH), embeddings are generated for cluster reps.
- `--casfinder-cmd`: path to CRISPRCasFinder.pl; required (or on PATH), arrays are ingested for canonical arrays.
- `--cronus-cmd`: path to CRONUS; required (or on PATH), arrays are ingested for noncanonical/imperfect arrays.

Integration hooks
-----------------
- FoldSeek: run on AlphaFold3 models (coarse reps) when available.
- ESM3: embed coarse reps to prioritize novelty.
- MMseqs2: already used for GPU clustering; reuse created DBs for further searches.
- AlphaFold3: model coarse cluster representatives; FoldSeek on resulting structures.

Optional post-processing recipes (manual steps)
-----------------------------------------------
- ESM3 embeddings (example; adapt to your install):
  ```bash
  esm3_embed --fasta cluster_reps_loose.faa --outdir esm3_embeddings --device cuda
  ```
- AlphaFold3 modeling: run externally on coarse reps; output PDBs to `af3_models/`.
- FoldSeek on AlphaFold3 models:
  ```bash
  foldseek createdb af3_models/ af3_db
  foldseek search af3_db af3_db af3_self aln_tmp --format-mode 4
  foldseek convertalis af3_db af3_db af3_self af3_self.tsv
  ```
  Swap the target DB in `foldseek search` to query against a structural reference database instead of self-search.
