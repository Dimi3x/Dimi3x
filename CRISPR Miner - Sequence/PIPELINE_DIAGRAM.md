Pipeline Diagram
================

ASCII flow
----------
```
MAG FASTA(s)
    |
    v
Prodigal (ORFs) ------------------------------+
    |                                         |
    v                                         |
CRONUS (noncanonical arrays)                  |
    |                                         |
    v                                         |
CRISPRCasFinder (canonical arrays)            |
    |                                         |
    v                                         |
crispr_operon_miner.py (arrays + ORFs) <-------+
    |
    v
operon_proteins.faa
    |
    v
MMseqs2 GPU (tight) -> clusters_tight.tsv + reps_tight.faa
    |
    v
MMseqs2 GPU (coarse) -> clusters_loose.tsv + reps_loose.faa
    |
    v
MAFFT + FastTree -> reps_loose.aln.faa + reps_loose.tree.nwk
    |
    +--> (optional) MMseqs2 search vs known DB -> known_hits.tsv
    +--> (optional) ESM3 embeddings -> esm3_embeddings/
    +--> (optional) AlphaFold3 models -> af3_models/ -> FoldSeek hits
```

Mermaid flow (optional)
-----------------------
```mermaid
flowchart TD
    A[MAG FASTA(s)] --> B[Prodigal ORF calling]
    B --> C[CRONUS arrays\n(noncanonical)]
    B --> D[CRISPRCasFinder arrays\n(canonical)]
    C --> E[crispr_operon_miner.py\n(assemble operons)]
    D --> E
    E --> F[operon_proteins.faa]
    F --> G[MMseqs2 GPU\nTight clustering]
    G --> H[cluster_reps_tight.faa]
    H --> I[MMseqs2 GPU\nCoarse clustering]
    I --> J[cluster_reps_loose.faa]
    J --> K[MAFFT + FastTree\nPhylogeny]
    I --> L[(optional) Known DB search]
    J --> M[(optional) ESM3 embeddings]
    J --> N[(optional) AlphaFold3 models]
    N --> O[(optional) FoldSeek search]
```
