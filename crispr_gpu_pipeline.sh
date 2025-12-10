#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <mag_fasta_dir_or_file> <output_dir> [--known-db /path/to/mmseqs_db] [--alphafold3-cmd /path/to/af3_runner] [--foldseek-target /path/to/struct_db] [--mafft-cmd /path/to/mafft] [--fasttree-cmd /path/to/fasttree] [--esm3-cmd /path/to/esm3_embed] [--casfinder-cmd /path/to/CRISPRCasFinder.pl] [--cronus-cmd /path/to/cronus]" >&2
  echo "Requires: prodigal, CRONUS, CRISPRCasFinder, mmseqs (GPU build), mafft, fasttree/FastTree, python3." >&2
  exit 1
fi

INPUT_PATH="$1"
OUT_DIR="$2"
shift 2 || true

KNOWN_DB=""
FOLDSEEK_TARGET=""
ALPHAFOLD3_CMD=""
MAFFT_CMD=""
FASTTREE_CMD=""
ESM3_CMD=""
CASFINDER_CMD=""
CRONUS_CMD=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --known-db)
      KNOWN_DB="$2"
      shift 2
      ;;
    --casfinder-cmd)
      CASFINDER_CMD="$2"
      shift 2
      ;;
    --esm3-cmd)
      ESM3_CMD="$2"
      shift 2
      ;;
    --mafft-cmd)
      MAFFT_CMD="$2"
      shift 2
      ;;
    --fasttree-cmd)
      FASTTREE_CMD="$2"
      shift 2
      ;;
    --foldseek-target)
      FOLDSEEK_TARGET="$2"
      shift 2
      ;;
    --alphafold3-cmd)
      ALPHAFOLD3_CMD="$2"
      shift 2
      ;;
    --cronus-cmd)
      CRONUS_CMD="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

if ! command -v mmseqs >/dev/null 2>&1; then
  echo "mmseqs not found in PATH" >&2
  exit 1
fi
if [[ -n "${FOLDSEEK_TARGET}" || -n "${ALPHAFOLD3_CMD}" ]]; then
  if ! command -v foldseek >/dev/null 2>&1; then
    echo "[WARN] foldseek not found; FoldSeek steps will be skipped." >&2
  fi
fi

# resolve MAFFT / FastTree if provided or on PATH (required for tree)
if [[ -z "${MAFFT_CMD}" ]]; then
  MAFFT_CMD="$(command -v mafft || true)"
fi
if [[ -z "${FASTTREE_CMD}" ]]; then
  FASTTREE_CMD="$(command -v fasttree || command -v FastTree || true)"
fi
if [[ -z "${MAFFT_CMD}" ]]; then
  echo "mafft is required for alignment/tree and was not found. Install or provide --mafft-cmd." >&2
  exit 1
fi
if [[ -z "${FASTTREE_CMD}" ]]; then
  echo "fasttree/FastTree is required for alignment/tree and was not found. Install or provide --fasttree-cmd." >&2
  exit 1
fi

# CRISPRCasFinder (required)
if [[ -z "${CASFINDER_CMD}" ]]; then
  CASFINDER_CMD="$(command -v CRISPRCasFinder.pl || true)"
fi
if [[ -z "${CASFINDER_CMD}" ]]; then
  echo "CRISPRCasFinder is required and not found on PATH. Install or provide --casfinder-cmd." >&2
  exit 1
fi

# ESM3 runner (optional novelty scoring)
if [[ -z "${ESM3_CMD}" ]]; then
  ESM3_CMD="$(command -v esm3_embed || true)"
fi
# CRONUS (required)
if [[ -z "${CRONUS_CMD}" ]]; then
  CRONUS_CMD="$(command -v cronus || true)"
fi
if [[ -z "${CRONUS_CMD}" ]]; then
  echo "CRONUS is required and not found on PATH. Install or provide --cronus-cmd." >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON=${PYTHON:-python3}

mkdir -p "${OUT_DIR}"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "${TMP_DIR}"

OPERON_JSON="${OUT_DIR}/operons.json"
OPERON_TSV="${OUT_DIR}/operons.tsv"
PROTEINS="${OUT_DIR}/operon_proteins.faa"

PRODIGAL_DIR="${TMP_DIR}/prodigal"
if ! command -v prodigal >/dev/null 2>&1; then
  echo "prodigal is required and not found on PATH." >&2
  exit 1
fi

# total steps for logging (prodigal + cronus + casfinder + tree mandatory)
TOTAL_STEPS=7  # prodigal, cronus, casfinder, miner, tight cluster, loose cluster, tree
[[ -n "${KNOWN_DB}" ]] && TOTAL_STEPS=$((TOTAL_STEPS+1))
[[ -n "${ESM3_CMD}" ]] && TOTAL_STEPS=$((TOTAL_STEPS+1))
[[ -n "${ALPHAFOLD3_CMD}" ]] && TOTAL_STEPS=$((TOTAL_STEPS+1))
STEP=1

mkdir -p "${PRODIGAL_DIR}"
echo "[${STEP}/${TOTAL_STEPS}] Running Prodigal on MAGs..."
mapfile -t FASTA_FILES < <(find "${INPUT_PATH}" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fas" \))
if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
  echo "No FASTA files found under ${INPUT_PATH}" >&2
  exit 1
fi
for f in "${FASTA_FILES[@]}"; do
  base="$(basename "$f")"
  stem="${base%.*}"
  prodigal -i "$f" -a "${PRODIGAL_DIR}/${stem}.faa" -f gff -o "${PRODIGAL_DIR}/${stem}.gff" -p meta
done
STEP=$((STEP+1))

echo "[${STEP}/${TOTAL_STEPS}] Detecting CRISPR arrays and operons with Python..."
MINER_ARGS=(
  --input "${INPUT_PATH}"
  --output-json "${OPERON_JSON}"
  --output-tsv "${OPERON_TSV}"
  --output-proteins "${PROTEINS}"
  --prodigal-dir "${PRODIGAL_DIR}"
)
CRONUS_DIR="${TMP_DIR}/cronus"
mkdir -p "${CRONUS_DIR}"
CRONUS_COMBINED="${CRONUS_DIR}/cronus_combined.tsv"
rm -f "${CRONUS_COMBINED}"
echo "[${STEP}/${TOTAL_STEPS}] Running CRONUS on MAGs..."
header_written=0
mapfile -t CF_FASTAS < <(find "${INPUT_PATH}" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" -o -name "*.fasta" -o -name "*.fas" \))
if [[ ${#CF_FASTAS[@]} -eq 0 ]]; then
  echo "No FASTA files found under ${INPUT_PATH}" >&2
  exit 1
fi
for f in "${CF_FASTAS[@]}"; do
  base="$(basename "$f")"
  stem="${base%.*}"
  out_tsv="${CRONUS_DIR}/${stem}.cronus.tsv"
  # Adjust CLI if your CRONUS invocation differs
  "${CRONUS_CMD}" -i "$f" -o "${out_tsv}" >/dev/null 2>&1 || {
    echo "[ERROR] CRONUS failed on $f" >&2
    exit 1
  }
  if [[ -f "${out_tsv}" ]]; then
    if [[ ${header_written} -eq 0 ]]; then
      head -n 1 "${out_tsv}" > "${CRONUS_COMBINED}"
      header_written=1
    fi
    tail -n +2 "${out_tsv}" >> "${CRONUS_COMBINED}"
  fi
done
if [[ ! -f "${CRONUS_COMBINED}" ]]; then
  echo "[ERROR] CRONUS did not produce arrays; aborting." >&2
  exit 1
fi
MINER_ARGS+=(--cronus-tsv "${CRONUS_COMBINED}")
STEP=$((STEP+1))
CASF_DIR="${TMP_DIR}/casfinder"
mkdir -p "${CASF_DIR}"
COMBINED_TSV="${CASF_DIR}/CRISPRs_REPORT.combined.tsv"
rm -f "${COMBINED_TSV}"
echo "[${STEP}/${TOTAL_STEPS}] Running CRISPRCasFinder on MAGs..."
header_written=0
for f in "${CF_FASTAS[@]}"; do
  base="$(basename "$f")"
  outdir="${CASF_DIR}/${base%.*}"
  mkdir -p "${outdir}"
  "${CASFINDER_CMD}" -in "$f" -out "${outdir}" >/dev/null 2>&1 || {
    echo "[ERROR] CRISPRCasFinder failed on $f" >&2
    exit 1
  }
  report="${outdir}/CRISPRs_REPORT.tsv"
  if [[ -f "${report}" ]]; then
    if [[ ${header_written} -eq 0 ]]; then
      head -n 1 "${report}" > "${COMBINED_TSV}"
      header_written=1
    fi
    tail -n +2 "${report}" >> "${COMBINED_TSV}"
  fi
done
if [[ ! -f "${COMBINED_TSV}" ]]; then
  echo "[ERROR] CRISPRCasFinder did not produce arrays; aborting." >&2
  exit 1
fi
MINER_ARGS+=(--casfinder-tsv "${COMBINED_TSV}")
STEP=$((STEP+1))
${PYTHON} "${SCRIPT_DIR}/crispr_operon_miner.py" "${MINER_ARGS[@]}"
STEP=$((STEP+1))

echo "[${STEP}/${TOTAL_STEPS}] MMseqs2 GPU linclust (tight) on operon proteins..."
PROT_DB="${TMP_DIR}/proteins.db"
CLUST_TIGHT="${TMP_DIR}/clusters_tight"
mmseqs createdb "${PROTEINS}" "${PROT_DB}"
mmseqs linclust "${PROT_DB}" "${CLUST_TIGHT}" "${TMP_DIR}/linclust_tight_tmp" \
  --min-seq-id 0.9 --cov-mode 1 -c 0.8 --cluster-mode 2 --gpu 1
mmseqs createtsv "${PROT_DB}" "${PROT_DB}" "${CLUST_TIGHT}" "${OUT_DIR}/clusters_tight.tsv"
STEP=$((STEP+1))

# Representatives from tight clustering
REP_DB_TIGHT="${TMP_DIR}/repseqs_tight"
REP_FASTA_TIGHT="${OUT_DIR}/cluster_reps_tight.faa"
mmseqs createrepseqdb "${PROT_DB}" "${CLUST_TIGHT}" "${REP_DB_TIGHT}"
mmseqs convert2fasta "${REP_DB_TIGHT}" "${REP_FASTA_TIGHT}"

# Coarse clustering on tight reps
echo "[${STEP}/${TOTAL_STEPS}] MMseqs2 GPU linclust (coarse) on reps..."
CLUST_LOOSE="${TMP_DIR}/clusters_loose"
REP_TIGHT_DB="${TMP_DIR}/rep_tight.db"
mmseqs createdb "${REP_FASTA_TIGHT}" "${REP_TIGHT_DB}"
mmseqs linclust "${REP_TIGHT_DB}" "${CLUST_LOOSE}" "${TMP_DIR}/linclust_loose_tmp" \
  --min-seq-id 0.5 --cov-mode 1 -c 0.33 --cluster-mode 2 --gpu 1
mmseqs createtsv "${REP_TIGHT_DB}" "${REP_TIGHT_DB}" "${CLUST_LOOSE}" "${OUT_DIR}/clusters_loose.tsv"
REP_DB_LOOSE="${TMP_DIR}/repseqs_loose"
REP_FASTA_LOOSE="${OUT_DIR}/cluster_reps_loose.faa"
mmseqs createrepseqdb "${REP_TIGHT_DB}" "${CLUST_LOOSE}" "${REP_DB_LOOSE}"
mmseqs convert2fasta "${REP_DB_LOOSE}" "${REP_FASTA_LOOSE}"
STEP=$((STEP+1))

# Optional: ESM3 embeddings on coarse cluster reps (novelty/scoring)
if [[ -n "${ESM3_CMD}" ]]; then
  ESM3_OUT="${OUT_DIR}/esm3_embeddings"
  mkdir -p "${ESM3_OUT}"
  echo "[${STEP}/${TOTAL_STEPS}] Running ESM3 embeddings on cluster reps..."
  "${ESM3_CMD}" --fasta "${REP_FASTA_LOOSE}" --outdir "${ESM3_OUT}" --device cuda || echo "[WARN] ESM3 embedding failed; continuing." >&2
  STEP=$((STEP+1))
else
  echo "[${STEP}/${TOTAL_STEPS}] Skipping ESM3 embeddings (not requested/found)." >&2
fi

# Phylogenetic tree on coarse cluster reps (MAFFT + FastTree required)
ALIGN_FA="${OUT_DIR}/cluster_reps_loose.aln.faa"
TREE_FILE="${OUT_DIR}/cluster_reps_loose.tree.nwk"
echo "[${STEP}/${TOTAL_STEPS}] Building alignment with MAFFT and tree with FastTree..."
"${MAFFT_CMD}" --auto "${REP_FASTA_LOOSE}" > "${ALIGN_FA}"
"${FASTTREE_CMD}" -quiet -out "${TREE_FILE}" "${ALIGN_FA}"
STEP=$((STEP+1))

if [[ -n "${KNOWN_DB}" ]]; then
  echo "[${STEP}/${TOTAL_STEPS}] GPU search against known CRISPR/Cas DB: ${KNOWN_DB}"
  KNOWN_RES="${TMP_DIR}/known_hits"
  mmseqs search "${PROT_DB}" "${KNOWN_DB}" "${KNOWN_RES}" "${TMP_DIR}/search_tmp" --gpu 1
  mmseqs createtsv "${PROT_DB}" "${KNOWN_DB}" "${KNOWN_RES}" "${OUT_DIR}/known_hits.tsv"
else
  echo "[${STEP}/${TOTAL_STEPS}] Skipping known DB search (no --known-db provided)."
fi
STEP=$((STEP+1))

# Optional: AlphaFold3 modeling on cluster reps, then FoldSeek search
if [[ -n "${ALPHAFOLD3_CMD}" ]]; then
  if [[ ! -x "${ALPHAFOLD3_CMD}" ]]; then
    echo "[WARN] AlphaFold3 command not executable: ${ALPHAFOLD3_CMD}. Skipping modeling." >&2
  else
    MODELS_DIR="${OUT_DIR}/af3_models"
    AF3_INPUT_DIR="${TMP_DIR}/af3_inputs"
    mkdir -p "${MODELS_DIR}" "${AF3_INPUT_DIR}"
    echo "[${STEP}/${TOTAL_STEPS}] Preparing AlphaFold3 inputs from cluster reps..."
    REP_FASTA="${REP_FASTA_LOOSE}" AF3_INPUT_DIR="${AF3_INPUT_DIR}" ${PYTHON} - <<'PY'
import os, sys
from pathlib import Path

rep_fasta = os.environ.get("REP_FASTA")
out_dir = os.environ.get("AF3_INPUT_DIR")
if not rep_fasta or not out_dir:
    sys.exit(0)

def read_fasta(path):
    name, seq = None, []
    with open(path, "r", encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            yield name, "".join(seq)

Path(out_dir).mkdir(parents=True, exist_ok=True)
for name, seq in read_fasta(rep_fasta):
    with open(os.path.join(out_dir, f"{name}.fa"), "w", encoding="ascii") as h:
        h.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            h.write(seq[i:i+80] + "\n")
PY
    echo "[${STEP}/${TOTAL_STEPS}] Running AlphaFold3 on cluster representatives..."
    shopt -s nullglob
    for f in "${AF3_INPUT_DIR}"/*.fa; do
      base="$(basename "$f" .fa)"
      "${ALPHAFOLD3_CMD}" --fasta_paths "$f" --output_dir "${MODELS_DIR}/${base}" || {
        echo "[WARN] AlphaFold3 failed on $f; continuing." >&2
      }
    done
    shopt -u nullglob

    if command -v foldseek >/dev/null 2>&1; then
      echo "[${STEP}/${TOTAL_STEPS}] Running FoldSeek on AlphaFold3 models..."
      AF3_DB="${TMP_DIR}/af3_db"
      foldseek createdb "${MODELS_DIR}" "${AF3_DB}"
      if [[ -n "${FOLDSEEK_TARGET}" ]]; then
        FS_RES="${TMP_DIR}/foldseek_hits"
        foldseek search "${AF3_DB}" "${FOLDSEEK_TARGET}" "${FS_RES}" "${TMP_DIR}/fs_tmp" --format-mode 4
        foldseek convertalis "${AF3_DB}" "${FOLDSEEK_TARGET}" "${FS_RES}" "${OUT_DIR}/foldseek_hits.tsv"
      else
        FS_RES="${TMP_DIR}/foldseek_self"
        foldseek search "${AF3_DB}" "${AF3_DB}" "${FS_RES}" "${TMP_DIR}/fs_tmp" --format-mode 4
        foldseek convertalis "${AF3_DB}" "${AF3_DB}" "${FS_RES}" "${OUT_DIR}/foldseek_self.tsv"
      fi
    else
      echo "[WARN] foldseek not available; skipping FoldSeek search." >&2
    fi
  fi
fi
STEP=$((STEP+1))

echo "Done. Outputs:"
echo "  - ${OPERON_JSON}"
echo "  - ${OPERON_TSV}"
echo "  - ${PROTEINS}"
echo "  - ${OUT_DIR}/clusters_tight.tsv"
echo "  - ${OUT_DIR}/cluster_reps_tight.faa"
echo "  - ${OUT_DIR}/clusters_loose.tsv"
echo "  - ${OUT_DIR}/cluster_reps_loose.faa"
[[ -n "${KNOWN_DB}" ]] && echo "  - ${OUT_DIR}/known_hits.tsv"
[[ -d "${OUT_DIR}/esm3_embeddings" ]] && echo "  - ${OUT_DIR}/esm3_embeddings/"
[[ -f "${OUT_DIR}/cluster_reps_loose.tree.nwk" ]] && echo "  - ${OUT_DIR}/cluster_reps_loose.tree.nwk (phylogenetic tree)"
[[ -f "${OUT_DIR}/cluster_reps_loose.aln.faa" ]] && echo "  - ${OUT_DIR}/cluster_reps_loose.aln.faa (alignment for tree)"
