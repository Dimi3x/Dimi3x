#!/usr/bin/env bash
# Minimal runner for the core pipeline steps listed in PIPELINE_SUMMARY.md.
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <mag_fasta_dir_or_file> <output_dir> [--known-db /path/to/mmseqs_db]" >&2
  exit 1
fi

INPUT_PATH="$1"
OUT_DIR="$2"
shift 2 || true

KNOWN_DB=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --known-db)
      KNOWN_DB="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

# Required tools
for cmd in prodigal cronus CRISPRCasFinder.pl mmseqs mafft fasttree; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required tool missing: $cmd" >&2
    exit 1
  fi
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON=${PYTHON:-python3}

mkdir -p "${OUT_DIR}"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "${TMP_DIR}"

PRODIGAL_DIR="${TMP_DIR}/prodigal"
CRONUS_DIR="${TMP_DIR}/cronus"
CASF_DIR="${TMP_DIR}/casfinder"
OPERON_JSON="${OUT_DIR}/operons.json"
OPERON_TSV="${OUT_DIR}/operons.tsv"
PROTEINS="${OUT_DIR}/operon_proteins.faa"

# 1) Prodigal ORFs
mkdir -p "${PRODIGAL_DIR}"
echo "[1/?] Running Prodigal..."
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

# 2) CRONUS arrays
mkdir -p "${CRONUS_DIR}"
CRONUS_COMBINED="${CRONUS_DIR}/cronus_combined.tsv"
rm -f "${CRONUS_COMBINED}"
echo "[2/?] Running CRONUS..."
header_written=0
for f in "${FASTA_FILES[@]}"; do
  base="$(basename "$f")"
  stem="${base%.*}"
  out_tsv="${CRONUS_DIR}/${stem}.cronus.tsv"
  cronus -i "$f" -o "${out_tsv}"
  if [[ -f "${out_tsv}" ]]; then
    if [[ ${header_written} -eq 0 ]]; then
      head -n 1 "${out_tsv}" > "${CRONUS_COMBINED}"
      header_written=1
    fi
    tail -n +2 "${out_tsv}" >> "${CRONUS_COMBINED}"
  fi
done
if [[ ! -f "${CRONUS_COMBINED}" ]]; then
  echo "CRONUS produced no arrays." >&2
  exit 1
fi

# 3) CRISPRCasFinder arrays
mkdir -p "${CASF_DIR}"
CASF_COMBINED="${CASF_DIR}/CRISPRs_REPORT.combined.tsv"
rm -f "${CASF_COMBINED}"
echo "[3/?] Running CRISPRCasFinder..."
header_written=0
for f in "${FASTA_FILES[@]}"; do
  base="$(basename "$f")"
  outdir="${CASF_DIR}/${base%.*}"
  mkdir -p "${outdir}"
  CRISPRCasFinder.pl -in "$f" -out "${outdir}"
  report="${outdir}/CRISPRs_REPORT.tsv"
  if [[ -f "${report}" ]]; then
    if [[ ${header_written} -eq 0 ]]; then
      head -n 1 "${report}" > "${CASF_COMBINED}"
      header_written=1
    fi
    tail -n +2 "${report}" >> "${CASF_COMBINED}"
  fi
done
if [[ ! -f "${CASF_COMBINED}" ]]; then
  echo "CRISPRCasFinder produced no arrays." >&2
  exit 1
fi

# 4) Miner: arrays -> operons
echo "[4/?] Running miner..."
${PYTHON} "${SCRIPT_DIR}/crispr_operon_miner.py" \
  --input "${INPUT_PATH}" \
  --output-json "${OPERON_JSON}" \
  --output-tsv "${OPERON_TSV}" \
  --output-proteins "${PROTEINS}" \
  --prodigal-dir "${PRODIGAL_DIR}" \
  --cronus-tsv "${CRONUS_COMBINED}" \
  --casfinder-tsv "${CASF_COMBINED}"

# 5) MMseqs2 clustering (tight, then coarse)
PROT_DB="${TMP_DIR}/proteins.db"
CLUST_TIGHT="${TMP_DIR}/clusters_tight"
mmseqs createdb "${PROTEINS}" "${PROT_DB}"
mmseqs linclust "${PROT_DB}" "${CLUST_TIGHT}" "${TMP_DIR}/linclust_tight_tmp" \
  --min-seq-id 0.9 --cov-mode 1 -c 0.8 --cluster-mode 2 --gpu 1
mmseqs createtsv "${PROT_DB}" "${PROT_DB}" "${CLUST_TIGHT}" "${OUT_DIR}/clusters_tight.tsv"
REP_DB_TIGHT="${TMP_DIR}/repseqs_tight"
REP_FASTA_TIGHT="${OUT_DIR}/cluster_reps_tight.faa"
mmseqs createrepseqdb "${PROT_DB}" "${CLUST_TIGHT}" "${REP_DB_TIGHT}"
mmseqs convert2fasta "${REP_DB_TIGHT}" "${REP_FASTA_TIGHT}"

echo "[5/?] Coarse clustering on tight reps..."
CLUST_LOOSE="${TMP_DIR}/clusters_loose"
REP_DB_LOOSE="${TMP_DIR}/repseqs_loose"
REP_FASTA_LOOSE="${OUT_DIR}/cluster_reps_loose.faa"
mmseqs createdb "${REP_FASTA_TIGHT}" "${TMP_DIR}/rep_tight.db"
mmseqs linclust "${TMP_DIR}/rep_tight.db" "${CLUST_LOOSE}" "${TMP_DIR}/linclust_loose_tmp" \
  --min-seq-id 0.5 --cov-mode 1 -c 0.33 --cluster-mode 2 --gpu 1
mmseqs createtsv "${TMP_DIR}/rep_tight.db" "${TMP_DIR}/rep_tight.db" "${CLUST_LOOSE}" "${OUT_DIR}/clusters_loose.tsv"
mmseqs createrepseqdb "${TMP_DIR}/rep_tight.db" "${CLUST_LOOSE}" "${REP_DB_LOOSE}"
mmseqs convert2fasta "${REP_DB_LOOSE}" "${REP_FASTA_LOOSE}"

# 6) Phylogeny (MAFFT + FastTree on coarse reps)
ALIGN_FA="${OUT_DIR}/cluster_reps_loose.aln.faa"
TREE_FILE="${OUT_DIR}/cluster_reps_loose.tree.nwk"
echo "[6/?] Building alignment and tree..."
mafft --auto "${REP_FASTA_LOOSE}" > "${ALIGN_FA}"
fasttree -quiet -out "${TREE_FILE}" "${ALIGN_FA}"

# 7) Optional known DB search
if [[ -n "${KNOWN_DB}" ]]; then
  echo "[7/?] MMseqs2 GPU search vs known DB..."
  KNOWN_RES="${TMP_DIR}/known_hits"
  mmseqs search "${PROT_DB}" "${KNOWN_DB}" "${KNOWN_RES}" "${TMP_DIR}/search_tmp" --gpu 1
  mmseqs createtsv "${PROT_DB}" "${KNOWN_DB}" "${KNOWN_RES}" "${OUT_DIR}/known_hits.tsv"
fi

echo "Done. Outputs in ${OUT_DIR}"
