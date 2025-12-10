#!/usr/bin/env python3
"""
Standalone CRISPR operon miner.

- Detects CRISPR arrays in assembled MAGs using simple repeat/spacer heuristics.
- Calls ORFs on both strands.
- Assembles operons (array + flanking ORFs), scores, and exports JSON/TSV/FASTA.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from collections import defaultdict, namedtuple
from glob import glob
from typing import Dict, Iterable, List, Tuple


ORF = namedtuple("ORF", ["id", "contig", "start", "end", "strand", "aa_seq"])


def read_fasta(path: str) -> Iterable[Tuple[str, str]]:
    header = None
    seq_parts: List[str] = []
    with open(path, "r", encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts).upper()
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts).upper()


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDN")
    return seq.translate(comp)[::-1]


def translate_dna(seq: str) -> str:
    table = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    aa = []
    for i in range(0, len(seq) - 2, 3):
        aa.append(table.get(seq[i:i+3], "X"))
    return "".join(aa)


def find_orfs(seq: str, contig: str, min_aa: int, max_aa: int) -> List[ORF]:
    start_codons = {"ATG", "GTG", "TTG"}
    stop_codons = {"TAA", "TAG", "TGA"}
    results: List[ORF] = []
    seq_len = len(seq)
    strands = [("+", seq), ("-", reverse_complement(seq))]

    for strand, nuc in strands:
        for frame in range(3):
            i = frame
            start_idx = None
            while i <= len(nuc) - 3:
                codon = nuc[i:i+3]
                if start_idx is None:
                    if codon in start_codons:
                        start_idx = i
                else:
                    if codon in stop_codons:
                        aa_len = (i + 3 - start_idx) // 3 - 1
                        if min_aa <= aa_len <= max_aa:
                            aa_seq = translate_dna(nuc[start_idx:i+3])[:-1]
                            if strand == "+":
                                start = start_idx
                                end = i + 3
                            else:
                                start = seq_len - (i + 3)
                                end = seq_len - start_idx
                            orf_id = f"{contig}|{strand}|{start+1}-{end}"
                            results.append(ORF(orf_id, contig, start, end, strand, aa_seq))
                        start_idx = None
                i += 3
    return results


def parse_prodigal_header(header: str) -> Tuple[str, int, int, str]:
    """
    Prodigal FASTA header format example:
    >contig_1 # 1 # 126 # 1 # ID=1_1;partial=00;...
    """
    parts = header.split("#")
    if len(parts) < 4:
        raise ValueError(f"Unexpected Prodigal header: {header}")
    contig = parts[0].split()[0].strip()
    start = int(parts[1].strip()) - 1  # convert to 0-based
    end = int(parts[2].strip())
    strand_raw = parts[3].strip()
    strand = "+" if strand_raw == "1" else "-"
    return contig, start, end, strand


def read_prodigal_orfs(faa_files: List[str], min_aa: int, max_aa: int) -> List[ORF]:
    orfs: List[ORF] = []
    for faa in faa_files:
        header = None
        seq_parts: List[str] = []
        with open(faa, "r", encoding="ascii") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if header:
                        contig, start, end, strand = parse_prodigal_header(header)
                        aa_seq = "".join(seq_parts)
                        if min_aa <= len(aa_seq) <= max_aa:
                            orf_id = f"{contig}|{strand}|{start+1}-{end}"
                            orfs.append(ORF(orf_id, contig, start, end, strand, aa_seq))
                    header = line[1:]
                    seq_parts = []
                else:
                    seq_parts.append(line)
            if header:
                contig, start, end, strand = parse_prodigal_header(header)
                aa_seq = "".join(seq_parts)
                if min_aa <= len(aa_seq) <= max_aa:
                    orf_id = f"{contig}|{strand}|{start+1}-{end}"
                    orfs.append(ORF(orf_id, contig, start, end, strand, aa_seq))
    return orfs


def read_cronus_tsv(path: str) -> List[Dict]:
    """
    Parse CRONUS TSV. Expected columns (flexible):
    contig/Sequence, start/Start, end/End, repeat_seq/DR, repeat_len/DR_Length, spacer_count/Spacer_nb.
    """
    arrays: List[Dict] = []
    with open(path, "r", encoding="ascii") as handle:
        header = None
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = [h.strip() for h in parts]
                continue
            row = dict(zip(header, parts))
            contig = row.get("contig") or row.get("Sequence") or row.get("seq")
            try:
                start = int(row.get("start", row.get("Start", "0"))) - 1
                end = int(row.get("end", row.get("End", "0")))
            except Exception:
                continue
            repeat_seq = row.get("repeat_seq") or row.get("DR") or row.get("DR_Consensus", "")
            try:
                repeat_len = int(row.get("repeat_len", row.get("DR_Length", len(repeat_seq))))
            except Exception:
                repeat_len = len(repeat_seq)
            try:
                spacer_count = int(row.get("spacer_count", row.get("Spacer_nb", "0")))
            except Exception:
                spacer_count = 0
            if not contig or end <= start:
                continue
            arrays.append(
                {
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "repeat_len": repeat_len,
                    "repeat_seq": repeat_seq,
                    "spacer_mean": 0.0,
                    "spacer_std": 0.0,
                    "spacer_count": spacer_count,
                    "repeat_count": spacer_count + 1 if spacer_count else 0,
                    "score": 0.0,
                }
            )
    return arrays


def _spacer_stats(spacers: List[int]) -> Tuple[float, float]:
    if not spacers:
        return 0.0, 0.0
    mean = sum(spacers) / len(spacers)
    var = sum((x - mean) ** 2 for x in spacers) / len(spacers)
    return mean, math.sqrt(var)


def detect_arrays(
    seq: str,
    contig: str,
    min_repeats: int,
    repeat_range: Tuple[int, int],
    spacer_range: Tuple[int, int],
) -> List[Dict]:
    arrays: List[Dict] = []
    repeat_min, repeat_max = repeat_range
    spacer_min, spacer_max = spacer_range
    candidate_lengths = list(range(repeat_min, repeat_max + 1, 4))

    for L in candidate_lengths:
        positions: Dict[str, List[int]] = defaultdict(list)
        for i in range(0, len(seq) - L + 1):
            repeat = seq[i:i+L]
            positions[repeat].append(i)

        for repeat_seq, pos_list in positions.items():
            if len(pos_list) < min_repeats:
                continue
            pos_list.sort()
            for idx in range(len(pos_list) - min_repeats + 1):
                anchors = [pos_list[idx]]
                spacers: List[int] = []
                last = pos_list[idx]
                for nxt in pos_list[idx + 1:]:
                    spacer = nxt - last - L
                    if spacer_min <= spacer <= spacer_max:
                        anchors.append(nxt)
                        spacers.append(spacer)
                        last = nxt
                    elif spacer > spacer_max:
                        break
                if len(anchors) >= min_repeats:
                    start = anchors[0]
                    end = anchors[-1] + L
                    mean_spacer, std_spacer = _spacer_stats(spacers)
                    base_score = len(anchors) * 1.2
                    length_bonus = 0.5 if 28 <= L <= 38 else 0.2
                    spacer_stability = max(0.1, 1.0 - (std_spacer / (mean_spacer + 1e-6)))
                    score = base_score + length_bonus + spacer_stability
                    arrays.append(
                        {
                            "contig": contig,
                            "start": start,
                            "end": end,
                            "repeat_len": L,
                            "repeat_seq": repeat_seq,
                            "spacer_mean": mean_spacer,
                            "spacer_std": std_spacer,
                            "spacer_count": len(spacers),
                            "repeat_count": len(anchors),
                            "score": round(score, 3),
                        }
                    )
    arrays.sort(key=lambda x: (x["contig"], x["start"]))
    filtered: List[Dict] = []
    for arr in arrays:
        if filtered and arr["contig"] == filtered[-1]["contig"] and abs(arr["start"] - filtered[-1]["start"]) < 20:
            if arr["score"] > filtered[-1]["score"]:
                filtered[-1] = arr
            continue
        filtered.append(arr)
    return filtered


def assemble_operons(
    arrays: List[Dict],
    orfs: List[ORF],
    flank_bp: int,
) -> List[Dict]:
    operons: List[Dict] = []
    orfs_by_contig: Dict[str, List[ORF]] = defaultdict(list)
    for orf in orfs:
        orfs_by_contig[orf.contig].append(orf)
    for contig in orfs_by_contig:
        orfs_by_contig[contig].sort(key=lambda o: o.start)

    for idx, arr in enumerate(arrays):
        contig_orfs = orfs_by_contig.get(arr["contig"], [])
        nearby: List[ORF] = []
        for orf in contig_orfs:
            if (arr["start"] - flank_bp) <= orf.start <= (arr["end"] + flank_bp):
                nearby.append(orf)
        prox_bonus = sum(math.exp(-abs((arr["start"] + arr["end"]) / 2 - (orf.start + orf.end) / 2) / 3000.0) for orf in nearby)
        density = len(nearby) / max(1.0, flank_bp / 1000.0)
        length_penalty = sum(0.5 for orf in nearby if len(orf.aa_seq) < 100 or len(orf.aa_seq) > 1800)
        score = arr["score"] * 1.5 + density + prox_bonus - length_penalty
        operons.append(
            {
                "id": f"operon_{idx+1}",
                "contig": arr["contig"],
                "array_start": arr["start"],
                "array_end": arr["end"],
                "array_score": arr["score"],
                "repeat_len": arr["repeat_len"],
                "repeat_count": arr["repeat_count"],
                "spacer_count": arr["spacer_count"],
                "orf_ids": [orf.id for orf in nearby],
                "orf_count": len(nearby),
                "score": round(score, 3),
            }
        )
    operons.sort(key=lambda x: x["score"], reverse=True)
    return operons


def write_proteins(orfs: List[ORF], path: str, max_per_contig: int | None = None) -> None:
    with open(path, "w", encoding="ascii") as handle:
        if max_per_contig is None:
            selected = orfs
        else:
            selected = []
            seen: Dict[str, int] = defaultdict(int)
            for orf in orfs:
                if seen[orf.contig] < max_per_contig:
                    selected.append(orf)
                    seen[orf.contig] += 1
        for orf in selected:
            handle.write(f">{orf.id}\n")
            seq = orf.aa_seq
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")


def write_operon_tsv(operons: List[Dict], path: str) -> None:
    with open(path, "w", encoding="ascii") as handle:
        header = [
            "operon_id",
            "contig",
            "array_start",
            "array_end",
            "repeat_len",
            "repeat_count",
            "spacer_count",
            "array_score",
            "orf_count",
            "score",
            "orf_ids",
        ]
        handle.write("\t".join(header) + "\n")
        for op in operons:
            row = [
                op["id"],
                op["contig"],
                str(op["array_start"] + 1),
                str(op["array_end"]),
                str(op["repeat_len"]),
                str(op["repeat_count"]),
                str(op["spacer_count"]),
                f"{op['array_score']:.3f}",
                str(op["orf_count"]),
                f"{op['score']:.3f}",
                ",".join(op["orf_ids"]),
            ]
            handle.write("\t".join(row) + "\n")


def discover_inputs(path: str) -> List[str]:
    if os.path.isfile(path):
        return [path]
    patterns = ["*.fa", "*.fasta", "*.fna", "*.fas"]
    files: List[str] = []
    for pat in patterns:
        files.extend(glob(os.path.join(path, pat)))
    return sorted(files)


def main() -> None:
    parser = argparse.ArgumentParser(description="Detect CRISPR arrays, assemble operons, emit proteins.")
    parser.add_argument("--input", required=True, help="FASTA file or directory with MAG FASTA files.")
    parser.add_argument("--output-json", required=True, help="Path to write operons.json.")
    parser.add_argument("--output-tsv", required=True, help="Path to write operons.tsv.")
    parser.add_argument("--output-proteins", required=True, help="Path to write operon_proteins.faa.")
    parser.add_argument("--prodigal-dir", help="Directory containing Prodigal .faa outputs (one per MAG).")
    parser.add_argument("--casfinder-tsv", help="Optional CRISPRCasFinder combined CRISPRs_REPORT.tsv; overrides built-in array detector.")
    parser.add_argument("--cronus-tsv", help="Optional CRONUS CRISPR array report TSV; highest priority for arrays.")
    parser.add_argument("--flank", type=int, default=10000, help="Flank size (bp) around arrays to collect ORFs.")
    parser.add_argument("--min-repeat", type=int, default=23, help="Minimum repeat length.")
    parser.add_argument("--max-repeat", type=int, default=47, help="Maximum repeat length.")
    parser.add_argument("--min-spacer", type=int, default=20, help="Minimum spacer length.")
    parser.add_argument("--max-spacer", type=int, default=60, help="Maximum spacer length.")
    parser.add_argument("--min-aa", type=int, default=60, help="Minimum ORF length (aa).")
    parser.add_argument("--max-aa", type=int, default=1800, help="Maximum ORF length (aa).")
    parser.add_argument("--max-proteins-per-contig", type=int, default=None, help="Optional cap on proteins exported per contig.")
    args = parser.parse_args()

    fasta_files = discover_inputs(args.input)
    if not fasta_files:
        raise SystemExit(f"No FASTA files found under {args.input}")

    arrays: List[Dict] = []
    all_orfs: List[ORF] = []

    prodigal_files: List[str] = []
    if args.prodigal_dir:
        for ext in ("*.faa", "*.fasta", "*.fa"):
            prodigal_files.extend(glob(os.path.join(args.prodigal_dir, ext)))
        prodigal_files = sorted(set(prodigal_files))
        if prodigal_files:
            all_orfs = read_prodigal_orfs(prodigal_files, args.min_aa, args.max_aa)

    cronus_arrays: List[Dict] = []
    casfinder_arrays: List[Dict] = []
    if args.cronus_tsv and os.path.isfile(args.cronus_tsv):
        cronus_arrays = read_cronus_tsv(args.cronus_tsv)
    if args.casfinder_tsv and os.path.isfile(args.casfinder_tsv):
        casfinder_arrays = read_casfinder_tsv(args.casfinder_tsv)

    for fasta in fasta_files:
        for contig, seq in read_fasta(fasta):
            contig_arrays: List[Dict] = []
            if cronus_arrays:
                contig_arrays.extend([a for a in cronus_arrays if a["contig"] == contig])
            if casfinder_arrays:
                contig_arrays.extend([a for a in casfinder_arrays if a["contig"] == contig])
            if not contig_arrays:
                contig_arrays = detect_arrays(
                    seq,
                    contig,
                    min_repeats=3,
                    repeat_range=(args.min_repeat, args.max_repeat),
                    spacer_range=(args.min_spacer, args.max_spacer),
                )
            arrays.extend(contig_arrays)
            if not args.prodigal_dir:
                orfs = find_orfs(seq, contig, args.min_aa, args.max_aa)
                all_orfs.extend(orfs)

    operons = assemble_operons(arrays, all_orfs, args.flank)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_json)), exist_ok=True)
    write_operon_tsv(operons, args.output_tsv)
    write_proteins(all_orfs, args.output_proteins, args.max_proteins_per_contig)

    payload = {
        "arrays": arrays,
        "orfs": [
            {
                "id": orf.id,
                "contig": orf.contig,
                "start": orf.start,
                "end": orf.end,
                "strand": orf.strand,
                "aa_len": len(orf.aa_seq),
            }
            for orf in all_orfs
        ],
        "operons": operons,
        "params": {
            "flank": args.flank,
            "repeat_range": [args.min_repeat, args.max_repeat],
            "spacer_range": [args.min_spacer, args.max_spacer],
            "aa_range": [args.min_aa, args.max_aa],
        },
    }
    with open(args.output_json, "w", encoding="ascii") as handle:
        json.dump(payload, handle, indent=2)

    print(f"Processed {len(fasta_files)} FASTA files.")
    print(f"Arrays: {len(arrays)} | ORFs: {len(all_orfs)} | Operons: {len(operons)}")
    print(f"Wrote: {args.output_json}, {args.output_tsv}, {args.output_proteins}")


if __name__ == "__main__":
    main()
def read_casfinder_tsv(path: str) -> List[Dict]:
    """
    Parse CRISPRCasFinder CRISPRs_REPORT.tsv.
    Expected columns (subset used): Sequence, Start, End, DR_Consensus, DR_Length, Spacer_nb.
    """
    arrays: List[Dict] = []
    with open(path, "r", encoding="ascii") as handle:
        header = None
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            row = dict(zip(header, parts))
            try:
                contig = row.get("Sequence") or row.get("seq")
                start = int(row.get("Start", "0")) - 1
                end = int(row.get("End", "0"))
                dr_len = int(row.get("DR_Length", "0"))
                dr_seq = row.get("DR_Consensus", "")
                spacer_nb = int(row.get("Spacer_nb", row.get("Spacer_count", "0")))
                if not contig or end <= start:
                    continue
            except Exception:
                continue
            arrays.append(
                {
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "repeat_len": dr_len,
                    "repeat_seq": dr_seq,
                    "spacer_mean": 0.0,
                    "spacer_std": 0.0,
                    "spacer_count": spacer_nb,
                    "repeat_count": spacer_nb + 1 if spacer_nb else 0,
                    "score": 0.0,  # CasFinder scores not parsed; miner will score via operon assembly step
                }
            )
    return arrays
