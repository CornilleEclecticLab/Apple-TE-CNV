#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple


def log(message: str) -> None:
    print(f"[INFO] {message}", file=sys.stderr)


def check_executable(name: str) -> None:
    result = subprocess.run(
        ["which", name],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Required executable not found in PATH: {name}")


def check_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"Input path is not a file: {path}")


def read_fasta(path: Path) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    current_id: str | None = None
    current_seq: List[str] = []

    with path.open("r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq).upper()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(current_seq).upper()

    return sequences


def write_fasta(sequences: Dict[str, str], output_path: Path) -> None:
    with output_path.open("w") as handle:
        for seq_id, seq in sequences.items():
            handle.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def sanitize_id(seq_id: str) -> str:
    return (
        seq_id.replace(" ", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace(":", "_")
        .replace(";", "_")
        .replace(",", "_")
        .replace("|", "_")
    )


def make_unique_ids(
    sequences: Dict[str, str],
    source_tag: str,
    seen_ids: set[str],
) -> Dict[str, str]:
    renamed: Dict[str, str] = {}

    for seq_id, seq in sequences.items():
        clean_id = sanitize_id(seq_id)
        candidate = clean_id
        counter = 2

        while candidate in seen_ids:
            candidate = f"{clean_id}__{source_tag}__{counter}"
            counter += 1

        seen_ids.add(candidate)
        renamed[candidate] = seq

    return renamed


def run_command(command: List[str]) -> None:
    log("Running: " + " ".join(command))
    subprocess.run(command, check=True)


def build_blast_database(fasta_path: Path, db_prefix: Path) -> None:
    run_command([
        "makeblastdb",
        "-in", str(fasta_path),
        "-dbtype", "nucl",
        "-out", str(db_prefix),
    ])


def run_blast_all_vs_all(
    fasta_path: Path,
    db_prefix: Path,
    output_path: Path,
    threads: int,
) -> None:
    outfmt = "6 qseqid sseqid pident length qlen slen bitscore evalue"

    run_command([
        "blastn",
        "-query", str(fasta_path),
        "-db", str(db_prefix),
        "-out", str(output_path),
        "-outfmt", outfmt,
        "-dust", "no",
        "-soft_masking", "false",
        "-evalue", "1e-20",
        "-num_threads", str(threads),
    ])


def parse_blast_and_deduplicate(
    blast_output: Path,
    sequences: Dict[str, str],
    identity_threshold: float,
    coverage_threshold: float,
) -> Tuple[Dict[str, str], List[Tuple[str, str, float, float, float, int, int]]]:
    removed: set[str] = set()
    summary_rows: List[Tuple[str, str, float, float, float, int, int]] = []

    with blast_output.open("r") as handle:
        for line in handle:
            fields = line.strip().split("\t")
            if len(fields) != 8:
                continue

            qseqid, sseqid, pident, aln_len, qlen, slen, bitscore, evalue = fields

            if qseqid == sseqid:
                continue
            if qseqid in removed or sseqid in removed:
                continue

            pident_f = float(pident)
            aln_len_i = int(aln_len)
            qlen_i = int(qlen)
            slen_i = int(slen)

            if qlen_i == 0 or slen_i == 0:
                continue

            qcov = 100.0 * aln_len_i / qlen_i
            scov = 100.0 * aln_len_i / slen_i

            if pident_f < identity_threshold:
                continue
            if qcov < coverage_threshold:
                continue
            if scov < coverage_threshold:
                continue

            if qlen_i > slen_i:
                keep_id, drop_id = qseqid, sseqid
            elif slen_i > qlen_i:
                keep_id, drop_id = sseqid, qseqid
            else:
                keep_id, drop_id = sorted([qseqid, sseqid])[0], sorted([qseqid, sseqid])[1]

            if drop_id not in removed:
                removed.add(drop_id)
                summary_rows.append((
                    drop_id,
                    keep_id,
                    pident_f,
                    qcov,
                    scov,
                    len(sequences[drop_id]),
                    len(sequences[keep_id]),
                ))

    final_sequences = {
        seq_id: seq
        for seq_id, seq in sequences.items()
        if seq_id not in removed
    }

    return final_sequences, summary_rows


def load_and_filter_fasta(
    fasta_path: Path,
    min_length: int,
    seen_ids: set[str],
) -> Dict[str, str]:
    source_tag = fasta_path.stem
    raw_sequences = read_fasta(fasta_path)

    filtered_sequences = {
        seq_id: seq
        for seq_id, seq in raw_sequences.items()
        if len(seq) >= min_length
    }

    return make_unique_ids(filtered_sequences, source_tag, seen_ids)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Merge TE FASTA files, remove sequences shorter than a minimum length, "
            "deduplicate using an 80/80/80 BLAST rule, and keep the longest sequence."
        )
    )

    parser.add_argument(
        "--repdenovo",
        required=True,
        help="REPdenovo FASTA file",
    )
    parser.add_argument(
        "--repeatexplorer",
        required=True,
        help="RepeatExplorer FASTA file",
    )
    parser.add_argument(
        "--gddh",
        required=True,
        help="TE_database_GDDH FASTA file",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output non-redundant FASTA",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=150,
        help="Minimum length to keep before deduplication (default: 150)",
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=80.0,
        help="Minimum percent identity (default: 80)",
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=80.0,
        help="Minimum coverage percentage for both query and subject (default: 80)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of BLAST threads (default: 1)",
    )

    args = parser.parse_args()

    repdenovo_path = Path(args.repdenovo)
    repeatexplorer_path = Path(args.repeatexplorer)
    gddh_path = Path(args.gddh)
    output_path = Path(args.output)

    check_executable("makeblastdb")
    check_executable("blastn")

    check_file(repdenovo_path)
    check_file(repeatexplorer_path)
    check_file(gddh_path)

    seen_ids: set[str] = set()

    log("Loading and filtering FASTA files...")
    repdenovo_sequences = load_and_filter_fasta(repdenovo_path, args.min_length, seen_ids)
    repeatexplorer_sequences = load_and_filter_fasta(repeatexplorer_path, args.min_length, seen_ids)
    gddh_sequences = load_and_filter_fasta(gddh_path, args.min_length, seen_ids)

    all_sequences: Dict[str, str] = {}
    all_sequences.update(repdenovo_sequences)
    all_sequences.update(repeatexplorer_sequences)
    all_sequences.update(gddh_sequences)

    if not all_sequences:
        raise RuntimeError("No sequences remained after filtering by minimum length.")

    log(f"REPdenovo sequences retained: {len(repdenovo_sequences)}")
    log(f"RepeatExplorer sequences retained: {len(repeatexplorer_sequences)}")
    log(f"GDDH database sequences retained: {len(gddh_sequences)}")
    log(f"Total sequences retained before deduplication: {len(all_sequences)}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="dedup_te_") as tmpdir_str:
        tmpdir = Path(tmpdir_str)

        combined_fasta = tmpdir / "combined_sequences.fasta"
        blast_db_prefix = tmpdir / "combined_db"
        blast_output = tmpdir / "all_vs_all.tsv"

        write_fasta(all_sequences, combined_fasta)
        build_blast_database(combined_fasta, blast_db_prefix)
        run_blast_all_vs_all(
            fasta_path=combined_fasta,
            db_prefix=blast_db_prefix,
            output_path=blast_output,
            threads=args.threads,
        )

        final_sequences, summary_rows = parse_blast_and_deduplicate(
            blast_output=blast_output,
            sequences=all_sequences,
            identity_threshold=args.identity,
            coverage_threshold=args.coverage,
        )

    write_fasta(final_sequences, output_path)
    summary_path = output_path.with_suffix(output_path.suffix + ".summary.tsv")
    with summary_path.open("w") as handle:
        handle.write(
            "removed_id\tkept_id\tpident\tqcov\tscov\tremoved_length\tkept_length\n"
        )
        for row in summary_rows:
            removed_id, kept_id, pident, qcov, scov, removed_len, kept_len = row
            handle.write(
                f"{removed_id}\t{kept_id}\t{pident:.2f}\t{qcov:.2f}\t{scov:.2f}\t{removed_len}\t{kept_len}\n"
            )

    log(f"Final non-redundant FASTA written to: {output_path}")
    log(f"Deduplication summary written to: {summary_path}")
    log(f"Sequences kept after deduplication: {len(final_sequences)}")
    log(f"Sequences removed as redundant: {len(summary_rows)}")


if __name__ == "__main__":
    main()