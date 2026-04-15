#!/usr/bin/env python3
import argparse
import csv
import os
import re
from pathlib import Path

REQUIRED_COLS = [
    "sample_id",
    "fastq_1",
    "fastq_2",
    "group",
    "mean_insert",
    "sd_insert",
    "genome_length",
    "read_length",
    "read_depth",
]

OPTIONAL_COLS = [
    "threads",
    "output_folder",
    "refiner_path",
    "contigs_merger_path",
    "memory_gb",
    "stage",  # optional: lets you label assembly/scaffolding choices; not required
    "extra_kv",  # optional: KEY=VALUE;KEY=VALUE overrides
]

DEFAULT_THREADS = 15
DEFAULT_MEMORY_GB = 90
DEFAULT_REFINER = "/usr/local/src/REPdenovo/TERefiner_1"
DEFAULT_MERGER = "/usr/local/src/REPdenovo/ContigsMerger"

SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

def parse_extra_kv(s: str) -> dict:
    """
    Parse "KEY=VALUE;KEY=VALUE" into dict.
    Empty/None -> {}.
    """
    if not s:
        return {}
    out = {}
    parts = [p.strip() for p in s.split(";") if p.strip()]
    for p in parts:
        if "=" not in p:
            raise ValueError(f"extra_kv entry must be KEY=VALUE (got: {p})")
        k, v = p.split("=", 1)
        out[k.strip()] = v.strip()
    return out

def load_template(path: Path) -> str:
    return path.read_text(encoding="utf-8")

def render_template(tpl: str, mapping: dict) -> str:
    out = tpl
    for k, v in mapping.items():
        out = out.replace(k, str(v))
    return out

def ensure_sample_id(sample_id: str):
    if not SAFE_ID.match(sample_id):
        raise ValueError(
            f"Invalid sample_id '{sample_id}'. Use letters/numbers/._- only, start with alnum."
        )

def main():
    ap = argparse.ArgumentParser(
        description="Generate REPdenovo per-sample configs + reads lists + HTCondor jobs.tsv"
    )
    ap.add_argument("--samples", required=True, help="TSV file (tab-separated) describing samples")
    ap.add_argument("--template", default="configs/template_repdenovo.conf", help="Path to template config")
    ap.add_argument("--out-configs", default="configs/generated", help="Output directory for configs")
    ap.add_argument("--out-reads", default="reads", help="Output directory for reads list files")
    ap.add_argument("--out-jobs", default="cluster/jobs.tsv", help="Output jobs.tsv for HTCondor queue")
    ap.add_argument("--default-output-root", default="outputs", help="Root output folder if output_folder not provided")
    args = ap.parse_args()

    samples_path = Path(args.samples)
    tpl_path = Path(args.template)

    out_configs = Path(args.out_configs)
    out_reads = Path(args.out_reads)
    out_jobs = Path(args.out_jobs)

    out_configs.mkdir(parents=True, exist_ok=True)
    out_reads.mkdir(parents=True, exist_ok=True)
    out_jobs.parent.mkdir(parents=True, exist_ok=True)

    tpl = load_template(tpl_path)

    with samples_path.open("r", encoding="utf-8", newline="") as f:
        # TSV
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in REQUIRED_COLS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(
                f"Missing required columns in {samples_path}: {missing}\n"
                f"Required: {REQUIRED_COLS}\nOptional: {OPTIONAL_COLS}"
            )

        jobs_rows = []
        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            ensure_sample_id(sample_id)

            fq1 = (row.get("fastq_1") or "").strip()
            fq2 = (row.get("fastq_2") or "").strip()

            group = (row.get("group") or "").strip()
            mean_insert = (row.get("mean_insert") or "").strip()
            sd_insert = (row.get("sd_insert") or "").strip()

            genome_length = (row.get("genome_length") or "").strip()
            read_length = (row.get("read_length") or "").strip()
            read_depth = (row.get("read_depth") or "").strip()

            threads = int((row.get("threads") or "").strip() or DEFAULT_THREADS)
            mem_gb = int((row.get("memory_gb") or "").strip() or DEFAULT_MEMORY_GB)

            out_folder = (row.get("output_folder") or "").strip()
            if not out_folder:
                out_folder = str(Path(args.default_output_root) / sample_id)

            refiner_path = (row.get("refiner_path") or "").strip() or DEFAULT_REFINER
            merger_path = (row.get("contigs_merger_path") or "").strip() or DEFAULT_MERGER

            extra = parse_extra_kv((row.get("extra_kv") or "").strip())

            # reads list file (like your RdN_Path_to_Ang_01)
            reads_txt = out_reads / f"{sample_id}.txt"
            if fq2:
                reads_txt.write_text(
                    f"{fq1} {group} {mean_insert} {sd_insert}\n{fq2} {group} {mean_insert} {sd_insert}\n",
                    encoding="utf-8",
                )
            else:
                # single-end
                reads_txt.write_text(f"{fq1} {group} -1 -1\n", encoding="utf-8")

            # config file
            conf_path = out_configs / f"{sample_id}.conf"
            mapping = {
                "__READ_LENGTH__": read_length,
                "__READ_DEPTH__": read_depth,
                "__GENOME_LENGTH__": genome_length,
                "__THREADS__": threads,
                "__OUTPUT_FOLDER__": out_folder,
                "__REFINER_PATH__": refiner_path,
                "__CONTIGS_MERGER_PATH__": merger_path,
            }
            conf_text = render_template(tpl, mapping)

            # Apply extra overrides: replace lines "KEY value" if present, else append.
            # This keeps it simple and lets you override any REPdenovo key.
            if extra:
                lines = conf_text.splitlines()
                keys_in_file = {ln.split()[0] for ln in lines if ln.strip() and not ln.strip().startswith("#")}
                # replace
                new_lines = []
                for ln in lines:
                    if not ln.strip() or ln.strip().startswith("#"):
                        new_lines.append(ln)
                        continue
                    k = ln.split()[0]
                    if k in extra:
                        new_lines.append(f"{k} {extra[k]}")
                    else:
                        new_lines.append(ln)
                # append keys not present
                for k, v in extra.items():
                    if k not in keys_in_file:
                        new_lines.append(f"{k} {v}")
                conf_text = "\n".join(new_lines) + "\n"
            else:
                conf_text = conf_text.rstrip("\n") + "\n"

            conf_path.write_text(conf_text, encoding="utf-8")

            jobs_rows.append(
                {
                    "sample_id": sample_id,
                    "config": str(conf_path),
                    "reads": str(reads_txt),
                    "cpus": str(threads),
                    "memory_gb": str(mem_gb),
                }
            )

    # Write cluster/jobs.tsv
    out_jobs.parent.mkdir(parents=True, exist_ok=True)
    with out_jobs.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["sample_id", "config", "reads", "cpus", "memory_gb"],
            delimiter="\t",
        )
        writer.writeheader()
        for r in jobs_rows:
            writer.writerow(r)

    print(f"OK: wrote {len(jobs_rows)} configs to {out_configs}/")
    print(f"OK: wrote reads lists to {out_reads}/")
    print(f"OK: wrote HTCondor queue file to {out_jobs}")

if __name__ == "__main__":
    main()