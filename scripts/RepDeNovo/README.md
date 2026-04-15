# REPdenovo reproducible runner (configs + HTCondor)

This repository provides:
- a REPdenovo config template,
- a generator that creates per-sample configs and reads-list files,
- a single HTCondor submit file that queues hundreds of genomes cleanly.

## Requirements
- REPdenovo installed somewhere accessible on your system/cluster.
- Python 3 for the generator script.
- HTCondor if you submit to a cluster.

## 1) Prepare `samples.tsv`

Create a **tab-separated** file with columns:

Required:
- `sample_id` : unique ID (letters/numbers/._-)
- `fastq_1` : path to R1 fastq
- `fastq_2` : path to R2 fastq (empty for single-end)
- `group` : integer group ID for the read pair (e.g. 1)
- `mean_insert` : integer
- `sd_insert` : integer
- `genome_length` : integer
- `read_length` : integer
- `read_depth` : integer

Optional:
- `threads` (default 15)
- `memory_gb` (default 90)
- `output_folder` (default outputs/<sample_id>)
- `refiner_path` (default /usr/local/src/REPdenovo/TERefiner_1)
- `contigs_merger_path` (default /usr/local/src/REPdenovo/ContigsMerger)
- `extra_kv` : overrides any key in config, format: `KEY=VALUE;KEY=VALUE`

## 2) Generate configs + reads lists + jobs.tsv

```bash
python3 scripts/build_inputs.py --samples samples.tsv




## s05_TEsorter (plants)

This repository includes a helper script to run [TEsorter](https://github.com/zhangrengang/TEsorter) on a single TE FASTA using the plant database (`rexdb-plant`).

### Install

Recommended (conda/mamba):

```bash
mamba install -c bioconda tesorter


## s06_1_s06_2_FASTQ list + SortMeRNA/Stampy jobs (2 scripts)

This repo contains 2 small scripts:

1) `scripts/make_fastq_list.pl`
   -> finds FASTQ files under a root folder and writes a list file.

2) `scripts/generate_sortmerna_stampy_jobs.pl`
   -> reads that list, pairs R1/R2, and generates runnable bash job scripts
      (SortMeRNA + Stampy).

Raw FASTQ data and generated outputs are NOT meant to be committed to GitHub.

## What you need installed

- Perl
- SortMeRNA (available in your PATH)
- Stampy (available in your PATH)

## Folder structure (recommended)

scripts/
  make_fastq_list.pl
  generate_sortmerna_stampy_jobs.pl
config/
  fastq_folders.txt

config/fastq_folders.txt contains one folder name per line (subfolders under --root).
Empty lines and lines starting with # are ignored.

## Step 1 — Create the FASTQ list

perl scripts/make_fastq_list.pl \
  --root /path/to/fastq_root \
  --folders-file config/fastq_folders.txt \
  --out fastq.gz.list

Output: fastq.gz.list (one FASTQ path per line)

Supported extensions: .fastq, .fq, .fastq.gz, .fq.gz

## Step 2 — Generate SortMeRNA + Stampy job scripts

perl scripts/generate_sortmerna_stampy_jobs.pl \
  --filelist fastq.gz.list \
  --jobs-dir jobs/sortmerna_stampy \
  --out-root results/sortmerna_stampy \
  --sortmerna-ref "/path/to/rRNA_db.fasta,/path/to/rRNA_db.idx" \
  --stampy-genome "Malus_db" \
  --threads 8 \
  --group-size 20

This creates job scripts in jobs/sortmerna_stampy/:
- job_001.sortmerna_stampy.sh
- job_002.sortmerna_stampy.sh
- etc.

Pairing rules (expected in filenames):
- _R1 / _R2  OR
- _1  / _2

## Step 3 — Run a generated job

bash jobs/sortmerna_stampy/job_001.sortmerna_stampy.sh

Outputs will be written under:
results/sortmerna_stampy/<sample_id>/

## Recommended .gitignore

Do not commit raw data, generated jobs, or outputs:

data/raw/
jobs/
results/
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.sam
*.bam
*.cram