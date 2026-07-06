# SEWAGE

**S**imulated **E**mulation of **W**astewater-**A**bundance **G**enome **E**nsembles

![SEWAGE](assets/sewage_logo.png)

SEWAGE simulates paired-end, wastewater-like FASTQ data for a defined mixture
of viral lineages. It is wired directly to a local clone of the
[Freyja-barcodes](https://github.com/andersen-lab/Freyja-barcodes) repository:
you supply a **pathogen name** and a table of **lineage proportions**, and
SEWAGE pulls the reference genome and the lineage/mutation barcode matrix from
the repo automatically, builds a full genome per lineage, and emits reads so
that each lineage contributes in proportion to its requested abundance.

This is useful for benchmarking wastewater relative-abundance callers such as
[Freyja](https://github.com/andersen-lab/Freyja), where you need FASTQ data with
a **known ground-truth mixture** of lineages.

---

## Features

- Pulls references/barcodes straight from a local `Freyja-barcodes` clone — no
  manual FASTA/VCF wrangling.
- Works for any pathogen present in the repo (DENV1–4, MPX, RSVa/RSVb, MEASLES,
  influenza segments, MTB, etc.).
- Builds each lineage's genome by applying **all** barcode mutations flagged for
  that lineage (IUPAC ambiguity codes resolved to a concrete alt base).
- Proportions supplied via a CSV/TSV file **or** auto-generated
  (`equal` / `dominant` / `random`).
- Coverage controlled by fold-**depth** (`--depth`) or absolute **read-pair
  count** (`--num-pairs`).
- Configurable read length, fragment size distribution, and per-base error rate.
- Fast, `numpy`-vectorized read generation; gzip output by default.
- Writes a manifest of the exact realized proportions and read-pair counts.

---

## Requirements

- Python **3.8+** (developed on 3.12)
- [`numpy`](https://numpy.org/) (only external dependency)
- A local clone of the **Freyja-barcodes** repository

Install the Python dependency:

```bash
pip install -r requirements.txt
```

Clone the barcode repository (anywhere you like):

```bash
git clone https://github.com/andersen-lab/Freyja-barcodes.git
```

SEWAGE finds the barcode repo, in order of precedence:

1. the `--repo /path/to/Freyja-barcodes` flag,
2. the `FREYJA_BARCODES` environment variable,
3. a `Freyja-barcodes/` folder sitting next to `sewage.py`.

```bash
export FREYJA_BARCODES=/path/to/Freyja-barcodes
```

---

## Quick start

```bash
# 1. See which pathogens are available in your barcode clone
./sewage.py --list

# 2. See which lineages exist for a pathogen
./sewage.py -p DENV1 --list-lineages

# 3. Simulate a 3-lineage mixture at 100x depth from a proportions file
./sewage.py -p DENV1 -i proportions.csv --depth 100 -o DENV1_sample

# 4. ...or let SEWAGE auto-generate a dominant-strain mixture of 4 lineages
./sewage.py -p DENV1 --generate-proportions --num-lineages 4 \
            --prop-mode dominant --depth 100 -o DENV1_sample --seed 42
```

### Proportions file format

A two-column CSV or TSV (comma or tab auto-detected; an optional header is
ignored). Column 1 is the lineage name (must exist in the barcode file),
column 2 is its proportion. Values are normalized automatically, so they need
not sum to exactly 1.

```csv
lineage,proportion
DENV1-1I-H,0.6
DENV1-III,0.3
DENV1-V,0.1
```

Every requested lineage must exist in the barcode file; otherwise SEWAGE exits
with an error naming the missing lineage(s).

---

## Output

For `-o sample` (gzip on by default):

| File | Description |
|------|-------------|
| `sample_R1.fastq.gz` | Forward reads |
| `sample_R2.fastq.gz` | Reverse reads (reverse-complemented mates) |
| `sample.proportions.tsv` | Manifest: requested proportion + realized read-pair count per lineage |

Read names follow an Illumina-like scheme encoding the source lineage, e.g.:

```
@DENV1:DENV1-1I-H:0 1:N:0:1
```

---

## Command-line options

| Flag | Description | Default |
|------|-------------|---------|
| `-p`, `--pathogen` | Pathogen name as it appears in the repo (e.g. `DENV4`, `MPX`, `RSVa`). | — |
| `--repo` | Path to the Freyja-barcodes repository. | `$FREYJA_BARCODES` or `./Freyja-barcodes` |
| `--version` | Barcode version subfolder (e.g. `latest` or a date). | `latest` |
| `--list` | List available pathogens and exit. | — |
| `--list-lineages` | List available lineages for `--pathogen` and exit. | — |
| `-i`, `--proportions` | CSV/TSV proportions file (`lineage,proportion`). | — |
| `--generate-proportions` | Auto-generate proportions instead of supplying a file. | — |
| `--num-lineages` | (generate) Number of lineages to include. | prompt |
| `--prop-mode` | (generate) `equal`, `dominant`, or `random`. | prompt |
| `--dominant-fraction` | (generate, dominant) Fraction for the dominant lineage. | random in [0.5, 0.8] |
| `--proportions-out` | Where to save the used/generated proportions table. | `<prefix>.proportions.tsv` |
| `--depth` | Target fold coverage for the whole sample. | — |
| `--num-pairs` | Total number of read **pairs** for the whole sample. | — |
| `--read-length` | Length of each mate. | `250` |
| `--fragment-mean` | Mean fragment (insert) size. | `500` |
| `--fragment-sd` | Std. dev. of fragment size. | `50` |
| `--error-rate` | Per-base sequencing error rate (`0` disables). | `0.005` |
| `-o`, `--output-prefix` | Prefix for output FASTQ files. | `sim_sample` |
| `--gzip` / `--no-gzip` | Gzip the FASTQ output (on by default). | gzip |
| `--seed` | Random seed for reproducibility. | — |

> **Note:** `--depth` and `--num-pairs` are mutually exclusive — provide exactly one.

Run `./sewage.py -h` for the full, authoritative help text.

---

## Using SEWAGE output in a pipeline

The reads are standard paired-end FASTQ and drop straight into tools like
`fastp`, `bwa`, `minimap2`, or nf-core-style workflows via a samplesheet:

```csv
sample,platform,fastq_1,fastq_2,primer_bed
DENV1-01-SEWAGE,illumina,/path/DENV1_sample_R1.fastq.gz,/path/DENV1_sample_R2.fastq.gz,
```

> **Gotcha:** if a downstream step reports `invalid gzip header`, the file that
> actually reached the tool was not gzip-compressed (e.g. an uncompressed copy
> got staged/uploaded under a `.fastq.gz` name). Confirm the real gzip files
> pass `gzip -t` and that whatever you upload/stage matches the extension in
> your samplesheet.

---

## How it works

1. Resolve `<repo>/<PATHOGEN>/<version>/{reference.fasta, barcode.csv}`.
   Barcode columns are mutation tokens like `A10019T` (ref `A` → alt `T` at
   1-based position `10019`); cells are `0.0` / `1.0`.
2. Read or generate the desired lineage proportions and normalize them.
3. Build a full genome per lineage by applying every mutation flagged for that
   lineage to the reference (there are no indels in these barcodes, so genomes
   keep the reference length).
4. Simulate whole-genome paired-end reads per lineage, with each lineage
   contributing read pairs in proportion to its abundance.
5. Write gzipped R1/R2 FASTQ plus the proportions manifest.

---

## Full help menu

Complete output of `./sewage.py -h`:

```
   .-==-.       _____ _______       _____   ____________         \ | /
  /::||::\     / ___// ____/ |     / /   | / ____/ ____/.     .--(#)--.
 |::-##-::|    \__ \/ __/..| |./| / / /| |/ / __/ __/....    -(  ###  )-
  \::||::/    ___/ / /___. | |/ |/ / ___ / /_/ / /___.        '--(#)--'
   '-==-'    /____/_____/. |__/|__/_/ .|_\____/_____/.           / | \
              ............  ..........  ..............
        Simulated Emulation of Wastewater-Abundance Genome Ensembles

usage: sewage [-h] [-p PATHOGEN] [--repo REPO] [--version VERSION] [--list]
              [--list-lineages] [-i PROPORTIONS] [--generate-proportions]
              [--num-lineages NUM_LINEAGES]
              [--prop-mode {equal,dominant,random}]
              [--dominant-fraction DOMINANT_FRACTION]
              [--proportions-out PROPORTIONS_OUT] [--depth DEPTH]
              [--num-pairs NUM_PAIRS] [--read-length READ_LENGTH]
              [--fragment-mean FRAGMENT_MEAN] [--fragment-sd FRAGMENT_SD]
              [--error-rate ERROR_RATE] [-o OUTPUT_PREFIX] [--gzip]
              [--no-gzip] [--gzip-level GZIP_LEVEL] [--seed SEED] [--timing]

SEWAGE (Simulated Emulation of Wastewater-Abundance
Genome Ensembles): simulate paired-end wastewater-like
FASTQ for a mixture of viral lineages, using
references/barcodes from a local Freyja-barcodes clone.

options:
  -h, --help            show this help message and exit
  -p PATHOGEN, --pathogen PATHOGEN
                        Pathogen name as it appears in the Freyja-barcodes
                        repo (e.g. DENV4, MPX, RSVa, MEASLESgenome). (default:
                        None)
  --repo REPO           Path to the Freyja-barcodes repository. (default:
                        /home/ofx5/GENOME_PROPORTION_GENERATOR/SEWAGE/Freyja-
                        barcodes)
  --version VERSION     Barcode version subfolder to use (e.g. latest or a
                        date like 2025-05-01). (default: latest)
  --list                List available pathogens in the repo and exit.
                        (default: False)
  --list-lineages       List available lineages for --pathogen and exit.
                        (default: False)
  -i PROPORTIONS, --proportions PROPORTIONS
                        CSV/TSV file: column 1 = lineage, column 2 =
                        proportion. Lineages must exist in the barcode file.
                        (default: None)
  --generate-proportions
                        Auto-generate the proportions table instead of
                        supplying one. Prompts if --num-lineages/--prop-mode
                        are not given. (default: False)
  --num-lineages NUM_LINEAGES
                        (generate) Number of lineages to include. (default:
                        None)
  --prop-mode {equal,dominant,random}
                        (generate) Proportion scheme. (default: None)
  --dominant-fraction DOMINANT_FRACTION
                        (generate, dominant mode) Fraction for the dominant
                        lineage (default: random in [0.5, 0.8]). (default:
                        None)
  --proportions-out PROPORTIONS_OUT
                        Where to save the (generated or used) proportions
                        table. Default: <output_prefix>.proportions.tsv
                        (default: None)
  --depth DEPTH         Target fold coverage depth for the whole sample. Read
                        pairs are derived from genome length and read length.
                        (default: None)
  --num-pairs NUM_PAIRS
                        Total number of read PAIRS for the whole sample.
                        (default: None)
  --read-length READ_LENGTH
                        Length of each mate. (default: 250)
  --fragment-mean FRAGMENT_MEAN
                        Mean fragment (insert) size. (default: 500.0)
  --fragment-sd FRAGMENT_SD
                        Std. dev. of fragment size. (default: 50.0)
  --error-rate ERROR_RATE
                        Per-base sequencing error rate (0 to disable).
                        (default: 0.005)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix for output FASTQ files (<prefix>_R1.fastq.gz /
                        <prefix>_R2.fastq.gz). (default: sim_sample)
  --gzip                Gzip the FASTQ output (default). (default: True)
  --no-gzip             Write plain (uncompressed) FASTQ. (default: True)
  --gzip-level GZIP_LEVEL
                        gzip compression level 1-9. Lower is much faster with
                        slightly larger files; 9 (Python's default) is the
                        slowest and a common bottleneck for big outputs.
                        (default: 6)
  --seed SEED           Random seed for reproducibility. (default: None)
  --timing              Print elapsed wall-time per phase (genome build, read
                        simulation) to expose bottlenecks. (default: False)
```

---

## License

Released under the [MIT License](LICENSE). Update the copyright holder in
`LICENSE` to your name/organization before publishing.
