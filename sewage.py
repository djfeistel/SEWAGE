#!/usr/bin/env python3
"""
SEWAGE - Simulated Emulation of Wastewater-Abundance Genome Ensembles.

Simulate paired-end wastewater-like FASTQ data for a mixture of viral lineages.

The tool is wired directly to a local clone of the Freyja-barcodes repository
(https://github.com/andersen-lab/Freyja-barcodes). The user only supplies the
*pathogen name* (as it appears in the repo, e.g. ``DENV4``, ``MPX``, ``RSVa``)
and a table of lineage proportions; the reference genome and the lineage/mutation
barcode matrix are pulled from the repo automatically.

Overview of the pipeline
------------------------
1. Resolve ``<repo>/<PATHOGEN>/<version>/{reference.fasta, barcode.csv}``.
   ``barcode.csv`` has lineages in rows, mutations in the header (``A10019T`` =
   ref ``A`` -> alt ``T`` at 1-based position ``10019``) and 0.0/1.0 cells.
2. Read the desired lineage proportions, either from a user CSV/TSV
   (``lineage<sep>proportion``) or auto-generate one (equal / dominant / random).
   Every requested lineage must exist in the barcode file, otherwise the tool
   exits with a clear error naming the missing lineage(s).
3. Build a full genome for each lineage by applying *all* mutations flagged 1.0
   for that lineage to the reference (positions not mutated keep the reference
   base). There are no indels in these barcodes, so every lineage genome keeps
   the reference length.
4. Simulate whole-genome paired-end reads per lineage such that each lineage
   contributes reads in proportion to its requested abundance. Coverage is set
   either as a fold-depth (``--depth``) or as an absolute read-pair count
   (``--num-pairs``).
5. Write gzipped R1/R2 FASTQ files plus a small manifest describing the exact
   proportions and read counts used.
"""

import argparse
import csv
import gzip
import os
import re
import sys
import time

import numpy as np


# ASCII wordmark shown at the top of the --help screen: a sewer manhole cover
# on the left, SEWAGE in tilted 3D (drop-shadowed) block letters, and a spiky
# virus particle on the right.
BANNER = r"""
   .-==-.       _____ _______       _____   ____________         \ | /
  /::||::\     / ___// ____/ |     / /   | / ____/ ____/.     .--(#)--.
 |::-##-::|    \__ \/ __/..| |./| / / /| |/ / __/ __/....    -(  ###  )-
  \::||::/    ___/ / /___. | |/ |/ / ___ / /_/ / /___.        '--(#)--'
   '-==-'    /____/_____/. |__/|__/_/ .|_\____/_____/.           / | \
              ............  ..........  ..............
        Simulated Emulation of Wastewater-Abundance Genome Ensembles
"""


# Default location of the Freyja-barcodes clone. Can be overridden with --repo
# or the FREYJA_BARCODES environment variable.
DEFAULT_REPO = os.environ.get(
    "FREYJA_BARCODES",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "Freyja-barcodes"),
)

# Integer encoding used internally for fast, vectorized read generation.
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)          # index -> base byte
BASE_TO_INT = {65: 0, 67: 1, 71: 2, 84: 3}              # A C G T -> 0..3
# Complement under the 0=A,1=C,2=G,3=T encoding is simply (3 - base).

# IUPAC ambiguity codes -> the concrete bases they can represent. Used when a
# barcode alt allele is ambiguous (a few tokens like ``T6541Y`` exist).
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

MUT_RE = re.compile(r"^([ACGT])(\d+)([ACGTRYSWKMBDHVN])$")


# --------------------------------------------------------------------------- #
# Repo / data loading
# --------------------------------------------------------------------------- #
def list_pathogens(repo):
    """Return the sorted list of pathogen folders that contain a ``latest`` dir."""
    out = []
    for name in sorted(os.listdir(repo)):
        path = os.path.join(repo, name)
        if os.path.isdir(path) and os.path.isdir(os.path.join(path, "latest")):
            out.append(name)
    return out


def resolve_pathogen_dir(repo, pathogen, version):
    """Resolve and validate the data directory for a pathogen + version."""
    if not os.path.isdir(repo):
        sys.exit(f"ERROR: Freyja-barcodes repo not found: {repo}\n"
                 f"       Pass --repo or set FREYJA_BARCODES.")
    pdir = os.path.join(repo, pathogen)
    if not os.path.isdir(pdir):
        avail = ", ".join(list_pathogens(repo))
        sys.exit(f"ERROR: pathogen '{pathogen}' not found in {repo}.\n"
                 f"       Available: {avail}")
    vdir = os.path.join(pdir, version)
    if not os.path.isdir(vdir):
        versions = sorted(d for d in os.listdir(pdir)
                          if os.path.isdir(os.path.join(pdir, d)))
        sys.exit(f"ERROR: version '{version}' not found for {pathogen}.\n"
                 f"       Available: {', '.join(versions)}")
    ref = os.path.join(vdir, "reference.fasta")
    bar = os.path.join(vdir, "barcode.csv")
    for f in (ref, bar):
        if not os.path.isfile(f):
            sys.exit(f"ERROR: expected file missing: {f}")
    return ref, bar


def load_reference(path):
    """Load a single-sequence FASTA. Returns (header, sequence_string)."""
    header = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break  # only the first record is used
                header = line[1:].split()[0]
            else:
                seq.append(line.strip().upper())
    if header is None or not seq:
        sys.exit(f"ERROR: no sequence found in {path}")
    return header, "".join(seq)


def load_barcode(path):
    """Load barcode.csv.

    Returns (mutations, lineage_muts) where ``mutations`` is the ordered list of
    header tokens and ``lineage_muts`` maps lineage name -> list of mutation
    tokens flagged 1.0 for that lineage.
    """
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        mutations = header[1:]
        lineage_muts = {}
        for row in reader:
            if not row or not row[0].strip():
                continue
            name = row[0].strip()
            muts = [mutations[i] for i, v in enumerate(row[1:])
                    if v.strip() not in ("", "0", "0.0")]
            lineage_muts[name] = muts
    return mutations, lineage_muts


# --------------------------------------------------------------------------- #
# Proportions handling
# --------------------------------------------------------------------------- #
def read_proportions(path):
    """Read a lineage/proportion CSV or TSV.

    Two columns: lineage, proportion. An optional header row (where column 2 is
    not a float) is skipped. Delimiter is auto-detected (comma or tab).
    """
    with open(path) as fh:
        sample = fh.read()
    if not sample.strip():
        sys.exit(f"ERROR: proportions file is empty: {path}")
    delim = "\t" if sample.count("\t") >= sample.count(",") else ","
    rows = [ln for ln in sample.splitlines() if ln.strip()]
    pairs = []
    for i, ln in enumerate(rows):
        parts = [p.strip() for p in ln.split(delim)]
        if len(parts) < 2:
            sys.exit(f"ERROR: line {i + 1} in {path} needs 2 columns "
                     f"(lineage, proportion): {ln!r}")
        name, val = parts[0], parts[1]
        try:
            prop = float(val)
        except ValueError:
            if i == 0:
                continue  # header row
            sys.exit(f"ERROR: non-numeric proportion on line {i + 1}: {val!r}")
        pairs.append((name, prop))
    if not pairs:
        sys.exit(f"ERROR: no data rows parsed from {path}")
    return pairs


def normalize(pairs):
    """Normalize proportions to sum to 1.0, preserving order and lineage names."""
    total = sum(p for _, p in pairs)
    if total <= 0:
        sys.exit("ERROR: proportions must sum to a positive value.")
    return [(n, p / total) for n, p in pairs]


def generate_proportions(available, num, mode, rng, dominant_fraction=None):
    """Auto-generate a lineage/proportion list from the available lineages.

    mode:
      - ``equal``    : every lineage gets 1/num.
      - ``dominant`` : one lineage dominates (dominant_fraction, default random
                       in [0.5, 0.8]); the rest share the remainder unequally.
      - ``random``   : fully random mix (Dirichlet).
    """
    if num < 1:
        sys.exit("ERROR: --num-lineages must be >= 1.")
    if num > len(available):
        sys.exit(f"ERROR: requested {num} lineages but only {len(available)} "
                 f"are available for this pathogen.")
    chosen = list(rng.choice(available, size=num, replace=False))

    if mode == "equal":
        props = np.full(num, 1.0 / num)
    elif mode == "dominant":
        if num == 1:
            props = np.array([1.0])
        else:
            frac = (dominant_fraction if dominant_fraction is not None
                    else float(rng.uniform(0.5, 0.8)))
            frac = min(max(frac, 1.0 / num), 0.99)
            rest = rng.dirichlet(np.ones(num - 1)) * (1.0 - frac)
            props = np.concatenate([[frac], rest])
    elif mode == "random":
        props = rng.dirichlet(np.ones(num))
    else:
        sys.exit(f"ERROR: unknown proportion mode: {mode}")

    return list(zip(chosen, (float(x) for x in props)))


def prompt_generation(available, rng):
    """Interactively ask the user how to auto-generate proportions."""
    print(f"\nAuto-generating proportions. "
          f"{len(available)} lineages available for this pathogen.")
    while True:
        try:
            num = int(input("1) How many lineages do you want in total? ").strip())
            break
        except ValueError:
            print("   Please enter an integer.")
    print("2) Proportion scheme:")
    print("   [1] equal proportions across all lineages")
    print("   [2] one dominant lineage, others lower and mixed (not equal)")
    print("   [3] completely mixed at random")
    choice = input("   Choose 1/2/3: ").strip()
    mode = {"1": "equal", "2": "dominant", "3": "random"}.get(choice)
    if mode is None:
        sys.exit("ERROR: invalid choice; expected 1, 2 or 3.")
    return generate_proportions(available, num, mode, rng)


def write_proportions_file(path, pairs):
    """Write the (lineage, proportion) table as TSV with a header."""
    with open(path, "w") as fh:
        fh.write("lineage\tproportion\n")
        for name, prop in pairs:
            fh.write(f"{name}\t{prop:.6f}\n")


# --------------------------------------------------------------------------- #
# Genome construction
# --------------------------------------------------------------------------- #
def build_lineage_genome(ref_ints, muts, rng, warn_prefix=""):
    """Return a copy of the reference (int-encoded) with the lineage's mutations.

    Applies every mutation token flagged for the lineage. Ambiguous IUPAC alt
    alleles are resolved to a concrete base (preferring one that differs from
    the reference). Malformed / out-of-range tokens are skipped with a warning.
    """
    genome = ref_ints.copy()
    n = len(genome)
    skipped = 0
    for tok in muts:
        m = MUT_RE.match(tok)
        if not m:
            skipped += 1
            continue
        _ref_base, pos_s, alt = m.group(1), m.group(2), m.group(3)
        pos = int(pos_s) - 1  # 1-based -> 0-based
        if pos < 0 or pos >= n:
            skipped += 1
            continue
        options = IUPAC.get(alt, "")
        if not options:
            skipped += 1
            continue
        # Prefer a concrete alt base that actually changes the reference.
        ref_here = genome[pos]
        alt_int = None
        for base in options:
            bi = BASE_TO_INT[ord(base)]
            if bi != ref_here:
                alt_int = bi
                break
        if alt_int is None:
            alt_int = BASE_TO_INT[ord(options[0])]
        genome[pos] = alt_int
    if skipped and warn_prefix:
        print(f"{warn_prefix}: skipped {skipped} unparseable/out-of-range "
              f"mutation token(s).", file=sys.stderr)
    return genome


# --------------------------------------------------------------------------- #
# Read simulation
# --------------------------------------------------------------------------- #
def phred_char_array(qs):
    """Convert an int quality array to a Phred+33 byte array."""
    return (np.clip(qs, 0, 40).astype(np.uint8) + 33)


def simulate_reads(genome, n_pairs, read_len, frag_mean, frag_sd, error_rate,
                   rng, name_prefix, r1_out, r2_out, chunk=100000):
    """Simulate ``n_pairs`` paired-end reads from one genome, writing FASTQ.

    Reads tile the whole genome uniformly. R1 is the forward strand at the
    fragment 5' end; R2 is the reverse-complement at the fragment 3' end.
    """
    glen = len(genome)
    if glen == 0 or n_pairs <= 0:
        return 0
    rl = min(read_len, glen)
    base_q = max(2, int(round(-10.0 * np.log10(max(error_rate, 1e-6)))))

    # Fixed zero-padded width for read ids so every FASTQ record has identical
    # length and can be assembled with pure-numpy (no per-read Python loop).
    id_width = max(1, len(str(n_pairs - 1)))

    written = 0
    offsets = np.arange(rl, dtype=np.int64)
    done = 0
    while done < n_pairs:
        m = min(chunk, n_pairs - done)

        # Fragment lengths, clamped to [rl, glen].
        frag = rng.normal(frag_mean, frag_sd, m)
        frag = np.clip(frag, rl, glen).astype(np.int64)
        max_start = glen - frag                       # inclusive upper bound
        start = (rng.random(m) * (max_start + 1)).astype(np.int64)

        # Forward (R1) and reverse-end (R2) base indices.
        r1_idx = start[:, None] + offsets              # m x rl
        r2_region_start = start + frag - rl
        r2_idx = r2_region_start[:, None] + offsets    # m x rl (forward region)

        r1 = genome[r1_idx]                            # int-encoded bases
        r2_fwd = genome[r2_idx]
        # Reverse-complement the R2 region: reverse along read, complement 3-b.
        r2 = 3 - r2_fwd[:, ::-1]

        # Introduce sequencing errors.
        q1 = np.full((m, rl), base_q, dtype=np.int64)
        q2 = np.full((m, rl), base_q, dtype=np.int64)
        if error_rate > 0:
            for arr, q in ((r1, q1), (r2, q2)):
                mask = rng.random(arr.shape) < error_rate
                if mask.any():
                    # Replace with a different random base (shift by 1..3 mod 4).
                    shift = rng.integers(1, 4, size=int(mask.sum()))
                    arr[mask] = (arr[mask] + shift) % 4
                    q[mask] = np.maximum(2, base_q // 2)

        # Encode to ASCII bytes.
        r1_seq = BASES[r1]                             # m x rl uint8
        r2_seq = BASES[r2]
        q1c = phred_char_array(q1)
        q2c = phred_char_array(q2)

        _write_fastq_chunk(r1_out, name_prefix, done, r1_seq, q1c, 1, id_width)
        _write_fastq_chunk(r2_out, name_prefix, done, r2_seq, q2c, 2, id_width)

        written += m
        done += m
    return written


def _write_fastq_chunk(out, name_prefix, base_idx, seq_bytes, qual_bytes, mate,
                       id_width):
    """Write a chunk of reads to an open FASTQ handle, fully vectorized.

    Every record has fixed byte length (read ids are zero-padded to
    ``id_width``), so the whole chunk is assembled as a single (m, L) uint8
    array and emitted with one ``tobytes()`` call - no per-read Python loop.
    """
    m, rl = seq_bytes.shape

    # Fixed byte segments surrounding the variable id / seq / qual.
    head = f"@{name_prefix}:".encode("ascii")          # before the id
    tail = f" {mate}:N:0:1\n".encode("ascii")          # after the id, ends line
    mid = b"\n+\n"                                      # between seq and qual
    end = b"\n"                                         # after qual

    head_a = np.frombuffer(head, dtype=np.uint8)
    tail_a = np.frombuffer(tail, dtype=np.uint8)
    mid_a = np.frombuffer(mid, dtype=np.uint8)
    end_a = np.frombuffer(end, dtype=np.uint8)

    lh, lt, lm, le = len(head_a), len(tail_a), len(mid_a), len(end_a)
    row_len = lh + id_width + lt + rl + lm + rl + le
    out_arr = np.empty((m, row_len), dtype=np.uint8)

    p = 0
    out_arr[:, p:p + lh] = head_a
    p += lh

    # Zero-padded ascii ids, computed digit-by-digit across the whole chunk.
    ids = np.arange(base_idx, base_idx + m, dtype=np.int64)
    powers = 10 ** np.arange(id_width - 1, -1, -1, dtype=np.int64)
    digits = (ids[:, None] // powers) % 10
    out_arr[:, p:p + id_width] = digits.astype(np.uint8) + ord("0")
    p += id_width

    out_arr[:, p:p + lt] = tail_a
    p += lt
    out_arr[:, p:p + rl] = seq_bytes
    p += rl
    out_arr[:, p:p + lm] = mid_a
    p += lm
    out_arr[:, p:p + rl] = qual_bytes
    p += rl
    out_arr[:, p:p + le] = end_a

    out.write(out_arr.tobytes())


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
class SewageHelpFormatter(argparse.RawDescriptionHelpFormatter,
                          argparse.ArgumentDefaultsHelpFormatter):
    """Keep the banner/description verbatim while still showing arg defaults."""


class SewageArgumentParser(argparse.ArgumentParser):
    """ArgumentParser that prints the manhole-cover banner above everything."""

    def format_help(self):
        return BANNER + "\n" + super().format_help()


def build_parser():
    p = SewageArgumentParser(
        prog="sewage",
        description="SEWAGE (Simulated Emulation of Wastewater-Abundance\n"
                    "Genome Ensembles): simulate paired-end wastewater-like\n"
                    "FASTQ for a mixture of viral lineages, using\n"
                    "references/barcodes from a local Freyja-barcodes clone.",
        formatter_class=SewageHelpFormatter,
    )
    p.add_argument("-p", "--pathogen",
                   help="Pathogen name as it appears in the Freyja-barcodes "
                        "repo (e.g. DENV4, MPX, RSVa, MEASLESgenome).")
    p.add_argument("--repo", default=DEFAULT_REPO,
                   help="Path to the Freyja-barcodes repository.")
    p.add_argument("--version", default="latest",
                   help="Barcode version subfolder to use (e.g. latest or a "
                        "date like 2025-05-01).")
    p.add_argument("--list", action="store_true",
                   help="List available pathogens in the repo and exit.")
    p.add_argument("--list-lineages", action="store_true",
                   help="List available lineages for --pathogen and exit.")

    # Proportions input.
    p.add_argument("-i", "--proportions",
                   help="CSV/TSV file: column 1 = lineage, column 2 = "
                        "proportion. Lineages must exist in the barcode file.")
    p.add_argument("--generate-proportions", action="store_true",
                   help="Auto-generate the proportions table instead of "
                        "supplying one. Prompts if --num-lineages/--prop-mode "
                        "are not given.")
    p.add_argument("--num-lineages", type=int,
                   help="(generate) Number of lineages to include.")
    p.add_argument("--prop-mode", choices=["equal", "dominant", "random"],
                   help="(generate) Proportion scheme.")
    p.add_argument("--dominant-fraction", type=float,
                   help="(generate, dominant mode) Fraction for the dominant "
                        "lineage (default: random in [0.5, 0.8]).")
    p.add_argument("--proportions-out",
                   help="Where to save the (generated or used) proportions "
                        "table. Default: <output_prefix>.proportions.tsv")

    # Coverage: exactly one of --depth / --num-pairs.
    p.add_argument("--depth", type=float,
                   help="Target fold coverage depth for the whole sample. Read "
                        "pairs are derived from genome length and read length.")
    p.add_argument("--num-pairs", type=int,
                   help="Total number of read PAIRS for the whole sample.")

    # Read geometry / errors.
    p.add_argument("--read-length", type=int, default=250,
                   help="Length of each mate.")
    p.add_argument("--fragment-mean", type=float, default=500.0,
                   help="Mean fragment (insert) size.")
    p.add_argument("--fragment-sd", type=float, default=50.0,
                   help="Std. dev. of fragment size.")
    p.add_argument("--error-rate", type=float, default=0.005,
                   help="Per-base sequencing error rate (0 to disable).")

    p.add_argument("-o", "--output-prefix", default="sim_sample",
                   help="Prefix for output FASTQ files "
                        "(<prefix>_R1.fastq.gz / <prefix>_R2.fastq.gz).")
    p.add_argument("--gzip", dest="gzip", action="store_true", default=True,
                   help="Gzip the FASTQ output (default).")
    p.add_argument("--no-gzip", dest="gzip", action="store_false",
                   help="Write plain (uncompressed) FASTQ.")
    p.add_argument("--gzip-level", type=int, default=6,
                   help="gzip compression level 1-9. Lower is much faster with "
                        "slightly larger files; 9 (Python's default) is the "
                        "slowest and a common bottleneck for big outputs.")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducibility.")
    p.add_argument("--timing", action="store_true",
                   help="Print elapsed wall-time per phase (genome build, "
                        "read simulation) to expose bottlenecks.")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    rng = np.random.default_rng(args.seed)

    # --- Informational modes ------------------------------------------------
    if args.list:
        print("\n".join(list_pathogens(args.repo)))
        return 0

    if not args.pathogen:
        sys.exit("ERROR: --pathogen is required (or use --list).")

    ref_path, bar_path = resolve_pathogen_dir(args.repo, args.pathogen,
                                              args.version)
    ref_name, ref_seq = load_reference(ref_path)
    mutations, lineage_muts = load_barcode(bar_path)
    available = sorted(lineage_muts)

    if args.list_lineages:
        print(f"# {args.pathogen} ({args.version}): {len(available)} lineages")
        print("\n".join(available))
        return 0

    # --- Coverage validation ------------------------------------------------
    if (args.depth is None) == (args.num_pairs is None):
        sys.exit("ERROR: specify exactly one of --depth or --num-pairs.")

    # --- Resolve proportions ------------------------------------------------
    if args.generate_proportions:
        if args.num_lineages is not None and args.prop_mode is not None:
            pairs = generate_proportions(
                available, args.num_lineages, args.prop_mode, rng,
                dominant_fraction=args.dominant_fraction)
        else:
            pairs = prompt_generation(available, rng)
    elif args.proportions:
        pairs = read_proportions(args.proportions)
    else:
        sys.exit("ERROR: supply -i/--proportions or --generate-proportions.")

    # Validate lineages exist in the barcode file.
    missing = [n for n, _ in pairs if n not in lineage_muts]
    if missing:
        sys.exit("ERROR: the following lineage(s) are not in the "
                 f"{args.pathogen} barcode file: {', '.join(missing)}\n"
                 f"       Use --list-lineages to see valid names.")

    pairs = normalize(pairs)

    # Save the proportions table used (for tracking / reuse).
    prop_out = (args.proportions_out
                or f"{args.output_prefix}.proportions.tsv")
    write_proportions_file(prop_out, pairs)

    # --- Build lineage genomes ---------------------------------------------
    ref_ints = np.frombuffer(ref_seq.encode("ascii"), dtype=np.uint8).copy()
    # Map ambiguous/other reference letters to A by default, real bases to 0..3.
    ref_encoded = np.zeros(len(ref_ints), dtype=np.int64)
    for byte, val in BASE_TO_INT.items():
        ref_encoded[ref_ints == byte] = val
    glen = len(ref_encoded)

    genomes = {}
    t_build = time.perf_counter()
    for name, _prop in pairs:
        genomes[name] = build_lineage_genome(
            ref_encoded, lineage_muts[name], rng, warn_prefix=f"WARN[{name}]")
    build_secs = time.perf_counter() - t_build

    # --- Compute per-lineage read-pair counts ------------------------------
    if args.depth is not None:
        # depth = total sample fold-coverage; lineage i gets prop_i of it.
        total_pairs = int(round(args.depth * glen / (2 * args.read_length)))
    else:
        total_pairs = int(args.num_pairs)
    if total_pairs <= 0:
        sys.exit("ERROR: computed 0 read pairs; increase --depth/--num-pairs.")

    counts = {name: int(round(prop * total_pairs)) for name, prop in pairs}

    # --- Simulate & write ---------------------------------------------------
    if args.gzip:
        def opener(path, mode):
            return gzip.open(path, mode, compresslevel=args.gzip_level)
    else:
        opener = open
    suffix = ".fastq.gz" if args.gzip else ".fastq"
    r1_path = f"{args.output_prefix}_R1{suffix}"
    r2_path = f"{args.output_prefix}_R2{suffix}"

    print("SEWAGE - Simulated Emulation of Wastewater-Abundance Genome Ensembles")
    print(f"Pathogen : {args.pathogen} ({args.version})")
    print(f"Reference: {ref_name}  ({glen} bp)")
    print(f"Coverage : {'depth ' + str(args.depth) + 'x' if args.depth else ''}"
          f"{'num-pairs ' + str(args.num_pairs) if args.num_pairs else ''} "
          f"-> {total_pairs} total read pairs")
    print(f"Read len : {args.read_length}  fragment: "
          f"{args.fragment_mean}±{args.fragment_sd}  err: {args.error_rate}")
    print(f"Lineages : {len(pairs)}")

    total_written = 0
    t_sim = time.perf_counter()
    with opener(r1_path, "wb") as r1, opener(r2_path, "wb") as r2:
        for name, prop in pairs:
            n = counts[name]
            prefix = f"{args.pathogen}:{name}"
            t_lin = time.perf_counter()
            w = simulate_reads(
                genomes[name], n, args.read_length, args.fragment_mean,
                args.fragment_sd, args.error_rate, rng, prefix, r1, r2)
            total_written += w
            lin_secs = time.perf_counter() - t_lin
            if args.timing:
                rate = w / lin_secs if lin_secs > 0 else 0.0
                print(f"  {name:<24} prop={prop:6.4f}  pairs={w}  "
                      f"[{lin_secs:6.2f}s, {rate:,.0f} pairs/s]")
            else:
                print(f"  {name:<24} prop={prop:6.4f}  pairs={w}")
    sim_secs = time.perf_counter() - t_sim

    # Update the manifest with realized read counts.
    with open(prop_out, "w") as fh:
        fh.write("lineage\tproportion\tread_pairs\n")
        for name, prop in pairs:
            fh.write(f"{name}\t{prop:.6f}\t{counts[name]}\n")

    print(f"\nWrote {total_written} read pairs")
    print(f"  R1: {r1_path}")
    print(f"  R2: {r2_path}")
    print(f"  proportions/manifest: {prop_out}")
    if args.timing:
        total_secs = build_secs + sim_secs
        rate = total_written / sim_secs if sim_secs > 0 else 0.0
        print("\nTiming:")
        print(f"  genome build : {build_secs:7.2f}s ({len(pairs)} lineages)")
        print(f"  read sim+write: {sim_secs:7.2f}s "
              f"({rate:,.0f} pairs/s, gzip level {args.gzip_level if args.gzip else 'off'})")
        print(f"  total        : {total_secs:7.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
