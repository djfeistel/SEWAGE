#!/usr/bin/env python3
"""
SEWAGE - Simulated Emulation of Wastewater-Abundance Genome Ensembles.

Simulate paired-end wastewater-like FASTQ data for a mixture of viral lineages.

The reference genomes and lineage/mutation barcode matrices come from the
Freyja-barcodes project (https://github.com/andersen-lab/Freyja-barcodes), but
that data is *bundled inside SEWAGE* (the ``data/`` folder) so users never have
to clone it themselves. Run ``sewage --update`` to download the latest upstream
data and rebuild ``data/`` (all pathogens, ``latest`` plus dated versions; only
``reference.fasta`` and ``barcode.csv`` are kept). The user only supplies the
*pathogen name* (as it appears in the data, e.g. ``DENV4``, ``MPX``, ``RSVa``)
and a table of lineage proportions.

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


# Bundled barcode data lives inside the SEWAGE repo itself (data/), so users
# do NOT have to clone Freyja-barcodes separately. Run ``sewage --update`` to
# download the latest upstream data and rebuild this folder. --repo or the
# FREYJA_BARCODES environment variable can point at a different data directory.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
DEFAULT_REPO = os.environ.get("FREYJA_BARCODES", DATA_DIR)

# Upstream source used by --update and the per-version files worth vendoring
# (everything else in the repo -- READMEs, HTML, auspice trees -- is skipped).
UPSTREAM_REPO = "andersen-lab/Freyja-barcodes"
UPSTREAM_BRANCH = "main"
KEEP_FILES = ("reference.fasta", "barcode.csv")

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
    if not os.path.isdir(repo):
        return []
    out = []
    for name in sorted(os.listdir(repo)):
        path = os.path.join(repo, name)
        if os.path.isdir(path) and os.path.isdir(os.path.join(path, "latest")):
            out.append(name)
    return out


def update_barcodes(data_dir, repo=UPSTREAM_REPO, branch=UPSTREAM_BRANCH):
    """Download the upstream Freyja-barcodes repo and (re)build ``data_dir``.

    Only per-version ``reference.fasta`` and ``barcode.csv`` files are kept, for
    every pathogen and every version (``latest`` plus dated folders). New
    pathogens/versions are added and ones removed upstream disappear locally,
    so the bundled data mirrors the upstream repo. The data is staged in a
    temporary directory and atomically swapped into place so an interrupted
    update never corrupts the existing data. Returns
    ``(n_pathogens, n_versions, n_files)``.
    """
    import io
    import shutil
    import tarfile
    import tempfile
    import urllib.request

    url = f"https://codeload.github.com/{repo}/tar.gz/refs/heads/{branch}"
    print(f"Downloading {repo}@{branch} ...")
    try:
        with urllib.request.urlopen(url, timeout=180) as resp:
            raw = resp.read()
    except Exception as exc:  # network / HTTP errors
        sys.exit(f"ERROR: failed to download barcode data from {url}\n"
                 f"       ({exc})")
    print(f"  downloaded {len(raw) / 1e6:.1f} MB; extracting pathogen files ...")

    data_dir = os.path.abspath(data_dir)
    parent = os.path.dirname(data_dir) or "."
    os.makedirs(parent, exist_ok=True)
    staging = tempfile.mkdtemp(prefix=".sewage_data_", dir=parent)

    pathogens, versions, n_files = set(), set(), 0
    try:
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:gz") as tar:
            for member in tar:
                if not member.isfile():
                    continue
                # Archive layout: "<repo>-<branch>/<PATHOGEN>/<version>/<file>".
                parts = member.name.split("/")
                if len(parts) != 4:
                    continue  # skip top-level and pathogen-level files
                _root, pathogen, version, fname = parts
                if fname not in KEEP_FILES:
                    continue
                src = tar.extractfile(member)
                if src is None:
                    continue
                dest_dir = os.path.join(staging, pathogen, version)
                os.makedirs(dest_dir, exist_ok=True)
                with open(os.path.join(dest_dir, fname), "wb") as out:
                    shutil.copyfileobj(src, out)
                n_files += 1
                pathogens.add(pathogen)
                versions.add((pathogen, version))
    except Exception as exc:
        shutil.rmtree(staging, ignore_errors=True)
        sys.exit(f"ERROR: failed to extract barcode archive ({exc})")

    if n_files == 0:
        shutil.rmtree(staging, ignore_errors=True)
        sys.exit("ERROR: no pathogen barcode files found in the archive; "
                 "nothing was updated.")

    # Atomic swap: move the current data aside, move staging into place.
    backup = data_dir + ".old"
    if os.path.isdir(data_dir):
        shutil.rmtree(backup, ignore_errors=True)
        os.rename(data_dir, backup)
    os.rename(staging, data_dir)
    shutil.rmtree(backup, ignore_errors=True)

    return len(pathogens), len(versions), n_files


def resolve_pathogen_dir(repo, pathogen, version):
    """Resolve and validate the data directory for a pathogen + version."""
    if not os.path.isdir(repo):
        sys.exit(f"ERROR: barcode data not found: {repo}\n"
                 f"       Run 'sewage --update' to download it, or pass --repo "
                 f"/ set FREYJA_BARCODES.")
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
QMAX = 41  # maximum Phred quality emitted (Illumina 1.9-style cap)


def phred_char_array(qs):
    """Convert an int quality array to a Phred+33 byte array."""
    return (np.clip(qs, 0, QMAX).astype(np.uint8) + 33)


def make_quality(m, rl, profile, qparams, base_q, rng):
    """Return an (m, rl) int64 array of per-base Phred qualities.

    ``profile`` == "flat" gives a constant quality of ``base_q`` for every base
    (the original behavior). ``profile`` == "illumina" draws position-dependent
    qualities: the mean declines and the spread (std. dev.) widens from the 5'
    end toward the 3' end of the read, mimicking real Illumina data.
    ``qparams`` is (q_start_mean, q_end_mean, sd_start, sd_end).
    """
    if profile != "illumina":
        return np.full((m, rl), base_q, dtype=np.int64)
    q_start, q_end, sd_start, sd_end = qparams
    frac = np.linspace(0.0, 1.0, rl) if rl > 1 else np.zeros(1)
    mean_q = q_start + (q_end - q_start) * frac        # (rl,) declines to 3'
    sd_q = sd_start + (sd_end - sd_start) * frac       # (rl,) widens to 3'
    q = mean_q + sd_q * rng.standard_normal((m, rl))   # broadcast over rows
    return np.clip(np.round(q), 2, QMAX).astype(np.int64)


def simulate_reads(genome, n_pairs, read_len, frag_mean, frag_sd, error_rate,
                   rng, name_prefix, r1_out, r2_out, chunk=100000,
                   cov_diff=None, rlen_counts=None, quality_profile="flat",
                   quality_params=None, qual_hist=None):
    """Simulate ``n_pairs`` paired-end reads from one genome, writing FASTQ.

    Reads tile the whole genome uniformly. R1 is the forward strand at the
    fragment 5' end; R2 is the reverse-complement at the fragment 3' end.

    ``quality_profile`` selects the base-quality model ("flat" or "illumina",
    see :func:`make_quality`). If ``cov_diff`` (length ``glen+1`` int64 diff
    array), ``rlen_counts`` (read-length -> count dict), and/or ``qual_hist``
    (an (rl, QMAX+1) int64 array of per-position quality counts) are supplied,
    the corresponding QC statistics are accumulated for plotting.
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

        # Accumulate QC stats before mutating arrays (cheap, O(reads) via a
        # coverage difference array rather than O(reads*read_len)).
        if cov_diff is not None:
            np.add.at(cov_diff, start, 1)
            np.add.at(cov_diff, start + rl, -1)
            np.add.at(cov_diff, r2_region_start, 1)
            np.add.at(cov_diff, r2_region_start + rl, -1)
        if rlen_counts is not None:
            rlen_counts[rl] = rlen_counts.get(rl, 0) + 2 * m  # both mates

        r1 = genome[r1_idx]                            # int-encoded bases
        r2_fwd = genome[r2_idx]
        # Reverse-complement the R2 region: reverse along read, complement 3-b.
        r2 = 3 - r2_fwd[:, ::-1]

        # Per-base qualities (position-dependent for the illumina profile).
        q1 = make_quality(m, rl, quality_profile, quality_params, base_q, rng)
        q2 = make_quality(m, rl, quality_profile, quality_params, base_q, rng)

        # Introduce sequencing errors. For the illumina profile the per-base
        # error probability is derived from the drawn quality (low-quality 3'
        # bases are more error-prone); otherwise a flat --error-rate is used.
        if quality_profile == "illumina":
            for arr, q in ((r1, q1), (r2, q2)):
                mask = rng.random(arr.shape) < 10.0 ** (-q / 10.0)
                if mask.any():
                    shift = rng.integers(1, 4, size=int(mask.sum()))
                    arr[mask] = (arr[mask] + shift) % 4
        elif error_rate > 0:
            for arr, q in ((r1, q1), (r2, q2)):
                mask = rng.random(arr.shape) < error_rate
                if mask.any():
                    # Replace with a different random base (shift by 1..3 mod 4).
                    shift = rng.integers(1, 4, size=int(mask.sum()))
                    arr[mask] = (arr[mask] + shift) % 4
                    q[mask] = np.maximum(2, base_q // 2)

        # Accumulate per-position quality counts for the FastQC-style boxplot.
        if qual_hist is not None:
            w = qual_hist.shape[1]
            cols = np.broadcast_to(np.arange(rl), (m, rl))
            for q in (q1, q2):
                idx = cols.ravel() * w + np.clip(q, 0, w - 1).ravel()
                qual_hist += np.bincount(idx, minlength=rl * w).reshape(rl, w)

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
# QC plots
# --------------------------------------------------------------------------- #
def write_qc_plots(qc_dir, prefix, coverage, rlen_counts, ref_name,
                   qual_hist=None):
    """Write QC PNGs (read-length histogram, per-base coverage, and — when
    ``qual_hist`` is supplied — a FastQC-style per-base quality boxplot) into
    ``qc_dir``.

    ``coverage`` is a length-``glen`` int array of per-base depth summed over
    all lineages; ``rlen_counts`` maps read length -> number of reads;
    ``qual_hist`` is an (read_len, QMAX+1) int64 array of per-position quality
    counts. Returns the list of written file paths.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")  # headless / no display required
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - only hit without matplotlib
        sys.exit("ERROR: --qc-plots requires matplotlib. Install it with "
                 f"'pip install matplotlib'.\n       ({exc})")

    os.makedirs(qc_dir, exist_ok=True)
    written = []

    # 1) Read-length distribution histogram.
    lengths = np.array(sorted(rlen_counts), dtype=np.int64)
    counts = np.array([rlen_counts[length] for length in lengths], dtype=np.int64)
    rl_path = os.path.join(qc_dir, f"{prefix}.read_length_hist.png")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(lengths, counts, width=max(1, 0.8), color="#2b6cb0", edgecolor="black")
    ax.set_xlabel("Read length (bp)")
    ax.set_ylabel("Number of reads")
    ax.set_title(f"Read length distribution — {prefix}")
    if lengths.size:
        pad = max(1, int(0.5 * (lengths.max() - lengths.min() + 1)))
        ax.set_xlim(lengths.min() - pad, lengths.max() + pad)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(rl_path, dpi=150)
    plt.close(fig)
    written.append(rl_path)

    # 2) Per-base depth of coverage across the genome (line only, no markers).
    cov_path = os.path.join(qc_dir, f"{prefix}.coverage.png")
    positions = np.arange(1, coverage.size + 1)
    mean_cov = float(coverage.mean()) if coverage.size else 0.0
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(positions, coverage, linewidth=0.8, color="#2f855a")
    ax.axhline(mean_cov, color="#c53030", linewidth=1.0, linestyle="--",
               label=f"mean = {mean_cov:.1f}x")
    ax.set_xlabel("Genome position (bp)")
    ax.set_ylabel("Depth of coverage")
    ax.set_title(f"Depth of coverage across {ref_name} — {prefix}")
    ax.set_xlim(1, coverage.size)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(cov_path, dpi=150)
    plt.close(fig)
    written.append(cov_path)

    # 3) FastQC-style per-base sequence quality box-and-whisker plot.
    if qual_hist is not None and qual_hist.sum() > 0:
        qb_path = os.path.join(qc_dir, f"{prefix}.per_base_quality.png")
        rlp, nq = qual_hist.shape
        qvals = np.arange(nq)
        stats = []
        for pos in range(rlp):
            col = qual_hist[pos]
            total = int(col.sum())
            if total == 0:
                stats.append(dict(label="", med=0, q1=0, q3=0,
                                  whislo=0, whishi=0, fliers=[]))
                continue
            cdf = np.cumsum(col)

            def pct(p):
                idx = int(np.searchsorted(cdf, p / 100.0 * total, side="left"))
                return int(qvals[min(idx, nq - 1)])

            stats.append(dict(label="", med=pct(50), q1=pct(25), q3=pct(75),
                              whislo=pct(10), whishi=pct(90), fliers=[]))

        fig, ax = plt.subplots(figsize=(min(24, max(10, rlp * 0.06)), 6))
        # FastQC-style background quality bands (poor / OK / good).
        ax.axhspan(0, 20, facecolor="#f2d5d5", zorder=0)   # red
        ax.axhspan(20, 28, facecolor="#f2eccf", zorder=0)  # orange
        ax.axhspan(28, nq, facecolor="#d5f2d5", zorder=0)  # green
        ax.bxp(stats, showfliers=False, patch_artist=True, widths=0.7,
               boxprops=dict(facecolor="#fff29a", edgecolor="black",
                             linewidth=0.5, zorder=3),
               medianprops=dict(color="#c53030", linewidth=1.0, zorder=4),
               whiskerprops=dict(color="black", linewidth=0.5, zorder=3),
               capprops=dict(color="black", linewidth=0.5, zorder=3))
        step = max(1, rlp // 25)
        ticks = list(range(1, rlp + 1, step))
        ax.set_xticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks])
        ax.set_xlim(0.5, rlp + 0.5)
        ax.set_ylim(0, nq)
        ax.set_xlabel("Position in read (bp)")
        ax.set_ylabel("Phred quality score")
        ax.set_title(f"Per-base sequence quality — {prefix}")
        ax.grid(axis="y", alpha=0.3, zorder=1)
        fig.tight_layout()
        fig.savefig(qb_path, dpi=150)
        plt.close(fig)
        written.append(qb_path)

    return written


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
                    "FASTQ for a mixture of viral lineages, using bundled\n"
                    "references/barcodes (refresh with --update).",
        formatter_class=SewageHelpFormatter,
    )
    # --- Reference / pathogen selection ------------------------------------
    g_ref = p.add_argument_group(
        "reference / pathogen selection",
        "Choose the pathogen and the Freyja-barcodes source to build from.")
    g_ref.add_argument("-p", "--pathogen",
                   help="Pathogen name as it appears in the Freyja-barcodes "
                        "repo (e.g. DENV4, MPX, RSVa, MEASLESgenome).")
    g_ref.add_argument("--repo", default=DEFAULT_REPO,
                   help="Path to the bundled barcode data directory. Defaults "
                        "to the data/ folder shipped with SEWAGE (populated by "
                        "--update); override or set FREYJA_BARCODES to use a "
                        "different location.")
    g_ref.add_argument("--version", default="latest",
                   help="Barcode version subfolder to use (e.g. latest or a "
                        "date like 2025-05-01).")
    g_ref.add_argument("--list", action="store_true",
                   help="List available pathogens in the repo and exit.")
    g_ref.add_argument("--list-lineages", action="store_true",
                   help="List available lineages for --pathogen and exit.")

    # --- Data management ---------------------------------------------------
    g_data = p.add_argument_group(
        "data management",
        "Manage the bundled Freyja-barcodes data (references + barcodes).")
    g_data.add_argument("--update", action="store_true",
                   help="Download the latest upstream barcode repo and rebuild "
                        "the local bundled data (all pathogens, latest + dated "
                        "versions), then exit. Only reference.fasta and "
                        "barcode.csv files are kept; new pathogens are added.")
    g_data.add_argument("--update-repo", default=UPSTREAM_REPO,
                   help="GitHub 'owner/name' to pull barcode data from.")
    g_data.add_argument("--update-branch", default=UPSTREAM_BRANCH,
                   help="Branch of the barcode repo to pull for --update.")

    # --- Lineage proportions -----------------------------------------------
    g_prop = p.add_argument_group(
        "lineage proportions",
        "Supply a proportions table or auto-generate one.")
    g_prop.add_argument("-i", "--proportions",
                   help="CSV/TSV file: column 1 = lineage, column 2 = "
                        "proportion. Lineages must exist in the barcode file.")
    g_prop.add_argument("--generate-proportions", action="store_true",
                   help="Auto-generate the proportions table instead of "
                        "supplying one. Prompts if --num-lineages/--prop-mode "
                        "are not given.")
    g_prop.add_argument("--num-lineages", type=int,
                   help="(generate) Number of lineages to include.")
    g_prop.add_argument("--prop-mode", choices=["equal", "dominant", "random"],
                   help="(generate) Proportion scheme.")
    g_prop.add_argument("--dominant-fraction", type=float,
                   help="(generate, dominant mode) Fraction for the dominant "
                        "lineage (default: random in [0.5, 0.8]).")
    g_prop.add_argument("--proportions-out",
                   help="Where to save the (generated or used) proportions "
                        "table. Default: <output_prefix>.proportions.tsv")

    # --- Sequencing depth --------------------------------------------------
    g_depth = p.add_argument_group(
        "sequencing depth",
        "Specify exactly one of --depth or --num-pairs.")
    g_depth.add_argument("--depth", type=float,
                   help="Target fold coverage depth for the whole sample. Read "
                        "pairs are derived from genome length and read length.")
    g_depth.add_argument("--num-pairs", type=int,
                   help="Total number of read PAIRS for the whole sample.")

    # --- Read geometry & errors --------------------------------------------
    g_read = p.add_argument_group(
        "read geometry & errors",
        "Read/fragment dimensions and the sequencing error rate.")
    g_read.add_argument("--read-length", type=int, default=250,
                   help="Length of each mate.")
    g_read.add_argument("--fragment-mean", type=float, default=500.0,
                   help="Mean fragment (insert) size.")
    g_read.add_argument("--fragment-sd", type=float, default=50.0,
                   help="Std. dev. of fragment size.")
    g_read.add_argument("--error-rate", type=float, default=0.005,
                   help="Per-base sequencing error rate (0 to disable). Used by "
                        "the 'flat' quality profile; ignored by 'illumina', "
                        "which derives errors from the per-base quality.")

    # --- Base-quality model ------------------------------------------------
    g_qual = p.add_argument_group(
        "base-quality model",
        "How per-base Phred quality scores are generated.")
    g_qual.add_argument("--quality-profile", choices=["flat", "illumina"],
                   default="flat",
                   help="Base-quality model. 'flat' (default) emits a constant "
                        "quality derived from --error-rate. 'illumina' emits "
                        "position-dependent qualities whose mean declines and "
                        "spread widens toward the 3' end of each read.")
    g_qual.add_argument("--quality-start", type=float, default=38.0,
                   help="(illumina) Mean Phred quality at the 5' end of reads.")
    g_qual.add_argument("--quality-end", type=float, default=30.0,
                   help="(illumina) Mean Phred quality at the 3' end of reads.")
    g_qual.add_argument("--quality-sd-start", type=float, default=1.0,
                   help="(illumina) Quality std. dev. at the 5' end (narrow).")
    g_qual.add_argument("--quality-sd-end", type=float, default=8.0,
                   help="(illumina) Quality std. dev. at the 3' end (wide).")

    # --- Output ------------------------------------------------------------
    g_out = p.add_argument_group(
        "output",
        "FASTQ output location, compression, and reproducibility.")
    g_out.add_argument("-o", "--output-prefix", default="sim_sample",
                   help="Prefix for output FASTQ files "
                        "(<prefix>_R1.fastq.gz / <prefix>_R2.fastq.gz).")
    g_out.add_argument("--gzip", dest="gzip", action="store_true", default=True,
                   help="Gzip the FASTQ output (default).")
    g_out.add_argument("--no-gzip", dest="gzip", action="store_false",
                   help="Write plain (uncompressed) FASTQ.")
    g_out.add_argument("--gzip-level", type=int, default=6,
                   help="gzip compression level 1-9. Lower is much faster with "
                        "slightly larger files; 9 (Python's default) is the "
                        "slowest and a common bottleneck for big outputs.")
    g_out.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducibility.")

    # --- QC & diagnostics --------------------------------------------------
    g_qc = p.add_argument_group(
        "QC & diagnostics",
        "Optional quality-control plots and timing output.")
    g_qc.add_argument("--qc-plots", action="store_true",
                   help="Generate QC PNGs (read-length histogram, per-base "
                        "depth-of-coverage line plot, and a FastQC-style "
                        "per-base quality boxplot) into a folder. Requires "
                        "matplotlib.")
    g_qc.add_argument("--qc-dir",
                   help="Folder for the QC plots when --qc-plots is set. "
                        "Default: <output_prefix>_qc")
    g_qc.add_argument("--timing", action="store_true",
                   help="Print elapsed wall-time per phase (genome build, "
                        "read simulation) to expose bottlenecks.")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    rng = np.random.default_rng(args.seed)

    # --- Informational modes ------------------------------------------------
    if args.update:
        n_p, n_v, n_f = update_barcodes(args.repo, args.update_repo,
                                        args.update_branch)
        print(f"Updated bundled barcode data: {args.repo}")
        print(f"  {n_p} pathogens, {n_v} pathogen/versions, {n_f} files")
        return 0

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
    if args.quality_profile == "illumina":
        print(f"Quality  : illumina  mean {args.quality_start}->"
              f"{args.quality_end}  sd {args.quality_sd_start}->"
              f"{args.quality_sd_end}")
    else:
        print("Quality  : flat")
    print(f"Lineages : {len(pairs)}")

    # QC accumulators (only allocated when --qc-plots is requested). A length
    # glen+1 difference array gives per-base coverage in O(reads) via cumsum;
    # qual_hist holds per-position quality counts for the FastQC-style boxplot.
    cov_diff = np.zeros(glen + 1, dtype=np.int64) if args.qc_plots else None
    rlen_counts = {} if args.qc_plots else None
    quality_params = (args.quality_start, args.quality_end,
                      args.quality_sd_start, args.quality_sd_end)
    if args.qc_plots:
        rl_qc = min(args.read_length, glen)
        qual_hist = np.zeros((rl_qc, QMAX + 1), dtype=np.int64)
    else:
        qual_hist = None

    total_written = 0
    t_sim = time.perf_counter()
    with opener(r1_path, "wb") as r1, opener(r2_path, "wb") as r2:
        for name, prop in pairs:
            n = counts[name]
            prefix = f"{args.pathogen}:{name}"
            t_lin = time.perf_counter()
            w = simulate_reads(
                genomes[name], n, args.read_length, args.fragment_mean,
                args.fragment_sd, args.error_rate, rng, prefix, r1, r2,
                cov_diff=cov_diff, rlen_counts=rlen_counts,
                quality_profile=args.quality_profile,
                quality_params=quality_params, qual_hist=qual_hist)
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

    # --- QC plots -----------------------------------------------------------
    if args.qc_plots:
        qc_dir = args.qc_dir or f"{args.output_prefix}_qc"
        coverage = np.cumsum(cov_diff)[:glen]  # diff array -> per-base depth
        qc_prefix = os.path.basename(args.output_prefix)
        qc_files = write_qc_plots(qc_dir, qc_prefix, coverage, rlen_counts,
                                  ref_name, qual_hist=qual_hist)
        print(f"  QC plots ({qc_dir}):")
        for f in qc_files:
            print(f"    {f}")

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
