"""Sequence scanning for the prespecified RNA-binding motifs."""

from __future__ import annotations

import itertools
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


def prepare_matrix(path: Path) -> np.ndarray:
    """Apply the DEWSeq-style rounding and low-information flank trimming."""
    table = pd.read_csv(path, sep="\t")
    pyrimidine = "U" if "U" in table else "T"
    columns = ["A", "C", "G", pyrimidine]
    if not set(columns).issubset(table):
        raise ValueError("motif matrix needs A, C, G, and U or T columns")
    matrix = table[columns].to_numpy(dtype=float)
    if matrix.ndim != 2 or not len(matrix) or not np.isfinite(matrix).all():
        raise ValueError("motif matrix must contain finite base probabilities")
    if np.any(matrix < 0) or not np.allclose(matrix.sum(axis=1), 1, atol=1e-5):
        raise ValueError("each motif row must contain nonnegative probabilities summing to one")
    matrix[matrix <= 0.025] = 0
    matrix /= matrix.sum(axis=1, keepdims=True)
    nonzero = np.where(matrix > 0, matrix * np.log2(np.maximum(matrix, 1e-300)), 0)
    information = 2 + nonzero.sum(axis=1)
    left, right = 0, len(matrix)
    while left < right and information[left] <= 0.1:
        left += 1
    while right > left and information[right - 1] <= 0.1:
        right -= 1
    matrix = matrix[left:right]
    if len(matrix) < 6:
        raise ValueError("trimmed motif must be at least six bases long")
    return matrix


def training_background(prepared_data: Path) -> np.ndarray:
    """Estimate zero-order A/C/G/T frequencies from training windows only."""
    with np.load(prepared_data, allow_pickle=False) as data:
        selected = np.asarray(data["split"]).astype(str) == "train"
        sequence = np.asarray(data["sequence"])
        if sequence.ndim != 3 or sequence.shape[1] != 4 or not selected.any():
            raise ValueError("training data needs four sequence channels and train windows")
        counts = sequence[selected].sum(axis=(0, 2), dtype=float) + 0.1
    return counts / counts.sum()


def _meme_payload(motif_id: str, matrix: np.ndarray,
                  background: np.ndarray) -> bytes:
    rows = [
        "MEME version 4", "", "ALPHABET= ACGT", "", "strands: + -", "",
        "Background letter frequencies",
        " ".join(f"{base} {value:.12g}" for base, value in zip("ACGT", background)),
        "", f"MOTIF {motif_id} {motif_id}",
        f"letter-probability matrix: alength= 4 w= {len(matrix)} nsites= 20 E= 0",
    ]
    rows.extend(" ".join(f"{value:.12g}" for value in row) for row in matrix)
    return ("\n".join(rows) + "\n").encode()


def qualifying_words(motif_id: str, matrix: np.ndarray,
                     background_frequencies: np.ndarray,
                     threshold: float = 0.001) -> dict[str, float]:
    """Enumerate DNA words accepted by the one-strand FIMO score threshold."""
    from pymemesuite.common import Array, Background, MotifFile, Sequence
    from pymemesuite.fimo import FIMO

    if len(matrix) > 10:
        raise ValueError("word enumeration is limited to motifs of at most 10 bases")
    words = ["".join(word) for word in itertools.product("ACGT", repeat=len(matrix))]
    with MotifFile(BytesIO(_meme_payload(motif_id, matrix, background_frequencies)),
                   pseudocount=0.1) as source:
        motif = source.read()
        background = Background(source.alphabet, Array(background_frequencies.tolist()))
        # N padding allows a hit at either boundary of a sequence of motif width.
        probes = [Sequence(f"N{word}N", name=str(i).encode())
                  for i, word in enumerate(words)]
        matches = FIMO(threshold=threshold, both_strands=False,
                       max_stored_scores=max(100_000, 2 * len(words))).score_motif(
                           motif, probes, background)
        pvalues: dict[str, float] = {}
        for match in matches.matched_elements:
            word = str(match.sequence)
            pvalues[word] = min(pvalues.get(word, 1.0), float(match.pvalue))
    return pvalues


def scan_intervals(table: pd.DataFrame, reference: Path, matrix_path: Path,
                   training_data: Path, motif_id: str,
                   threshold: float = 0.001) -> pd.DataFrame:
    """Count motif hits in transcript-oriented genomic intervals."""
    required = {"chrom", "start", "end", "strand"}
    if not required.issubset(table):
        raise ValueError(f"interval table needs columns {sorted(required)}")
    matrix = prepare_matrix(matrix_path)
    words = qualifying_words(motif_id, matrix, training_background(training_data), threshold)
    width = len(matrix)
    out = table.copy()
    hit_counts = []
    lengths = []
    complement = str.maketrans("ACGTN", "TGCAN")
    with pysam.FastaFile(str(reference)) as genome:
        for row in out.itertuples(index=False):
            start, end = int(row.start), int(row.end)
            if start < 0 or end <= start or row.strand not in {"+", "-"}:
                raise ValueError("interval coordinates or strand are invalid")
            sequence = genome.fetch(str(row.chrom), start, end).upper()
            if row.strand == "-":
                sequence = sequence.translate(complement)[::-1]
            hits = sum(sequence[i:i + width] in words
                       for i in range(max(0, len(sequence) - width + 1)))
            hit_counts.append(hits)
            lengths.append(len(sequence))
    out["motif_hits"] = hit_counts
    out["motif_present"] = (np.asarray(hit_counts) > 0).astype(int)
    out["scan_nt"] = lengths
    return out


def motif_block_counts(table: pd.DataFrame,
                       seeds: tuple[int, ...] = (42, 1337, 2026)) -> pd.DataFrame:
    """Aggregate scanned original and seed-specific model intervals by 1 Mb block."""
    required = {"rbp", "condition", "chrom", "start", "approach", "seed",
                "motif_hits", "scan_nt"}
    if not required.issubset(table):
        raise ValueError(f"scanned interval table needs {sorted(required)}")
    if not table["approach"].isin(["clipper", "model"]).all():
        raise ValueError("approach must be clipper or model")
    if not np.isfinite(table[["motif_hits", "scan_nt"]].to_numpy(dtype=float)).all():
        raise ValueError("motif counts and scanned lengths must be finite")
    if (table[["motif_hits", "scan_nt"]] < 0).any().any():
        raise ValueError("motif counts and scanned lengths must be nonnegative")
    frame = table.copy()
    frame["block"] = frame["chrom"].astype(str) + ":" + (
        frame["start"].astype(int) // 1_000_000).astype(str)
    rows = []
    for (rbp, condition, block), section in frame.groupby(
            ["rbp", "condition", "block"], sort=False):
        clipper = section[section["approach"] == "clipper"]
        if not clipper.empty and clipper["seed"].nunique() != 1:
            raise ValueError("CLIPper intervals must appear once, with one seed label")
        row = {"rbp": rbp, "condition": condition, "block": block,
               "clipper_hits": int(clipper["motif_hits"].sum()),
               "clipper_nt": int(clipper["scan_nt"].sum())}
        for seed in seeds:
            model = section[(section["approach"] == "model") &
                            (section["seed"] == seed)]
            row[f"model_hits_{seed}"] = int(model["motif_hits"].sum())
            row[f"model_nt_{seed}"] = int(model["scan_nt"].sum())
        rows.append(row)
    return pd.DataFrame(rows)
