from __future__ import annotations

import numpy as np
import pandas as pd
import pysam
import pytest

from eclip_thesis_analysis.motifs import motif_block_counts, scan_intervals


def test_motif_scanning_respects_transcript_strand(tmp_path):
    pytest.importorskip("pymemesuite")
    reference = tmp_path / "genome.fa"
    reference.write_text(">chr1\nTTTAAAAAATTTTTTTTT\n")
    pysam.faidx(str(reference))
    matrix = tmp_path / "motif.tsv"
    matrix.write_text("A\tC\tG\tU\n" + "0.997\t0.001\t0.001\t0.001\n" * 6)
    training = tmp_path / "training.npz"
    sequence = np.zeros((1, 4, 300), dtype=np.float32)
    sequence[0, np.arange(300) % 4, np.arange(300)] = 1
    np.savez_compressed(training, sequence=sequence, split=np.array(["train"]))
    intervals = pd.DataFrame({
        "chrom": ["chr1", "chr1"], "start": [0, 12], "end": [12, 18],
        "strand": ["+", "-"],
    })
    scanned = scan_intervals(intervals, reference, matrix, training, "AAAAAA")
    assert scanned["motif_hits"].tolist() == [1, 1]
    assert scanned["motif_present"].tolist() == [1, 1]
    assert scanned["scan_nt"].tolist() == [12, 6]


def test_motif_blocks_keep_seed_counts_separate():
    rows = pd.DataFrame({
        "rbp": ["TARDBP"] * 4,
        "condition": ["matched"] * 4,
        "chrom": ["chr3"] * 4,
        "start": [10, 20, 30, 40],
        "approach": ["clipper", "model", "model", "model"],
        "seed": [0, 42, 1337, 2026],
        "motif_hits": [2, 1, 3, 4],
        "scan_nt": [100, 70, 80, 90],
    })
    result = motif_block_counts(rows).iloc[0]
    assert result["clipper_hits"] == 2
    assert result["model_hits_42"] == 1
    assert result["model_hits_1337"] == 3
    assert result["model_hits_2026"] == 4
