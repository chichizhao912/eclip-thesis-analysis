from __future__ import annotations

from collections import Counter

import numpy as np
import pytest


def test_dinucleotide_shuffle_retains_pairs():
    pytest.importorskip("eclip_reconstruct")
    from eclip_thesis_analysis.controls import dinucleotide_shuffle

    source = "ACGTTGCAACGT" * 25
    shuffled = dinucleotide_shuffle(source, np.random.default_rng(42))
    assert len(shuffled) == len(source)
    assert Counter(zip(shuffled[:-1], shuffled[1:])) == Counter(
        zip(source[:-1], source[1:]))


def test_all_input_controls_write_standard_predictions(tmp_path):
    pytest.importorskip("eclip_reconstruct")
    from eclip_thesis_analysis.controls import run_controls

    dataset = tmp_path / "windows.npz"
    n = 6
    sequence = np.zeros((n, 4, 300), dtype=np.float32)
    sequence[:, np.arange(300) % 4, np.arange(300)] = 1
    baseline = np.linspace(0.1, 1.0, 300, dtype=np.float32)
    profiles = np.stack([baseline + i * 0.01 for i in range(n)])
    np.savez_compressed(dataset, sequence=sequence,
                        ip1=profiles, ip2=profiles * 1.02,
                        sminput=profiles * 0.2,
                        split=np.array(["train"] * 4 + ["validation"] * 2),
                        chrom=np.array(["chr1"] * 4 + ["chr2"] * 2),
                        start=np.arange(n) * 300, end=np.arange(n) * 300 + 300,
                        strand=np.array(["+"] * n),
                        window_id=np.array([f"w{i}" for i in range(n)]),
                        window_class=np.array(["peak"] * n))
    output = tmp_path / "controls"
    run_controls(dataset, output, seeds=(42,), maximum_epochs=1)
    for configuration in ("coverage", "full", "shuffled", "no_sminput"):
        with np.load(output / f"{configuration}_seed42_predictions.npz") as prediction:
            assert prediction["p1"].shape == (n, 300)
            assert np.isfinite(prediction["p1"]).all()
