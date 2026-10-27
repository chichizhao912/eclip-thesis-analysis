from __future__ import annotations

import numpy as np
import pandas as pd

from eclip_thesis_analysis.analyses import (average_precision, jaccard_differences,
                                            reconstruction_differences,
                                            outside_peak_motifs, paired_tests)
from eclip_thesis_analysis.intervals import jaccard_block_table
from eclip_thesis_analysis.ranking import rank_windows
from eclip_thesis_analysis.statistics import bootstrap_indices, holm, sign_flip_p


def test_repeated_blocks_retain_all_contributions():
    blocks = np.array(["a", "a", "b"])
    rng = np.random.default_rng(2)
    observed = [bootstrap_indices(blocks, rng) for _ in range(30)]
    assert any(len(indices) == 4 for indices in observed)
    assert any(np.array_equal(indices, np.array([0, 1, 0, 1])) for indices in observed)


def test_holm_and_sign_flip():
    assert np.allclose(holm(np.array([0.03, 0.01, 0.04])), [0.06, 0.03, 0.06])
    estimate, p = sign_flip_p(np.ones(8), np.array([f"b{i}" for i in range(8)]),
                              "mean", 1000, 42)
    assert estimate == 1.0
    assert p < 0.05


def test_paired_family_and_ap():
    table = pd.DataFrame({
        "family": ["metric_a"] * 8 + ["metric_b"] * 8,
        "comparison": ["x"] * 8 + ["y"] * 8,
        "chrom": ["chr1"] * 16,
        "start": np.arange(16) * 1_000_000,
        "difference": [1.0] * 8 + [-1.0] * 8,
    })
    result = paired_tests(table, method="sign-flip", statistic="mean",
                          bootstrap_repetitions=200, test_repetitions=1000, seed=3)
    assert len(result) == 2
    assert np.allclose(result["estimate"], [1.0, -1.0])
    assert np.allclose(result["p_holm"], result["p_raw"])
    assert np.isclose(average_precision([True, False, True], [0.9, 0.8, 0.7]), 5 / 6)


def test_outside_peak_and_jaccard_summaries():
    rows = []
    for rbp in ("A", "B", "C"):
        for i in range(40):
            rows.append({"rbp": rbp, "chrom": "chr3", "start": i * 1_000_000,
                         "selected": i < 20, "motif_present": i < 15 or 20 <= i < 25})
    result = outside_peak_motifs(pd.DataFrame(rows), repetitions=200, seed=4)
    assert len(result) == 3
    assert np.allclose(result["difference_pp"], 50.0)
    assert result["p_holm"].between(0, 1).all()

    counts = pd.DataFrame({
        "rbp": ["A", "A", "B", "B"], "seed": [42, 42, 42, 42],
        "block": ["chr3:0", "chr3:1", "chr3:0", "chr3:1"],
        "model_intersection": [4, 4, 4, 4], "model_union": [10, 10, 10, 10],
        "original_intersection": [2, 2, 2, 2], "original_union": [10, 10, 10, 10],
    })
    values, sign_p = jaccard_differences(counts, repetitions=100, seed=4)
    assert np.allclose(values["difference"], 0.2)
    assert sign_p == 0.5


def test_jaccard_base_counts_across_block_boundary(tmp_path):
    paths = [tmp_path / name for name in ("o1.bed", "o2.bed", "m1.bed", "m2.bed")]
    lines = [
        "chr3\t999990\t1000010\ta\t0\t+\n",
        "chr3\t999995\t1000005\tb\t0\t+\n",
        "chr3\t999990\t1000010\tc\t0\t+\n",
        "chr3\t999990\t1000010\td\t0\t+\n",
    ]
    for path, line in zip(paths, lines):
        path.write_text(line)
    counts = jaccard_block_table("RBP", 42, *paths)
    assert len(counts) == 2
    assert counts["original_intersection"].sum() == 10
    assert counts["original_union"].sum() == 20
    assert counts["model_intersection"].sum() == 20


def test_reconstruction_comparisons_average_seeds_within_window():
    rows = []
    for config in ("full", "baseline"):
        for seed in (42, 1337):
            rows.append({"configuration": config, "seed": seed,
                         "window_id": "w1", "chrom": "chr2", "start": 100,
                         "reconstruction_agreement": (0.8 if config == "full" else 0.6)
                         + (0.02 if seed == 1337 else 0),
                         "target_agreement": 0.55 if config == "full" else 0.50})
    result = reconstruction_differences(pd.DataFrame(rows), "full", ["baseline"],
                                        "input_ablation")
    assert len(result) == 2
    assert np.allclose(sorted(result["difference"]), [0.05, 0.20])


def test_rank_before_excluding_original_peaks(tmp_path):
    data = tmp_path / "windows.npz"
    prediction = tmp_path / "prediction.npz"
    peak = tmp_path / "original.bed"
    peak.write_text("chr3\t0\t300\tpeak\t0\t+\n")
    np.savez_compressed(
        data, window_id=np.array(["first", "second", "third"]),
        chrom=np.array(["chr3"] * 3), start=np.array([0, 300, 600]),
        end=np.array([300, 600, 900]), strand=np.array(["+"] * 3),
        ip1=np.ones((3, 300), dtype=float), ip2=np.ones((3, 300), dtype=float),
    )
    signal = np.array([3, 2, 1], dtype=float)[:, None] * np.ones((3, 300))
    np.savez_compressed(prediction, window_id=np.array(["first", "second", "third"]),
                        consensus=signal, sminput=np.zeros((3, 300)))
    result = rank_windows(data, prediction, [peak], "RBP", top_n=2)
    assert result["window_id"].tolist() == ["second", "third"]
    assert result["rank"].tolist() == [2, 3]
    assert result["selected"].tolist() == [1, 0]

