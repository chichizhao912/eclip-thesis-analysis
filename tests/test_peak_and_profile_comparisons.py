from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pyBigWig

from eclip_thesis_analysis.analyses import coverage_presence
from eclip_thesis_analysis.caller_comparison import caller_score_comparison
from eclip_thesis_analysis.candidate_lengths import candidate_lengths
from eclip_thesis_analysis.peak_diagnostics import peak_diagnostics, peak_shape_tests
from eclip_thesis_analysis.profile_comparison import (
    dominant_fwhm, gaussian_comparison_rows, gaussian_differences, gaussian_summary,
)


def _bigwig(path, values):
    with pyBigWig.open(str(path), "w") as handle:
        handle.addHeader([("chr1", len(values))])
        handle.addEntries("chr1", 0, values=values.tolist(), span=1, step=1)
    return path


def test_original_peak_overlap_shape_and_partner_support(tmp_path):
    first = tmp_path / "ip1.bed"
    second = tmp_path / "ip2.bed"
    first.write_text("chr1\t170\t180\ta\t0\t+\nchr1\t200\t210\tb\t0\t+\n")
    second.write_text("chr1\t172\t182\tc\t0\t+\nchr1\t230\t240\td\t0\t+\n")
    ip1 = np.zeros(500)
    ip2 = np.zeros(500)
    ip1[170:180], ip1[200:210] = 2, 3
    ip2[172:182], ip2[230:240] = 2, 3
    input_track = np.full(500, 0.001)
    tracks = {}
    for name, values in (("ip1", ip1), ("ip2", ip2), ("sminput", input_track)):
        tracks[name] = (_bigwig(tmp_path / f"{name}.plus.bw", values),
                        _bigwig(tmp_path / f"{name}.minus.bw", np.zeros(500)))
    result = peak_diagnostics("TEST", first, second, tracks, resamples=30)
    summary = result["summary"]
    assert summary.n_shared.tolist() == [1, 1]
    assert np.allclose(summary.base_pair_jaccard, 0.25)
    assert summary.median_summit_displacement_nt.tolist() == [2.0, 2.0]
    shape = result["shape"]
    assert len(shape) == 4
    assert shape.half_height_span_nt.tolist() == [10, 10, 10, 10]
    support = result["partner_support"]
    assert len(support) == 4
    assert set(support[support.interval == "peak"].support) == {"none"}
    assert set(support[support.interval == "window"].support) == {"above_sminput"}
    tests = peak_shape_tests(shape, resamples=30)
    assert len(tests) == 6
    assert set(tests.metric) == {"width_nt", "half_height_span_nt", "concentration"}


def test_gaussian_comparison_keeps_seeds_within_windows(tmp_path):
    x = np.arange(300)
    profiles = np.stack([np.exp(-0.5 * ((x - centre) / 10) ** 2)
                         for centre in (120, 140)])
    data = tmp_path / "data.npz"
    predictions = tmp_path / "predictions.npz"
    preprocessing = tmp_path / "preprocessing.json"
    np.savez(data, window_id=np.array(["a", "b"]), chrom=np.array(["chr1", "chr1"]),
             start=np.array([0, 300]), ip1=profiles, ip2=profiles * 0.9)
    transformed = np.log1p(profiles)
    np.savez(predictions, window_id=np.array(["a", "b"]),
             p1=transformed, p2=transformed * 0.9)
    preprocessing.write_text(json.dumps({"clip_value": 10}))
    presence = coverage_presence(data, "TEST")
    assert presence.both.tolist() == [2]
    assert presence.both_percent.tolist() == [100.0]
    first = gaussian_comparison_rows(data, predictions, preprocessing, 42, sigma=3)
    assert len(first) == 2
    assert np.isfinite(first.autoencoder_fwhm_nt).all()
    assert np.isclose(dominant_fwhm(np.array([0., 1., 2., 1., 0.])), 2)
    second = first.copy()
    second["seed"] = 1337
    broad = np.exp(-0.5 * ((x - 120) / 20) ** 2)
    second.loc[second.window_id == "a", "autoencoder_consensus_profile"] = json.dumps(broad.tolist())
    combined = pd.concat([first, second], ignore_index=True)
    summary = gaussian_summary(combined)
    assert summary.n_windows.tolist() == [2, 2]
    assert summary.n_width_windows.tolist() == [2, 2]
    expected = dominant_fwhm((np.asarray(json.loads(first.iloc[0].autoencoder_consensus_profile))
                              + broad) / 2)
    widths = gaussian_differences(combined, "width")
    actual = widths[(widths.window_id == "a") &
                    (widths.comparison == "autoencoder_minus_observed")].difference.iloc[0]
    assert np.isclose(actual, expected - first.iloc[0].observed_fwhm_nt)
    assert len(gaussian_differences(combined, "correlation")) == 4
    assert len(widths) == 4


def test_caller_scores_get_separate_background_thresholds(tmp_path):
    path = tmp_path / "calls.npz"
    levels = np.array([0.8, 0.9, 1.0, 1.1, 1.4, 1.5])
    consensus = np.repeat(levels[:, None], 300, axis=1)
    sminput = np.repeat(np.array([0.0, 0.0, 0.0, 0.0, 0.8, 0.8])[:, None], 300, axis=1)
    classes = np.array(["hard_background"] * 4 + ["peak_derived"] * 2)
    np.savez(path, window_id=np.array([str(i) for i in range(6)]),
             chrom=np.array(["chr1"] * 6), start=np.arange(6) * 300,
             end=(np.arange(6) + 1) * 300, strand=np.array(["+"] * 6),
             split=np.array(["validation"] * 6), window_class=classes,
             consensus=consensus, sminput=sminput)
    result = caller_score_comparison(path, 42)
    assert set(result.score_scheme) == {"consensus_only", "signal_minus_input"}
    assert result.hard_background_called_percent.tolist() == [25.0, 25.0]
    assert result.peak_recovery_percent.tolist() == [100.0, 100.0]


def test_candidate_lengths_keep_out_of_range_regions(tmp_path):
    original = tmp_path / "original.bed"
    refined = tmp_path / "refined.bed"
    original.write_text("chr1\t0\t235\ta\t0\t+\nchr1\t0\t301\tb\t0\t+\n")
    refined.write_text("chr1\t0\t52\ta\t0\t+\n")
    lengths, summary = candidate_lengths(original, refined)
    assert lengths.length_nt.tolist() == [235, 301, 52]
    assert summary.n_over_300_nt.tolist() == [1, 0]
