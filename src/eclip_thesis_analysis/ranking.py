"""Exploratory ranking of windows outside original CLIPper calls."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


def _merged_peaks(paths: list[Path]) -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    intervals = defaultdict(list)
    for path in paths:
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or fields[5] not in {"+", "-"}:
                    raise ValueError(f"invalid stranded BED record in {path}")
                intervals[(fields[0], fields[5])].append((int(fields[1]), int(fields[2])))
    result = {}
    for key, values in intervals.items():
        merged = []
        for start, end in sorted(values):
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        result[key] = (np.asarray([item[0] for item in merged]),
                       np.asarray([item[1] for item in merged]))
    return result


def rank_windows(data_path: Path, predictions_path: Path,
                 original_peak_paths: list[Path], rbp: str,
                 top_n: int = 1000) -> pd.DataFrame:
    """Rank all coverage-eligible windows, then exclude original-peak overlaps."""
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    with np.load(data_path, allow_pickle=False) as source, \
         np.load(predictions_path, allow_pickle=False) as predicted:
        if not np.array_equal(source["window_id"], predicted["window_id"]):
            raise ValueError("window IDs in data and predictions differ")
        eligible = (source["ip1"] + source["ip2"]).sum(axis=1) >= 1
        if not eligible.any():
            raise ValueError("no windows pass the measured-IP coverage filter")
        net = predicted["consensus"] - 0.25 * predicted["sminput"]
        starts = np.arange(0, 201, 5)
        prefix = np.pad(np.cumsum(net, axis=1), ((0, 0), (1, 0)))
        score = ((prefix[:, starts + 100] - prefix[:, starts]) / 100).max(axis=1)
        if not np.isfinite(score[eligible]).all():
            raise ValueError("eligible windows must have finite candidate scores")
        order = np.flatnonzero(eligible)[np.argsort(-score[eligible], kind="stable")]
        peaks = _merged_peaks(original_peak_paths)
        rows = []
        for rank, index in enumerate(order, 1):
            chrom, strand = str(source["chrom"][index]), str(source["strand"][index])
            start, end = int(source["start"][index]), int(source["end"][index])
            peak_starts, peak_ends = peaks.get((chrom, strand),
                                               (np.array([], dtype=int), np.array([], dtype=int)))
            position = int(np.searchsorted(peak_starts, end, side="left") - 1)
            overlap = position >= 0 and peak_ends[position] > start
            if overlap:
                continue
            rows.append({"rbp": rbp, "window_id": str(source["window_id"][index]),
                         "chrom": chrom, "start": start, "end": end,
                         "strand": strand, "rank": rank, "score": float(score[index]),
                         "selected": int(rank <= top_n)})
    return pd.DataFrame(rows)

