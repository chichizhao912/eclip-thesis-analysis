"""Original CLIPper peak overlap, partner support, and retained signal shape."""

from __future__ import annotations

from contextlib import ExitStack
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig

from .intervals import _intersection, _merge
from .statistics import block_ids, bootstrap, centred_bootstrap_p, holm, percentile_interval

CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


def _peaks(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2, 5],
                        compression="infer", dtype={0: str, 5: str})
    table.columns = ["chrom", "start", "end", "strand"]
    table[["start", "end"]] = table[["start", "end"]].apply(pd.to_numeric, errors="coerce")
    valid = (table.chrom.isin(CANONICAL) & table.strand.isin(("+", "-")) &
             table.start.notna() & table.end.notna() & (table.end > table.start))
    table = table.loc[valid].copy().reset_index(drop=True)
    table[["start", "end"]] = table[["start", "end"]].astype(np.int64)
    table["peak_id"] = np.arange(len(table))
    return table


def _match(source: pd.DataFrame, partner: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    shared = np.zeros(len(source), dtype=bool)
    partner_id = np.full(len(source), -1, dtype=np.int64)
    groups = {(chrom, strand): group.sort_values(["start", "end", "peak_id"])
              for (chrom, strand), group in partner.groupby(["chrom", "strand"], sort=False)}
    for key, rows in source.groupby(["chrom", "strand"], sort=False):
        other = groups.get(key)
        if other is None:
            continue
        starts, ends = other.start.to_numpy(), other.end.to_numpy()
        ids = other.peak_id.to_numpy()
        longest = int(np.max(ends - starts))
        for row in rows.itertuples(index=False):
            lo = np.searchsorted(starts, row.start - longest, side="right")
            hi = np.searchsorted(starts, row.end, side="left")
            candidates = np.flatnonzero(ends[lo:hi] > row.start) + lo
            if not len(candidates):
                continue
            overlap = np.minimum(row.end, ends[candidates]) - np.maximum(row.start, starts[candidates])
            best = candidates[int(np.argmax(overlap))]
            shared[row.peak_id] = True
            partner_id[row.peak_id] = ids[best]
    return shared, partner_id


def _profile(pair: tuple, chrom: str, start: int, end: int, strand: str) -> np.ndarray | None:
    try:
        values = pair[0 if strand == "+" else 1].values(chrom, start, end, numpy=True)
    except (RuntimeError, ValueError, KeyError):
        return None
    if values is None or len(values) != end - start:
        return None
    values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(values if strand == "+" else -values[::-1], 0.0)


def _summits(peaks: pd.DataFrame, pair: tuple) -> np.ndarray:
    result = np.full(len(peaks), -1, dtype=np.int64)
    for row in peaks.itertuples(index=False):
        signal = _profile(pair, row.chrom, row.start, row.end, row.strand)
        if signal is None or not np.any(signal > 0):
            continue
        genomic = signal if row.strand == "+" else signal[::-1]
        result[row.peak_id] = row.start + int(np.argmax(genomic))
    return result


def _half_height_span(signal: np.ndarray) -> int:
    maximum = float(signal.max())
    centre = int(np.argmax(signal))
    left = right = centre
    while left > 0 and signal[left - 1] >= maximum / 2:
        left -= 1
    while right + 1 < len(signal) and signal[right + 1] >= maximum / 2:
        right += 1
    return right - left + 1


def _jaccard_blocks(first: pd.DataFrame, second: pd.DataFrame) -> pd.DataFrame:
    grouped = []
    for table in (first, second):
        blocks = {}
        for (chrom, strand), rows in table.groupby(["chrom", "strand"], sort=False):
            for start, end in _merge(list(rows[["start", "end"]].itertuples(index=False, name=None))):
                while start < end:
                    index = start // 1_000_000
                    stop = min(end, (index + 1) * 1_000_000)
                    blocks.setdefault((chrom, strand, index), []).append((start, stop))
                    start = stop
        grouped.append(blocks)
    rows = []
    for chrom, strand, index in grouped[0].keys() | grouped[1].keys():
        a = _merge(grouped[0].get((chrom, strand, index), []))
        b = _merge(grouped[1].get((chrom, strand, index), []))
        intersection = _intersection(a, b)
        union = sum(end - start for start, end in a + b) - intersection
        rows.append({"chrom": chrom, "block": f"{chrom}:{index}",
                     "intersection": intersection, "union": union})
    if not rows:
        return pd.DataFrame(columns=["chrom", "block", "intersection", "union"])
    return pd.DataFrame(rows).groupby(["chrom", "block"], as_index=False).sum()


def peak_diagnostics(rbp: str, ip1: Path, ip2: Path, tracks: dict[str, tuple[Path, Path]],
                     *, shape_sample: int = 5_000, seed: int = 42,
                     resamples: int = 2_000) -> dict[str, pd.DataFrame]:
    """Calculate the main original-call diagnostics for one RBP.

    Shape sampling precedes exclusion of peaks with no positive calling-IP coverage.
    """
    if shape_sample < 1:
        raise ValueError("shape_sample must be positive")
    first, second = _peaks(ip1), _peaks(ip2)
    for source, partner in ((first, second), (second, first)):
        source["shared"], source["partner_peak_id"] = _match(source, partner)
    blocks = _jaccard_blocks(first, second)
    if blocks.empty or not blocks.union.sum():
        raise ValueError("no called bases in the peak sets")
    inter, union = blocks.intersection.to_numpy(), blocks.union.to_numpy()
    observed, draws = bootstrap(lambda idx: float(inter[idx].sum() / union[idx].sum()),
                                blocks.block.to_numpy(), resamples, seed,
                                chromosome=blocks.chrom.to_numpy())
    ci_low, ci_high = percentile_interval(draws)
    summaries, shapes, support, displacements = [], [], [], []
    with ExitStack() as stack:
        handles = {}
        for name, paths in tracks.items():
            pair = []
            for path in paths:
                handle = pyBigWig.open(str(path))
                stack.callback(handle.close)
                pair.append(handle)
            handles[name] = tuple(pair)
        for name, source, partner in (("IP1", first, second), ("IP2", second, first)):
            source["summit"] = _summits(source, handles[name.lower()])
        for name, source, partner, other in (("IP1", first, second, "IP2"),
                                              ("IP2", second, first, "IP1")):
            rng = np.random.default_rng(seed)
            for row in source.loc[source.shared].itertuples(index=False):
                partner_summit = int(partner.iloc[row.partner_peak_id].summit)
                if row.summit >= 0 and partner_summit >= 0:
                    displacements.append({"rbp": rbp, "calling_ip": name, "chrom": row.chrom,
                                          "start": row.start, "peak_id": row.peak_id,
                                          "partner_peak_id": row.partner_peak_id,
                                          "displacement_nt": abs(row.summit - partner_summit)})
            for shared in (True, False):
                group = source.loc[source.shared == shared]
                if len(group) > shape_sample:
                    group = group.iloc[rng.choice(len(group), shape_sample, replace=False)]
                for row in group.itertuples(index=False):
                    signal = _profile(handles[name.lower()], row.chrom, row.start, row.end, row.strand)
                    if signal is None or not np.any(signal > 0):
                        continue
                    shapes.append({"rbp": rbp, "calling_ip": name, "chrom": row.chrom,
                                   "start": (row.start + row.end) // 2,
                                   "peak_id": row.peak_id,
                                   "peak_class": "shared" if shared else "specific",
                                   "width_nt": row.end - row.start,
                                   "half_height_span_nt": _half_height_span(signal),
                                   "concentration": float(signal.max() / signal.mean())})
            for row in source.loc[~source.shared & (source.summit >= 0)].itertuples(index=False):
                window_start, window_end = row.summit - 150, row.summit + 150
                if window_start < 0 or window_end > handles["sminput"][0].chroms().get(row.chrom, 0):
                    continue
                peak_rows = []
                for span, lo, hi in (("peak", row.start, row.end),
                                     ("window", window_start, window_end)):
                    values = [_profile(handles[key], row.chrom, lo, hi, row.strand)
                              for key in (name.lower(), other.lower(), "sminput")]
                    if any(value is None for value in values):
                        break
                    caller, partner_signal, input_signal = (float(value.sum()) for value in values)
                    category = ("none" if partner_signal == 0 else
                                "above_sminput" if partner_signal > input_signal else
                                "positive_at_or_below_sminput")
                    peak_rows.append({"rbp": rbp, "calling_ip": name, "chrom": row.chrom,
                                    "start": row.start, "peak_id": row.peak_id,
                                    "interval": span, "caller_total": caller,
                                    "partner_total": partner_signal, "sminput_total": input_signal,
                                    "support": category})
                if len(peak_rows) == 2:
                    support.extend(peak_rows)
            displacement = [x["displacement_nt"] for x in displacements
                            if x["calling_ip"] == name]
            summaries.append({"rbp": rbp, "calling_ip": name, "n_peaks": len(source),
                              "n_shared": int(source.shared.sum()),
                              "shared_fraction": float(source.shared.mean()),
                              "base_pair_jaccard": observed,
                              "jaccard_ci_low": ci_low, "jaccard_ci_high": ci_high,
                              "n_summit_pairs": len(displacement),
                              "median_summit_displacement_nt":
                              float(np.median(displacement)) if displacement else np.nan,
                              "summits_within_50_nt":
                              float(np.mean(np.asarray(displacement) <= 50)) if displacement else np.nan})
    return {"summary": pd.DataFrame(summaries), "shape": pd.DataFrame(shapes),
            "partner_support": pd.DataFrame(support),
            "summit_pairs": pd.DataFrame(displacements)}


def peak_shape_tests(table: pd.DataFrame, *, resamples: int = 2_000,
                     seed: int = 42) -> pd.DataFrame:
    """Test specific-minus-shared median differences by metric and RBP/IP."""
    needed = {"rbp", "calling_ip", "chrom", "start", "peak_class", "width_nt",
              "half_height_span_nt", "concentration"}
    if not needed.issubset(table):
        raise ValueError(f"shape table requires {sorted(needed)}")
    rows = []
    for metric in ("width_nt", "half_height_span_nt", "concentration"):
        for (rbp, calling_ip), group in table.groupby(["rbp", "calling_ip"], sort=False):
            group = group.dropna(subset=[metric])
            values = group[metric].to_numpy(float)
            specific = group.peak_class.to_numpy() == "specific"
            if not specific.any() or specific.all():
                raise ValueError(f"both peak classes are required for {rbp}/{calling_ip}")
            blocks = block_ids(group.chrom.to_numpy(), group.start.to_numpy())
            def difference(index: np.ndarray) -> float:
                mask = specific[index]
                return float(np.median(values[index][mask]) -
                             np.median(values[index][~mask])) if mask.any() and (~mask).any() else np.nan
            observed, draws = bootstrap(difference, blocks, resamples, seed,
                                        chromosome=group.chrom.to_numpy())
            draws = draws[np.isfinite(draws)]
            if not len(draws):
                raise ValueError(f"no valid bootstrap draws for {rbp}/{calling_ip}/{metric}")
            low, high = percentile_interval(draws)
            rows.append({"metric": metric, "rbp": rbp, "calling_ip": calling_ip,
                         "n_shared": int((~specific).sum()), "n_specific": int(specific.sum()),
                         "median_shared": float(np.median(values[~specific])),
                         "median_specific": float(np.median(values[specific])),
                         "difference": observed, "ci_low": low, "ci_high": high,
                         "p_raw": centred_bootstrap_p(observed, draws),
                         "valid_resamples": len(draws)})
    result = pd.DataFrame(rows)
    result["p_holm"] = np.nan
    for _, indices in result.groupby("metric", sort=False).groups.items():
        result.loc[indices, "p_holm"] = holm(result.loc[indices, "p_raw"].to_numpy())
    return result
