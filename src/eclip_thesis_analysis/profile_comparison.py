"""Variance-matched Gaussian comparison and dominant-feature measurements."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .analyses import pearson


def smooth_profiles(values: np.ndarray, sigma: float) -> np.ndarray:
    """Apply a truncated Gaussian kernel with reflected boundary padding."""
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    radius = int(np.ceil(4 * sigma))
    positions = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (positions / sigma) ** 2)
    kernel /= kernel.sum()
    padded = np.pad(np.asarray(values, dtype=float), ((0, 0), (radius, radius)), mode="reflect")
    return np.stack([np.convolve(row, kernel, mode="valid") for row in padded])


def dominant_fwhm(profile: np.ndarray) -> float:
    """Interpolate crossings around the dominant internal feature's half height."""
    y = np.asarray(profile, dtype=float)
    if y.ndim != 1 or len(y) < 3 or not np.isfinite(y).all():
        return np.nan
    peak = int(np.argmax(y))
    if peak in (0, len(y) - 1):
        return np.nan
    baseline = float(y.min())
    if y[peak] - baseline <= 1e-8:
        return np.nan
    half = baseline + (y[peak] - baseline) / 2
    left = peak
    while left > 0 and y[left] >= half:
        left -= 1
    right = peak
    while right < len(y) - 1 and y[right] >= half:
        right += 1
    if left == 0 and y[left] >= half or right == len(y) - 1 and y[right] >= half:
        return np.nan
    left_cross = left + (half - y[left]) / (y[left + 1] - y[left])
    right_cross = right - 1 + (y[right - 1] - half) / (y[right - 1] - y[right])
    width = right_cross - left_cross
    return float(width) if np.isfinite(width) and width > 0 else np.nan


def gaussian_comparison_rows(data_path: Path, prediction_path: Path,
                             preprocessing_path: Path, seed: int,
                             sigma: float = 26.19) -> pd.DataFrame:
    """Measure model, crossed Gaussian, and observed profiles in each parent."""
    cap = float(json.loads(preprocessing_path.read_text())["clip_value"])
    with np.load(data_path, allow_pickle=False) as source, \
         np.load(prediction_path, allow_pickle=False) as predicted:
        if not np.array_equal(source["window_id"], predicted["window_id"]):
            raise ValueError("data and prediction window IDs differ")
        ip1 = np.minimum(np.log1p(np.maximum(source["ip1"], 0)), cap)
        ip2 = np.minimum(np.log1p(np.maximum(source["ip2"], 0)), cap)
        # The smoothed opposite source is paired with each measured target,
        # preserving the directional comparison used for the autoencoder.
        g1, g2 = smooth_profiles(ip2, sigma), smooth_profiles(ip1, sigma)
        p1, p2 = predicted["p1"], predicted["p2"]
        rows = []
        for i, window_id in enumerate(source["window_id"]):
            row = {"window_id": str(window_id), "chrom": str(source["chrom"][i]),
                   "start": int(source["start"][i]), "seed": seed, "sigma_nt": sigma,
                   "raw_agreement": pearson(ip1[i], ip2[i]),
                   "observed_fwhm_nt": dominant_fwhm((ip1[i] + ip2[i]) / 2),
                   "target_var_p1": float(np.var(ip1[i])),
                   "target_var_p2": float(np.var(ip2[i]))}
            for name, first, second in (("autoencoder", p1[i], p2[i]),
                                        ("gaussian", g1[i], g2[i])):
                row[f"{name}_agreement"] = pearson(first, second)
                targets = (pearson(first, ip1[i]), pearson(second, ip2[i]))
                row[f"{name}_target_agreement"] = (float(np.mean(targets))
                                                   if np.isfinite(targets).all() else np.nan)
                row[f"{name}_fwhm_nt"] = dominant_fwhm((first + second) / 2)
                if name == "autoencoder":
                    row["autoencoder_consensus_profile"] = json.dumps(
                        ((first + second) / 2).tolist(), separators=(",", ":"))
                row[f"{name}_pred_var_p1"] = float(np.var(first))
                row[f"{name}_pred_var_p2"] = float(np.var(second))
                row[f"{name}_log_total_error"] = float(np.mean((
                    abs(np.log1p(first.sum()) - np.log1p(ip1[i].sum())),
                    abs(np.log1p(second.sum()) - np.log1p(ip2[i].sum())))))
            rows.append(row)
    return pd.DataFrame(rows)


def gaussian_summary(table: pd.DataFrame, *, seeds: tuple[int, ...] = (42, 1337)) -> pd.DataFrame:
    """Aggregate correlation, width, variance, and total-signal measures."""
    selected = _complete_seeds(table, seeds)
    correlation_columns = ["raw_agreement", "autoencoder_agreement", "gaussian_agreement",
                           "autoencoder_target_agreement", "gaussian_target_agreement"]
    correlation_valid = selected.groupby("window_id")[correlation_columns].apply(
        lambda group: np.isfinite(group.to_numpy(float)).all())
    correlation = selected[selected.window_id.isin(correlation_valid[correlation_valid].index)]
    widths = _ensemble_widths(selected).dropna()
    rows = []
    for method in ("autoencoder", "gaussian"):
        by_window = correlation.groupby("window_id")[[f"{method}_agreement",
                                                        f"{method}_target_agreement"]].mean()
        width_by_window = widths[f"{method}_fwhm_nt"]
        retention = []
        for seed, group in selected.groupby("seed"):
            directional = []
            for direction in ("p1", "p2"):
                denominator = group[f"target_var_{direction}"].sum()
                directional.append(group[f"{method}_pred_var_{direction}"].sum() / denominator)
            retention.append(float(np.mean(directional)))
        rows.append({"method": method, "n_windows": selected.window_id.nunique(),
                     "n_correlation_windows": len(by_window),
                     "n_width_windows": len(width_by_window),
                     "mean_agreement": float(by_window[f"{method}_agreement"].mean()),
                     "mean_target_agreement": float(by_window[f"{method}_target_agreement"].mean()),
                     "median_fwhm_nt": float(width_by_window.median()),
                     "variance_retention": float(np.mean(retention)),
                     "mean_log_total_error": float(selected[f"{method}_log_total_error"].mean())})
    return pd.DataFrame(rows)


def _complete_seeds(table: pd.DataFrame, seeds: tuple[int, ...]) -> pd.DataFrame:
    needed = {"window_id", "chrom", "start", "seed", "autoencoder_agreement",
              "gaussian_agreement", "autoencoder_target_agreement",
              "gaussian_target_agreement", "observed_fwhm_nt",
              "autoencoder_fwhm_nt", "gaussian_fwhm_nt",
              "autoencoder_consensus_profile"}
    if not needed.issubset(table):
        raise ValueError(f"comparison table requires {sorted(needed)}")
    selected = table[table.seed.isin(seeds)].copy()
    if selected.duplicated(["window_id", "seed"]).any():
        raise ValueError("window/seed rows must be unique")
    complete = selected.groupby("window_id").seed.agg(set)
    return selected[selected.window_id.isin(complete[complete == set(seeds)].index)]


def _ensemble_widths(selected: pd.DataFrame) -> pd.DataFrame:
    """Measure width after averaging the seed consensus profiles per window."""
    rows = []
    for window_id, group in selected.groupby("window_id", sort=False):
        profiles = [np.asarray(json.loads(value), dtype=float)
                    for value in group.autoencoder_consensus_profile]
        if len({profile.shape for profile in profiles}) != 1:
            raise ValueError(f"inconsistent consensus profile lengths for {window_id}")
        rows.append({"window_id": window_id,
                     "autoencoder_fwhm_nt": dominant_fwhm(np.mean(profiles, axis=0)),
                     "observed_fwhm_nt": group.observed_fwhm_nt.iloc[0],
                     "gaussian_fwhm_nt": group.gaussian_fwhm_nt.iloc[0]})
    return pd.DataFrame(rows).set_index("window_id")


def gaussian_differences(table: pd.DataFrame, comparison: str,
                         *, seeds: tuple[int, ...] = (42, 1337)) -> pd.DataFrame:
    """Produce paired differences for mean-correlation or median-width tests."""
    if comparison not in {"correlation", "width"}:
        raise ValueError("comparison must be correlation or width")
    selected = _complete_seeds(table, seeds)
    if comparison == "correlation":
        columns = ["autoencoder_agreement", "gaussian_agreement",
                   "autoencoder_target_agreement", "gaussian_target_agreement"]
        valid = selected.groupby("window_id")[columns].apply(
            lambda group: np.isfinite(group.to_numpy(float)).all())
        mean = selected[selected.window_id.isin(valid[valid].index)].groupby("window_id")[columns].mean()
    else:
        mean = _ensemble_widths(selected).dropna()
    metadata = selected.drop_duplicates("window_id").set_index("window_id")
    pairs = ([ ("agreement", "autoencoder_agreement", "gaussian_agreement"),
               ("target_agreement", "autoencoder_target_agreement", "gaussian_target_agreement")]
             if comparison == "correlation" else
             [("autoencoder_minus_observed", "autoencoder_fwhm_nt", "observed_fwhm_nt"),
              ("gaussian_minus_autoencoder", "gaussian_fwhm_nt", "autoencoder_fwhm_nt")])
    rows = []
    for window_id, values in mean.iterrows():
        for label, first, second in pairs:
            rows.append({"family": f"gaussian_{comparison}", "comparison": label,
                         "window_id": window_id, "chrom": metadata.loc[window_id, "chrom"],
                         "start": int(metadata.loc[window_id, "start"]),
                         "difference": float(values[first] - values[second])})
    return pd.DataFrame(rows)
