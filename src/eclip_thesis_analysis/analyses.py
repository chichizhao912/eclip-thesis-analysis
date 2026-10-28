"""Selected thesis evaluations from tidy, study-specific tables."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from .statistics import (block_ids, bootstrap, centred_bootstrap_p, holm,
                         percentile_interval, sign_flip_p)


def pearson(first: np.ndarray, second: np.ndarray) -> float:
    if not np.isfinite(first).all() or not np.isfinite(second).all():
        return float("nan")
    if np.std(first) == 0 or np.std(second) == 0:
        return float("nan")
    return float(np.corrcoef(first, second)[0, 1])


def emd(first: np.ndarray, second: np.ndarray) -> float:
    """One-dimensional Wasserstein distance between positive coverage profiles."""
    if first.sum() <= 0 or second.sum() <= 0:
        return float("nan")
    left = np.cumsum(first / first.sum())
    right = np.cumsum(second / second.sum())
    return float(np.abs(left - right).sum())


def coverage_presence(data_path: Path, rbp: str) -> pd.DataFrame:
    """Count measurable IP coverage in each replicate of the supplied windows."""
    with np.load(data_path, allow_pickle=False) as source:
        first = np.sum(np.maximum(source["ip1"], 0), axis=1) > 0
        second = np.sum(np.maximum(source["ip2"], 0), axis=1) > 0
    counts = {"both": int(np.sum(first & second)),
              "ip1_only": int(np.sum(first & ~second)),
              "ip2_only": int(np.sum(~first & second)),
              "neither": int(np.sum(~first & ~second))}
    total = len(first)
    return pd.DataFrame([{"rbp": rbp, "n_windows": total, **counts,
                          "both_percent": 100 * counts["both"] / total if total else np.nan}])


def diagnostic_contrasts(data_path: Path, rbp: str) -> pd.DataFrame:
    """Within-window Pearson and EMD contrasts against SMInput."""
    with np.load(data_path, allow_pickle=False) as source:
        ip1, ip2, smi = (source[name] for name in ("ip1", "ip2", "sminput"))
        chrom, start, ids = (source[name] for name in ("chrom", "start", "window_id"))
        rows = []
        for i in range(len(ids)):
            tracks = (ip1[i], ip2[i], smi[i])
            if any(np.sum(track) <= 0 or np.std(track) == 0 for track in tracks):
                continue
            log1, log2, logsmi = (np.log1p(track) for track in tracks)
            ip_pearson = pearson(log1, log2)
            background_pearson = (pearson(log1, logsmi) + pearson(log2, logsmi)) / 2
            ip_emd = emd(ip1[i], ip2[i])
            background_emd = (emd(ip1[i], smi[i]) + emd(ip2[i], smi[i])) / 2
            for family, difference in (("pearson", ip_pearson - background_pearson),
                                       ("emd", background_emd - ip_emd)):
                rows.append({"rbp": rbp, "window_id": str(ids[i]),
                             "family": family, "comparison": rbp,
                             "chrom": str(chrom[i]), "start": int(start[i]),
                             "difference": float(difference)})
    return pd.DataFrame(rows)


def reconstruction_rows(data_path: Path, prediction_path: Path,
                        preprocessing_path: Path, configuration: str,
                        seed: int) -> pd.DataFrame:
    """Compute the two window-level correlations used in Figure 3.11."""
    cap = float(json.loads(preprocessing_path.read_text())["clip_value"])
    with np.load(data_path, allow_pickle=False) as source, \
         np.load(prediction_path, allow_pickle=False) as predicted:
        ids = source["window_id"]
        if not np.array_equal(ids, predicted["window_id"]):
            raise ValueError("data and prediction window IDs differ")
        ip1 = np.minimum(np.log1p(np.maximum(source["ip1"], 0)), cap)
        ip2 = np.minimum(np.log1p(np.maximum(source["ip2"], 0)), cap)
        p1, p2 = predicted["p1"], predicted["p2"]
        chrom, start = source["chrom"], source["start"]
        rows = []
        for i in range(len(ids)):
            agreement = pearson(p1[i], p2[i])
            target_1 = pearson(p1[i], ip1[i])
            target_2 = pearson(p2[i], ip2[i])
            rows.append({
                "configuration": configuration, "seed": seed,
                "window_id": str(ids[i]), "chrom": str(chrom[i]),
                "start": int(start[i]),
                "block": f"{chrom[i]}:{int(start[i]) // 1_000_000}",
                "reconstruction_agreement": agreement,
                "target_agreement": float(np.mean((target_1, target_2))),
            })
    return pd.DataFrame(rows)


def reconstruction_differences(table: pd.DataFrame, reference: str,
                               alternatives: list[str], family: str,
                               seeds: tuple[int, ...] = (42, 1337)) -> pd.DataFrame:
    """Average seed-level metrics within windows before paired comparisons."""
    metrics = ("reconstruction_agreement", "target_agreement")
    needed = {"configuration", "seed", "window_id", "chrom", "start", *metrics}
    if not needed.issubset(table):
        raise ValueError(f"table requires {sorted(needed)}")
    if not alternatives or reference in alternatives:
        raise ValueError("provide distinct reference and alternative configurations")
    selected = table[table["configuration"].isin([reference, *alternatives])].copy()
    if selected.duplicated(["configuration", "seed", "window_id"]).any():
        raise ValueError("configuration/seed/window rows must be unique")
    expected = set(seeds)
    counts = selected.groupby(["window_id", "configuration"])["seed"].agg(set)
    complete_seed_pairs = counts[counts.map(lambda found: found == expected)].index
    selected = selected.set_index(["window_id", "configuration"]).loc[
        selected.set_index(["window_id", "configuration"]).index.isin(complete_seed_pairs)
    ].reset_index()
    mean = selected.groupby(["window_id", "configuration"])[list(metrics)].mean().unstack()
    required = [reference, *alternatives]
    for metric in metrics:
        for config in required:
            if (metric, config) not in mean.columns:
                raise ValueError(f"missing {metric} measurements for {config}")
    complete = mean.dropna(subset=[(metric, config) for metric in metrics
                                   for config in required])
    metadata = selected.drop_duplicates("window_id").set_index("window_id")
    rows = []
    for window_id in complete.index:
        chrom, start = str(metadata.loc[window_id, "chrom"]), int(metadata.loc[window_id, "start"])
        for alternative in alternatives:
            for metric in metrics:
                rows.append({"family": family,
                             "comparison": f"{metric}: {reference} vs {alternative}",
                             "window_id": window_id, "chrom": chrom, "start": start,
                             "difference": float(complete.loc[window_id, (metric, reference)]
                                                 - complete.loc[window_id, (metric, alternative)])})
    return pd.DataFrame(rows)


def paired_tests(table: pd.DataFrame, *, method: str, statistic: str,
                 bootstrap_repetitions: int, test_repetitions: int, seed: int,
                 confidence: float = 0.95,
                 stratify_chromosome: bool = False) -> pd.DataFrame:
    """Test comparisons in a predeclared family, preserving genomic blocks."""
    needed = {"family", "comparison", "chrom", "start", "difference"}
    if not needed.issubset(table):
        raise ValueError(f"table requires {sorted(needed)}")
    rows = []
    for (family, comparison), group in table.groupby(["family", "comparison"], sort=False):
        clean = group.dropna(subset=["difference"]).copy()
        values = clean["difference"].to_numpy(float)
        if not len(values):
            raise ValueError(f"no finite differences for {family}/{comparison}")
        blocks = block_ids(clean["chrom"].to_numpy(), clean["start"].to_numpy())
        summary = np.mean if statistic == "mean" else np.median
        chromosomes = clean["chrom"].to_numpy() if stratify_chromosome else None
        observed, estimates = bootstrap(lambda idx: float(summary(values[idx])),
                                        blocks, bootstrap_repetitions, seed,
                                        chromosome=chromosomes)
        low, high = percentile_interval(estimates, confidence)
        if method == "sign-flip":
            _, p = sign_flip_p(values, blocks, statistic, test_repetitions, seed)
        elif method == "centred-bootstrap":
            p = centred_bootstrap_p(observed, estimates)
        else:
            raise ValueError("method must be sign-flip or centred-bootstrap")
        rows.append({"family": family, "comparison": comparison,
                     "n_windows": len(values), "n_blocks": len(np.unique(blocks)),
                     "estimate": observed, "ci_low": low, "ci_high": high,
                     "p_raw": p})
    result = pd.DataFrame(rows)
    result["p_holm"] = np.nan
    for _, indices in result.groupby("family", sort=False).groups.items():
        result.loc[indices, "p_holm"] = holm(result.loc[indices, "p_raw"].to_numpy())
    return result


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    if len(labels) != len(scores) or not labels.any() or not np.isfinite(scores).all():
        raise ValueError("finite scores and at least one positive label are required")
    ordered = labels[np.argsort(-scores, kind="stable")]
    precision = np.cumsum(ordered) / np.arange(1, len(ordered) + 1)
    return float(precision[ordered].mean())


def depth_scores(data_path: Path, prediction_path: Path, preprocessing_path: Path,
                 rbp: str, depth: float, sample: int) -> pd.DataFrame:
    """Score matched raw and reconstructed coverage from one thinning sample."""
    cap = float(json.loads(preprocessing_path.read_text())["clip_value"])
    with np.load(data_path, allow_pickle=False) as source, \
         np.load(prediction_path, allow_pickle=False) as predicted:
        if not np.array_equal(source["window_id"], predicted["window_id"]):
            raise ValueError("thinned data and prediction window IDs differ")
        selected = source["split"] == "test"
        if not selected.any():
            raise ValueError("no test windows in thinned dataset")
        ip1 = np.minimum(np.log1p(np.maximum(source["ip1"][selected], 0)), cap)
        ip2 = np.minimum(np.log1p(np.maximum(source["ip2"][selected], 0)), cap)
        sminput = np.minimum(np.log1p(np.maximum(source["sminput"][selected], 0)), cap)
        raw = (ip1 + ip2) / 2
        reconstructed = predicted["consensus"][selected]
        starts = np.arange(0, 201, 5)
        def parent_score(profiles):
            net = profiles - 0.25 * sminput
            prefix = np.pad(np.cumsum(net, axis=1), ((0, 0), (1, 0)))
            return ((prefix[:, starts + 100] - prefix[:, starts]) / 100).max(axis=1)
        ids = source["window_id"][selected]
        labels = source["window_class"][selected] == "peak"
        rows = []
        for approach, scores in (("original", parent_score(raw)),
                                 ("reconstructed", parent_score(reconstructed))):
            rows.extend({"rbp": rbp, "depth": depth, "sample": sample,
                         "approach": approach, "window_id": str(ids[i]),
                         "label": int(labels[i]), "score": float(scores[i])}
                        for i in range(len(ids)))
    return pd.DataFrame(rows)


def depth_ap(table: pd.DataFrame) -> pd.DataFrame:
    needed = {"rbp", "depth", "sample", "approach", "window_id", "label", "score"}
    if not needed.issubset(table):
        raise ValueError(f"table requires {sorted(needed)}")
    rows = []
    for keys, group in table.groupby(["rbp", "depth", "sample", "approach"], sort=False):
        rows.append(dict(zip(("rbp", "depth", "sample", "approach"), keys),
                         n_windows=len(group),
                         average_precision=average_precision(group["label"], group["score"])))
    result = pd.DataFrame(rows)
    return result.sort_values(["rbp", "depth", "sample", "approach"])


def outside_peak_motifs(table: pd.DataFrame, *, repetitions: int = 50_000,
                        seed: int = 2026) -> pd.DataFrame:
    """Compare motif-positive proportions for selected and other windows."""
    needed = {"rbp", "chrom", "start", "selected", "motif_present"}
    if not needed.issubset(table):
        raise ValueError(f"table requires {sorted(needed)}")
    rows = []
    for rbp, group in table.groupby("rbp", sort=False):
        for field in ("selected", "motif_present"):
            if not group[field].isin([0, 1]).all():
                raise ValueError(f"{field} must contain only 0 and 1")
        chosen = group["selected"].to_numpy(dtype=bool)
        motif = group["motif_present"].to_numpy(dtype=bool)
        if not chosen.any() or chosen.all():
            raise ValueError(f"{rbp} must have selected and comparison windows")
        blocks = block_ids(group["chrom"].to_numpy(), group["start"].to_numpy())
        def difference(idx):
            selected = chosen[idx]
            if not selected.any() or selected.all():
                return np.nan
            return float(motif[idx][selected].mean() - motif[idx][~selected].mean())
        observed, estimates = bootstrap(difference, blocks, repetitions, seed)
        valid = estimates[np.isfinite(estimates)]
        if len(valid) < 0.95 * repetitions:
            raise ValueError(f"too many invalid block resamples for {rbp}")
        low, high = percentile_interval(valid)
        rows.append({"rbp": rbp, "selected_n": int(chosen.sum()),
                     "comparison_n": int((~chosen).sum()),
                     "selected_percent": float(100 * motif[chosen].mean()),
                     "comparison_percent": float(100 * motif[~chosen].mean()),
                     "difference_pp": float(100 * observed),
                     "ci_low_pp": float(100 * low), "ci_high_pp": float(100 * high),
                     "p_raw": centred_bootstrap_p(observed, valid)})
    result = pd.DataFrame(rows)
    result["p_holm"] = holm(result["p_raw"].to_numpy())
    return result


def motif_rate_ratios(table: pd.DataFrame, *, seeds: tuple[int, ...] = (42, 1337, 2026),
                      repetitions: int = 50_000, random_seed: int = 2026) -> pd.DataFrame:
    """Resample block-level motif counts and scanned lengths together."""
    common = {"rbp", "condition", "block", "clipper_hits", "clipper_nt"}
    model_cols = {f"model_hits_{seed}" for seed in seeds}
    model_cols |= {f"model_nt_{seed}" for seed in seeds}
    if not (common | model_cols).issubset(table):
        raise ValueError(f"table requires {sorted(common | model_cols)}")
    rows = []
    for (rbp, condition), group in table.groupby(["rbp", "condition"], sort=False):
        if group["block"].duplicated().any():
            raise ValueError(f"duplicate genomic block in {rbp}/{condition}")
        counts = group[["clipper_hits", "clipper_nt",
                        *(f"model_hits_{seed}" for seed in seeds),
                        *(f"model_nt_{seed}" for seed in seeds)]].to_numpy(float)
        if not np.isfinite(counts).all() or np.any(counts < 0):
            raise ValueError(f"motif counts and lengths must be finite and nonnegative for {rbp}/{condition}")
        def ratio(idx):
            totals = counts[idx].sum(axis=0)
            clipper_rate = totals[0] / totals[1] if totals[1] > 0 else np.nan
            seed_ratios = [
                (totals[2+i] / totals[2+len(seeds)+i]) / clipper_rate
                for i in range(len(seeds))
            ]
            return float(np.mean(seed_ratios))
        blocks = group["block"].to_numpy()
        observed, estimates = bootstrap(ratio, blocks, repetitions, random_seed)
        valid = estimates[np.isfinite(estimates) & (estimates > 0)]
        if not np.isfinite(observed) or observed <= 0 or len(valid) < 0.95 * repetitions:
            raise ValueError(f"invalid motif-rate bootstrap for {rbp}/{condition}")
        low, high = percentile_interval(valid)
        rows.append({"rbp": rbp, "condition": condition,
                     "n_blocks": len(blocks), "mean_ratio": observed,
                     "ci_low": low, "ci_high": high,
                     "p_raw": centred_bootstrap_p(np.log(observed), np.log(valid))})
    result = pd.DataFrame(rows)
    result["p_holm"] = holm(result["p_raw"].to_numpy())
    return result


def jaccard_differences(table: pd.DataFrame, *, repetitions: int = 2_000,
                        seed: int = 2026) -> tuple[pd.DataFrame, float]:
    """Bootstrap seed-specific call-set differences from per-block base counts."""
    needed = {"rbp", "seed", "block", "model_intersection", "model_union",
              "original_intersection", "original_union"}
    if not needed.issubset(table):
        raise ValueError(f"table requires {sorted(needed)}")
    rows = []
    for (rbp, model_seed), group in table.groupby(["rbp", "seed"], sort=False):
        if group["block"].duplicated().any():
            raise ValueError(f"duplicate block for {rbp}/{model_seed}")
        values = group[["model_intersection", "model_union",
                        "original_intersection", "original_union"]].to_numpy(float)
        def delta(idx):
            total = values[idx].sum(axis=0)
            return float(total[0] / total[1] - total[2] / total[3])
        observed, estimates = bootstrap(delta, group["block"].to_numpy(),
                                        repetitions, seed)
        low, high = percentile_interval(estimates)
        rows.append({"rbp": rbp, "seed": model_seed, "difference": observed,
                     "ci_low": low, "ci_high": high, "n_blocks": len(group)})
    result = pd.DataFrame(rows)
    protein_means = result.groupby("rbp")["difference"].mean()
    positive = int((protein_means > 0).sum())
    sign_p = float(binomtest(positive, len(protein_means), 0.5,
                            alternative="two-sided").pvalue)
    return result, sign_p

