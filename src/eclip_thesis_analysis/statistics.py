"""Genomic-block resampling and multiplicity adjustment."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


def block_ids(chrom: np.ndarray, start: np.ndarray, block_size: int = 1_000_000) -> np.ndarray:
    """Identify occupied genomic blocks without merging chromosomes."""
    return np.asarray([f"{c}:{int(s) // block_size}" for c, s in zip(chrom, start)])


def bootstrap_indices(blocks: np.ndarray, rng: np.random.Generator,
                      chromosome: np.ndarray | None = None) -> np.ndarray:
    """Sample occupied blocks with replacement, retaining repeated draws."""
    blocks = np.asarray(blocks)
    if chromosome is None:
        strata = [np.arange(len(blocks))]
    else:
        chromosome = np.asarray(chromosome)
        strata = [np.flatnonzero(chromosome == name) for name in np.unique(chromosome)]
    parts = []
    for stratum in strata:
        names = np.unique(blocks[stratum])
        chosen = rng.choice(names, size=len(names), replace=True)
        lookup = {name: stratum[blocks[stratum] == name] for name in names}
        parts.extend(lookup[name] for name in chosen)
    return np.concatenate(parts)


def bootstrap(statistic: Callable[[np.ndarray], float], blocks: np.ndarray,
              repetitions: int, seed: int, *,
              chromosome: np.ndarray | None = None) -> tuple[float, np.ndarray]:
    if repetitions <= 0:
        raise ValueError("repetitions must be positive")
    observed = float(statistic(np.arange(len(blocks))))
    rng = np.random.default_rng(seed)
    estimates = np.asarray([
        statistic(bootstrap_indices(blocks, rng, chromosome))
        for _ in range(repetitions)
    ], dtype=float)
    return observed, estimates


def percentile_interval(estimates: np.ndarray, confidence: float = 0.95) -> tuple[float, float]:
    tail = (1 - confidence) / 2
    low, high = np.quantile(estimates, [tail, 1 - tail])
    return float(low), float(high)


def centred_bootstrap_p(observed: float, estimates: np.ndarray,
                         null: float = 0.0) -> float:
    """Two-sided Monte Carlo test using bootstrap errors centred on the estimate."""
    errors = np.asarray(estimates) - observed
    extreme = np.count_nonzero(np.abs(errors) >= abs(observed - null))
    return float((extreme + 1) / (len(errors) + 1))


def sign_flip_p(differences: np.ndarray, blocks: np.ndarray,
                statistic: str, repetitions: int, seed: int) -> tuple[float, float]:
    """Two-sided block sign-flip test for paired window differences."""
    values = np.asarray(differences, dtype=float)
    blocks = np.asarray(blocks)
    if len(values) != len(blocks) or not len(values) or not np.isfinite(values).all():
        raise ValueError("finite paired differences and matching blocks are required")
    if statistic not in {"mean", "median"}:
        raise ValueError("statistic must be mean or median")
    summary = np.mean if statistic == "mean" else np.median
    observed = float(summary(values))
    names, inverse = np.unique(blocks, return_inverse=True)
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(repetitions):
        signs = rng.choice(np.array([-1, 1]), size=len(names))
        simulated = float(summary(values * signs[inverse]))
        extreme += abs(simulated) >= abs(observed)
    return observed, float((extreme + 1) / (repetitions + 1))


def holm(p_values: np.ndarray) -> np.ndarray:
    """Holm-adjusted p-values in their original order."""
    values = np.asarray(p_values, dtype=float)
    if np.any((values < 0) | (values > 1) | ~np.isfinite(values)):
        raise ValueError("p-values must be finite and within [0, 1]")
    order = np.argsort(values, kind="stable")
    adjusted = np.empty(len(values))
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(values) - rank) * values[index])
        adjusted[index] = min(running, 1.0)
    return adjusted

