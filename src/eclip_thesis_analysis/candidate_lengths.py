"""Summarise candidate intervals before and after boundary refinement."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def candidate_lengths(original: Path, refined: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read all BED intervals and retain every reported length."""
    tables = []
    for label, path in (("original", original), ("refined", refined)):
        table = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[1, 2],
                            compression="infer", names=["start", "end"])
        lengths = pd.to_numeric(table.end, errors="raise") - pd.to_numeric(table.start, errors="raise")
        if (lengths <= 0).any():
            raise ValueError(f"{path} contains a nonpositive interval length")
        tables.append(pd.DataFrame({"boundary": label, "length_nt": lengths.astype(np.int64)}))
    all_lengths = pd.concat(tables, ignore_index=True)
    summary = all_lengths.groupby("boundary", sort=False).length_nt.agg(
        n_candidates="size", median_length_nt="median", minimum_nt="min", maximum_nt="max",
        n_over_300_nt=lambda values: int((values > 300).sum())).reset_index()
    return all_lengths, summary
