"""Compare candidate detector scores at matched hard-background calibration."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from eclip_reconstruct.caller import load_predictions, local_scores

from .analyses import average_precision


def caller_score_comparison(predictions: Path, seed: int, *,
                            peak_class: str = "peak_derived",
                            zero_class: str = "zero_ip") -> pd.DataFrame:
    """Assess consensus-only and SMInput-adjusted scores on validation parents."""
    data = load_predictions(predictions)
    valid = data["split"] == "validation"
    classes = data["window_class"][valid].astype(str)
    hard = classes == "hard_background"
    peak = classes == peak_class
    zero = classes == zero_class
    if not hard.any() or not peak.any():
        raise ValueError("validation windows need hard-background and peak-derived classes")
    consensus = data["consensus"][valid]
    sminput = data["sminput"][valid]
    schemes = {
        "consensus_only": local_scores(consensus, np.zeros_like(sminput))[0],
        "signal_minus_input": local_scores(consensus, sminput)[0],
    }
    rows = []
    for scheme, scans in schemes.items():
        maximum = scans.max(axis=1)
        threshold = float(np.quantile(maximum[hard], 0.95, method="higher"))
        called = maximum >= threshold
        rows.append({"seed": seed, "score_scheme": scheme,
                     "threshold": threshold, "n_validation": len(classes),
                     "n_peak_derived": int(peak.sum()),
                     "n_hard_background": int(hard.sum()),
                     "n_zero_ip": int(zero.sum()),
                     "peak_recovery_percent": float(100 * called[peak].mean()),
                     "hard_background_called_percent": float(100 * called[hard].mean()),
                     "zero_ip_called_percent": float(100 * called[zero].mean()) if zero.any() else np.nan,
                     "average_precision": average_precision(peak.astype(int), maximum)})
    return pd.DataFrame(rows)
