"""Base-pair call agreement from strand-matched BED intervals."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import pandas as pd


def _read_bed(path: Path) -> dict[tuple[str, str, int], list[tuple[int, int]]]:
    groups = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"BED file requires at least six columns: {path}")
            chrom, start, end, strand = fields[0], int(fields[1]), int(fields[2]), fields[5]
            if not 0 <= start < end or strand not in {"+", "-"}:
                raise ValueError(f"invalid BED interval in {path}")
            position = start
            while position < end:
                block = position // 1_000_000
                right = min(end, (block + 1) * 1_000_000)
                groups[(chrom, strand, block)].append((position, right))
                position = right
    return groups


def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def _intersection(left: list[tuple[int, int]], right: list[tuple[int, int]]) -> int:
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def pair_counts(first: Path, second: Path) -> dict[tuple[str, str, int], tuple[int, int]]:
    left, right = _read_bed(first), _read_bed(second)
    result = {}
    for key in left.keys() | right.keys():
        a, b = _merge(left.get(key, [])), _merge(right.get(key, []))
        intersection = _intersection(a, b)
        union = sum(end - start for start, end in a)
        union += sum(end - start for start, end in b) - intersection
        result[key] = (intersection, union)
    return result


def jaccard_block_table(rbp: str, seed: int, original_ip1: Path,
                        original_ip2: Path, model_p1: Path,
                        model_p2: Path) -> pd.DataFrame:
    """Count intersecting and union bases in each occupied 1 Mb block."""
    original = pair_counts(original_ip1, original_ip2)
    model = pair_counts(model_p1, model_p2)
    rows = []
    for chrom, strand, index in sorted(original.keys() | model.keys()):
        orig_inter, orig_union = original.get((chrom, strand, index), (0, 0))
        model_inter, model_union = model.get((chrom, strand, index), (0, 0))
        rows.append({"rbp": rbp, "seed": seed,
                     "block": f"{chrom}:{index}", "strand": strand,
                     "model_intersection": model_inter, "model_union": model_union,
                     "original_intersection": orig_inter, "original_union": orig_union})
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("no intervals in either original or reconstructed call sets")
    return frame.drop(columns="strand").groupby(
        ["rbp", "seed", "block"], as_index=False
    ).sum()

