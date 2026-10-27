"""Window-level read thinning for the reduced-depth evaluation."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pysam


def _window_reads(bam: pysam.AlignmentFile, chrom: str, start: int,
                  end: int, strand: str) -> list[np.ndarray]:
    """Collect clipped blocks from the transcript-orienting second mate."""
    blocks = []
    for read in bam.fetch(chrom, start, end):
        if (read.is_unmapped or read.is_duplicate or read.is_secondary
                or read.is_supplementary or not read.is_read2):
            continue
        if ("-" if read.is_reverse else "+") != strand:
            continue
        clipped = []
        for left, right in read.get_blocks():
            left, right = max(left, start), min(right, end)
            if right > left:
                if strand == "+":
                    clipped.append((left - start, right - start))
                else:
                    clipped.append((end - right, end - left))
        if clipped:
            blocks.append(np.asarray(clipped, dtype=np.int32))
    return blocks


def _coverage(reads: list[np.ndarray], length: int,
              keep: np.ndarray | None = None) -> np.ndarray:
    profile = np.zeros(length, dtype=np.float64)
    for index, intervals in enumerate(reads):
        if keep is not None and not keep[index]:
            continue
        for left, right in intervals:
            profile[left:right] += 1
    return profile


def simulate(data_path: Path, bam_paths: dict[str, Path], output: Path, *,
             fractions: tuple[float, ...] = (1.0, 0.75, 0.50, 0.30, 0.20, 0.10, 0.05),
             samples: int = 5, seed: int = 2026) -> None:
    """Recreate bigWig-scaled coverage, then thin reads independently per window."""
    if set(bam_paths) != {"ip1", "ip2", "sminput"}:
        raise ValueError("provide IP1, IP2, and SMInput BAMs")
    if samples <= 0 or any(not 0 < fraction <= 1 for fraction in fractions):
        raise ValueError("samples and read-retention fractions must be positive")
    with np.load(data_path, allow_pickle=False) as archive:
        data = {name: archive[name] for name in archive.files}
    needed = {"ip1", "ip2", "sminput", "chrom", "start", "end", "strand", "sequence"}
    if not needed.issubset(data):
        raise ValueError(f"prepared data is missing {sorted(needed - set(data))}")
    count, length = data["ip1"].shape
    if length != 300:
        raise ValueError("thinning expects 300 nt parent windows")
    output.mkdir(parents=True, exist_ok=True)
    metadata = {}
    read_sets = {}
    for track_index, name in enumerate(("ip1", "ip2", "sminput")):
        with pysam.AlignmentFile(str(bam_paths[name]), "rb") as bam:
            windows = [_window_reads(bam, str(data["chrom"][i]),
                                     int(data["start"][i]), int(data["end"][i]),
                                     str(data["strand"][i])) for i in range(count)]
        full = np.stack([_coverage(reads, length) for reads in windows])
        if full.sum() <= 0:
            raise ValueError(f"{name}: no contributing reads; check mate and strand conventions")
        scale = float(data[name].sum() / full.sum())
        scaled = full * scale
        if np.std(scaled) == 0 or np.std(data[name]) == 0:
            raise ValueError(f"{name}: full-depth coverage is constant")
        correlation = float(np.corrcoef(scaled.ravel(), data[name].ravel())[0, 1])
        if correlation < 0.999:
            raise ValueError(f"{name}: BAM-to-bigWig audit failed (Pearson r={correlation:.6f})")
        metadata[name] = {"scale_factor": scale, "full_depth_pearson": correlation,
                          "mate_convention": "read2"}
        read_sets[name] = windows
    (output / "read_coverage_audit.json").write_text(json.dumps(metadata, indent=2) + "\n")

    for fraction_index, fraction in enumerate(fractions):
        for sample in range(1, samples + 1):
            thinned = {}
            for track_index, name in enumerate(("ip1", "ip2", "sminput")):
                rng = np.random.default_rng(seed + 1_000_000 * track_index
                                            + 10_000 * sample + fraction_index)
                windows = read_sets[name]
                profiles = []
                for reads in windows:
                    keep = (rng.random(len(reads)) < fraction
                            if fraction < 1 else None)
                    profiles.append((_coverage(reads, length, keep) / fraction
                                     * metadata[name]["scale_factor"]).astype(np.float32))
                thinned[name] = np.stack(profiles)
            path = output / f"depth_{fraction:.2f}_sample{sample}.npz"
            np.savez_compressed(path, **{**data, **thinned})

