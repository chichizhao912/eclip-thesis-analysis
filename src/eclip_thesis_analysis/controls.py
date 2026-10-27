"""The four prespecified TARDBP input-feature controls."""

from __future__ import annotations

import csv
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from eclip_reconstruct.data import fit_cap, load_data, transform
from eclip_reconstruct.model import CoverageAutoencoder


CONFIGURATIONS = ("coverage", "full", "shuffled", "no_sminput")
CHANNELS = {"coverage": 2, "full": 6, "shuffled": 6, "no_sminput": 5}
SHUFFLE_SEED = 20260827


def dinucleotide_shuffle(sequence: str, rng: np.random.Generator) -> str:
    """Randomised Eulerian traversal that retains every adjacent-base count."""
    adjacent: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(sequence[:-1], sequence[1:]):
        adjacent[left].append(right)
    for edges in adjacent.values():
        rng.shuffle(edges)
    stack, trail = [sequence[0]], []
    while stack:
        base = stack[-1]
        if adjacent[base]:
            stack.append(adjacent[base].pop())
        else:
            trail.append(stack.pop())
    shuffled = "".join(reversed(trail))
    if len(shuffled) != len(sequence) or Counter(zip(sequence[:-1], sequence[1:])) != Counter(
            zip(shuffled[:-1], shuffled[1:])):
        raise RuntimeError("dinucleotide shuffle failed its composition check")
    return shuffled


def shuffled_one_hot(sequence: np.ndarray, window_ids: np.ndarray) -> np.ndarray:
    """Use one reproducible shuffled sequence for both training seeds."""
    output = np.zeros_like(sequence, dtype=np.float32)
    for index, (encoded, window_id) in enumerate(zip(sequence, window_ids)):
        letters = np.array(list("ACGT"))
        original = "".join(letters[np.argmax(encoded, axis=0)][position]
                           if encoded[:, position].any() else "N"
                           for position in range(encoded.shape[1]))
        best, best_identity = original, 1.0
        for attempt in range(64):
            payload = f"{SHUFFLE_SEED}|{window_id}|{attempt}".encode()
            seed = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
            candidate = dinucleotide_shuffle(original, np.random.default_rng(seed))
            identity = sum(a == b for a, b in zip(candidate, original)) / len(original)
            if identity < best_identity:
                best, best_identity = candidate, identity
            if identity < 0.95:
                break
        for channel, base in enumerate("ACGT"):
            output[index, channel] = np.fromiter((letter == base for letter in best),
                                                  dtype=np.float32)
    return output


def inputs(configuration: str, sequence: torch.Tensor,
           source_ip: torch.Tensor, sminput: torch.Tensor) -> torch.Tensor:
    if configuration == "coverage":
        return torch.stack((source_ip, sminput), dim=1)
    if configuration == "no_sminput":
        return torch.cat((sequence, source_ip[:, None, :]), dim=1)
    if configuration in {"full", "shuffled"}:
        return torch.cat((sequence, source_ip[:, None, :], sminput[:, None, :]), dim=1)
    raise ValueError(f"unknown configuration: {configuration}")


def _device() -> torch.device:
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def run_controls(data_path: Path, output: Path,
                 configurations: tuple[str, ...] = CONFIGURATIONS,
                 seeds: tuple[int, ...] = (42, 1337),
                 maximum_epochs: int = 8, patience: int = 2) -> None:
    """Train controls and write predictions in the tool's standard NPZ format."""
    if not configurations or any(name not in CONFIGURATIONS for name in configurations):
        raise ValueError(f"configurations must come from {CONFIGURATIONS}")
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError("provide unique training seeds")
    data = load_data(data_path)
    train_idx = np.flatnonzero(data["split"] == "train")
    validation_idx = np.flatnonzero(data["split"] == "validation")
    if not len(train_idx) or not len(validation_idx):
        raise ValueError("train and validation windows are required")
    if set(data["chrom"][train_idx]) & set(data["chrom"][validation_idx]):
        raise ValueError("train and validation chromosomes must be disjoint")
    cap = fit_cap(data)
    observed = {name: transform(data[name], cap)
                for name in ("ip1", "ip2", "sminput")}
    sequence = data["sequence"].astype(np.float32)
    shuffled = (shuffled_one_hot(sequence, data["window_id"])
                if "shuffled" in configurations else None)
    output.mkdir(parents=True, exist_ok=True)
    (output / "preprocessing.json").write_text(json.dumps({
        "transform": "min(log1p(max(coverage, 0)), cap)",
        "fit_partition": "train", "clip_percentile": 99.9,
        "clip_value": cap,
    }, indent=2) + "\n")
    device = _device()
    history = []
    for configuration in configurations:
        chosen_sequence = shuffled if configuration == "shuffled" else sequence
        tensors = [torch.from_numpy(observed[name])
                   for name in ("ip1", "ip2", "sminput")]
        tensors.append(torch.from_numpy(chosen_sequence))
        dataset = TensorDataset(*tensors)
        for seed in seeds:
            random.seed(seed)
            np.random.seed(seed)
            torch.manual_seed(seed)
            model = CoverageAutoencoder(input_channels=CHANNELS[configuration]).to(device)
            optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
            best, stale = float("inf"), 0
            checkpoint_path = output / f"{configuration}_seed{seed}.pt"
            for epoch in range(1, maximum_epochs + 1):
                losses = {}
                for phase, indices in (("train", train_idx),
                                       ("validation", validation_idx)):
                    model.train(phase == "train")
                    loader = DataLoader(torch.utils.data.Subset(dataset, indices),
                                        batch_size=64, shuffle=phase == "train")
                    values = []
                    for batch in loader:
                        ip1, ip2, smi, seq = (item.to(device) for item in batch)
                        with torch.set_grad_enabled(phase == "train"):
                            p2, z1 = model(inputs(configuration, seq, ip1, smi))
                            p1, z2 = model(inputs(configuration, seq, ip2, smi))
                            reconstruction = (nn.functional.huber_loss(p1[:, 0], ip1)
                                              + nn.functional.huber_loss(p2[:, 0], ip2))
                            loss = reconstruction + 0.05 * nn.functional.mse_loss(z1, z2)
                            if phase == "train":
                                optimizer.zero_grad()
                                loss.backward()
                                optimizer.step()
                        values.append(float(loss.detach().cpu()))
                    losses[phase] = float(np.mean(values))
                history.append({"configuration": configuration, "seed": seed,
                                "epoch": epoch, "train_loss": losses["train"],
                                "validation_loss": losses["validation"]})
                if losses["validation"] < best:
                    best, stale = losses["validation"], 0
                    torch.save({"model_state": model.state_dict(),
                                "configuration": configuration, "seed": seed,
                                "clip_value": cap, "input_channels": CHANNELS[configuration],
                                "best_epoch": epoch, "validation_loss": best}, checkpoint_path)
                else:
                    stale += 1
                    if stale >= patience:
                        break
            saved = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
            model.load_state_dict(saved["model_state"])
            model.eval()
            p1_parts, p2_parts = [], []
            with torch.no_grad():
                for first in range(0, len(dataset), 64):
                    ip1, ip2, smi, seq = (item.to(device) for item in
                                           dataset[first:first + 64])
                    p1, _ = model(inputs(configuration, seq, ip2, smi))
                    p2, _ = model(inputs(configuration, seq, ip1, smi))
                    p1_parts.append(p1[:, 0].cpu().numpy())
                    p2_parts.append(p2[:, 0].cpu().numpy())
            p1 = np.concatenate(p1_parts).astype(np.float32)
            p2 = np.concatenate(p2_parts).astype(np.float32)
            np.savez_compressed(
                output / f"{configuration}_seed{seed}_predictions.npz",
                p1=p1, p2=p2, consensus=(p1 + p2) / 2,
                sminput=observed["sminput"],
                **{name: data[name] for name in
                   ("window_id", "chrom", "start", "end", "strand", "split", "window_class")},
            )
    with (output / "training_history.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=history[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(history)
