"""Small command-line entry points for the principal thesis comparisons."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .analyses import (depth_ap, depth_scores, diagnostic_contrasts,
                       jaccard_differences, motif_rate_ratios,
                       outside_peak_motifs, paired_tests, reconstruction_differences,
                       reconstruction_rows)
from .intervals import jaccard_block_table
from .motifs import motif_block_counts, scan_intervals
from .ranking import rank_windows
from .read_depth import simulate


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(prog="eclip-thesis-analysis")
    sub = command.add_subparsers(dest="command", required=True)

    collect = sub.add_parser("reconstruction-metrics")
    collect.add_argument("--data", type=Path, required=True)
    collect.add_argument("--predictions", type=Path, required=True)
    collect.add_argument("--preprocessing", type=Path, required=True)
    collect.add_argument("--configuration", required=True)
    collect.add_argument("--seed", type=int, required=True)
    collect.add_argument("--output", type=Path, required=True)

    compare = sub.add_parser("reconstruction-differences")
    compare.add_argument("--input", type=Path, required=True)
    compare.add_argument("--reference", required=True)
    compare.add_argument("--alternatives", nargs="+", required=True)
    compare.add_argument("--family", required=True)
    compare.add_argument("--seeds", nargs="+", type=int, default=[42, 1337])
    compare.add_argument("--output", type=Path, required=True)

    controls = sub.add_parser("train-controls")
    controls.add_argument("--data", type=Path, required=True)
    controls.add_argument("--output", type=Path, required=True)
    controls.add_argument("--configurations", nargs="+",
                          choices=("coverage", "full", "shuffled", "no_sminput"),
                          default=["coverage", "full", "shuffled", "no_sminput"])
    controls.add_argument("--seeds", nargs="+", type=int, default=[42, 1337])

    diagnostics = sub.add_parser("diagnostic-contrasts")
    diagnostics.add_argument("--data", type=Path, required=True)
    diagnostics.add_argument("--rbp", required=True)
    diagnostics.add_argument("--output", type=Path, required=True)

    tests = sub.add_parser("paired-tests")
    tests.add_argument("--input", type=Path, required=True)
    tests.add_argument("--output", type=Path, required=True)
    tests.add_argument("--method", choices=("sign-flip", "centred-bootstrap"), required=True)
    tests.add_argument("--statistic", choices=("mean", "median"), required=True)
    tests.add_argument("--bootstrap-resamples", type=int, default=3_000)
    tests.add_argument("--test-resamples", type=int, default=20_000)
    tests.add_argument("--stratify-chromosome", action="store_true")
    tests.add_argument("--seed", type=int, default=2026)

    depth = sub.add_parser("depth-ap")
    depth.add_argument("--input", type=Path, required=True)
    depth.add_argument("--output", type=Path, required=True)

    thinning = sub.add_parser("thin-reads")
    thinning.add_argument("--data", type=Path, required=True)
    thinning.add_argument("--ip1-bam", type=Path, required=True)
    thinning.add_argument("--ip2-bam", type=Path, required=True)
    thinning.add_argument("--sminput-bam", type=Path, required=True)
    thinning.add_argument("--output", type=Path, required=True)
    thinning.add_argument("--samples", type=int, default=5)
    thinning.add_argument("--seed", type=int, default=2026)

    ranking = sub.add_parser("depth-scores")
    ranking.add_argument("--data", type=Path, required=True)
    ranking.add_argument("--predictions", type=Path, required=True)
    ranking.add_argument("--preprocessing", type=Path, required=True)
    ranking.add_argument("--rbp", required=True)
    ranking.add_argument("--depth", type=float, required=True)
    ranking.add_argument("--sample", type=int, required=True)
    ranking.add_argument("--output", type=Path, required=True)

    motifs = sub.add_parser("outside-peak-motifs")
    motifs.add_argument("--input", type=Path, required=True)
    motifs.add_argument("--output", type=Path, required=True)
    motifs.add_argument("--resamples", type=int, default=50_000)
    motifs.add_argument("--seed", type=int, default=2026)

    motif_scan = sub.add_parser("motif-scan")
    motif_scan.add_argument("--input", type=Path, required=True)
    motif_scan.add_argument("--reference", type=Path, required=True)
    motif_scan.add_argument("--matrix", type=Path, required=True)
    motif_scan.add_argument("--training-data", type=Path, required=True)
    motif_scan.add_argument("--motif-id", required=True)
    motif_scan.add_argument("--output", type=Path, required=True)

    motif_blocks = sub.add_parser("motif-blocks")
    motif_blocks.add_argument("--input", type=Path, required=True)
    motif_blocks.add_argument("--output", type=Path, required=True)

    rates = sub.add_parser("motif-rates")
    rates.add_argument("--input", type=Path, required=True)
    rates.add_argument("--output", type=Path, required=True)
    rates.add_argument("--resamples", type=int, default=50_000)
    rates.add_argument("--seed", type=int, default=2026)

    jaccard = sub.add_parser("call-jaccard")
    jaccard.add_argument("--input", type=Path, required=True)
    jaccard.add_argument("--output", type=Path, required=True)
    jaccard.add_argument("--resamples", type=int, default=2_000)
    jaccard.add_argument("--seed", type=int, default=2026)

    blocks = sub.add_parser("jaccard-blocks")
    blocks.add_argument("--rbp", required=True)
    blocks.add_argument("--seed", type=int, required=True)
    blocks.add_argument("--original-ip1", type=Path, required=True)
    blocks.add_argument("--original-ip2", type=Path, required=True)
    blocks.add_argument("--model-p1", type=Path, required=True)
    blocks.add_argument("--model-p2", type=Path, required=True)
    blocks.add_argument("--output", type=Path, required=True)

    ranking_outside = sub.add_parser("rank-outside-peaks")
    ranking_outside.add_argument("--data", type=Path, required=True)
    ranking_outside.add_argument("--predictions", type=Path, required=True)
    ranking_outside.add_argument("--original-ip1", type=Path, required=True)
    ranking_outside.add_argument("--original-ip2", type=Path, required=True)
    ranking_outside.add_argument("--rbp", required=True)
    ranking_outside.add_argument("--top", type=int, default=1000)
    ranking_outside.add_argument("--output", type=Path, required=True)
    return command


def main(argv: list[str] | None = None) -> None:
    args = parser().parse_args(argv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.command == "reconstruction-metrics":
        result = reconstruction_rows(args.data, args.predictions, args.preprocessing,
                                     args.configuration, args.seed)
    elif args.command == "train-controls":
        from .controls import run_controls

        run_controls(args.data, args.output, tuple(args.configurations), tuple(args.seeds))
        print(f"Control checkpoints and predictions: {args.output}")
        return
    elif args.command == "diagnostic-contrasts":
        result = diagnostic_contrasts(args.data, args.rbp)
    elif args.command == "reconstruction-differences":
        table = pd.read_csv(args.input, sep="\t")
        result = reconstruction_differences(table, args.reference,
                                            args.alternatives, args.family,
                                            tuple(args.seeds))
    elif args.command == "thin-reads":
        simulate(args.data, {"ip1": args.ip1_bam, "ip2": args.ip2_bam,
                             "sminput": args.sminput_bam}, args.output,
                 samples=args.samples, seed=args.seed)
        print(f"Thinned profiles: {args.output}")
        return
    elif args.command == "depth-scores":
        result = depth_scores(args.data, args.predictions, args.preprocessing,
                              args.rbp, args.depth, args.sample)
    elif args.command == "jaccard-blocks":
        result = jaccard_block_table(args.rbp, args.seed, args.original_ip1,
                                     args.original_ip2, args.model_p1,
                                     args.model_p2)
    elif args.command == "rank-outside-peaks":
        result = rank_windows(args.data, args.predictions,
                              [args.original_ip1, args.original_ip2],
                              args.rbp, args.top)
    elif args.command == "motif-scan":
        result = scan_intervals(pd.read_csv(args.input, sep="\t"), args.reference,
                                args.matrix, args.training_data, args.motif_id)
    else:
        table = pd.read_csv(args.input, sep="\t")
        if args.command == "paired-tests":
            result = paired_tests(table, method=args.method, statistic=args.statistic,
                                  bootstrap_repetitions=args.bootstrap_resamples,
                                  test_repetitions=args.test_resamples,
                                  stratify_chromosome=args.stratify_chromosome,
                                  seed=args.seed)
        elif args.command == "depth-ap":
            result = depth_ap(table)
        elif args.command == "outside-peak-motifs":
            result = outside_peak_motifs(table, repetitions=args.resamples,
                                         seed=args.seed)
        elif args.command == "motif-rates":
            result = motif_rate_ratios(table, repetitions=args.resamples,
                                       random_seed=args.seed)
        elif args.command == "motif-blocks":
            result = motif_block_counts(table)
        else:
            result, sign_p = jaccard_differences(table, repetitions=args.resamples,
                                                 seed=args.seed)
            args.output.with_suffix(".sign_test.json").write_text(
                json.dumps({"two_sided_exact_sign_test_p": sign_p}, indent=2) + "\n"
            )
    result.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(result)} rows to {args.output}")


if __name__ == "__main__":
    main()

