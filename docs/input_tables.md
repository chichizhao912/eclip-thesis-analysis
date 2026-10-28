# Input table formats

All tables are UTF-8, tab-delimited, and have a header row. Genomic coordinates are zero-based and half-open. `chrom` values must use the same assembly and naming convention across inputs. A genomic block is `chrom` plus the integer quotient of `start / 1,000,000`.

## `coverage-presence` and `diagnostic-contrasts`

Both commands read the prepared-window NPZ from `eclip-reconstruct prepare`, with aligned `ip1`, `ip2`, and `sminput` arrays and `window_id`, `chrom`, and `start` metadata. Supply the intended peak-derived evaluation population. `coverage-presence` counts positive measured IP coverage per window. `diagnostic-contrasts` emits one Pearson and one EMD contrast per eligible window; the Pearson calculation uses log-transformed coverage, while EMD uses raw positive coverage.

## `peak-diagnostics` and `peak-shape-tests`

`peak-diagnostics` takes the two complete, replicate-specific published CLIPper narrowPeak files and six strand-specific bigWigs: positive- and negative-strand tracks for IP1, IP2, and SMInput. The peak files need columns 1–3 for genomic coordinates and column 6 for strand. Negative-strand bigWig signal is interpreted as negative-valued genomic signal and reversed into transcript orientation. The command uses canonical GRCh38 chromosomes and reports only valid peak intervals.

The per-IP summary contains peak counts, same-strand overlap fraction, a shared base-pair Jaccard estimate with a pointwise 95% chromosome-stratified 1 Mb block-bootstrap interval, and matched summit displacement. For each peak called by only one IP, `partner_support` classifies the partner IP's summed coverage relative to SMInput within both the exact published interval and a 300 nt window centred on the calling-IP summit. `shape` contains sampled shared and replicate-specific peak measurements: reported width, contiguous half-height span around the calling-IP maximum, and maximum-to-mean coverage ratio. Peak-shape sampling happens before removing peaks with no positive calling-IP coverage.

Combine the six RBP shape tables into one table before running `peak-shape-tests`; each must retain `rbp`, `calling_ip`, `chrom`, `start`, `peak_class`, `width_nt`, `half_height_span_nt`, and `concentration`. The command tests the difference between replicate-specific and shared medians with a two-sided centred chromosome-stratified genomic-block bootstrap. It adjusts the 12 RBP-by-IP p-values separately for width, half-height span, and concentration.

## `paired-tests`

Required columns: `family`, `comparison`, `chrom`, `start`, `difference`. Each row is one window-level paired contrast. `family` defines which p-values are adjusted together with Holm's method. For the six-RBP diagnostic tests use separate `pearson` and `emd` families. For the input ablation, use one family containing all six comparisons shown in Figure 3.11. Do not put unrelated exploratory tests into the same file merely because they share a command.

For sign-flip tests, use `--bootstrap-resamples 3000 --test-resamples 20000`; the former controls the pointwise percentile interval and the latter controls Monte Carlo p-value resolution. For diagnostic centred-bootstrap tests, use `--bootstrap-resamples 2000 --stratify-chromosome`. The `--test-resamples` option is ignored for a centred-bootstrap test because its bootstrap estimates also form the null distribution.

## `reconstruction-metrics`

`train-controls` uses the prepared-window NPZ created by the tool repository. It trains the IP+SMInput, full sequence, dinucleotide-shuffled sequence, and sequence-without-SMInput configurations on identical train and validation partitions. The same fixed shuffle is reused for both model seeds. Each output prediction NPZ can be passed directly to `reconstruction-metrics` with the shared preprocessing JSON file.

`--data` is the prepared NPZ from `eclip-reconstruct prepare`. `--predictions` is its prediction NPZ for one configuration and one model seed. `--preprocessing` is the training output `preprocessing.json`. The command writes one row per window with `reconstruction_agreement = r(P1, P2)` and `target_agreement = [r(P1, IP1) + r(P2, IP2)] / 2`. Correlations are over the 300 nucleotide positions. Undefined correlations are recorded as missing.

`reconstruction-differences` combines these rows across configurations and seeds. It requires complete results for every requested configuration, both reported metrics, and both named training seeds. It averages the seed-specific correlations within each window before subtracting each alternative from the reference.

## `gaussian-comparison`, `gaussian-summary`, and `gaussian-differences`

`gaussian-comparison` accepts the same prepared NPZ, prediction NPZ, and preprocessing JSON as `reconstruction-metrics`, plus a seed identifier. Run it for both model seeds on the same chromosome 1 windows and combine the two output tables. Its default Gaussian standard deviation is the chromosome 2 variance-matched value of 26.19 nt. It reports the per-window P1–P2 and prediction–target correlations for the model and crossed Gaussian comparator, each direction's positional variance, absolute log-total signal error, and dominant-feature width. The output also retains each model consensus profile so that width can be measured after the two seed profiles are averaged.

`gaussian-summary` uses complete two-seed windows for correlations and reports their mean across windows. Positional-variance retention is the ratio of summed predicted to summed measured within-window variances, calculated per direction and seed, then averaged. Log-total signal error is the mean absolute difference between `log1p` predicted and target profile sums. Width uses the dominant internal feature's half-height crossings after averaging each window's two seed consensus profiles; it excludes flat or edge-truncated features.

`gaussian-differences --comparison correlation` generates two model-minus-Gaussian paired mean-correlation contrasts. `--comparison width` generates two paired median-width contrasts: autoencoder minus observed and Gaussian minus autoencoder. Each output declares its own Holm family for `paired-tests`.

## `caller-score-comparison`

Supply a prediction NPZ produced by `eclip-reconstruct predict` for one model seed. It must include `consensus`, transformed `sminput`, `split`, and `window_class` alongside the standard window metadata. The command selects validation windows and expects `hard_background`, `peak_derived`, and optionally `zero_ip` class labels by default; use `--peak-class` or `--zero-class` if the prepared manifest uses different names.

The command calculates the maximum of 41 local 100 nt scores per 300 nt parent. It separately calibrates consensus-only and SMInput-adjusted scores against the higher empirical 95th percentile of the validation hard-background maxima. Reported calling percentages refer to parents with a maximum score at or above the corresponding threshold. Average precision treats peak-derived parents as the positive internal reference and all other validation parents as the comparison class. Combine the two seed outputs and average corresponding scheme-level estimates for the thesis development summary.

## `candidate-lengths`

Supply the saved BED call sets from the original merged-interval boundary procedure and the refined boundary procedure applied to the same development sample. Each file needs chromosome, start, and end fields; any additional BED fields are ignored. The command writes one length per interval and a summary with counts, medians, minima, maxima, and the number of lengths exceeding 300 nt. It does not truncate the lengths to the plotted histogram range.

## `jaccard-blocks` and `call-jaccard`

Each BED file must have at least six columns: chromosome, start, end, name, score, and strand. Input BED call sets must already be restricted to the same evaluation regions. `jaccard-blocks` splits intervals at 1 Mb boundaries, merges overlaps within a call set, and counts intersecting and union bases separately by strand before combining both strands within a block. It writes `rbp`, `seed`, `block`, `model_intersection`, `model_union`, `original_intersection`, and `original_union`. Combine the tables for all RBP/seed combinations before running `call-jaccard`.

The Jaccard bootstrap keeps repeated block draws. The separate sign test uses the RBP-level mean across training seeds as its six observations. Seed-specific intervals are pointwise and are not a family-wise multiple-testing correction.

## `depth-scores` and `depth-ap`

`depth-scores` uses the raw NPZ generated by `thin-reads` and a prediction NPZ generated from that same sample. It writes `rbp`, `depth`, `sample`, `approach`, `window_id`, `label`, and `score` for test windows. `label=1` identifies original peak-derived parents; constructed hard-background and zero-IP parents have `label=0`. `depth-ap` groups this table by RBP, depth, sample, and approach and uses stable descending score order for ties.

## `motif-rates`

Use `motif-scan` on strand-aware interval tables to obtain `motif_hits`, `motif_present`, and `scan_nt`. The command needs a tab-delimited motif matrix with `A`, `C`, `G`, and `U` or `T` columns, an indexed reference FASTA, and the matching prepared training NPZ for the nucleotide background. Its default threshold is FIMO p ≤ 0.001 on transcript-oriented sequence. Run `motif-blocks` on combined scanned intervals to produce the required block table. The combined table must label each interval with `rbp`, `condition`, `approach`, and `seed`; use `approach=clipper, seed=0` for original calls, each included only once, and `approach=model` with seed 42, 1337, or 2026 for reconstructed calls. Overlapping intervals contribute their motif hits and lengths separately, as in the thesis.

One row represents one occupied genomic block for one RBP and one scanning condition. Required columns are `rbp`, `condition`, `block`, `clipper_hits`, `clipper_nt`, `model_hits_42`, `model_nt_42`, `model_hits_1337`, `model_nt_1337`, `model_hits_2026`, and `model_nt_2026`. Include zero-hit blocks with their scanned lengths. The command sums hits and lengths within each resample, calculates each seed's model-to-CLIPper rate ratio, and then averages the three ratios. It uses the logarithm of that mean ratio for the centred-bootstrap null test. Its percentile interval is on the mean-ratio scale.

The thesis conditions were unequal 5-prime extension, matched 50 nt extension, and matched total scanned span. For the matched-span condition, counts must come from actual threshold-selected call sets and their scanned sequences. Do not supply interpolated counts.

## `outside-peak-motifs`

Required columns are `rbp`, `chrom`, `start`, `selected`, and `motif_present`. Boolean fields can be `0` or `1`. Include only eligible windows outside original CLIPper calls: top-ranked selected windows and the remaining comparison windows. The command calculates the selected-minus-comparison difference in motif-positive proportions within each RBP and resamples the occupied 1 Mb blocks. Both groups retain their within-block relationship. The three RBP-specific p-values are adjusted together with Holm's method.

`rank-outside-peaks` produces the window groups before motif scanning. Its prepared NPZ must cover the designated chromosome 3 region on both strands. The two original CLIPper inputs are strand-aware BED files for the same RBP. Ranking occurs before the original-peak overlap filter, so `rank` is the position among *all* measured-IP-eligible windows, including windows later removed for overlapping a published peak.

