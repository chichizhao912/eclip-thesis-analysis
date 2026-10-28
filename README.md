# eCLIP reconstruction analyses

This companion package provides analyses for paired-replicate eCLIP reconstruction. It quantifies measured replicate disagreement, evaluates reconstructed profiles and candidate calls, examines ranking under reduced read depth, and measures motif support for prioritised regions.

Its commands turn strand-oriented coverage windows, model predictions, BED calls, and motif matrices into comparison tables with genomic-block uncertainty estimates. Published CLIPper peaks provide an internal reference call set for technical comparisons. The reconstruction model and candidate caller are available in the companion [eclip-reconstruct](https://github.com/chichizhao912/eclip-reconstruct) package.

## Install

Use Python 3.10 or newer. From this repository, install the analysis package and the published reconstruction tool in the same environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m pip install 'git+https://github.com/chichizhao912/eclip-reconstruct.git'
# Add the optional FIMO-compatible motif scanner when reproducing motif analyses:
python -m pip install -e '.[motifs]'
eclip-thesis-analysis --help
```

The tool repository prepares strand-oriented window NPZ files and generates model predictions and calls. This repository accepts those files, together with study-specific tables described in [Input tables](docs/input_tables.md). Outputs are tab-delimited text unless otherwise noted.

## Main analyses

The commands follow the Results workflow, from disagreement in the measured data through model development and evaluation. A *window* is an aligned 300 nt region; a *call* is a reported genomic interval. These answer different questions, even when they come from the same experiment.

| Results section | Question | Commands | Measurements |
| --- | --- | --- | --- |
| 3.1, measured IP coverage | How often do both replicates show coverage, and do their profiles agree more closely with each other than with SMInput? | `coverage-presence`, `diagnostic-contrasts`, `paired-tests` | Coverage presence, within-window Pearson correlation and EMD contrasts, genomic-block tests |
| 3.1, published CLIPper peaks | How much do the original IP1 and IP2 call sets overlap, and what signal remains in peaks called by only one replicate? | `peak-diagnostics`, `peak-shape-tests` | Same-strand peak overlap and base-pair Jaccard, summit displacement, partner-IP support, peak width, half-height span, and coverage concentration |
| 3.2, input-feature controls | How do SMInput and correctly ordered sequence affect reconstruction? | `train-controls`, `reconstruction-metrics`, `reconstruction-differences`, `paired-tests` | Four input configurations, P1–P2 and prediction–target Pearson correlations, paired block sign-flip tests |
| 3.2, Gaussian comparator and local shape | Does reconstruction retain the measured profile structure better than variance-matched smoothing? | `gaussian-comparison`, `gaussian-summary`, `gaussian-differences`, `paired-tests` | Correlations, positional-variance retention, log-total signal error, and dominant-feature width at half maximum |
| 3.2, candidate-calling development | Does the SMInput penalty improve detection at a matched background-calling rate, and how do refined boundaries change interval lengths? | `caller-score-comparison`, `candidate-lengths`, plus `calibrate` and `call` in [eclip-reconstruct](https://github.com/chichizhao912/eclip-reconstruct) | Separately calibrated score comparisons, refined candidate boundaries, and length distributions |
| 3.3, technical evaluation | Do directional calls agree more closely, and how does candidate ranking change at reduced read depth? | `jaccard-blocks`, `call-jaccard`, `thin-reads`, `depth-scores`, `depth-ap` | Base-pair Jaccard differences, matched read-thinning profiles, and average precision |
| 3.3, biological support | Do calls recover expected motifs, and do ranked windows outside published peaks contain more motifs? | `motif-scan`, `motif-blocks`, `motif-rates`, `rank-outside-peaks`, `outside-peak-motifs` | Motif hits and scanned length, block-resampled rate ratios, and motif-positive window proportions |

Examples below show the command structure. The complete field definitions and required inputs are in [Input tables](docs/input_tables.md).

### Diagnostic coverage contrasts

Use the same prepared peak-derived window population for each RBP when counting measurable coverage and measuring profile agreement. `coverage-presence` counts windows with positive coverage in both IPs, in just one IP, or in neither IP.

```bash
eclip-thesis-analysis coverage-presence --data outputs/tardbp_windows.npz \
  --rbp TARDBP --output outputs/tardbp_coverage_presence.tsv
```

```bash
eclip-thesis-analysis diagnostic-contrasts --data outputs/tardbp_windows.npz \
  --rbp TARDBP --output outputs/tardbp_diagnostic.tsv
```

Combine the six RBP tables without changing their `family` labels. The prespecified median contrasts are then tested within separate six-test Pearson and EMD families:

```bash
eclip-thesis-analysis paired-tests --input outputs/all_diagnostic.tsv \
  --method centred-bootstrap --statistic median \
  --bootstrap-resamples 2000 --stratify-chromosome \
  --output outputs/diagnostic_tests.tsv
```

### Published peak calls and local signal

`peak-diagnostics` reads the complete replicate-specific narrowPeak/BED call sets and the corresponding strand-specific IP and SMInput bigWigs. It writes four tables: a per-IP overlap summary, a sample of peak-shape measurements, partner-replicate support for replicate-specific peaks, and matched summit displacements. The shape sample uses up to 5,000 shared and 5,000 replicate-specific peaks per IP, sampled before excluding peaks with no positive calling-IP coverage.

```bash
eclip-thesis-analysis peak-diagnostics --rbp TARDBP \
  --ip1-peaks ip1.narrowPeak --ip2-peaks ip2.narrowPeak \
  --ip1-plus ip1.plus.bw --ip1-minus ip1.minus.bw \
  --ip2-plus ip2.plus.bw --ip2-minus ip2.minus.bw \
  --sminput-plus sminput.plus.bw --sminput-minus sminput.minus.bw \
  --output-dir outputs/peak_diagnostics
```

Combine the shape tables from all six RBPs before running `peak-shape-tests`. It compares replicate-specific with shared peaks for each IP, then applies Holm correction separately to the 12 RBP-by-IP comparisons for each of the three measurements.

```bash
eclip-thesis-analysis peak-shape-tests --input outputs/all_peak_shapes.tsv \
  --output outputs/peak_shape_tests.tsv
```

### Reconstruction comparisons

`train-controls` fits the four prespecified TARDBP configurations on the same prepared windows and two seeds. The `shuffled` configuration uses a reproducible within-window Eulerian dinucleotide shuffle, held fixed across seeds. It writes checkpoints and prediction NPZ files for each configuration. This command requires `eclip-reconstruct`, installed as shown above.

```bash
eclip-thesis-analysis train-controls --data outputs/tardbp_windows.npz \
  --output outputs/input_controls
```

For each configuration and training seed, collect correlations from its prediction NPZ and the matching prepared-window NPZ. The prepared-window data must be the same for all compared configurations.

```bash
eclip-thesis-analysis reconstruction-metrics --data outputs/windows.npz \
  --predictions outputs/input_controls/full_seed42_predictions.npz \
  --preprocessing outputs/input_controls/preprocessing.json \
  --configuration full --seed 42 --output outputs/full_seed42_metrics.tsv
```

Combine all configuration/seed metrics tables, then form full-minus-alternative paired differences. The command averages the two seed-specific measurements **within each window** and keeps windows with complete results for every requested configuration and both metrics:

```bash
eclip-thesis-analysis reconstruction-differences \
  --input outputs/all_configuration_metrics.tsv \
  --reference full --alternatives coverage shuffled no_sminput \
  --family input_ablation --output outputs/ablation_differences.tsv
```

The resulting six comparisons form one declared family and are tested with:

```bash
eclip-thesis-analysis paired-tests --input outputs/ablation_differences.tsv \
  --method sign-flip --statistic mean \
  --bootstrap-resamples 3000 --test-resamples 20000 \
  --output outputs/ablation_tests.tsv
```

The sign-flip test treats each occupied 1 Mb block as one sign unit. Holm adjustment is applied within the `family` column. The same command can assess other prespecified paired comparisons by supplying their window-level differences and declaring the correct family.

### Gaussian comparator and local shape

For the chromosome 1 historical comparison, use the fixed Gaussian standard deviation of 26.19 nt selected on chromosome 2. Run `gaussian-comparison` on the same prepared windows for each of the two model seeds, then combine its tables. The Gaussian prediction for IP1 smooths measured IP2, and the prediction for IP2 smooths measured IP1, preserving the model's cross-replicate comparison.

```bash
eclip-thesis-analysis gaussian-comparison --data outputs/chr1_windows.npz \
  --predictions outputs/chr1_seed42_predictions.npz \
  --preprocessing outputs/preprocessing.json --seed 42 \
  --output outputs/chr1_seed42_gaussian.tsv

eclip-thesis-analysis gaussian-summary --input outputs/chr1_both_seeds_gaussian.tsv \
  --output outputs/chr1_gaussian_summary.tsv
```

The summary gives P1–P2 and prediction–target correlations, positional-variance retention, log-total signal error, and the dominant-feature width of observed, smoothed, and reconstructed consensus profiles. For width, the two model-seed consensus profiles are averaged **before** measuring width. `gaussian-differences` prepares paired window-level comparisons for `paired-tests`: use `--comparison correlation --statistic mean` for the two correlation tests and `--comparison width --statistic median` for the two width tests. Run these as separate Holm families, with `--method sign-flip --bootstrap-resamples 3000 --test-resamples 20000`.

### Candidate score comparison

`caller-score-comparison` scans the reconstructed consensus in 100 nt intervals at 5 nt steps. It compares the unadjusted consensus mean with consensus mean minus 0.25 times mean transformed SMInput. Each score gets its **own** threshold at the higher empirical 95th percentile of validation hard-background parent maxima. The output gives peak-derived recovery, hard-background and zero-IP calling percentages, and average precision for one training seed. Run it for both seeds and average the two rows for each score scheme when reproducing the development comparison.

```bash
eclip-thesis-analysis caller-score-comparison \
  --predictions outputs/chr2_seed42_predictions.npz --seed 42 \
  --output outputs/chr2_seed42_score_comparison.tsv
```

Use the selected SMInput-adjusted score with the tool repository's `calibrate` and `call` commands to produce refined intervals.

`candidate-lengths` compares the saved original and refined BED call sets from the same boundary-development run. It reports every interval length, the candidate count, and the median for each procedure; lengths above 300 nt remain in the output and are counted explicitly.

```bash
eclip-thesis-analysis candidate-lengths --original outputs/original_calls.bed \
  --refined outputs/refined_calls.bed --output-dir outputs/boundary_lengths
```

### Directional call agreement

Run the tool caller separately on `p1` and `p2` for each seed. Restrict both original CLIPper and model call BED files to the same evaluation population before constructing block counts.

```bash
eclip-thesis-analysis jaccard-blocks --rbp TARDBP --seed 42 \
  --original-ip1 original_ip1.bed --original-ip2 original_ip2.bed \
  --model-p1 model_p1.bed --model-p2 model_p2.bed \
  --output outputs/tardbp_seed42_blocks.tsv
```

Combine all protein and seed block tables, then run `call-jaccard`. It writes seed-specific pointwise 95% genomic-block intervals and a separate two-sided exact sign-test result over the six protein-level mean differences.

### Read thinning and ranking

`thin-reads` uses aligned blocks of the transcript-orienting second read in each paired-end BAM. Before generating any depth samples, it checks that full-depth reconstructed coverage matches the supplied bigWig-derived NPZ profiles with pooled Pearson correlation at least 0.999 for each track. This gate must pass for the input libraries. Reads are sampled independently in each parent window and sampled coverage is divided by the retained fraction before applying a fixed track scaling factor.

```bash
eclip-thesis-analysis thin-reads --data outputs/windows.npz \
  --ip1-bam ip1.bam --ip2-bam ip2.bam --sminput-bam sminput.bam \
  --output outputs/thinning
```

The default fractions are 1.00, 0.75, 0.50, 0.30, 0.20, 0.10, and 0.05, with five samples per fraction. Predict each generated NPZ using the relevant model checkpoint. For a matched thinned-data/prediction pair, `depth-scores` calculates the highest 100 nt score for the mean original IP coverage and the reconstructed consensus. Combine its outputs and run `depth-ap` to calculate AP by protein, depth, sample, and approach. AP describes ranking against the fixed CLIPper-derived versus constructed-background labels.

### Motif comparisons

For the exploratory tiled search, prepare non-overlapping 300 nt windows on both strands in the prespecified chromosome 3 region and predict with the three-seed ensemble. `rank-outside-peaks` applies the original measured-IP coverage filter, ranks *all* eligible windows by the highest internal 100 nt score, selects the top 1,000, and only then removes overlaps with either original same-strand CLIPper peak set. It labels the remaining CLIPper-negative top-ranked and comparison windows for motif scanning.

```bash
eclip-thesis-analysis rank-outside-peaks --data outputs/tiled_windows.npz \
  --predictions outputs/tiled_ensemble.npz \
  --original-ip1 original_ip1.bed --original-ip2 original_ip2.bed \
  --rbp TARDBP --output outputs/tardbp_ranked_windows.tsv
```

`motif-scan` takes an interval table, indexed reference FASTA, the RBP's prespecified position-frequency matrix, and the prepared training windows. It applies matrix rounding and flank trimming, estimates the zero-order nucleotide background from training sequences, and counts transcript-oriented hits passing the one-strand FIMO threshold of 0.001. The matrix table needs `A`, `C`, `G`, and `U` (or `T`) columns. For the motif-rate analysis, label each interval with `rbp`, `condition`, `approach` (`clipper` or `model`), and `seed` (use `0` for CLIPper). Scan the intervals, combine the scanned tables, then run `motif-blocks` and `motif-rates`. Include the original CLIPper intervals once, rather than repeating them for each model seed. For the outside-peak analysis, scan the `rank-outside-peaks` output and pass its `selected` and `motif_present` columns to `outside-peak-motifs`.

```bash
eclip-thesis-analysis motif-scan --input outputs/tardbp_ranked_windows.tsv \
  --reference GRCh38.fa --matrix tardbp_motif.tsv \
  --training-data outputs/tardbp_windows.npz --motif-id M00096_2.00 \
  --output outputs/tardbp_ranked_motifs.tsv
```

The motif-rate analysis keeps hits and scanned lengths together within genomic blocks and retains the three model seeds as repeated fits. `motif-rates` adjusts the 18 planned protein-by-scanning-condition p-values with Holm's method. `outside-peak-motifs` adjusts the three planned protein comparisons with Holm's method. Both default to 50,000 block resamples. Motif presence supports biological plausibility but does not independently verify binding, since the model also received reference sequence.

## Reproducibility

The resampling code samples occupied 1 Mb blocks with replacement and retains a repeated block's observations *each time it is drawn*. Paired measurements, model seeds, and group labels stay together within a block. Model-development comparisons use two-sided block sign flips of within-window differences. Diagnostic and motif analyses use two-sided centred block-bootstrap tests. The two tests construct their null distributions differently and are kept separate in the code.

The package focuses on numerical comparisons and table generation from traceable input files. The required table formats and analysis populations are documented in [Input tables](docs/input_tables.md), allowing the same workflows to be applied to appropriately prepared eCLIP datasets.

