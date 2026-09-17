# Reproducing the generated-instance statistical analyses

This directory reproduces the Skillings--Mack analyses reported for the generated instances:

- the algorithm comparisons in Table 5; and
- the factor-combination difficulty analysis underlying Figure 3.

## Requirements

The script was tested with R 4.5.2 and `NSM3` 1.20 and requires the CRAN package:

```r
install.packages("NSM3")
```

No package is installed automatically. The script stops with an explicit message if `NSM3` is unavailable.

## Run

From the repository root:

```bash
Rscript analysis/reproduce.R
```

An optional input CSV and output directory may be specified:

```bash
Rscript analysis/reproduce.R data/instance_generator/all_data.csv analysis/output
```

The input has 10,500 run-level observations. The analysis only ranks and aggregates these observations; it does not run any WSP solver or generate instances.

## Input data

`data/instance_generator/all_data.csv` contains the run-level results. In addition to the instance identifier, it exposes the eight experimental factors (`grid`, `slope`, `wind`, `delay`, `num_resources`, `num_decision_points`, `first_release_time`, and `last_release_time`), the algorithm, replication seed, objective value (`objv`), and lower bound (`lb`).

The default factor configuration occurs in more than one of the four parameter groups, so the aggregate file includes repeated `(instance, algorithm, seed)` records. The script retains the first such record, yielding the 9,750 observations used in the paper's normalized analysis. The experiments use ten algorithm replications for each generated instance. The time limits, platform, and other experimental settings are documented in the accompanying paper. Elapsed running times are not used in these statistical analyses and are not included in this cooked release.

## Outputs

The command creates these CSV files in the selected output directory:

- `skillings_mack_summary.csv`: omnibus-test p-values, post-hoc thresholds, and design sizes;
- `table5_scores.csv` and `table5_comparisons.csv`: algorithm scores and pairwise post-hoc comparisons for Table 5;
- `difficulty_scores.csv` and `difficulty_comparisons.csv`: factor-combination scores and pairwise comparisons underlying Figure 3.

Lower scores are better. A pairwise comparison is significant when its score difference exceeds `delta`, using the family-wise significance level \(\alpha=0.001\) specified in the paper.
