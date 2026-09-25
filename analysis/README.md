# Reproducing the generated-instance statistical analyses

This directory reproduces the generated-instance analyses reported in:

- the algorithm comparisons in Table 5;
- the factor-combination difficulty analysis underlying Figure 6; and
- the performance profiles in Figure 7.

## Requirements

The script was tested with R 4.5.2 and `NSM3` 1.20 and requires the CRAN package:

```r
install.packages("NSM3")
```

The Figure 7 reproduction additionally requires Python 3, pandas, Matplotlib, and a LaTeX installation. No package is installed automatically; the script stops with an explicit message when a dependency is unavailable.

## Run

From the repository root:

```bash
Rscript analysis/reproduce.R
```

An optional input CSV and output directory may be specified:

```bash
Rscript analysis/reproduce.R data/instance_generator/all_data.csv analysis/output
```

The input describes 10,500 algorithm runs. The analysis only ranks and aggregates these observations; it does not run any WSP solver or generate instances.

## Input data

`data/instance_generator/all_data.csv` contains the results of the individual algorithm runs. In addition to the instance identifier and parameter group, it exposes the eight experimental factors (`grid`, `slope`, `wind`, `delay`, `num_resources`, `num_decision_points`, `first_release_time`, and `last_release_time`), the algorithm, replication seed, objective value (`objv`), and lower bound (`lb`).

The default factor configuration occurs in each of the four parameter groups. The `parameter_group` column distinguishes these separately executed batches even though they share factor values and instance identifiers. The Table 5 and Figure 6 analyses use one copy of the shared default configuration, whereas the Figure 7 profiles preserve the group-specific batches used to produce the published panels. The experiments use ten algorithm replications for each generated instance. The time limits, platform, and other experimental settings are documented in the accompanying paper. Elapsed running times are not used in these statistical analyses and are not included in this cooked release.

## Outputs

The command creates these files in the selected output directory:

- `skillings_mack_summary.csv`: omnibus-test p-values, post-hoc thresholds, and design sizes;
- `table5_scores.csv` and `table5_comparisons.csv`: algorithm scores and pairwise post-hoc comparisons for Table 5;
- `difficulty_scores.csv` and `difficulty_comparisons.csv`: factor-combination scores and pairwise comparisons underlying Figure 6;
- `performance_ratios.csv` and `performance_profiles.csv`: the median objectives, performance ratios, and plotted profile points underlying Figure 7; and
- `pp_instance_size.pdf`, `pp_suppression_capacity.pdf`, `pp_environmental_factors.pdf`, and `pp_release_window.pdf`: the four panels of Figure 7.

Lower scores are better. A pairwise comparison is significant when its score difference exceeds `delta`, using the family-wise significance level \(\alpha=0.001\) specified in the paper.
