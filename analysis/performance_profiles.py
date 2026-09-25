#!/usr/bin/env python3
"""Reproduce the performance profiles reported in Figure 7."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ALGORITHMS = ["IBS", "MIP", "ILS", "LBBD", "RS"]
FACTOR_COLUMNS = [
    "grid",
    "slope",
    "wind",
    "delay",
    "num_resources",
    "num_decision_points",
    "first_release_time",
    "last_release_time",
]
GROUPS = [
    "instance_size",
    "suppression_capacity",
    "environmental_factors",
    "release_window",
]

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 12,
    }
)


def ecdf_points(values: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    taus, counts = np.unique(np.sort(values.to_numpy()), return_counts=True)
    return taus, np.cumsum(counts) / len(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    data = pd.read_csv(args.input)
    required = {"instance", "parameter_group", "algorithm", "seed", "objv", *FACTOR_COLUMNS}
    missing = required - set(data.columns)
    if missing:
        parser.error(f"input is missing required columns: {', '.join(sorted(missing))}")

    args.output.mkdir(parents=True, exist_ok=True)

    ratios = []
    profiles = []
    for group in GROUPS:
        selected = data[data["parameter_group"] == group].copy()
        medians = (
            selected.groupby(["instance", "algorithm"], as_index=False)["objv"]
            .median()
            .rename(columns={"objv": "median_objv"})
        )
        best = (
            medians.groupby("instance", as_index=False)["median_objv"]
            .min()
            .rename(columns={"median_objv": "best_median_objv"})
        )
        group_ratios = medians.merge(best, on="instance", validate="many_to_one")
        group_ratios["ratio"] = group_ratios["median_objv"] / group_ratios["best_median_objv"]
        group_ratios.insert(0, "group", group)
        ratios.append(group_ratios)

        plt.figure()
        for algorithm in ALGORITHMS:
            values = group_ratios.loc[group_ratios["algorithm"] == algorithm, "ratio"]
            if values.empty:
                raise ValueError(f"no observations for {algorithm} in {group}")
            taus, probabilities = ecdf_points(values)
            profiles.append(
                pd.DataFrame(
                    {
                        "group": group,
                        "algorithm": algorithm,
                        "tau": taus,
                        "proportion": probabilities,
                    }
                )
            )
            line, = plt.plot(taus, probabilities, label=algorithm)
            plt.plot(
                [taus[-1]],
                [probabilities[-1]],
                marker="*",
                markersize=8,
                color=line.get_color(),
                markeredgecolor="black",
                markeredgewidth=0.5,
                zorder=3,
            )

        if group == "instance_size":
            plt.legend(fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel(r"$\tau$", fontsize=15)
        plt.ylabel(r"$P_a(\tau)$", fontsize=15)
        plt.xlim(left=1.0)
        plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)
        plt.tight_layout()
        plt.savefig(args.output / f"pp_{group}.pdf")
        plt.close()

    pd.concat(ratios, ignore_index=True).to_csv(args.output / "performance_ratios.csv", index=False)
    pd.concat(profiles, ignore_index=True).to_csv(args.output / "performance_profiles.csv", index=False)


if __name__ == "__main__":
    main()
